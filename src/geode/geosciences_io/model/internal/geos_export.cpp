/*
 * Copyright (c) 2019 - 2026 Geode-solutions
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>

#include <pugixml.hpp>

#include <absl/algorithm/container.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_join.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/aabb.hpp>
#include <geode/geometry/distance.hpp>

#include <geode/geosciences_io/model/common.hpp>

#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/polyhedral_solid.hpp>
#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/surface_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/helpers/aabb_solid_helpers.hpp>

#include <geode/mesh/io/hybrid_solid_output.hpp>
#include <geode/mesh/io/point_set_output.hpp>
#include <geode/mesh/io/polyhedral_solid_output.hpp>
#include <geode/mesh/io/tetrahedral_solid_output.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences_io/model/internal/geos_export.hpp>
#include <geode/model/representation/core/brep.hpp>

namespace geode
{
    namespace internal
    {
        constexpr auto REGION_ID_ATTRIBUTE_NAME = "attribute";
        constexpr auto PERMEABILITY_GEOS_NAME = "rockPerm_permeability";
        constexpr auto POROSITY_GEOS_NAME = "rockPorosity_referencePorosity";

        template < typename T, typename Model >
        void transfer_block_attribute_by_id( const Model& model,
            SolidMesh3D& solid,
            const ModelToMeshMappings& model2solid,
            std::string_view attribute_name,
            const uuid& attribute_id )
        {
            AttributeValues< T > solid_attribute_values;
            solid_attribute_values.default_value = T{};
            solid_attribute_values.no_value = T{};
            AttributeProperties solid_attribute_properties;
            solid_attribute_properties.assignable = false;
            solid_attribute_properties.interpolable = false;
            solid_attribute_properties.transferable = true;
            const auto solid_property_id =
                solid.polyhedron_attribute_manager()
                    .template create_attribute< VariableAttribute, T >(
                        attribute_name, solid_attribute_values,
                        solid_attribute_properties );
            auto solid_property =
                solid.polyhedron_attribute_manager()
                    .template find_attribute< VariableAttribute, T >(
                        solid_property_id );
            const auto& polyhedra_mapping = model2solid.solid_polyhedra_mapping;
            for( const auto& block : model.blocks() )
            {
                const auto& block_mesh = block.mesh();
                const auto block_property =
                    block_mesh.polyhedron_attribute_manager()
                        .template find_read_only_attribute< T >( attribute_id );
                for( const auto polyhedron_id :
                    Range( block_mesh.nb_polyhedra() ) )
                {
                    const MeshElement block_polyhedron{ block.id(),
                        polyhedron_id };
                    if( !polyhedra_mapping.has_mapping_input(
                            block_polyhedron ) )
                    {
                        continue;
                    }
                    const auto& value = block_property->value( polyhedron_id );
                    for( const auto solid_polyhedron_id :
                        polyhedra_mapping.in2out( block_polyhedron ) )
                    {
                        solid_property->set_value( solid_polyhedron_id, value );
                    }
                }
            }
        }

        template < typename Model >
        std::string transfer_block_attribute_by_id( const Model& model,
            SolidMesh3D& solid,
            const ModelToMeshMappings& model2solid,
            const PhysicalProperties::Info& property_info )
        {
            const auto attribute_id = property_info.attribute_id;
            OpenGeodeGeosciencesIOModelException::check_exception(
                property_info.component_type
                    == Block3D::component_type_static(),
                nullptr, OpenGeodeException::TYPE::data,
                "[GeosExporter] Physical property attribute ",
                attribute_id.string(),
                " must be defined on Blocks to be exported." );
            std::optional< std::string > attribute_name;
            local_index_t nb_items{ 0 };
            for( const auto& block : model.blocks() )
            {
                const auto& manager =
                    block.mesh().polyhedron_attribute_manager();
                OpenGeodeGeosciencesIOModelException::check_exception(
                    manager.attribute_exists( attribute_id ), nullptr,
                    OpenGeodeException::TYPE::data,
                    "[GeosExporter] Physical property attribute ",
                    attribute_id.string(), " is missing on Block ",
                    block.id().string(), "." );
                if( !attribute_name )
                {
                    const auto attribute =
                        manager.find_generic_attribute( attribute_id );
                    attribute_name = attribute->name();
                    nb_items = attribute->nb_items();
                }
            }
            OpenGeodeGeosciencesIOModelException::check_exception(
                attribute_name.has_value(), nullptr,
                OpenGeodeException::TYPE::data,
                "[GeosExporter] Physical property attribute ",
                attribute_id.string(), " has no name." );
            if( nb_items == 1 )
            {
                transfer_block_attribute_by_id< double >( model, solid,
                    model2solid, attribute_name.value(), attribute_id );
            }
            else if( nb_items == 3 )
            {
                transfer_block_attribute_by_id< std::array< double, 3 > >(
                    model, solid, model2solid, attribute_name.value(),
                    attribute_id );
            }
            return attribute_name.value();
        }

        template < typename Model >
        GeosExporterImpl< Model >::GeosExporterImpl(
            std::string_view files_directory, const Model& model )
            : model_( model ),
              files_directory_{
                  std::filesystem::path{ to_string( files_directory ) }.string()
              },
              prefix_{ filename_without_extension( files_directory ).string() }
        {
            auto curve_conversion = convert_brep_into_curve( model_ );
            model_curve_ = std::move( std::get< 0 >( curve_conversion ) );
            auto surface_conversion = convert_brep_into_surface( model_ );
            model_surface_ = std::move( std::get< 0 >( surface_conversion ) );
            std::tie( model_solid_, model2solid_ ) =
                convert_brep_into_solid( model_ );
            AttributeValues< index_t > region_attribute_values;
            region_attribute_values.default_value = NO_ID;
            region_attribute_values.no_value = NO_ID;
            AttributeProperties region_attribute_properties;
            region_attribute_properties.assignable = false;
            region_attribute_properties.interpolable = false;
            region_attribute_properties.transferable = true;
            const auto region_attribute_id =
                model_solid_->polyhedron_attribute_manager()
                    .template create_attribute< VariableAttribute, index_t >(
                        REGION_ID_ATTRIBUTE_NAME, region_attribute_values,
                        region_attribute_properties );
            region_attribute_ =
                model_solid_->polyhedron_attribute_manager()
                    .find_attribute< VariableAttribute, index_t >(
                        region_attribute_id );
            if( std::filesystem::path{ to_string( files_directory ) }
                    .is_relative() )
            {
                std::filesystem::create_directory(
                    std::filesystem::current_path() / files_directory_ );
            }
            else
            {
                std::filesystem::create_directory( files_directory_ );
            }
        }

        template < typename Model >
        void GeosExporterImpl< Model >::write_files() const
        {
            pugi::xml_document doc_xml;
            auto pb_node = doc_xml.append_child( "Problem" );

            auto mesh_node = pb_node.append_child( "Mesh" );
            write_mesh_files( mesh_node );

            if( !well_perforations_.empty() )
            {
                auto geometry_node = pb_node.append_child( "Geometry" );
                write_well_perforations_boxes( geometry_node );
                write_well_perforation_file();
            }

            const auto filename_xml = absl::StrCat(
                files_directory(), "/", prefix(), "_simulation.xml" );
            doc_xml.save_file( filename_xml.c_str(), PUGIXML_TEXT( "    " ),
                pugi::format_indent_attributes );
        }

        template < typename Model >
        void GeosExporterImpl< Model >::add_well_perforations(
            const PointSet3D& perforations )
        {
            add_well_perforations( perforations,
                absl::StrCat( "well_", well_perforations_.size() ) );
        }

        template < typename Model >
        void GeosExporterImpl< Model >::add_well_perforations(
            const PointSet3D& perforations, std::string_view name )
        {
            OpenGeodeGeosciencesIOModelException::check_exception(
                absl::c_find( well_names_, name ) == well_names_.end(), nullptr,
                OpenGeodeException::TYPE::data,
                "[GeosExporter] Well perforations named ", name,
                " already added." );
            well_perforations_.push_back( perforations.clone() );
            well_names_.emplace_back( to_string( name ) );
        }

        template < typename Model >
        void GeosExporterImpl< Model >::prepare_export()
        {
            initialize_solid_region_attribute();
            transfer_physical_properties();
        }

        template < typename Model >
        std::string_view GeosExporterImpl< Model >::files_directory() const
        {
            return files_directory_;
        }
        template < typename Model >
        std::string_view GeosExporterImpl< Model >::prefix() const
        {
            return prefix_;
        }

        template < typename Model >
        index_t GeosExporterImpl< Model >::initialize_solid_region_attribute()
        {
            auto region_map_id = create_region_attribute_map( model_ );
            for( const auto polyhedron_id :
                Range( model_solid_->nb_polyhedra() ) )
            {
                region_attribute_->set_value( polyhedron_id,
                    region_map_id
                        .find( model2solid_.solid_polyhedra_mapping
                                .out2in( polyhedron_id )
                                .front()
                                .mesh_id )
                        ->second );
            }
            return region_map_id.size();
        }

        template < typename Model >
        void GeosExporterImpl< Model >::write_well_perforations_boxes(
            pugi::xml_node& root ) const
        {
            auto aabb = create_aabb_tree( *model_solid_ );
            for( const auto well_id : Indices{ well_perforations_ } )
            {
                const auto& well = well_perforations_[well_id];
                BoundingBox3D perf_box;
                for( const auto point : Range( well->nb_vertices() ) )
                {
                    const auto neigh_cells =
                        aabb.containing_boxes( well->point( point ) );
                    if( neigh_cells.empty() )
                    {
                        continue;
                    }
                    double distance_to_nearest_cell_center{
                        std::numeric_limits< double >::max()
                    };
                    index_t selected_cell_id{ NO_ID };
                    for( const auto cell_id : neigh_cells )
                    {
                        const auto tmp_dist = point_point_distance< 3 >(
                            well->point( point ),
                            model_solid_->polyhedron_barycenter( cell_id ) );
                        if( distance_to_nearest_cell_center < tmp_dist )
                        {
                            continue;
                        }
                        distance_to_nearest_cell_center = tmp_dist;
                        selected_cell_id = cell_id;
                    }
                    for( auto& vertex_id :
                        model_solid_->polyhedron_vertices( selected_cell_id ) )
                    {
                        perf_box.add_point( model_solid_->point( vertex_id ) );
                    }
                }
                auto box_node = root.append_child( "Box" );
                box_node.append_attribute( "name" ).set_value(
                    well_names_[well_id].c_str() );
                static constexpr auto SAFETY_OFFSET = 100. * GLOBAL_EPSILON;
                box_node.append_attribute( "xMin" ).set_value( absl::StrCat(
                    "{", perf_box.min().value( 0 ) - SAFETY_OFFSET, ", ",
                    perf_box.min().value( 1 ) - SAFETY_OFFSET, ", ",
                    perf_box.min().value( 2 ) - SAFETY_OFFSET, "}" )
                        .c_str() );
                box_node.append_attribute( "xMax" ).set_value( absl::StrCat(
                    "{", perf_box.max().value( 0 ) + SAFETY_OFFSET, ", ",
                    perf_box.max().value( 1 ) + SAFETY_OFFSET, ", ",
                    perf_box.max().value( 2 ) + SAFETY_OFFSET, "}" )
                        .c_str() );
            }
        }

        template < typename Model >
        void GeosExporterImpl< Model >::write_mesh_files(
            pugi::xml_node& root ) const
        {
            const auto file_vtu = write_solid_file();

            auto vtk_mesh_node = root.append_child( "VTKMesh" );
            vtk_mesh_node.append_attribute( "name" ).set_value(
                to_string( prefix() ).c_str() );
            vtk_mesh_node.append_attribute( "file" ).set_value(
                absl::StrCat( "./", file_vtu ).c_str() );
            if( imported_fields_.empty() )
            {
                return;
            }
            const auto join_fields = [this]( bool vtu_names ) {
                return absl::StrCat( "{ ",
                    absl::StrJoin( imported_fields_, ", ",
                        [vtu_names]( std::string* out, const auto& field ) {
                            absl::StrAppend(
                                out, vtu_names ? field.first : field.second );
                        } ),
                    " }" );
            };
            vtk_mesh_node.append_attribute( "fieldsToImport" )
                .set_value( join_fields( true ).c_str() );
            vtk_mesh_node.append_attribute( "fieldNamesInGEOS" )
                .set_value( join_fields( false ).c_str() );
        }

        template < typename Model >
        bool GeosExporterImpl< Model >::check_property_name(
            std::string_view property_name ) const
        {
            for( const auto& block : model_.blocks() )
            {
                if( !block.mesh()
                        .polyhedron_attribute_manager()
                        .attribute_ids_matching_name( property_name )
                        .has_value() )
                {
                    Logger::info( "The property ", property_name,
                        " will not be exported because it is not defined on "
                        "every block of the model." );
                    return false;
                }
            }
            return true;
        }

        template < typename Model >
        void GeosExporterImpl< Model >::transfer_physical_properties()
        {
            if( model_.has_physical_property(
                    PHYSICAL_PROPERTY_NAME::permeability ) )
            {
                auto attribute_name = transfer_block_attribute_by_id( model_,
                    *model_solid_, model2solid_,
                    model_.physical_property_info(
                        PHYSICAL_PROPERTY_NAME::permeability ) );
                imported_fields_.emplace_back(
                    std::move( attribute_name ), PERMEABILITY_GEOS_NAME );
            }
            if( model_.has_physical_property(
                    PHYSICAL_PROPERTY_NAME::porosity ) )
            {
                auto attribute_name = transfer_block_attribute_by_id( model_,
                    *model_solid_, model2solid_,
                    model_.physical_property_info(
                        PHYSICAL_PROPERTY_NAME::porosity ) );
                imported_fields_.emplace_back(
                    std::move( attribute_name ), POROSITY_GEOS_NAME );
            }
        }

        template < typename Model >
        std::string GeosExporterImpl< Model >::write_solid_file() const
        {
            const auto filename = absl::StrCat( prefix(), ".vtu" );
            const auto file_vtu =
                absl::StrCat( files_directory(), "/", filename );
            if( const auto* tetra = dynamic_cast< const TetrahedralSolid3D* >(
                    model_solid_.get() ) )
            {
                save_tetrahedral_solid( *tetra, file_vtu );
            }
            else if( const auto* hybrid = dynamic_cast< const HybridSolid3D* >(
                         model_solid_.get() ) )
            {
                save_hybrid_solid( *hybrid, file_vtu );
            }
            else if( const auto* poly =
                         dynamic_cast< const PolyhedralSolid3D* >(
                             model_solid_.get() ) )
            {
                save_polyhedral_solid( *poly, file_vtu );
            }
            else
            {
                throw OpenGeodeGeosciencesIOModelException{ nullptr,
                    OpenGeodeException::TYPE::data,
                    "[Blocks::save_geos] Cannot find the explicit "
                    "SolidMesh type" };
            }
            return filename;
        }
        template < typename Model >
        void GeosExporterImpl< Model >::write_well_perforation_file() const
        {
            index_t well_id{ 0 };
            for( const auto& well : well_perforations_ )
            {
                const auto file = absl::StrCat( files_directory(), "/",
                    prefix(), "_well", well_id, ".vtp" );
                save_point_set( *well, file );
                well_id++;
            }
        }
        template class opengeode_geosciencesio_model_api
            GeosExporterImpl< BRep >;
        template class opengeode_geosciencesio_model_api
            GeosExporterImpl< StructuralModel >;
    } // namespace internal
} // namespace geode