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

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

#include <pugixml.hpp>

#include <absl/algorithm/container.h>
#include <absl/container/linked_hash_map.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_join.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/mapping.hpp>
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

#include <geode/mesh/io/point_set_output.hpp>

#include <geode/model/helpers/component_mesh_polygons.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>
#include <geode/model/mixin/core/surface.hpp>

#include <geode/io/mesh/detail/vtu_solid_output_impl.hpp>

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences_io/model/internal/geos_export.hpp>
#include <geode/model/representation/core/brep.hpp>

namespace
{
    // Writes the solid polyhedra followed by the surface polygons in a single
    // UnstructuredGrid, all the cells sharing the solid vertices
    template < typename SolidOutputImpl >
    class GeosVTUOutputImpl : public SolidOutputImpl
    {
    public:
        template < typename Solid >
        GeosVTUOutputImpl( std::string_view filename,
            const Solid& solid,
            absl::Span< const geode::PolygonVertices > surface_cells,
            const geode::AttributeManager& cell_attributes )
            : SolidOutputImpl{ filename, solid },
              surface_cells_{ surface_cells },
              cell_attributes_( cell_attributes )
        {
        }

    private:
        [[nodiscard]] geode::index_t nb_additional_polygons() const override
        {
            return surface_cells_.size();
        }

        [[nodiscard]] absl::Span< const geode::index_t >
            additional_polygon_vertices(
                geode::index_t polygon_id ) const override
        {
            return surface_cells_[polygon_id];
        }

        [[nodiscard]] const geode::AttributeManager&
            cell_attribute_manager() const override
        {
            return cell_attributes_;
        }

        absl::Span< const geode::PolygonVertices > surface_cells_;
        const geode::AttributeManager& cell_attributes_;
    };

    template < typename SolidOutputImpl, typename Solid >
    void save_geos_vtu( std::string_view filename,
        const Solid& solid,
        absl::Span< const geode::PolygonVertices > surface_cells,
        const geode::AttributeManager& cell_attributes )
    {
        GeosVTUOutputImpl< SolidOutputImpl > writer{ filename, solid,
            surface_cells, cell_attributes };
        writer.write_file();
    }
} // namespace

namespace geode
{
    namespace internal
    {
        constexpr auto REGION_ID_ATTRIBUTE_NAME = "attribute";

        const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME, std::string_view >
            PHYSICAL_PROPERTY_GEOS_NAMES{
                { PHYSICAL_PROPERTY_NAME::permeability,
                    "rockPerm_permeability" },
                { PHYSICAL_PROPERTY_NAME::porosity,
                    "rockPorosity_referencePorosity" },
            };

        template < typename Model >
        std::string transfer_block_attribute( const Model& model,
            SolidMesh3D& solid,
            const ModelToMeshMappings& model2solid,
            const PhysicalProperties::Info& property_info )
        {
            const auto& attribute_id = property_info.attribute_id;
            OpenGeodeGeosciencesIOModelException::check_exception(
                property_info.component_type
                    == Block3D::component_type_static(),
                nullptr, OpenGeodeException::TYPE::data,
                "[GeosExporter] Physical property attribute ",
                attribute_id.string(),
                " must be defined on Blocks to be exported." );
            auto& solid_manager = solid.polyhedron_attribute_manager();
            const auto& polyhedra_mapping = model2solid.solid_polyhedra_mapping;
            for( const auto& block : model.blocks() )
            {
                const auto& block_mesh = block.mesh();
                const auto& block_manager =
                    block_mesh.polyhedron_attribute_manager();
                OpenGeodeGeosciencesIOModelException::check_exception(
                    block_manager.attribute_exists( attribute_id ), nullptr,
                    OpenGeodeException::TYPE::data,
                    "[GeosExporter] Physical property attribute ",
                    attribute_id.string(), " is missing on Block ",
                    block.id().string(), "." );
                GenericMapping< index_t > block2solid;
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
                    for( const auto solid_polyhedron_id :
                        polyhedra_mapping.in2out( block_polyhedron ) )
                    {
                        block2solid.map( polyhedron_id, solid_polyhedron_id );
                    }
                }
                solid_manager.import(
                    block_manager, block2solid, attribute_id );
            }
            const auto attribute =
                solid_manager.find_generic_attribute( attribute_id );
            OpenGeodeGeosciencesIOModelException::check_exception(
                attribute && attribute->name().has_value(), nullptr,
                OpenGeodeException::TYPE::data,
                "[GeosExporter] Physical property attribute ",
                attribute_id.string(), " has no name." );
            return attribute->name().value();
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
            std::tie( model_solid_, model2solid_ ) =
                convert_brep_into_solid( model_ );
            AttributeValues< index_t > region_attribute_values;
            region_attribute_values.default_value = NO_ID;
            region_attribute_values.no_value = NO_ID;
            AttributeProperties region_attribute_properties;
            region_attribute_properties.assignable = false;
            region_attribute_properties.interpolable = false;
            region_attribute_properties.transferable = true;
            region_attribute_id_ =
                model_solid_->polyhedron_attribute_manager()
                    .template create_attribute< VariableAttribute, index_t >(
                        REGION_ID_ATTRIBUTE_NAME, region_attribute_values,
                        region_attribute_properties );
            region_attribute_ =
                model_solid_->polyhedron_attribute_manager()
                    .find_attribute< VariableAttribute, index_t >(
                        region_attribute_id_ );
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
            const auto nb_solid_regions = initialize_solid_region_attribute();
            initialize_surface_cells( nb_solid_regions );
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
            index_t nb_regions{ 0 };
            for( const auto& [block_id, region_id] : region_map_id )
            {
                nb_regions = std::max( nb_regions, region_id + 1 );
            }
            return nb_regions;
        }

        template < typename Model >
        void GeosExporterImpl< Model >::initialize_surface_cells(
            index_t first_surface_region_id )
        {
            auto region_id = first_surface_region_id;
            for( const auto& surface : model_.surfaces() )
            {
                for( const auto polygon_id :
                    Range{ surface.mesh().nb_polygons() } )
                {
                    auto vertices =
                        polygon_unique_vertices( model_, surface, polygon_id );
                    for( auto& vertex : vertices )
                    {
                        vertex = model2solid_.unique_vertices_mapping.in2out(
                            vertex );
                    }
                    surface_cells_.emplace_back( std::move( vertices ) );
                    surface_cells_region_.push_back( region_id );
                }
                region_id++;
            }
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
            for( const auto& [property, geos_name] :
                PHYSICAL_PROPERTY_GEOS_NAMES )
            {
                if( !model_.has_physical_property( property ) )
                {
                    continue;
                }
                auto attribute_name = transfer_block_attribute( model_,
                    *model_solid_, model2solid_,
                    model_.physical_property_info( property ) );
                imported_fields_.emplace_back(
                    std::move( attribute_name ), to_string( geos_name ) );
            }
        }

        template < typename Model >
        std::string GeosExporterImpl< Model >::write_solid_file() const
        {
            const auto filename = absl::StrCat( prefix(), ".vtu" );
            const auto file_vtu =
                absl::StrCat( files_directory(), "/", filename );
            const auto nb_polyhedra = model_solid_->nb_polyhedra();
            AttributeManager cell_attributes;
            cell_attributes.copy(
                model_solid_->polyhedron_attribute_manager() );
            cell_attributes.resize( nb_polyhedra + surface_cells_.size() );
            auto region_attribute =
                cell_attributes.find_attribute< VariableAttribute, index_t >(
                    region_attribute_id_ );
            for( const auto surface_cell : Indices{ surface_cells_region_ } )
            {
                region_attribute->set_value( nb_polyhedra + surface_cell,
                    surface_cells_region_[surface_cell] );
            }
            if( const auto* tetra = dynamic_cast< const TetrahedralSolid3D* >(
                    model_solid_.get() ) )
            {
                save_geos_vtu< detail::VTUTetrahedralOutputImpl >(
                    file_vtu, *tetra, surface_cells_, cell_attributes );
            }
            else if( const auto* hybrid = dynamic_cast< const HybridSolid3D* >(
                         model_solid_.get() ) )
            {
                save_geos_vtu< detail::VTUHybridOutputImpl >(
                    file_vtu, *hybrid, surface_cells_, cell_attributes );
            }
            else if( const auto* poly =
                         dynamic_cast< const PolyhedralSolid3D* >(
                             model_solid_.get() ) )
            {
                save_geos_vtu< detail::VTUPolyhedralOutputImpl >(
                    file_vtu, *poly, surface_cells_, cell_attributes );
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