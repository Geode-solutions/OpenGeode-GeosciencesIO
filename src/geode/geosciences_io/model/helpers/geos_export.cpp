/*
 * Copyright (c) 2019 - 2024 Geode-solutions
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

#pragma once

#include <fstream>
#include <memory>
#include <string>

#include <absl/strings/string_view.h>

#include <ghc/filesystem.hpp>
#include <pugixml.hpp>

#include <geode/basic/attribute.h>
#include <geode/basic/attribute_manager.h>
#include <geode/basic/filename.h>
#include <geode/basic/logger.h>

#include <geode/geometry/aabb.h>
#include <geode/geometry/distance.h>

#include <geode/geosciences_io/model/common.h>

#include <geode/mesh/core/hybrid_solid.h>
#include <geode/mesh/core/point_set.h>
#include <geode/mesh/core/polyhedral_solid.h>
#include <geode/mesh/core/solid_mesh.h>
#include <geode/mesh/core/surface_mesh.h>
#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/mesh/helpers/aabb_solid_helpers.h>

#include <geode/mesh/io/hybrid_solid_output.h>
#include <geode/mesh/io/point_set_output.h>
#include <geode/mesh/io/polyhedral_solid_output.h>
#include <geode/mesh/io/tetrahedral_solid_output.h>

#include <geode/model/helpers/convert_to_mesh.h>
#include <geode/model/mixin/core/block.h>

#include <geode/geosciences/explicit/representation/core/structural_model.h>
#include <geode/geosciences_io/model/helpers/geos_export.h>
#include <geode/model/representation/core/brep.h>

namespace geode
{
    constexpr auto REGION_ID_ATTRIBUTE_NAME = "attribute";

    template < typename Model >
    GeosExporterImpl< Model >::GeosExporterImpl(
        absl::string_view files_directory, const Model& model )
        : model_( model ),
          files_directory_{
              ghc::filesystem::path{ to_string( files_directory ) }.string()
          },
          prefix_{ filename_without_extension( files_directory ) }
    {
        std::tie( model_curve_, model_surface_, model_solid_ ) =
            convert_brep_into_curve_and_surface_and_solid( model_ );
        region_attribute_ =
            model_solid_->polyhedron_attribute_manager()
                .template find_or_create_attribute< VariableAttribute,
                    index_t >( REGION_ID_ATTRIBUTE_NAME, NO_ID );
        if( ghc::filesystem::path{ to_string( files_directory ) }
                .is_relative() )
        {
            ghc::filesystem::create_directory(
                ghc::filesystem::current_path() / files_directory_ );
        }
        else
        {
            ghc::filesystem::create_directory( files_directory_ );
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

        const auto filename_xml =
            absl::StrCat( files_directory(), "/", prefix(), "_simulation.xml" );
        doc_xml.save_file( filename_xml.c_str(), PUGIXML_TEXT( "    " ),
            pugi::format_indent_attributes );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::add_well_perforations(
        const PointSet3D& perforations )
    {
        well_perforations_.push_back( perforations.clone() );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::add_cell_property1d(
        absl::string_view property_name )
    {
        if( check_property_name( property_name ) )
        {
            cell_1Dproperty_names_.push_back( to_string( property_name ) );
        }
    }
    template < typename Model >
    void GeosExporterImpl< Model >::add_cell_property2d(
        absl::string_view property_name )
    {
        if( check_property_name( property_name ) )
        {
            cell_2Dproperty_names_.push_back( to_string( property_name ) );
        }
    }
    template < typename Model >
    void GeosExporterImpl< Model >::add_cell_property3d(
        absl::string_view property_name )
    {
        if( check_property_name( property_name ) )
        {
            cell_3Dproperty_names_.push_back( to_string( property_name ) );
        }
    }
    template < typename Model >
    void GeosExporterImpl< Model >::prepare_export()
    {
        initialize_solid_region_attribute();
        transfer_cell_properties();
        delete_mapping_attributes();
    }

    template < typename Model >
    absl::string_view GeosExporterImpl< Model >::files_directory() const
    {
        return files_directory_;
    }
    template < typename Model >
    absl::string_view GeosExporterImpl< Model >::prefix() const
    {
        return prefix_;
    }

    template < typename Model >
    index_t GeosExporterImpl< Model >::initialize_solid_region_attribute()
    {
        auto region_map_id = create_region_attribute_map( model_ );
        const auto brep_component_uuid_attribute =
            model_solid_->polyhedron_attribute_manager()
                .find_attribute< uuid_from_conversion_attribute_type >(
                    uuid_from_conversion_attribute_name );
        for( const auto polyhedron_id : Range( model_solid_->nb_polyhedra() ) )
        {
            region_attribute_->set_value( polyhedron_id,
                region_map_id
                    .find(
                        brep_component_uuid_attribute->value( polyhedron_id ) )
                    ->second );
        }
        return region_map_id.size();
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_well_perforations_boxes(
        pugi::xml_node& root ) const
    {
        auto aabb = create_aabb_tree( *model_solid_ );
        index_t well_id{ 0 };
        for( const auto& well : well_perforations_ )
        {
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
                    const auto tmp_dist =
                        point_point_distance< 3 >( well->point( point ),
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
                absl::StrCat( "well_", well_id++ ).c_str() );
            box_node.append_attribute( "xMin" ).set_value(
                absl::StrCat( "{", perf_box.min().value( 0 ) - global_epsilon,
                    ", ", perf_box.min().value( 1 ) - global_epsilon, ", ",
                    perf_box.min().value( 2 ) - global_epsilon, "}" )
                    .c_str() );
            box_node.append_attribute( "xMax" ).set_value(
                absl::StrCat( "{", perf_box.max().value( 0 ) + global_epsilon,
                    ", ", perf_box.max().value( 1 ) + global_epsilon, ", ",
                    perf_box.max().value( 2 ) + global_epsilon, "}" )
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
        std::string property_names{ "{" };
        auto first = true;
        for( const auto& cell_prop_name : cell_1Dproperty_names_ )
        {
            if( !first )
            {
                absl::StrAppend( &property_names, "," );
            }
            else
            {
                first = false;
            }
            absl::StrAppend( &property_names, cell_prop_name );
        }
        for( const auto& cell_prop_name : cell_2Dproperty_names_ )
        {
            if( !first )
            {
                absl::StrAppend( &property_names, "," );
            }
            else
            {
                first = false;
            }
            absl::StrAppend( &property_names, cell_prop_name );
        }
        for( const auto& cell_prop_name : cell_3Dproperty_names_ )

        {
            if( !first )
            {
                absl::StrAppend( &property_names, "," );
            }
            else
            {
                first = false;
            }
            absl::StrAppend( &property_names, cell_prop_name );
        }
        absl::StrAppend( &property_names, "}" );
        vtk_mesh_node.append_attribute( "fieldsToImport" )
            .set_value( property_names.c_str() );
        vtk_mesh_node.append_attribute( "fieldNamesInGEOSX" )
            .set_value( "{ please enter property name in Geos}" );
    }

    template < typename Model >
    bool GeosExporterImpl< Model >::check_property_name(
        absl::string_view property_name ) const
    {
        for( const auto& block : model_.blocks() )
        {
            if( !block.mesh().polyhedron_attribute_manager().attribute_exists(
                    property_name ) )

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
    void GeosExporterImpl< Model >::transfer_cell_properties()
    {
        const auto brep_mesh_elements =
            model_solid_->polyhedron_attribute_manager()
                .template find_attribute< mesh_elements_attribute_type >(
                    MESH_ELEMENT_ATTRIBUTE_NAME );
        for( const auto& property_name : cell_1Dproperty_names_ )
        {
            auto solid_property =
                model_solid_->polyhedron_attribute_manager()
                    .template find_or_create_attribute< VariableAttribute,
                        double >( property_name, 0. );
            for( const auto polyhedron_id :
                Range( model_solid_->nb_polyhedra() ) )
            {
                auto polygon_mesh_element =
                    brep_mesh_elements->value( polyhedron_id );

                const auto model_property =
                    model_.block( polygon_mesh_element.mesh_id )
                        .mesh()
                        .polyhedron_attribute_manager()
                        .template find_attribute< double >( property_name );
                auto value =
                    model_property->value( polygon_mesh_element.element_id );
                solid_property->set_value( polyhedron_id, value );
            }
        }
        for( const auto& property_name : cell_2Dproperty_names_ )
        {
            auto solid_property =
                model_solid_->polyhedron_attribute_manager()
                    .template find_or_create_attribute< VariableAttribute,
                        std::array< double, 2 > >( property_name, { 0., 0. } );
            for( const auto polyhedron_id :
                Range( model_solid_->nb_polyhedra() ) )
            {
                auto polygon_mesh_element =
                    brep_mesh_elements->value( polyhedron_id );

                const auto model_property =
                    model_.block( polygon_mesh_element.mesh_id )
                        .mesh()
                        .polyhedron_attribute_manager()
                        .template find_attribute< std::array< double, 2 > >(
                            property_name );
                auto value =
                    model_property->value( polygon_mesh_element.element_id );
                solid_property->set_value( polyhedron_id, value );
            }
        }
        for( const auto& property_name : cell_3Dproperty_names_ )
        {
            auto solid_property =
                model_solid_->polyhedron_attribute_manager()
                    .template find_or_create_attribute< VariableAttribute,
                        std::array< double, 3 > >(
                        property_name, { 0., 0., 0. } );
            for( const auto polyhedron_id :
                Range( model_solid_->nb_polyhedra() ) )
            {
                auto polygon_mesh_element =
                    brep_mesh_elements->value( polyhedron_id );

                const auto model_property =
                    model_.block( polygon_mesh_element.mesh_id )
                        .mesh()
                        .polyhedron_attribute_manager()
                        .template find_attribute< std::array< double, 3 > >(
                            property_name );
                auto value =
                    model_property->value( polygon_mesh_element.element_id );
                solid_property->set_value( polyhedron_id, value );
            }
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::delete_mapping_attributes()
    {
        model_solid_->vertex_attribute_manager().delete_attribute(
            unique_vertex_from_conversion_attribute_name );
        model_solid_->polyhedron_attribute_manager().delete_attribute(
            MESH_ELEMENT_ATTRIBUTE_NAME );
        model_solid_->polyhedron_attribute_manager().delete_attribute(
            uuid_from_conversion_attribute_name );
    }

    template < typename Model >
    std::string GeosExporterImpl< Model >::write_solid_file() const
    {
        const auto filename = absl::StrCat( prefix(), ".vtu" );
        const auto file_vtu = absl::StrCat( files_directory(), "/", filename );
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
        else if( const auto* poly = dynamic_cast< const PolyhedralSolid3D* >(
                     model_solid_.get() ) )
        {
            save_polyhedral_solid( *poly, file_vtu );
        }
        else
        {
            throw OpenGeodeException(
                "[Blocks::save_geos] Cannot find the explicit "
                "SolidMesh type" );
        }
        return filename;
    }
    template < typename Model >
    void GeosExporterImpl< Model >::write_well_perforation_file() const
    {
        index_t well_id{ 0 };
        for( const auto& well : well_perforations_ )
        {
            const auto file = absl::StrCat(
                files_directory(), "/", prefix(), "_well", well_id, ".vtp" );
            save_point_set( *well, file );
            well_id++;
        }
    }
    template class opengeode_geosciencesio_model_api GeosExporterImpl< BRep >;
    template class opengeode_geosciencesio_model_api
        GeosExporterImpl< StructuralModel >;

} // namespace geode