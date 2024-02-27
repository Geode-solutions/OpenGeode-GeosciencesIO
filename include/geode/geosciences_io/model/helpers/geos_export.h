/*
 * Copyright (c) 2019 - 2023 Geode-solutions
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

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <absl/container/flat_hash_map.h>
#include <absl/strings/string_view.h>

#include <ghc/filesystem.hpp>

#include <geode/basic/attribute.h>
#include <geode/basic/attribute_manager.h>
#include <geode/basic/filename.h>

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
#include <geode/model/representation/core/brep.h>

namespace geode
{
    template < typename Model >
    class GeosExporterImpl
    {
    public:
        GeosExporterImpl(
            absl::string_view files_directory, const Model& model )
            : model_( model ),
              files_directory_{ ghc::filesystem::path{
                  geode::to_string( files_directory ) }
                                    .string() },
              prefix_{ geode::filename_without_extension( files_directory ) }
        {
            std::tie( model_curve_, model_surface_, model_solid_ ) =
                geode::convert_brep_into_curve_and_surface_and_solid( model_ );
            static constexpr auto REGION_ID_ATTRIBUTE_NAME{ "attribute" };
            region_attribute_ =
                model_solid_->polyhedron_attribute_manager()
                    .template find_or_create_attribute<
                        geode::VariableAttribute, geode::index_t >(
                        REGION_ID_ATTRIBUTE_NAME, geode::NO_ID );
            if( ghc::filesystem::path{ geode::to_string( files_directory ) }
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

        virtual ~GeosExporterImpl() = default;

        void write_files() const
        {
            auto geometry_file = write_well_perforations_boxes();
            write_well_perforation_file();
            auto mesh_file = write_mesh_files();

            const auto simu_xml = absl::StrCat(
                files_directory(), "/", prefix(), "_simulation.xml" );
            std::ofstream file( simu_xml, std::ofstream::out );
            file << "<?xml version=\" 1.0 \"?>" << std::endl;
            file << "<Problem>" << std::endl;
            file << "<!-- include mesh xml file -->" << std::endl;
            file << "    <Included>" << std::endl;
            file << "        <File name = \"./" << mesh_file << "\"/>"
                 << std::endl;
            file << "    </Included>" << std::endl;
            file << "    <!-- include geometry xml file -->" << std::endl;
            file << "    <Included>" << std::endl;
            file << "        <File name = \"./" << geometry_file << "\"/>"
                 << std::endl;
            file << "    </Included>" << std::endl;
            file << "    <!-- please complete simulation file here -->"
                 << std::endl;
            file << "</Problem>" << std::endl;
        }

        void add_well_perforations( const PointSet3D& perforations )
        {
            well_perforations_.push_back( perforations.clone() );
        }

        void add_cell_property1D( absl::string_view property_name )
        {
            if( check_property_name( property_name ) )
            {
                cell_1Dproperty_names_.push_back(
                    geode::to_string( property_name ) );
            }
        }
        void add_cell_property2D( absl::string_view property_name )
        {
            if( check_property_name( property_name ) )
            {
                cell_2Dproperty_names_.push_back(
                    geode::to_string( property_name ) );
            }
        }
        void add_cell_property3D( absl::string_view property_name )
        {
            if( check_property_name( property_name ) )
            {
                cell_3Dproperty_names_.push_back(
                    geode::to_string( property_name ) );
            }
        }
        void prepare_export()
        {
            init_solid_region_attribute();
            transfer_cell_properties();
            delete_mapping_attributes();
        }

    protected:
        absl::string_view files_directory() const
        {
            return files_directory_;
        }

        absl::string_view prefix() const
        {
            return prefix_;
        }

        index_t init_solid_region_attribute()
        {
            auto region_map_id = create_region_attribute_map( model_ );
            const auto brep_component_uuid_attribute =
                model_solid_->polyhedron_attribute_manager()
                    .find_attribute<
                        geode::uuid_from_conversion_attribute_type >(
                        geode::uuid_from_conversion_attribute_name );
            for( const auto polyhedron_id :
                geode::Range( model_solid_->nb_polyhedra() ) )
            {
                region_attribute_->set_value( polyhedron_id,
                    region_map_id
                        .find( brep_component_uuid_attribute->value(
                            polyhedron_id ) )
                        ->second );
            }
            return region_map_id.size();
        }

        virtual absl::flat_hash_map< geode::uuid, geode::index_t >
            create_region_attribute_map( const Model& model ) const = 0;

        std::string write_well_perforations_boxes() const
        {
            auto aabb = create_aabb_tree( *model_solid_ );
            const auto filename = absl::StrCat( prefix(), "_wellbox.xml" );
            const auto file_xml = absl::StrCat(
                files_directory(), "/", prefix(), "_wellbox.xml" );
            std::ofstream file( file_xml, std::ofstream::out );
            file << "<?xml version=\" 1.0 \"?>" << std::endl;
            file << "<Problem>" << std::endl;
            file << "    <Geometry>" << std::endl;
            index_t well_id{ 0 };
            for( const auto& well : well_perforations_ )
            {
                BoundingBox3D perf_box;
                for( const auto point : geode::Range( well->nb_vertices() ) )
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
                        DEBUG_CONST auto tmp_dist = point_point_distance< 3 >(
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
                auto safety_lag = perf_box.diagonal().length() * 0.001;
                const auto boxstr = absl::StrCat( "        <Box name=\"well_",
                    well_id++, "\" xMin = \"{",
                    perf_box.min().value( 0 ) - safety_lag, ", ",
                    perf_box.min().value( 1 ) - safety_lag, ", ",
                    perf_box.min().value( 2 ) - safety_lag, "}\" xMax = \"{",
                    perf_box.max().value( 0 ) + safety_lag, ", ",
                    perf_box.max().value( 1 ) + safety_lag, ", ",
                    perf_box.max().value( 2 ) + safety_lag, "}\"  />" );
                file << boxstr << std::endl;
            }
            file << "    </Geometry>" << std::endl;
            file << "</Problem>" << std::endl;
            return filename;
        }

        std::string write_mesh_files() const
        {
            const auto file_vtu = write_solid_file();

            const auto filename = absl::StrCat( prefix(), "_mesh.xml" );
            const auto file_xml =
                absl::StrCat( files_directory(), "/", filename );
            std::ofstream file( file_xml, std::ofstream::out );
            file << "<?xml version=\" 1.0 \"?>" << std::endl;
            file << "<Problem>" << std::endl;
            file << "    <Mesh>" << std::endl;
            file << "        <VTKMesh" << std::endl;
            file << "            name = \"" << prefix() << "\"" << std::endl;
            file << "            file = \"./" << file_vtu << "\"" << std::endl;
            file << "            fieldsToImport = \"{";
            auto first = true;
            for( const auto& cell_prop_name : cell_1Dproperty_names_ )
            {
                if( !first )
                {
                    file << ",";
                }
                else
                {
                    first = false;
                }
                file << cell_prop_name;
            }
            for( const auto& cell_prop_name : cell_2Dproperty_names_ )
            {
                if( !first )
                {
                    file << ",";
                }
                else
                {
                    first = false;
                }
                file << cell_prop_name;
            }
            for( const auto& cell_prop_name : cell_3Dproperty_names_ )
            {
                if( !first )
                {
                    file << ",";
                }
                else
                {
                    first = false;
                }
                file << cell_prop_name;
            }
            file << "}\"" << std::endl;
            file << "            fieldNamesInGEOSX = \"{ please enter property "
                    "name in Geos}\"/>"
                 << std::endl;
            file << "    </Mesh>" << std::endl;
            file << "</Problem>" << std::endl;
            return filename;
        }

        bool check_property_name( absl::string_view property_name ) const
        {
            for( const auto& block : model_.blocks() )
            {
                if( !block.mesh()
                         .polyhedron_attribute_manager()
                         .attribute_exists( property_name ) )

                {
                    Logger::info( "The property ", property_name,
                        " will not be exported because it is not defined on "
                        "every block of the model." );
                    return false;
                }
            }
            return true;
        }

        void transfer_cell_properties()
        {
            const auto brep_mesh_elements =
                model_solid_->polyhedron_attribute_manager()
                    .template find_attribute<
                        geode::mesh_elements_attribute_type >(
                        geode::MESH_ELEMENT_ATTRIBUTE_NAME );
            for( const auto& property_name : cell_1Dproperty_names_ )
            {
                auto solid_property =
                    model_solid_->polyhedron_attribute_manager()
                        .template find_or_create_attribute<
                            geode::VariableAttribute, double >(
                            property_name, 0. );
                for( const auto polyhedron_id :
                    geode::Range( model_solid_->nb_polyhedra() ) )
                {
                    auto polygon_mesh_element =
                        brep_mesh_elements->value( polyhedron_id );

                    const auto model_property =
                        model_.block( polygon_mesh_element.mesh_id )
                            .mesh()
                            .polyhedron_attribute_manager()
                            .template find_attribute< double >( property_name );
                    auto value = model_property->value(
                        polygon_mesh_element.element_id );
                    solid_property->set_value( polyhedron_id, value );
                }
            }
            for( const auto& property_name : cell_2Dproperty_names_ )
            {
                auto solid_property =
                    model_solid_->polyhedron_attribute_manager()
                        .template find_or_create_attribute<
                            geode::VariableAttribute, std::array< double, 2 > >(
                            property_name, { 0., 0. } );
                for( const auto polyhedron_id :
                    geode::Range( model_solid_->nb_polyhedra() ) )
                {
                    auto polygon_mesh_element =
                        brep_mesh_elements->value( polyhedron_id );

                    const auto model_property =
                        model_.block( polygon_mesh_element.mesh_id )
                            .mesh()
                            .polyhedron_attribute_manager()
                            .template find_attribute< std::array< double, 2 > >(
                                property_name );
                    auto value = model_property->value(
                        polygon_mesh_element.element_id );
                    solid_property->set_value( polyhedron_id, value );
                }
            }
            for( const auto& property_name : cell_3Dproperty_names_ )
            {
                auto solid_property =
                    model_solid_->polyhedron_attribute_manager()
                        .template find_or_create_attribute<
                            geode::VariableAttribute, std::array< double, 3 > >(
                            property_name, { 0., 0., 0. } );
                for( const auto polyhedron_id :
                    geode::Range( model_solid_->nb_polyhedra() ) )
                {
                    auto polygon_mesh_element =
                        brep_mesh_elements->value( polyhedron_id );

                    const auto model_property =
                        model_.block( polygon_mesh_element.mesh_id )
                            .mesh()
                            .polyhedron_attribute_manager()
                            .template find_attribute< std::array< double, 3 > >(
                                property_name );
                    auto value = model_property->value(
                        polygon_mesh_element.element_id );
                    solid_property->set_value( polyhedron_id, value );
                }
            }
        }

        void delete_mapping_attributes()
        {
            model_solid_->vertex_attribute_manager().delete_attribute(
                geode::unique_vertex_from_conversion_attribute_name );
            model_solid_->polyhedron_attribute_manager().delete_attribute(
                geode::MESH_ELEMENT_ATTRIBUTE_NAME );
            model_solid_->polyhedron_attribute_manager().delete_attribute(
                geode::uuid_from_conversion_attribute_name );
        }

        std::string write_solid_file() const
        {
            const auto filename = absl::StrCat( prefix(), ".vtu" );
            const auto file_vtu =
                absl::StrCat( files_directory(), "/", filename );
            if( const auto* tetra =
                    dynamic_cast< const geode::TetrahedralSolid3D* >(
                        model_solid_.get() ) )
            {
                geode::save_tetrahedral_solid( *tetra, file_vtu );
            }
            else if( const auto* hybrid =
                         dynamic_cast< const geode::HybridSolid3D* >(
                             model_solid_.get() ) )
            {
                geode::save_hybrid_solid( *hybrid, file_vtu );
            }
            else if( const auto* poly =
                         dynamic_cast< const geode::PolyhedralSolid3D* >(
                             model_solid_.get() ) )
            {
                geode::save_polyhedral_solid( *poly, file_vtu );
            }
            else
            {
                throw geode::OpenGeodeException(
                    "[Blocks::save_geos] Cannot find the explicit "
                    "SolidMesh type" );
            }
            return filename;
        }

        void write_well_perforation_file() const
        {
            auto id{ 0 };
            for( const auto& well : well_perforations_ )
            {
                const auto file = absl::StrCat(
                    files_directory(), "/", prefix(), "_well", id, ".vtp" );
                geode::save_point_set( *well, file );
                id++;
            }
        }

    private:
        const Model& model_;
        std::unique_ptr< geode::EdgedCurve3D > model_curve_;
        std::unique_ptr< geode::SurfaceMesh3D > model_surface_;
        std::unique_ptr< geode::SolidMesh3D > model_solid_;

        std::shared_ptr< geode::VariableAttribute< geode::index_t > >
            region_attribute_;

        DEBUG_CONST std::string files_directory_;
        DEBUG_CONST std::string prefix_;

        std::vector< std::string > cell_1Dproperty_names_;
        std::vector< std::string > cell_2Dproperty_names_;
        std::vector< std::string > cell_3Dproperty_names_;

        std::vector< std::unique_ptr< PointSet3D > > well_perforations_;
    };
} // namespace geode