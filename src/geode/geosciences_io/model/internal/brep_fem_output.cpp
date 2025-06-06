/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

#include <geode/geosciences_io/model/internal/brep_fem_output.hpp>

#include <fstream>
#include <string_view>

#include <geode/basic/attribute.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/io.hpp>
#include <geode/basic/types.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/solid_edges.hpp>
#include <geode/mesh/core/solid_facets.hpp>
#include <geode/mesh/core/surface_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/io/tetrahedral_solid_output.hpp>

#include <geode/model/helpers/convert_model_meshes.hpp>
#include <geode/model/helpers/convert_to_mesh.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/line.hpp>

namespace
{

    class BRepFemOutputImpl
    {
        static constexpr const char APERTURE_ATTRIBUTE_NAME[]{
            "diagres_discontinuity_aperture"
        };
        static constexpr const char CONDUIT_AREA_ATTRIBUTE_NAME[]{
            "diagres_conduit_area"
        };
        static constexpr const char CONDUCTIVITY_ATTRIBUTE_NAME[]{
            "diagres_conductivity"
        };
        static constexpr const char SURFACE_NAME_ATTRIBUTE[]{ "surface_name" };
        static constexpr const char LINE_NAME_ATTRIBUTE[]{ "line_name" };

    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };
        BRepFemOutputImpl( std::string_view filename, const geode::BRep& brep )
            : file_{ geode::to_string( filename ) },
              file_str_view_{ filename },
              brep_( brep )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void write_file()
        {
            geode::Logger::info(
                "[BRepSolidFemOutput::write] Writing fem file." );
            const auto directories =
                geode::filepath_without_filename( "filename" );
            if( !directories.empty() )
            {
                std::filesystem::create_directories( directories );
            }
            auto [solid, model_to_mesh_mapping] =
                geode::convert_brep_into_solid( brep_ );
            auto& tet_solid =
                static_cast< const geode::TetrahedralSolid3D& >( *solid );
            copy_tetrahedra_attributes_in_solid(
                tet_solid, model_to_mesh_mapping );
            paint_polyhedron_block_id_attributes(
                tet_solid, model_to_mesh_mapping );
            paint_vertices_block_id_attributes(
                tet_solid, model_to_mesh_mapping );
            paint_surface_infos( tet_solid, model_to_mesh_mapping );
            paint_line_infos( tet_solid, model_to_mesh_mapping );
            geode::save_tetrahedral_solid( tet_solid, file_str_view_ );
        }

    private:
        void paint_facets_attribute( const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping,
            std::string_view attribute_name,
            double default_value ) const
        {
            auto solid_attribute =
                tet_solid.facets()
                    .facet_attribute_manager()
                    .find_or_create_attribute< geode::VariableAttribute,
                        double >( attribute_name, default_value );
            auto solid_surface_name_attribute =
                tet_solid.facets()
                    .facet_attribute_manager()
                    .find_or_create_attribute< geode::VariableAttribute,
                        std::string >( SURFACE_NAME_ATTRIBUTE, "No_name" );
            for( const auto& surface : brep_.surfaces() )
            {
                const auto& mesh = surface.mesh();
                if( !mesh.polygon_attribute_manager().attribute_exists(
                        attribute_name ) )
                {
                    continue;
                }
                const auto& surface_attribute =
                    mesh.polygon_attribute_manager().find_generic_attribute(
                        attribute_name );
                for( const auto polygon : geode::Range{ mesh.nb_polygons() } )
                {
                    const auto facets_in_solid =
                        model_to_mesh_mapping.surface_polygons_mapping.in2out(
                            { surface.id(), polygon } );
                    const auto surface_attribute_value =
                        surface_attribute->generic_value( polygon );
                    for( const auto facet_in_solid : facets_in_solid )
                    {
                        solid_attribute->set_value(
                            facet_in_solid, surface_attribute_value );
                        std::string surface_name{ surface.name() };
                        solid_surface_name_attribute->set_value(
                            facet_in_solid, surface_name );
                    }
                }
            }
        }

        void paint_surface_infos( const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            paint_facets_attribute( tet_solid, model_to_mesh_mapping,
                APERTURE_ATTRIBUTE_NAME, -1.0 );
            paint_facets_attribute( tet_solid, model_to_mesh_mapping,
                CONDUCTIVITY_ATTRIBUTE_NAME, -1.0 );
        }

        void paint_line_attribute( const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping,
            std::string_view attribute_name,
            double default_value ) const
        {
            const auto& solid_attribute =
                tet_solid.edges()
                    .edge_attribute_manager()
                    .find_or_create_attribute< geode::VariableAttribute,
                        double >( attribute_name, default_value );
            auto solid_line_name_attribute =
                tet_solid.edges()
                    .edge_attribute_manager()
                    .find_or_create_attribute< geode::VariableAttribute,
                        std::string >( LINE_NAME_ATTRIBUTE, "No_name" );
            for( const auto& line : brep_.lines() )
            {
                const auto& mesh = line.mesh();
                if( !mesh.edge_attribute_manager().attribute_exists(
                        attribute_name ) )
                {
                    continue;
                }
                const auto& line_attribute =
                    mesh.edge_attribute_manager().find_generic_attribute(
                        attribute_name );
                for( const auto edge : geode::Range{ mesh.nb_edges() } )
                {
                    const auto edges_in_solid =
                        model_to_mesh_mapping.line_edges_mapping.in2out(
                            { line.component_id().id(), edge } );
                    const auto line_attribute_value =
                        line_attribute->generic_value( edge );
                    for( const auto edge_in_solid : edges_in_solid )
                    {
                        solid_attribute->set_value(
                            edge_in_solid, line_attribute_value );
                        std::string line_name{ line.name() };
                        solid_line_name_attribute->set_value(
                            edge_in_solid, line_name );
                    }
                }
            }
        }

        void paint_line_infos( const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            paint_line_attribute( tet_solid, model_to_mesh_mapping,
                CONDUIT_AREA_ATTRIBUTE_NAME, -1.0 );
            paint_line_attribute( tet_solid, model_to_mesh_mapping,
                CONDUCTIVITY_ATTRIBUTE_NAME, -1.0 );
        }

        void paint_polyhedron_block_id_attributes(
            const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            auto attribute_p =
                tet_solid.polyhedron_attribute_manager()
                    .find_or_create_attribute< geode::VariableAttribute,
                        std::string_view >( "Block_ID_polyhedron", "No_name" );
            for( const auto& block : brep_.blocks() )
            {
                for( const auto polyhedron_id :
                    geode::Range{ block.mesh().nb_polyhedra() } )
                {
                    for( const auto polyhedron_out :
                        model_to_mesh_mapping.solid_polyhedra_mapping.in2out(
                            { block.id(), polyhedron_id } ) )
                    {
                        attribute_p->set_value( polyhedron_out, block.name() );
                    }
                }
            }
        }

        void paint_vertices_block_id_attributes(
            const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            auto attribute_v =
                tet_solid.vertex_attribute_manager()
                    .find_or_create_attribute< geode::VariableAttribute,
                        std::vector< std::string_view > >(
                        "Block_ID_vertex", {} );
            for( const auto& block : brep_.blocks() )
            {
                for( const auto vertex_id :
                    geode::Range{ block.mesh().nb_vertices() } )
                {
                    const auto unique_vertex = brep_.unique_vertex(
                        { block.component_id(), vertex_id } );
                    const auto vertex_out =
                        model_to_mesh_mapping.unique_vertices_mapping.in2out(
                            unique_vertex );
                    auto vertex_blocks = attribute_v->value( vertex_out );
                    vertex_blocks.push_back( block.name() );
                    attribute_v->set_value( vertex_out, vertex_blocks );
                }
            }
        }

        void copy_tetrahedra_attributes_in_solid(
            const geode::TetrahedralSolid3D& solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            for( const auto& block : brep_.blocks() )
            {
                const auto& mesh = block.mesh();
                for( const auto& attribute_name :
                    mesh.polyhedron_attribute_manager().attribute_names() )
                {
                    const auto& block_polyhedron_attribute =
                        mesh.polyhedron_attribute_manager()
                            .find_generic_attribute( attribute_name );
                    if( !block_polyhedron_attribute )
                    {
                        continue;
                    }
                    if( !block_polyhedron_attribute->properties().transferable )
                    {
                        continue;
                    }
                    if( block_polyhedron_attribute->type()
                            != typeid( double ).name()
                        && block_polyhedron_attribute->type()
                               != typeid( float ).name() )
                    {
                        continue;
                    }
                    auto solid_polyhedron_attribute =
                        solid.polyhedron_attribute_manager()
                            .find_or_create_attribute< geode::VariableAttribute,
                                double >( attribute_name,
                                block_polyhedron_attribute->generic_value(
                                    0 ) );
                    for( const auto polyhedron :
                        geode::Range{ mesh.nb_polyhedra() } )
                    {
                        const auto polyhedron_out =
                            model_to_mesh_mapping.solid_polyhedra_mapping
                                .in2out( { block.id(), polyhedron } );
                        solid_polyhedron_attribute->set_value(
                            polyhedron_out.front(),
                            block_polyhedron_attribute->generic_value(
                                polyhedron ) );
                    }
                }
            }
        }

    private:
        std::ofstream file_;
        std::string_view file_str_view_;
        const geode::BRep& brep_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::vector< std::string > BRepFemOutput::write(
            const BRep& brep ) const
        {
            BRepFemOutputImpl impl{ filename(), brep };
            impl.write_file();
            return { to_string( filename() ) };
        }

        bool BRepFemOutput::is_saveable( const BRep& brep ) const
        {
            for( const auto unique_vertex : Range{ brep.nb_unique_vertices() } )
            {
                if( !brep.has_component_mesh_vertices(
                        unique_vertex, Block3D::component_type_static() ) )
                {
                    return false;
                }
            }
            return true;
        }
    } // namespace internal
} // namespace geode
