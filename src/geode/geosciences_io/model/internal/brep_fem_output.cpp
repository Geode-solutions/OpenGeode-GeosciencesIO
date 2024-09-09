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

#include <geode/geosciences_io/model/internal/brep_fem_output.hpp>

#include <fstream>
#include <string_view>

#include <geode/basic/attribute.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/io.hpp>
#include <geode/basic/types.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/io/tetrahedral_solid_output.hpp>

#include <geode/model/helpers/convert_model_meshes.hpp>
#include <geode/model/helpers/convert_to_mesh.hpp>
#include <geode/model/mixin/core/block.hpp>

namespace
{

    class BRepFemOutputImpl
    {
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
            auto attribute_p =
                tet_solid.polyhedron_attribute_manager()
                    .find_or_create_attribute< geode::VariableAttribute,
                        std::string_view >( "Block_ID_polyhedron", "No name" );
            auto attribute_v =
                tet_solid.vertex_attribute_manager()
                    .find_or_create_attribute< geode::VariableAttribute,
                        std::string_view >( "Block_ID_vertex", "No name" );
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
                for( const auto vertex_id :
                    geode::Range{ block.mesh().nb_vertices() } )
                {
                    const auto vertex_out =
                        model_to_mesh_mapping.unique_vertices_mapping.in2out(
                            vertex_id );

                    attribute_v->set_value( vertex_out, block.name() );
                }
            }
            geode::save_tetrahedral_solid( tet_solid, file_str_view_ );
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
    } // namespace internal
} // namespace geode
