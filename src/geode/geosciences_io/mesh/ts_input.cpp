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

#include <geode/geosciences_io/mesh/internal/ts_input.hpp>

#include <fstream>

#include <geode/basic/file.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/triangulated_surface_builder.hpp>
#include <geode/mesh/core/triangulated_surface.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

namespace
{
    class TSInputImpl
    {
    public:
        TSInputImpl(
            std::string_view filename, geode::TriangulatedSurface3D& surface )
            : file_{ geode::to_string( filename ),std::ios::binary },
              surface_( surface ),
              builder_(
                  geode::TriangulatedSurfaceBuilder< 3 >::create( surface ) )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void read_file()
        {
            while( const auto tsurf = geode::internal::read_tsurf( file_ ) )
            {
                build_surface( tsurf.value() );
            }
            builder_->compute_polygon_adjacencies();
        }

    private:
        void build_surface( const geode::internal::TSurfData& tsurf )
        {
            const auto offset = surface_.nb_vertices();
            if( tsurf.header.name )
            {
                builder_->set_name( tsurf.header.name.value() );
            }
            for( const auto& point : tsurf.points )
            {
                builder_->create_point( point );
            }
            if( offset == 0 )
            {
                for( const auto& triangle : tsurf.triangles )
                {
                    builder_->create_triangle( triangle );
                }
            }
            else
            {
                for( auto triangle : tsurf.triangles )
                {
                    for( auto& vertex : triangle )
                    {
                        vertex += offset;
                    }
                    builder_->create_triangle( triangle );
                }
            }
            std::vector< geode::index_t > inverse_vertex_mapping(
                tsurf.points.size() );
            absl::c_iota( inverse_vertex_mapping, 0 );
            geode::internal::create_attributes(
                tsurf.vertices_properties_header,
                tsurf.vertices_attribute_values,
                surface_.vertex_attribute_manager(), tsurf.points.size(),
                inverse_vertex_mapping );
        }

    private:
        std::ifstream file_;
        geode::TriangulatedSurface3D& surface_;
        std::unique_ptr< geode::TriangulatedSurfaceBuilder3D > builder_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::unique_ptr< TriangulatedSurface3D > TSInput::read(
            const MeshImpl& impl )
        {
            auto surface = TriangulatedSurface3D::create( impl );
            TSInputImpl reader{ this->filename(), *surface };
            reader.read_file();
            return surface;
        }

        Percentage TSInput::is_loadable() const
        {
            std::ifstream file{ to_string( this->filename() ) };
            if( goto_keyword_if_it_exists( file, "GOCAD TSurf" ) )
            {
                return Percentage{ 1 };
            }
            return Percentage{ 0 };
        }
    } // namespace internal
} // namespace geode
