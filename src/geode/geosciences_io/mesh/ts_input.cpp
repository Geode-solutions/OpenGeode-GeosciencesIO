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

#include <geode/geosciences_io/mesh/private/ts_input.h>

#include <fstream>

#include <geode/geometry/point.h>

#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/core/triangulated_surface.h>

#include <geode/geosciences_io/mesh/private/gocad_common.h>

namespace
{
    class TSInputImpl
    {
    public:
        TSInputImpl(
            absl::string_view filename, geode::TriangulatedSurface3D& surface )
            : file_{ geode::to_string( filename ) },
              surface_( surface ),
              builder_(
                  geode::TriangulatedSurfaceBuilder< 3 >::create( surface ) )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void read_file()
        {
            while( const auto tsurf = geode::detail::read_tsurf( file_ ) )
            {
                build_surface( tsurf.value() );
            }
            builder_->compute_polygon_adjacencies();
        }

    private:
        void build_surface( const geode::detail::TSurfData& tsurf )
        {
            const auto offset = surface_.nb_vertices();
            builder_->set_name( tsurf.header.name );
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
            for( const auto v : geode::Indices{ inverse_vertex_mapping } )
            {
                inverse_vertex_mapping[v] = v;
            }
            geode::detail::create_attributes( tsurf.vertices_properties_header,
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
    namespace detail
    {
        std::unique_ptr< TriangulatedSurface3D > TSInput::read(
            const MeshImpl& impl )
        {
            auto surface = TriangulatedSurface3D::create( impl );
            TSInputImpl reader{ this->filename(), *surface };
            reader.read_file();
            return surface;
        }
    } // namespace detail
} // namespace geode
