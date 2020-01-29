/*
 * Copyright (c) 2019 Geode-solutions
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

#include <geode/geosciences/private/ts_input.h>

#include <fstream>

#include <geode/geometry/point.h>

#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/core/triangulated_surface.h>

#include <geode/geosciences/private/gocad_common.h>

namespace
{
    class TSInputImpl
    {
    public:
        TSInputImpl( absl::string_view filename,
            geode::TriangulatedSurfaceBuilder3D& builder )
            : file_( filename.data() ), builder_( builder )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void read_file()
        {
            build_surface( geode::detail::read_tsurf( file_ ) );
        }

    private:
        void build_surface( const geode::detail::TSurfData& tsurf )
        {
            for( const auto& point : tsurf.points )
            {
                builder_.create_point( point );
            }

            for( const auto& triangle : tsurf.triangles )
            {
                builder_.create_triangle( triangle );
            }

            builder_.compute_polygon_adjacencies();
        }

    private:
        std::ifstream file_;
        geode::TriangulatedSurfaceBuilder3D& builder_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        void TSInput::do_read()
        {
            auto builder = TriangulatedSurfaceBuilder< 3 >::create(
                this->triangulated_surface() );
            TSInputImpl impl{ this->filename(), *builder };
            impl.read_file();
        }
    } // namespace detail
} // namespace geode
