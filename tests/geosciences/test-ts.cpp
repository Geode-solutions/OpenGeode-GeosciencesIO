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

#include <geode/tests_config.h>

#include <geode/basic/assert.h>
#include <geode/basic/logger.h>

#include <geode/mesh/core/geode_triangulated_surface.h>
#include <geode/mesh/core/triangulated_surface.h>
#include <geode/mesh/io/triangulated_surface_input.h>
#include <geode/mesh/io/triangulated_surface_output.h>

#include <geode/model/mixin/core/surface.h>

#include <geode/geosciences/private/ts_input.h>

void test_tsurf_3d()
{
    auto surface = geode::TriangulatedSurface3D::create();

    // Load file
    load_triangulated_surface(
        *surface, absl::StrCat( geode::data_path, "/surf2d.",
                      geode::detail::TSInput::extension() ) );

    OPENGEODE_EXCEPTION( surface->nb_vertices() == 46,
        "Number of vertices in the loaded TSurf 3D is not correct" );
    OPENGEODE_EXCEPTION( surface->nb_polygons() == 46,
        "Number of polygons in the loaded TSurf 3D is not correct" );

    // Save triangulated tsurf
    const auto output_file_native =
        absl::StrCat( "surf3d.", surface->native_extension() );
    save_triangulated_surface( *surface, output_file_native );

    // Reload
    auto reloaded_surface = geode::TriangulatedSurface3D::create();

    // Load file
    load_triangulated_surface( *reloaded_surface, output_file_native );
    OPENGEODE_EXCEPTION( reloaded_surface->nb_vertices() == 46,
        "Number of vertices in the reloaded TSurf 3D is not correct" );
    OPENGEODE_EXCEPTION( reloaded_surface->nb_polygons() == 46,
        "Number of polygons in the reloaded TSurf 3D is not correct" );
}

int main()
{
    using namespace geode;

    try
    {
        detail::initialize_geosciences_io();
        test_tsurf_3d();

        Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
