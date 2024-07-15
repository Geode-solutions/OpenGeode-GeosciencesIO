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

#include <geode/tests_config.h>

#include <geode/basic/assert.h>
#include <geode/basic/logger.h>

#include <geode/mesh/core/geode/geode_triangulated_surface.h>
#include <geode/mesh/core/triangulated_surface.h>
#include <geode/mesh/io/triangulated_surface_input.h>
#include <geode/mesh/io/triangulated_surface_output.h>

#include <geode/model/mixin/core/surface.h>

#include <geode/geosciences_io/mesh/internal/ts_input.h>

void check_surface( const geode::SurfaceMesh3D& surface,
    geode::index_t nb_vertices,
    geode::index_t nb_polygons )
{
    OPENGEODE_EXCEPTION( surface.nb_vertices() == nb_vertices,
        "Number of vertices in the TSurf 3D is not correct" );
    OPENGEODE_EXCEPTION( surface.nb_polygons() == nb_polygons,
        "Number of polygons in the TSurf 3D is not correct" );
}

void check_file(
    std::string file, geode::index_t nb_vertices, geode::index_t nb_polygons )
{
    // Load file
    auto surface = geode::load_triangulated_surface< 3 >( file );
    check_surface( *surface, nb_vertices, nb_polygons );

    // Save triangulated tsurf
    const auto output_file_native =
        absl::StrCat( "test_output.", surface->native_extension() );
    geode::save_triangulated_surface( *surface, output_file_native );
    const auto output_file_ts = "test_output.ts";
    geode::save_triangulated_surface( *surface, output_file_ts );

    // Load file
    auto reloaded_surface =
        geode::load_triangulated_surface< 3 >( output_file_native );
    check_surface( *reloaded_surface, nb_vertices, nb_polygons );
    auto reloaded_surface_ts =
        geode::load_triangulated_surface< 3 >( output_file_ts );
    check_surface( *reloaded_surface_ts, nb_vertices, nb_polygons );
}

int main()
{
    try
    {
        geode::GeosciencesIOMeshLibrary::initialize();
        check_file( absl::StrCat( geode::DATA_PATH, "/surf2d_multi.",
                        geode::internal::TSInput::extension() ),
            92, 92 );
        check_file( absl::StrCat( geode::DATA_PATH, "/surf2d.",
                        geode::internal::TSInput::extension() ),
            46, 46 );
        check_file( absl::StrCat( geode::DATA_PATH, "/2triangles.",
                        geode::internal::TSInput::extension() ),
            4, 2 );
        check_file( absl::StrCat( geode::DATA_PATH, "/sgrid_tsurf.",
                        geode::internal::TSInput::extension() ),
            4, 2 );
        check_file( absl::StrCat( geode::DATA_PATH, "Fault_without_crs.",
                        geode::internal::TSInput::extension() ),
            189, 324 );
        check_file( absl::StrCat( geode::DATA_PATH, "ts-2props.",
                        geode::internal::TSInput::extension() ),
            4, 2 );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
