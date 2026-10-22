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

#include <geode/tests_config.hpp>

#include <fstream>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>

#include <geode/mesh/core/geode/geode_triangulated_surface.hpp>
#include <geode/mesh/core/triangulated_surface.hpp>
#include <geode/mesh/io/triangulated_surface_input.hpp>
#include <geode/mesh/io/triangulated_surface_output.hpp>

#include <geode/geosciences_io/mesh/internal/ts_input.hpp>

namespace
{
    std::unique_ptr< geode::TriangulatedSurface3D > load_test_file(
        std::string_view filename, std::string_view content )
    {
        {
            std::ofstream file{ geode::to_string( filename ),
                std::ios::binary };
            file << content;
        }
        return geode::load_triangulated_surface< 3 >( filename );
    }

    void test_blank_lines()
    {
        const auto surface = load_test_file( "test_ts_blank_lines.ts",
            "GOCAD TSurf 1\nHEADER {\nname: blank\n}\n"
            "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n\nNAME Default\nZPOSITIVE\n"
            "END_ORIGINAL_COORDINATE_SYSTEM\nTFACE\n\n# comment\n"
            "VRTX 1 0 0 0\nVRTX 2 1 0 0\n   \nVRTX 3 0 1 0\nTRGL 1 2\n"
            "TRGL 1 2 3\nBSTONE\nBORDER 1 2\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 3,
            "Wrong number of vertices with blank lines" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_polygons() == 1,
            "Wrong number of triangles with blank lines" );
    }
} // namespace

void check_surface( const geode::SurfaceMesh3D& surface,
    geode::index_t nb_vertices,
    geode::index_t nb_polygons,
    std::string name )
{
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface.nb_vertices() == nb_vertices,
        "Number of vertices in the TSurf 3D is not correct" );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface.nb_polygons() == nb_polygons,
        "Number of polygons in the TSurf 3D is not correct" );
    geode::OpenGeodeGeosciencesIOMeshException::test( surface.name() == name,
        "Wrong name: ", surface.name().value(), " should be ", name );
}

void check_file( std::string file,
    geode::index_t nb_vertices,
    geode::index_t nb_polygons,
    std::string name )
{
    // Load file
    auto surface = geode::load_triangulated_surface< 3 >( file );
    check_surface( *surface, nb_vertices, nb_polygons, name );

    // Save triangulated tsurf
    const auto output_file_native =
        absl::StrCat( "test_output.", surface->native_extension() );
    geode::save_triangulated_surface( *surface, output_file_native );
    const auto output_file_ts = "test_output.ts";
    geode::save_triangulated_surface( *surface, output_file_ts );

    // Load file
    auto reloaded_surface =
        geode::load_triangulated_surface< 3 >( output_file_native );
    check_surface( *reloaded_surface, nb_vertices, nb_polygons, name );
    auto reloaded_surface_ts =
        geode::load_triangulated_surface< 3 >( output_file_ts );
    check_surface( *reloaded_surface_ts, nb_vertices, nb_polygons, name );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::trace );
        check_file( absl::StrCat( geode::DATA_PATH, "/surf2d_multi.",
                        geode::internal::TSInput::extension() ),
            92, 92, "section1" );
        check_file( absl::StrCat( geode::DATA_PATH, "/surf2d.",
                        geode::internal::TSInput::extension() ),
            46, 46, "section1" );
        check_file( absl::StrCat( geode::DATA_PATH, "/2triangles.",
                        geode::internal::TSInput::extension() ),
            4, 2, "2triangles.ts" );
        check_file( absl::StrCat( geode::DATA_PATH, "/sgrid_tsurf.",
                        geode::internal::TSInput::extension() ),
            4, 2, "FractureSet_N20_1fract_visualisation" );
        check_file( absl::StrCat( geode::DATA_PATH, "Fault_without_crs.",
                        geode::internal::TSInput::extension() ),
            189, 324, "Fault" );
        check_file( absl::StrCat( geode::DATA_PATH, "ts-2props.",
                        geode::internal::TSInput::extension() ),
            4, 2, "test" );
        test_blank_lines();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
