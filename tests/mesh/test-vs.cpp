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

#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/io/point_set_input.hpp>
#include <geode/mesh/io/point_set_output.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>
#include <geode/geosciences_io/mesh/internal/vs_input.hpp>

namespace
{
    std::string write_test_file(
        std::string_view filename, std::string_view content )
    {
        std::ofstream file{ geode::to_string( filename ), std::ios::binary };
        file << content;
        return geode::to_string( filename );
    }

    std::unique_ptr< geode::PointSet3D > load_test_file(
        std::string_view filename, std::string_view content )
    {
        return geode::load_point_set< 3 >(
            write_test_file( filename, content ) );
    }

    void test_blank_lines()
    {
        const auto pointset = load_test_file( "test_vs_blank_lines.vs",
            "GOCAD VSet 1\nHEADER {\nname: a\n}\nVRTX 1 0 0 0\n\n"
            "# comment\nVRTX 2 1 1 1\n\t\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            pointset->nb_vertices() == 2,
            "Wrong number of vertices with blank lines" );
    }

    void test_quotes()
    {
        const auto filename = write_test_file( "test_vs_quotes.txt",
            "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
            "NAME \"SKUA Local\"\n"
            "DATUM \" a  b \"\n"
            "PROJECTION \"\n"
            "END_ORIGINAL_COORDINATE_SYSTEM\n" );
        std::ifstream file{ filename, std::ios::binary };
        const auto crs = geode::internal::read_CRS( file );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            crs.name == "SKUA Local", "Wrong CRS name: ", crs.name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            crs.datum == "a b", "Wrong CRS datum: ", crs.datum );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            crs.projection.empty(), "Wrong CRS projection: ", crs.projection );

        const auto pointset = load_test_file( "test_vs_unclosed_quote.vs",
            "GOCAD VSet 1\nHEADER {\nname: a\n}\nPROPERTIES \"\n"
            "PVRTX 1 0 0 0 4\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            pointset->nb_vertices() == 1
                && pointset->vertex_attribute_manager()
                    .attribute_ids_matching_name( "property_0" )
                    .has_value(),
            "Wrong point set with unclosed quote" );
    }
} // namespace

void check_pointset(
    const geode::PointSet3D& pointset, geode::index_t nb_vertices )
{
    geode::OpenGeodeGeosciencesIOMeshException::test(
        pointset.nb_vertices() == nb_vertices,
        "Number of vertices in the VSet 3D is not correct" );
}

void check_file( std::string file, geode::index_t nb_vertices )
{
    // Load file
    auto pointset = geode::load_point_set< 3 >( file );
    check_pointset( *pointset, nb_vertices );

    // Save point set
    const auto output_file_native =
        absl::StrCat( "test_output.", pointset->native_extension() );
    geode::save_point_set( *pointset, output_file_native );
    const auto output_file_vs = "test_output.vs";
    geode::save_point_set( *pointset, output_file_vs );

    // Load file
    auto reloaded_pointset = geode::load_point_set< 3 >( output_file_native );
    check_pointset( *reloaded_pointset, nb_vertices );
    auto reloaded_pointset_vs = geode::load_point_set< 3 >( output_file_vs );
    check_pointset( *reloaded_pointset_vs, nb_vertices );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        check_file( absl::StrCat( geode::DATA_PATH, "points.",
                        geode::internal::VSInput::extension() ),
            6 );
        test_quotes();
        test_blank_lines();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
