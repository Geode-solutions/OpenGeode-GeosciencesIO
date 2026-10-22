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

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/io/edged_curve_input.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

namespace
{
    std::unique_ptr< geode::EdgedCurve3D > load_test_file(
        std::string_view filename, std::string_view content )
    {
        {
            std::ofstream file{ geode::to_string( filename ),
                std::ios::binary };
            file << content;
        }
        return geode::load_edged_curve< 3 >( filename );
    }

    void check_curve( const geode::EdgedCurve3D& curve,
        geode::index_t nb_vertices,
        geode::index_t nb_edges,
        std::string_view test_name )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.nb_vertices() == nb_vertices,
            "Wrong number of vertices: ", test_name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.nb_edges() == nb_edges,
            "Wrong number of edges: ", test_name );
    }

    void test_few_points()
    {
        check_curve(
            *load_test_file( "test_well_empty.txt", "" ), 0, 0, "empty well" );
        check_curve( *load_test_file( "test_well_one_point.txt", "1 2 3\n" ), 1,
            0, "one point well" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        auto curve = geode::load_edged_curve< 3 >(
            absl::StrCat( geode::DATA_PATH, "test_well.txt" ) );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 11, "Wrong number of vertices" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_edges() == 10, "Wrong number of edges" );
        test_few_points();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
