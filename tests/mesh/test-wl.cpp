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

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/io/edged_curve_input.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

namespace
{
    void test_first_path()
    {
        {
            std::ofstream file{ "test_wl_first_path.wl", std::ios::binary };
            file << "GOCAD Well 1\nHEADER {\nname: w\n}\n"
                    "WREF 10 20 0\nPATH 5 5 1 2\nPATH 10 10 1 2\nEND\n";
        }
        const auto curve =
            geode::load_edged_curve< 3 >( "test_wl_first_path.wl" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 3 && curve->nb_edges() == 2,
            "First PATH point should not be skipped" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->point( 1 ) == geode::Point3D{ { 11, 22, 5 } },
            "Wrong first PATH point" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        auto curve = geode::load_edged_curve< 3 >(
            absl::StrCat( geode::DATA_PATH, "test.wl" ) );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 55, "Wrong number of vertices" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_edges() == 54, "Wrong number of edges" );
        test_first_path();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
