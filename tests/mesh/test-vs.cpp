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

#include <geode/tests_config.hpp>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>

#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/io/point_set_input.hpp>
#include <geode/mesh/io/point_set_output.hpp>

#include <geode/geosciences_io/mesh/internal/vs_input.hpp>

void check_pointset(
    const geode::PointSet3D& pointset, geode::index_t nb_vertices )
{
    OPENGEODE_EXCEPTION( pointset.nb_vertices() == nb_vertices,
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
        geode::GeosciencesIOMeshLibrary::initialize();
        check_file( absl::StrCat( geode::DATA_PATH, "points.",
                        geode::internal::VSInput::extension() ),
            6 );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
