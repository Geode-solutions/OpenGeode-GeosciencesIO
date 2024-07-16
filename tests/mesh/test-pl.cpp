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
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/geode/geode_edged_curve_builder.hpp>
#include <geode/mesh/core/geode/geode_edged_curve.hpp>
#include <geode/mesh/io/edged_curve_input.hpp>
#include <geode/mesh/io/edged_curve_output.hpp>

#include <geode/geosciences_io/mesh/internal/pl_input.hpp>
#include <geode/geosciences_io/mesh/internal/pl_output.hpp>

namespace
{
    void check_curve( const geode::EdgedCurve3D& curve,
        geode::index_t nb_vertices,
        geode::index_t nb_edges )
    {
        OPENGEODE_EXCEPTION( curve.nb_vertices() == nb_vertices,
            "Number of vertices in the ECurve 3D is not correct" );
        OPENGEODE_EXCEPTION( curve.nb_edges() == nb_edges,
            "Number of edges in the ECurve 3D is not correct" );
    }

    void check_file( std::string file,
        geode::index_t nb_vertices,
        geode::index_t nb_edges,
        std::string output_file )
    {
        // Load file
        auto curve = geode::load_edged_curve< 3 >( file );
        check_curve( *curve, nb_vertices, nb_edges );

        // Save Edged Curve
        const auto output_file_native =
            absl::StrCat( output_file, curve->native_extension() );
        geode::save_edged_curve( *curve, output_file_native );
        const auto output_file_pl = "test_output.pl";
        geode::save_edged_curve( *curve, output_file_pl );

        // Load file
        auto reloaded_curve =
            geode::load_edged_curve< 3 >( output_file_native );
        check_curve( *reloaded_curve, nb_vertices, nb_edges );
        auto reloaded_curve_ec = geode::load_edged_curve< 3 >( output_file_pl );
        check_curve( *reloaded_curve_ec, nb_vertices, nb_edges );
        save_edged_curve(
            *reloaded_curve, absl::StrCat( "reloaded_", output_file_native ) );
    }
} // namespace

int main()
{
    try
    {
        geode::GeosciencesIOMeshLibrary::initialize();
        // NOLINTBEGIN(*-magic-numbers)
        check_file( absl::StrCat( geode::DATA_PATH, "/normal_lines.",
                        geode::internal::PLInput::extension() ),
            11391, 11374, "normal_lines." );
        check_file( absl::StrCat( geode::DATA_PATH, "/closed_lines.",
                        geode::internal::PLInput::extension() ),
            9395, 9395, "closed_lines." );
        // NOLINTEND(*-magic-numbers)
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}