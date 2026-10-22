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
        geode::index_t nb_edges );

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

    void test_blank_lines()
    {
        const auto curve = load_test_file( "test_pl_blank_lines.pl",
            "GOCAD PLine 1\nHEADER {\nname: a\n}\nILINE\n\n# comment\n"
            "VRTX 1 0 0 0\n\nVRTX 2 1 0 0\nSEG 1\nSEG 1 2\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 2,
            "Wrong number of vertices with blank lines" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_edges() == 1, "Wrong number of edges with blank lines" );
    }

    void test_invalid_segments()
    {
        const auto curve = load_test_file( "test_pl_invalid_segments.pl",
            "GOCAD PLine 1\nHEADER {\nname: a\n}\nILINE\n"
            "VRTX 1 0 0 0\nVRTX 2 1 0 0\nVRTX 3 2 0 0\nSEG 1 5\nSEG 0 1\n"
            "SEG 2 2\nSEG 1 2\nSEG 2 3\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 3,
            "Wrong number of vertices with invalid segments" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_edges() == 2,
            "Wrong number of edges with invalid segments" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->edge_vertices( 1 )[0] == 1
                && curve->edge_vertices( 1 )[1] == 2,
            "Wrong edge vertices with invalid segments" );
    }

    void test_output_precision()
    {
        auto curve = geode::EdgedCurve3D::create();
        auto builder = geode::EdgedCurveBuilder3D::create( *curve );
        const geode::Point3D point0{ { 1056340.123456789, 0.1, -2.5e-7 } };
        const geode::Point3D point1{ { 1.0 / 3.0, 7654321.987654321, 3 } };
        builder->create_point( point0 );
        builder->create_point( point1 );
        builder->create_edge( 0, 1 );
        geode::save_edged_curve( *curve, "test_pl_precision.pl" );
        const auto reloaded =
            geode::load_edged_curve< 3 >( "test_pl_precision.pl" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            reloaded->nb_vertices() == 2 && reloaded->point( 0 ) == point0
                && reloaded->point( 1 ) == point1,
            "Coordinates are not saved with full precision" );
    }

    void test_non_contiguous_ids()
    {
        const auto curve = load_test_file( "test_pl_sparse_ids.pl",
            "GOCAD PLine 1\nHEADER {\nname: a\n}\nILINE\nVRTX 5 0 0 0\n"
            "VRTX 3 1 0 0\nVRTX 8 2 0 0\nSEG 5 3\nSEG 3 8\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 3 && curve->nb_edges() == 2
                && curve->edge_vertices( 1 )[0] == 1
                && curve->edge_vertices( 1 )[1] == 2,
            "Wrong curve with non contiguous ids" );
    }

    double vertex_value( const geode::EdgedCurve3D& curve,
        std::string_view name,
        geode::index_t vertex )
    {
        const auto& manager = curve.vertex_attribute_manager();
        const auto ids = manager.attribute_ids_matching_name( name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            ids.has_value() && ids->size() == 1,
            "Attribute should exist once: ", name );
        return manager.find_read_only_attribute< double >( ids->front() )
            ->value( vertex );
    }

    void test_properties_and_atoms()
    {
        const auto curve = load_test_file( "test_pl_properties.pl",
            "GOCAD PLine 1\nHEADER {\nname: a\n}\nPROPERTIES p\nESIZES 1\n"
            "ILINE\nPVRTX 1 0 0 0 1.5\nPVRTX 2 1 0 0 2.5\nSEG 1 2\nILINE\n"
            "ATOM 3 2\nPVRTX 4 2 0 0 3.5\nSEG 3 4\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 3 && curve->nb_edges() == 2
                && curve->edge_vertices( 1 )[0] == 1
                && curve->edge_vertices( 1 )[1] == 2,
            "Wrong curve with ATOM" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_value( *curve, "p", 0 ) == 1.5
                && vertex_value( *curve, "p", 2 ) == 3.5,
            "Wrong property values in PLine" );

        geode::save_edged_curve( *curve, "test_pl_properties_saved.pl" );
        const auto reloaded =
            geode::load_edged_curve< 3 >( "test_pl_properties_saved.pl" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_value( *reloaded, "p", 1 ) == 2.5,
            "Wrong property values in saved PLine" );
    }

    std::vector< std::array< geode::Point3D, 2 > > sorted_edge_points(
        const geode::EdgedCurve3D& curve )
    {
        std::vector< std::array< geode::Point3D, 2 > > edges;
        for( const auto e : geode::Range{ curve.nb_edges() } )
        {
            std::array< geode::Point3D, 2 > points{
                curve.point( curve.edge_vertices( e )[0] ),
                curve.point( curve.edge_vertices( e )[1] )
            };
            if( points[1] < points[0] )
            {
                std::swap( points[0], points[1] );
            }
            edges.push_back( points );
        }
        absl::c_sort( edges, []( const auto& lhs, const auto& rhs ) {
            return lhs[0] < rhs[0] || ( lhs[0] == rhs[0] && lhs[1] < rhs[1] );
        } );
        return edges;
    }

    void test_closed_loop_output()
    {
        auto curve = geode::EdgedCurve3D::create();
        auto builder = geode::EdgedCurveBuilder3D::create( *curve );
        builder->create_point( geode::Point3D{ { 0, 0, 0 } } );
        builder->create_point( geode::Point3D{ { 1, 0, 0 } } );
        builder->create_point( geode::Point3D{ { 0, 1, 0 } } );
        builder->create_point( geode::Point3D{ { 5, 5, 5 } } );
        builder->create_point( geode::Point3D{ { 6, 5, 5 } } );
        builder->create_edge( 0, 1 );
        builder->create_edge( 1, 2 );
        builder->create_edge( 2, 0 );
        builder->create_edge( 3, 4 );
        geode::save_edged_curve( *curve, "test_pl_closed_loop.pl" );
        const auto reloaded =
            geode::load_edged_curve< 3 >( "test_pl_closed_loop.pl" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            sorted_edge_points( *reloaded ) == sorted_edge_points( *curve ),
            "Closed loop is not saved correctly" );
    }

    void test_shared_and_isolated_vertices_output()
    {
        auto curve = geode::EdgedCurve3D::create();
        auto builder = geode::EdgedCurveBuilder3D::create( *curve );
        builder->create_point( geode::Point3D{ { 0, 0, 0 } } );
        builder->create_point( geode::Point3D{ { 1, 0, 0 } } );
        builder->create_point( geode::Point3D{ { 2, 0, 0 } } );
        builder->create_point( geode::Point3D{ { 2, 1, 0 } } );
        builder->create_point( geode::Point3D{ { 9, 9, 9 } } );
        builder->create_edge( 0, 1 );
        builder->create_edge( 1, 2 );
        builder->create_edge( 2, 3 );
        builder->create_edge( 3, 1 );
        builder->set_name( "my lollipop" );
        geode::save_edged_curve( *curve, "test_pl_lollipop.pl" );
        const auto reloaded =
            geode::load_edged_curve< 3 >( "test_pl_lollipop.pl" );
        check_curve( *reloaded, 5, 4 );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            sorted_edge_points( *reloaded ) == sorted_edge_points( *curve ),
            "Shared vertices are not saved correctly" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            reloaded->edges_around_vertex( 1 ).size() == 3,
            "Shared vertex is disconnected after reload" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            reloaded->name() == "my lollipop", "Wrong saved PLine name" );
    }

    void check_curve( const geode::EdgedCurve3D& curve,
        geode::index_t nb_vertices,
        geode::index_t nb_edges )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.nb_vertices() == nb_vertices,
            "Number of vertices in the ECurve 3D is not correct" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.nb_edges() == nb_edges,
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
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        // NOLINTBEGIN(*-magic-numbers)
        check_file( absl::StrCat( geode::DATA_PATH, "/normal_lines.",
                        geode::internal::PLInput::extension() ),
            11391, 11374, "normal_lines." );
        check_file( absl::StrCat( geode::DATA_PATH, "/closed_lines.",
                        geode::internal::PLInput::extension() ),
            9395, 9395, "closed_lines." );
        // NOLINTEND(*-magic-numbers)
        test_blank_lines();
        test_invalid_segments();
        test_output_precision();
        test_closed_loop_output();
        test_non_contiguous_ids();
        test_properties_and_atoms();
        test_shared_and_isolated_vertices_output();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}