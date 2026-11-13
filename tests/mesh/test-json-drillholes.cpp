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

#include <cmath>
#include <fstream>

#include <geode/basic/assert.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/coordinate_reference_system_manager.hpp>
#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/io/edged_curve_input.hpp>

#include <geode/geosciences/explicit/geometry/geographic_coordinate_system.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

namespace
{
    using Curves = std::vector< std::unique_ptr< geode::EdgedCurve3D > >;

    std::shared_ptr< geode::ReadOnlyAttribute< double > > find_attribute(
        const geode::EdgedCurve3D& curve, std::string_view name )
    {
        const auto ids =
            curve.vertex_attribute_manager().attribute_ids_matching_name(
                name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            ids.has_value() && ids->size() == 1, "Attribute ", name,
            " should exist once" );
        auto attribute =
            curve.vertex_attribute_manager().find_read_only_attribute< double >(
                ids->front() );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            attribute != nullptr, "Attribute ", name, " has wrong type" );
        return attribute;
    }

    void check_no_attribute(
        const geode::EdgedCurve3D& curve, std::string_view name )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !curve.vertex_attribute_manager()
                .attribute_ids_matching_name( name )
                .has_value(),
            "Attribute ", name, " should not exist" );
    }

    void check_curve( const geode::EdgedCurve3D& curve,
        std::string_view name,
        geode::index_t nb_vertices )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.name() == name, "Wrong curve name for ", name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.nb_vertices() == nb_vertices, "Wrong number of vertices for ",
            name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.nb_edges() == nb_vertices - 1, "Wrong number of edges for ",
            name );
        for( const auto edge : geode::Range{ curve.nb_edges() } )
        {
            geode::OpenGeodeGeosciencesIOMeshException::test(
                curve.edge_vertex( { edge, 0 } ) == edge
                    && curve.edge_vertex( { edge, 1 } ) == edge + 1,
                "Wrong edge ", edge, " for ", name );
        }
    }

    void check_coordinate_system(
        const geode::EdgedCurve3D& curve, std::string_view authority_code )
    {
        const auto& crs_manager =
            curve.main_coordinate_reference_system_manager();
        geode::OpenGeodeGeosciencesIOMeshException::test(
            crs_manager.nb_coordinate_reference_systems() == 1
                && crs_manager.active_coordinate_reference_system_name()
                       == authority_code,
            "Wrong curve CRS name" );
        const auto* crs =
            dynamic_cast< const geode::GeographicCoordinateSystem3D* >(
                &crs_manager.active_coordinate_reference_system() );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            crs != nullptr && crs->info().authority_code() == authority_code
                && !crs->info().name.empty(),
            "Wrong curve geographic CRS" );
    }

    void check_no_coordinate_system( const geode::EdgedCurve3D& curve )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.main_coordinate_reference_system_manager()
                    .active_coordinate_reference_system_name()
                == "points",
            "Curve without coordinate system should keep the default CRS" );
    }

    void check_point( const geode::EdgedCurve3D& curve,
        geode::index_t vertex,
        const geode::Point3D& point )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve.point( vertex ).inexact_equal( point ),
            "Wrong point for vertex ", vertex, " of ",
            curve.name().value_or( "" ) );
    }

    void write_file( std::string_view filename, std::string_view content )
    {
        std::ofstream file{ geode::to_string( filename ) };
        file << content;
    }

    void test_data_file()
    {
        const auto filename = absl::StrCat( geode::DATA_PATH, "data.json" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            geode::is_edged_curve_loadable< 3 >( filename ).value() == 1,
            "data.json should be loadable as EdgedCurve3D" );
        const auto curves = geode::load_edged_curve< 3 >( filename );
        constexpr std::array< geode::index_t, 8 > nb_vertices{ 6, 5, 14, 8, 10,
            6, 6, 7 };
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curves.size() == nb_vertices.size(), "Wrong number of drillholes" );
        for( const auto i : geode::Indices{ nb_vertices } )
        {
            const auto& curve = *curves[i];
            check_curve(
                curve, absl::StrCat( "drillhole_", i ), nb_vertices[i] );
            check_coordinate_system( curve, "EPSG:3949" );
            for( const auto name :
                { "x", "y", "z", "anomaly", "zone", "drillhole" } )
            {
                check_no_attribute( curve, name );
            }
            find_attribute( curve, "pollutant" );
        }
        const auto& curve = *curves.front();
        check_point( curve, 0,
            geode::Point3D{ { 59.15448930548192, 52.53876457274344, 18. } } );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            std::fabs( curve.point( 5 ).value( 2 ) - 15.5 ) < 1e-9,
            "Wrong z for last vertex of drillhole_0" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            find_attribute( curve, "pollutant" )->value( 0 )
                == 19.611011754760533,
            "Wrong pollutant value" );
    }

    void test_synthetic_file()
    {
        constexpr auto filename = "drillholes_synthetic.json";
        write_file( filename, R"({
            "drillholes_df": [
                { "hole": "B", "x": 10, "y": 0, "z": 0, "v": 2, "tag": "t" },
                { "hole": "A", "x": 0, "y": 0, "z": 0, "v": 1, "tag": "t" },
                { "hole": "B", "x": 10, "y": 0, "z": -1 },
                { "hole": "A", "x": 0, "y": 0, "v": 99 },
                { "x": 5, "y": 5, "z": 5, "v": 98 },
                { "hole": "A", "x": 0, "y": 0, "z": -1, "v": 3 },
                { "hole": 7, "x": 6, "y": 6, "z": 6, "v": 97 },
                { "hole": "B", "x": 10, "y": 0, "z": -2, "v": 6 },
                { "hole": "A", "x": 0, "y": 0, "z": -2, "v": 4 }
            ],
            "drillhole_col_name": "hole"
        })" );
        const auto curves = geode::load_edged_curve< 3 >( filename );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curves.size() == 2, "Wrong number of synthetic drillholes" );
        const auto& b = *curves[0];
        const auto& a = *curves[1];
        check_curve( b, "B", 3 );
        check_curve( a, "A", 3 );
        check_point( b, 0, geode::Point3D{ { 10, 0, 0 } } );
        check_point( b, 1, geode::Point3D{ { 10, 0, -1 } } );
        check_point( b, 2, geode::Point3D{ { 10, 0, -2 } } );
        check_point( a, 0, geode::Point3D{ { 0, 0, 0 } } );
        check_point( a, 1, geode::Point3D{ { 0, 0, -1 } } );
        check_point( a, 2, geode::Point3D{ { 0, 0, -2 } } );
        const auto v_b = find_attribute( b, "v" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            v_b->value( 0 ) == 2 && std::isnan( v_b->value( 1 ) )
                && v_b->value( 2 ) == 6,
            "Wrong v values for B" );
        const auto v_a = find_attribute( a, "v" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            v_a->value( 0 ) == 1 && v_a->value( 1 ) == 3
                && v_a->value( 2 ) == 4,
            "Wrong v values for A" );
        check_no_attribute( a, "tag" );
        check_no_attribute( a, "hole" );
        check_no_coordinate_system( a );
        check_no_coordinate_system( b );
    }

    void test_default_column_and_invalid_crs()
    {
        constexpr auto filename = "drillholes_default_column.json";
        write_file( filename, R"({
            "drillholes_df": [
                { "drillhole": "A", "x": 0, "y": 0, "z": 0 },
                { "drillhole": "A", "x": 0, "y": 0, "z": -1 }
            ],
            "coordinate_system": "99999"
        })" );
        const auto curves = geode::load_edged_curve< 3 >( filename );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curves.size() == 1, "Wrong number of default column drillholes" );
        check_curve( *curves.front(), "A", 2 );
        check_no_coordinate_system( *curves.front() );
    }

    void test_truncated_after_drillholes()
    {
        constexpr auto filename = "drillholes_truncated.json";
        write_file( filename, R"({
            "drillholes_df": [
                { "drillhole": "A", "x": 0, "y": 0, "z": 0 },
                { "drillhole": "A", "x": 0, "y": 0, "z": -1 }
            ],
            "coordinate_system": )" );
        const auto curves = geode::load_edged_curve< 3 >( filename );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curves.size() == 1, "Drillholes before the parse error are lost" );
    }

    void test_not_loadable(
        std::string_view filename, std::string_view content )
    {
        write_file( filename, content );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            geode::is_edged_curve_loadable< 3 >( filename ).value() == 0,
            filename, " should not be loadable as EdgedCurve3D" );
        try
        {
            [[maybe_unused]] const auto curves =
                geode::load_edged_curve< 3 >( filename );
        }
        catch( const geode::OpenGeodeException& )
        {
            return;
        }
        geode::OpenGeodeGeosciencesIOMeshException::test(
            false, "Invalid file ", filename, " should not be loaded" );
    }

    void test_invalid_files()
    {
        test_not_loadable( "drillholes_missing.json", R"({
            "grid_df": [ { "x": 0, "y": 0, "z": 0 } ]
        })" );
        test_not_loadable(
            "drillholes_empty.json", R"({ "drillholes_df": [] })" );
        test_not_loadable( "drillholes_no_valid_record.json", R"({
            "drillholes_df": [
                { "drillhole": "A", "x": 0, "y": 0 },
                { "x": 0, "y": 0, "z": 0 }
            ]
        })" );
        test_not_loadable( "drillholes_not_an_object.json", "[ 1, 2 ]" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::trace );

        test_data_file();
        test_synthetic_file();
        test_default_column_and_invalid_crs();
        test_truncated_after_drillholes();
        test_invalid_files();

        geode::Logger::info( "[TEST SUCCESS]" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
