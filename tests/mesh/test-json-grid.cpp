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
#include <geode/basic/range.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/light_regular_grid.hpp>
#include <geode/mesh/io/light_regular_grid_input.hpp>
#include <geode/mesh/io/light_regular_grid_output.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

template < typename T >
std::shared_ptr< geode::ReadOnlyAttribute< T > > find_attribute(
    const geode::LightRegularGrid3D& grid, std::string_view name )
{
    const auto ids =
        grid.grid_vertex_attribute_manager().attribute_ids_matching_name(
            name );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        ids.has_value() && ids->size() == 1, "Attribute ", name,
        " should exist once" );
    auto attribute =
        grid.grid_vertex_attribute_manager().find_read_only_attribute< T >(
            ids->front() );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        attribute != nullptr, "Attribute ", name, " has wrong type" );
    return attribute;
}

void check_vertex( const geode::LightRegularGrid3D& grid,
    const geode::LightRegularGrid3D::VertexIndices& indices,
    const geode::Point3D& point,
    bool topo,
    bool site,
    double layer_id,
    double pollutant )
{
    const auto vertex = grid.vertex_index( indices );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        grid.point( vertex ).inexact_equal( point ), "Wrong point for vertex ",
        vertex );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        find_attribute< bool >( grid, "topo" )->value( vertex ) == topo,
        "Wrong topo value for vertex ", vertex );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        find_attribute< bool >( grid, "site" )->value( vertex ) == site,
        "Wrong site value for vertex ", vertex );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        find_attribute< double >( grid, "layer_id" )->value( vertex )
            == layer_id,
        "Wrong layer_id value for vertex ", vertex );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        find_attribute< double >( grid, "pollutant" )->value( vertex )
            == pollutant,
        "Wrong pollutant value for vertex ", vertex );
}

void check_grid( const geode::LightRegularGrid3D& grid )
{
    geode::OpenGeodeGeosciencesIOMeshException::test(
        grid.nb_cells_in_direction( 0 ) == 19
            && grid.nb_cells_in_direction( 1 ) == 19
            && grid.nb_cells_in_direction( 2 ) == 13,
        "Wrong number of cells" );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        grid.nb_grid_vertices() == 5600, "Wrong number of vertices" );
    for( const auto& name : { "x", "y", "z" } )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !grid.grid_vertex_attribute_manager()
                .attribute_ids_matching_name( name )
                .has_value(),
            "Coordinate ", name, " should not be an attribute" );
    }
    check_vertex( grid, { 0, 0, 0 }, geode::Point3D{ { 2.5, 2.5, 10 } }, true,
        false, 1, 4.2169650342859424e-132 );
    check_vertex( grid, { 0, 0, 8 }, geode::Point3D{ { 2.5, 2.5, 18 } }, false,
        false, 0, 4.2169650342859424e-132 );
    check_vertex( grid, { 9, 1, 0 }, geode::Point3D{ { 47.5, 7.5, 10 } }, true,
        true, 1, 4.216965034285796e-51 );
    check_vertex( grid, { 5, 3, 7 }, geode::Point3D{ { 27.5, 17.5, 17 } }, true,
        false, 4, 1.3335214321633148e-43 );
}

void write_file( std::string_view filename, std::string_view content )
{
    std::ofstream file{ geode::to_string( filename ) };
    file << content;
}

void test_unordered_records()
{
    constexpr auto filename = "unordered_grid.json";
    write_file( filename, R"({
        "grid_df": [
            { "x": 2, "y": 1, "z": 6, "value": 112 },
            { "x": 0, "y": 0, "z": 5, "value": 0 },
            { "x": 0, "y": 1, "z": 6, "value": 110 },
            { "x": 2, "y": 0, "z": 5, "value": 2 },
            { "x": 0, "y": 1, "z": 5, "value": 10 },
            { "x": 2, "y": 0, "z": 6, "value": 102 },
            { "x": 2, "y": 1, "z": 5, "value": 12 },
            { "x": 0, "y": 0, "z": 6, "value": 100 }
        ],
        "name": "unordered"
    })" );
    const auto grid = geode::load_light_regular_grid< 3 >( filename );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        grid.nb_cells_in_direction( 0 ) == 1
            && grid.nb_cells_in_direction( 1 ) == 1
            && grid.nb_cells_in_direction( 2 ) == 1,
        "Wrong number of cells for unordered grid" );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        grid.name() == "unordered", "Wrong name for unordered grid" );
    const auto attribute = find_attribute< double >( grid, "value" );
    for( const auto vertex : geode::Range{ grid.nb_grid_vertices() } )
    {
        const auto& point = grid.point( vertex );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            attribute->value( vertex )
                == point.value( 0 ) + 10. * point.value( 1 )
                       + 100. * ( point.value( 2 ) - 5. ),
            "Wrong value for unordered vertex ", vertex );
    }
}

void test_invalid_grid( std::string_view filename, std::string_view content )
{
    write_file( filename, content );
    try
    {
        [[maybe_unused]] const auto grid =
            geode::load_light_regular_grid< 3 >( filename );
    }
    catch( const geode::OpenGeodeException& )
    {
        return;
    }
    geode::OpenGeodeGeosciencesIOMeshException::test(
        false, "Invalid grid ", filename, " should not be loaded" );
}

void test_invalid_grids()
{
    test_invalid_grid( "duplicated_vertex_grid.json", R"({
        "grid_df": [
            { "x": 0, "y": 0, "z": 0 },
            { "x": 1, "y": 0, "z": 0 },
            { "x": 0, "y": 1, "z": 0 },
            { "x": 1, "y": 1, "z": 0 },
            { "x": 0, "y": 0, "z": 1 },
            { "x": 1, "y": 0, "z": 1 },
            { "x": 0, "y": 1, "z": 1 },
            { "x": 0, "y": 1, "z": 1 }
        ]
    })" );
    test_invalid_grid( "irregular_grid.json", R"({
        "grid_df": [
            { "x": 0, "y": 0, "z": 0 },
            { "x": 1, "y": 0, "z": 0 },
            { "x": 3, "y": 0, "z": 0 },
            { "x": 0, "y": 1, "z": 0 },
            { "x": 1, "y": 1, "z": 0 },
            { "x": 3, "y": 1, "z": 0 },
            { "x": 0, "y": 0, "z": 1 },
            { "x": 1, "y": 0, "z": 1 },
            { "x": 3, "y": 0, "z": 1 },
            { "x": 0, "y": 1, "z": 1 },
            { "x": 1, "y": 1, "z": 1 },
            { "x": 3, "y": 1, "z": 1 }
        ]
    })" );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::trace );

        const auto grid = geode::load_light_regular_grid< 3 >(
            absl::StrCat( geode::DATA_PATH, "data.json" ) );
        check_grid( grid );

        const auto output =
            absl::StrCat( "json_grid.", grid.native_extension() );
        geode::save_light_regular_grid( grid, output );
        check_grid( geode::load_light_regular_grid< 3 >( output ) );

        test_unordered_records();
        test_invalid_grids();

        geode::Logger::info( "[TEST SUCCESS]" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
