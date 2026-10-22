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
// #include <iostream>
// #include <stdexcept>
// #include <string.h>

#include <geode/tests_config.hpp>

#include <fstream>

#include <geode/basic/assert.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/regular_grid_solid.hpp>
#include <geode/mesh/io/regular_grid_input.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

void test_grid_input()
{
    const auto grid = geode::load_regular_grid< 3 >(
        absl::StrCat( geode::DATA_PATH, "test.vo" ) );
    auto attribute_ids =
        grid->cell_attribute_manager().attribute_ids_matching_name( "random" );
    auto attribute_to_test =
        grid->cell_attribute_manager().find_read_only_attribute< double >(
            attribute_ids.value().front() );

    const auto first_value =
        attribute_to_test->value( grid->cell_index( { 0, 0, 0 } ) );
    geode::OpenGeodeGeosciencesIOMeshException::test( first_value == 6.48414,
        "[TEST] Error in grid attributes, value for attribute 'random' at "
        "cell [0,0,0] is ",
        first_value, " where it should be 6.48414" );

    const auto second_value =
        attribute_to_test->value( grid->cell_index( { 5, 0, 9 } ) );
    geode::OpenGeodeGeosciencesIOMeshException::test( second_value == 8.95907,
        "[TEST] Error in grid attributes, value for attribute 'random' at "
        "cell [5,0,9] is ",
        second_value, " where it should be 8.95907" );

    const auto third_value =
        attribute_to_test->value( grid->cell_index( { 9, 9, 9 } ) );
    geode::OpenGeodeGeosciencesIOMeshException::test( third_value == 7.21909,
        "[TEST] Error in grid attributes, value for attribute 'random' at "
        "cell [9,9,9] is ",
        third_value, " where it should be 7.21909" );
}

namespace
{
    std::unique_ptr< geode::RegularGrid3D > load_test_voxet(
        std::string_view name,
        std::string_view geometry,
        std::string_view data )
    {
        const auto data_name = absl::StrCat( name, "__ascii@@" );
        {
            std::ofstream file{ data_name, std::ios::binary };
            file << data;
        }
        const auto filename = absl::StrCat( name, ".vo" );
        {
            std::ofstream file{ filename, std::ios::binary };
            file << "GOCAD Voxet 1\nHEADER {\nname: " << name << "\n}\n"
                 << geometry << "ASCII_DATA_FILE " << data_name
                 << "\nPROPERTY 1 p\nEND\n";
        }
        return geode::load_regular_grid< 3 >( filename );
    }

    constexpr auto SIMPLE_GEOMETRY =
        "AXIS_O 0 0 0\nAXIS_U 1 0 0\nAXIS_V 0 1 0\nAXIS_W 0 0 1\n"
        "AXIS_MIN 0 0 0\nAXIS_MAX 1 1 1\nAXIS_N 2 2 2\n";

    double cell_value( const geode::RegularGrid3D& grid,
        geode::index_t i,
        geode::index_t j,
        geode::index_t k )
    {
        const auto& manager = grid.cell_attribute_manager();
        const auto attribute = manager.find_read_only_attribute< double >(
            manager.attribute_ids_matching_name( "p" ).value().front() );
        return attribute->value( grid.cell_index( { i, j, k } ) );
    }

    void test_cells_outside_grid()
    {
        const auto grid = load_test_voxet( "test_vo_outside", SIMPLE_GEOMETRY,
            "*\n* I J K p\n*\n0 0 0 1\n2 0 0 2\n0 5 0 3\n0 0 1000000 4\n"
            "1 1 1 5\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            cell_value( *grid, 0, 0, 0 ) == 1
                && cell_value( *grid, 1, 1, 1 ) == 5,
            "[TEST] Wrong values with cells outside of the grid" );
    }

    void test_data_file_format()
    {
        const auto grid = load_test_voxet( "test_vo_blank", SIMPLE_GEOMETRY,
            "*\n\n# comment\n* I J K p\n*\n0 0 0 1\n\n1 1 1 5\n0 1\n"
            "1 0 0 2\n\n\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            cell_value( *grid, 0, 0, 0 ) == 1
                && cell_value( *grid, 1, 1, 1 ) == 5
                && cell_value( *grid, 1, 0, 0 ) == 2,
            "[TEST] Wrong values with blank lines in data file" );

        auto has_thrown = false;
        try
        {
            load_test_voxet(
                "test_vo_no_columns", SIMPLE_GEOMETRY, "*\n*\n*\n0 0 0 1\n" );
        }
        catch( const geode::OpenGeodeException& exception )
        {
            geode::Logger::info( "Expected exception: ", exception.what() );
            has_thrown = true;
        }
        geode::OpenGeodeGeosciencesIOMeshException::test(
            has_thrown, "[TEST] Data file without column header should throw" );
    }

    void test_geometry()
    {
        const auto grid = geode::load_regular_grid< 3 >(
            absl::StrCat( geode::DATA_PATH, "test.vo" ) );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            grid->cell_barycenter( { 0, 0, 0 } )
                    .inexact_equal( geode::Point3D{ { 0, 0, 0 } } )
                && grid->cell_barycenter( { 9, 9, 9 } )
                    .inexact_equal( geode::Point3D{ { 900, 900, 900 } } ),
            "[TEST] Wrong test.vo geometry" );

        const auto rotated = load_test_voxet( "test_vo_rotated",
            "GOCAD_ORIGINAL_COORDINATE_SYSTEM\nZPOSITIVE Depth\n"
            "END_ORIGINAL_COORDINATE_SYSTEM\nAXIS_O 10 20 30\nAXIS_U 0 4 0\n"
            "AXIS_V -4 0 0\nAXIS_W 0 0 8\nAXIS_MIN 0 0 0.5\nAXIS_MAX 1 1 1\n"
            "AXIS_N 3 5 2\n",
            "*\n* I J K p\n*\n0 0 0 1\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            rotated->nb_cells_in_direction( 0 ) == 3
                && rotated->nb_cells_in_direction( 1 ) == 5
                && rotated->nb_cells_in_direction( 2 ) == 2,
            "[TEST] Wrong rotated voxet size" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            rotated->cell_barycenter( { 0, 0, 0 } )
                    .inexact_equal( geode::Point3D{ { 10, 20, -34 } } )
                && rotated->cell_barycenter( { 2, 4, 1 } )
                    .inexact_equal( geode::Point3D{ { 6, 24, -38 } } ),
            "[TEST] Wrong rotated voxet geometry" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::debug );
        test_grid_input();
        test_cells_outside_grid();
        test_geometry();
        test_data_file_format();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
