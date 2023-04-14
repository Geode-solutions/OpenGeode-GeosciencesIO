/*
 * Copyright (c) 2019 - 2023 Geode-solutions
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

#include <geode/tests_config.h>

#include <geode/basic/assert.h>
#include <geode/basic/attribute_manager.h>
#include <geode/basic/logger.h>

#include <geode/mesh/core/regular_grid_solid.h>
#include <geode/mesh/io/regular_grid_input.h>

#include <geode/geosciences_io/mesh/common.h>

void test_grid_input()
{
    const auto grid = geode::load_regular_grid< 3 >(
        absl::StrCat( geode::data_path, "test.vo" ) );

    auto attribute_to_test =
        grid->cell_attribute_manager().find_attribute< double >( "random" );

    const auto first_value =
        attribute_to_test->value( grid->cell_index( { 0, 0, 0 } ) );
    OPENGEODE_EXCEPTION( first_value == 6.48414,
        "[TEST] Error in grid attributes, value for attribute 'random' at "
        "cell [0,0,0] is ",
        first_value, " where it should be 6.48414" );

    const auto second_value =
        attribute_to_test->value( grid->cell_index( { 5, 0, 9 } ) );
    OPENGEODE_EXCEPTION( second_value == 8.95907,
        "[TEST] Error in grid attributes, value for attribute 'random' at "
        "cell [5,0,9] is ",
        second_value, " where it should be 8.95907" );

    const auto third_value =
        attribute_to_test->value( grid->cell_index( { 9, 9, 9 } ) );
    OPENGEODE_EXCEPTION( third_value == 7.21909,
        "[TEST] Error in grid attributes, value for attribute 'random' at "
        "cell [9,9,9] is ",
        third_value, " where it should be 7.21909" );
}

int main()
{
    try
    {
        geode::GeosciencesIOMeshLibrary::initialize();
        test_grid_input();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
