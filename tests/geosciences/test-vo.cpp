/*
 * Copyright (c) 2019 - 2022 Geode-solutions
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
#include <stdexcept>
#include <iostream>
#include <string.h>

#include <geode/tests_config.h>

#include <geode/basic/assert.h>
#include <geode/basic/logger.h>
#include <geode/basic/attribute_manager.h>

#include <geode/mesh/builder/regular_grid_solid_builder.h>
#include <geode/mesh/core/regular_grid_solid.h>
#include <geode/mesh/io/regular_grid_input.h>

#include <geode/io/geosciences/common.h>

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOGeosciences::initialize();
        auto grid = geode::load_regular_grid< 3 >(
            absl::StrCat( geode::data_path, "test.vo" ) );
       
        auto attribute_random = grid->cell_attribute_manager().find_attribute<double>("random");
        
        if(attribute_random->value(grid->cell_index({0,0,0})) != 6.48414) {
            throw std::runtime_error( "Attribute creation has failed somewhere, [0,0,0] value for 'random' is "
                + std::to_string(attribute_random->value(grid->cell_index({0,0,0})))
                + " and it should be 6.48414" );
        }
        if(attribute_random->value(grid->cell_index({5,0,9})) != 8.95907) {
            throw std::runtime_error( "Attribute creation has failed somewhere, [5,0,9] value for 'random' is "
                + std::to_string(attribute_random->value(grid->cell_index({5,0,9})))
                + " and it should be 8.95907" );
        }
        if(attribute_random->value(grid->cell_index({9,9,9})) != 7.21909) {
            throw std::runtime_error( "Attribute creation has failed somewhere, [9,9,9] value for 'random' is "
                + std::to_string(attribute_random->value(grid->cell_index({9,9,9})))
                + " and it should be 7.21909" );
        }

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
