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

#include <geode/tests_config.h>

#include <geode/basic/assert.h>
#include <geode/basic/filename.h>
#include <geode/basic/logger.h>

#include <geode/geosciences/explicit/mixin/core/horizon.h>
#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.h>
#include <geode/geosciences/implicit/representation/core/stratigraphic_units_stack.h>
#include <geode/geosciences/implicit/representation/io/stratigraphic_units_stack_output.h>
#include <geode/geosciences_io/model/private/su_stack_skua_input.h>

void test_file()
{
    const auto su_stack = geode::load_stratigraphic_units_stack< 3 >(
        absl::StrCat( geode::data_path, "test_skua_su_stack.xml" ) );

    OPENGEODE_EXCEPTION( su_stack.name() == "skua_model",
        "[TEST] StratigraphicUnitsStack should be named 'skua_model'" );

    OPENGEODE_EXCEPTION( su_stack.nb_horizons() == 4,
        "[TEST] Wrong number of horizons in the "
        "loaded StratigraphicUnitsStack." );
    OPENGEODE_EXCEPTION( su_stack.nb_stratigraphic_units() == 5,
        "[TEST] Wrong number of units in the loaded StratigraphicUnitsStack." );
    geode::uuid erosion_horizon_uuid;
    geode::uuid h1_horizon_uuid;
    for( const auto& horizon : su_stack.horizons() )
    {
        if( horizon.name() == "model_erosion" )
        {
            erosion_horizon_uuid = horizon.id();
        }
        if( horizon.name() == "model_horizon_h1" )
        {
            h1_horizon_uuid = horizon.id();
        }
    }
    geode::uuid eroded_unit_uuid;
    for( const auto& unit : su_stack.stratigraphic_units() )
    {
        if( unit.name() == "eroded_unit" )
        {
            eroded_unit_uuid = unit.id();
        }
    }
    OPENGEODE_EXCEPTION(
        su_stack.is_above( erosion_horizon_uuid, eroded_unit_uuid ),
        "[TEST] Horizon 'model_erosion' should be above unit 'eroded_unit'" );
    OPENGEODE_EXCEPTION( su_stack.is_above( eroded_unit_uuid, h1_horizon_uuid ),
        "[TEST] Unit 'eroded_unit' should be above horizon "
        "'model_horizon_h1'" );
    OPENGEODE_EXCEPTION(
        su_stack.is_eroded_by( eroded_unit_uuid, erosion_horizon_uuid ),
        "[TEST] Horizon 'model_erosion' should erode unit 'eroded_unit'" );

    geode::save_stratigraphic_units_stack(
        su_stack, "test_su_stack_import.og_sus3d" );
}

int main()
{
    try
    {
        geode::GeosciencesIOModelLibrary::initialize();

        geode::Logger::info( "Reading stratigraphic column file." );
        test_file();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
