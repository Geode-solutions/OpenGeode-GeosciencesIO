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

#include <geode/tests_config.h>

#include <geode/basic/assert.h>
#include <geode/basic/filename.h>
#include <geode/basic/logger.h>

#include <geode/geosciences/explicit/mixin/core/horizon.h>
#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.h>
#include <geode/geosciences/implicit/representation/core/horizons_stack.h>
#include <geode/geosciences/implicit/representation/io/horizons_stack_output.h>
#include <geode/geosciences_io/model/internal/horizons_stack_skua_input.h>

void test_file()
{
    const auto horizons_stack = geode::load_horizons_stack< 3 >(
        absl::StrCat( geode::DATA_PATH, "test_skua_horizons_stack.xml" ) );

    OPENGEODE_EXCEPTION( horizons_stack.name() == "skua_model",
        "[TEST] HorizonsStack should be named 'skua_model'" );

    OPENGEODE_EXCEPTION( horizons_stack.nb_horizons() == 4,
        "[TEST] Wrong number of horizons in the "
        "loaded HorizonsStack." );
    OPENGEODE_EXCEPTION( horizons_stack.nb_stratigraphic_units() == 5,
        "[TEST] Wrong number of units in the loaded HorizonsStack." );
    geode::uuid erosion_horizon_uuid;
    geode::uuid h1_horizon_uuid;
    for( const auto& horizon : horizons_stack.horizons() )
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
    for( const auto& unit : horizons_stack.stratigraphic_units() )
    {
        if( unit.name() == "eroded_unit" )
        {
            eroded_unit_uuid = unit.id();
        }
    }
    OPENGEODE_EXCEPTION(
        horizons_stack.is_above( erosion_horizon_uuid, eroded_unit_uuid ),
        "[TEST] Horizon 'model_erosion' should be above unit 'eroded_unit'" );
    OPENGEODE_EXCEPTION(
        horizons_stack.is_above( eroded_unit_uuid, h1_horizon_uuid ),
        "[TEST] Unit 'eroded_unit' should be above horizon "
        "'model_horizon_h1'" );
    OPENGEODE_EXCEPTION(
        horizons_stack.is_eroded_by(
            horizons_stack.stratigraphic_unit( eroded_unit_uuid ),
            horizons_stack.horizon( erosion_horizon_uuid ) ),
        "[TEST] Horizon 'model_erosion' should erode unit 'eroded_unit'" );

    geode::save_horizons_stack(
        horizons_stack, "test_horizons_stack_import.og_hst3d" );
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
