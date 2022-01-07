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

#include <geode/tests_config.h>

#include <geode/basic/assert.h>
#include <geode/basic/logger.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>

#include <geode/geosciences/private/lso_input.h>
#include <geode/geosciences/representation/core/structural_model.h>
#include <geode/geosciences/representation/io/structural_model_output.h>

void check_model( const geode::StructuralModel& model )
{
    OPENGEODE_EXCEPTION( model.nb_corners() == 22,
        "[Test] Number of Corners in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_lines() == 39,
        "[Test] Number of Lines in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_surfaces() == 23,
        "[Test] Number of Surfaces in the loaded StructuralModel is not "
        "correct" );
    OPENGEODE_EXCEPTION( model.nb_blocks() == 4,
        "[Test] Number of Blocks in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_horizons() == 4,
        "[Test] Number of Horizons in the loaded "
        "StructuralModel is not correct" );

    geode::index_t nb_block_internals{ 0 };
    for( const auto& block : model.blocks() )
    {
        nb_block_internals += model.nb_internals( block.id() );
    }
    OPENGEODE_EXCEPTION( nb_block_internals == 2,
        "[Test] Number of Block internals in the "
        "loaded StructuralModel is not correct" );
}

int main()
{
    try
    {
        geode::detail::initialize_geosciences_io();

        // Load structural model
        auto model = geode::load_structural_model( absl::StrCat(
            geode::data_path, "test.", geode::detail::LSOInput::extension() ) );
        check_model( model );

        geode::save_structural_model( model, "test.lso" );
        auto reload_model = geode::load_structural_model( "test.lso" );
        check_model( reload_model );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
