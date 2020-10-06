/*
 * Copyright (c) 2019 - 2020 Geode-solutions
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

#include <geode/geosciences/private/ml_input.h>
#include <geode/geosciences/private/ml_output.h>
#include <geode/geosciences/representation/core/structural_model.h>
#include <geode/geosciences/representation/io/structural_model_output.h>

void check_model( const geode::StructuralModel& model,
    geode::index_t nb_corners,
    geode::index_t nb_lines,
    geode::index_t nb_surfaces,
    geode::index_t nb_blocks,
    geode::index_t nb_faults,
    geode::index_t nb_horizons,
    geode::index_t nb_model_boundaries )
{
    OPENGEODE_EXCEPTION( model.nb_corners() == nb_corners,
        "[Test] Number of Corners in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_lines() == nb_lines,
        "[Test] Number of Lines in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_surfaces() == nb_surfaces,
        "[Test] Number of Surfaces in the loaded StructuralModel is not "
        "correct" );
    OPENGEODE_EXCEPTION( model.nb_blocks() == nb_blocks,
        "[Test] Number of Blocks in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_faults() == nb_faults,
        "[Test] Number of Faults in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_horizons() == nb_horizons,
        "[Test] Number of Horizons in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_model_boundaries() == nb_model_boundaries,
        "[Test] Number of ModelBoundary in the loaded StructuralModel is "
        "not correct" );
}

void test_modelA4()
{
    // Load structural model
    auto model = geode::load_structural_model( absl::StrCat(
        geode::data_path, "/modelA4.", geode::detail::MLInput::extension() ) );
    check_model( model, 52, 98, 55, 8, 2, 3, 6 );

    geode::index_t nb_block_internals{ 0 };
    for( const auto& block : model.blocks() )
    {
        const auto nb_internals = model.nb_internals( block.id() );
        if( nb_internals )
        {
            auto token = block.name().substr( block.name().size() - 3 );
            OPENGEODE_EXCEPTION(
                token == "b_2", "[Test] Block name should end by b_2" );
        }
        nb_block_internals += nb_internals;
    }
    OPENGEODE_EXCEPTION( nb_block_internals == 4,
        "[Test] Number of Block internals in the "
        "loaded StructuralModel is not correct" );

    geode::index_t nb_surface_internals{ 0 };
    for( const auto& surface : model.surfaces() )
    {
        const auto nb_internals = model.nb_internals( surface.id() );
        if( nb_internals )
        {
            for( const auto& collection : model.collections( surface.id() ) )
            {
                const auto& name =
                    model.model_boundary( collection.id() ).name();
                OPENGEODE_EXCEPTION(
                    name == "voi_top_boundary" || name == "voi_bottom_boundary",
                    "[Test] ModelBoundary name is not correct" );
            }
        }
        nb_surface_internals += model.nb_internals( surface.id() );
    }
    OPENGEODE_EXCEPTION( nb_surface_internals == 2,
        "[Test] Number of Surface internals in the loaded StructuralModel "
        "is not correct" );

    // Save structural model
    geode::save_structural_model(
        model, absl::StrCat( "modelA4.", model.native_extension() ) );

    geode::save_structural_model( model, "modelA4_saved.ml" );
    auto reload = geode::load_structural_model( "modelA4_saved.ml" );
    check_model( reload, 52, 98, 55, 8, 2, 3, 6 );
}

int main()
{
    using namespace geode;

    try
    {
        detail::initialize_geosciences_io();
        test_modelA4();

        Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
