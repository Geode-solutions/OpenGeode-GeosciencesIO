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

int main()
{
    using namespace geode;

    try
    {
        detail::initialize_geosciences_io();
        StructuralModel model;

        // Load structural model
        load_structural_model( model, absl::StrCat( data_path, "/modelA4.",
                                          detail::MLInput::extension() ) );

        OPENGEODE_EXCEPTION( model.nb_corners() == 52,
            "[Test] Number of Corners in the loaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_lines() == 98,
            "[Test] Number of Lines in the loaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_surfaces() == 55,
            "[Test] Number of Surfaces in the loaded StructuralModel is not "
            "correct" );
        OPENGEODE_EXCEPTION( model.nb_blocks() == 8,
            "[Test] Number of Blocks in the loaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_faults() == 2,
            "[Test] Number of Faults in the loaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_horizons() == 3,
            "[Test] Number of Horizons in the loaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_model_boundaries() == 6,
            "[Test] Number of ModelBoundary in the loaded StructuralModel is "
            "not correct" );

        index_t nb_block_internals{ 0 };
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

        index_t nb_surface_internals{ 0 };
        for( const auto& surface : model.surfaces() )
        {
            const auto nb_internals = model.nb_internals( surface.id() );
            if( nb_internals )
            {
                for( const auto& collection :
                    model.collections( surface.id() ) )
                {
                    const auto& name =
                        model.model_boundary( collection.id() ).name();
                    OPENGEODE_EXCEPTION( name == "voi_top_boundary"
                                             || name == "voi_bottom_boundary",
                        "[Test] ModelBoundary name is not correct" );
                }
            }
            nb_surface_internals += model.nb_internals( surface.id() );
        }
        OPENGEODE_EXCEPTION( nb_surface_internals == 2,
            "[Test] Number of Surface internals in the loaded StructuralModel "
            "is not correct" );

        // Save structural model
        save_structural_model(
            model, absl::StrCat( "modelA4.", model.native_extension() ) );

        save_structural_model( model, "modelA4_saved.ml" );
        StructuralModel reload;
        load_structural_model( reload, "modelA4_saved.ml" );

        OPENGEODE_EXCEPTION( reload.nb_corners() == 52,
            "[Test] Number of Corners in the reloaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( reload.nb_lines() == 98,
            "[Test] Number of Lines in the reloaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( reload.nb_surfaces() == 55,
            "[Test] Number of Surfaces in the reloaded StructuralModel is not "
            "correct" );
        OPENGEODE_EXCEPTION( reload.nb_blocks() == 8,
            "[Test] Number of Blocks in the reloaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( reload.nb_faults() == 2,
            "[Test] Number of Faults in the reloaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( reload.nb_horizons() == 3,
            "[Test] Number of Horizons in the reloaded "
            "StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( reload.nb_model_boundaries() == 6,
            "[Test] Number of ModelBoundary in the reloaded StructuralModel is "
            "not correct" );

        // StructuralModel model2;
        // load_brep( model2, "/home/anquez/Bureau/brep_mengsu.og_brep" );
        // save_structural_model( model2, "brep_mengsu.ml" );

        Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
