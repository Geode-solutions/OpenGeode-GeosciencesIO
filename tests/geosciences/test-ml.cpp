/*
 * Copyright (c) 2019 Geode-solutions
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

#include <geode/geosciences/detail/ml_input.h>
#include <geode/geosciences/representation/core/structural_model.h>
#include <geode/geosciences/representation/io/structural_model_output.h>

int main()
{
    using namespace geode;

    try
    {
        initialize_geosciences_io();
        StructuralModel model;

        // Load structural model
        load_structural_model( model,
            test_path + "geosciences/data/modelA2." + MLInput::extension() );

        OPENGEODE_EXCEPTION( model.nb_corners() == 52,
            "Number of Corners in the loaded StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_lines() == 100,
            "Number of Lines in the loaded StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_surfaces() == 61,
            "Number of Surfaces in the loaded StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_blocks() == 12,
            "Number of Blocks in the loaded StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_faults() == 2,
            "Number of Faults in the loaded StructuralModel is not correct" );
        OPENGEODE_EXCEPTION( model.nb_horizons() == 3,
            "Number of Horizons in the loaded StructuralModel is not correct" );

        // Save structural model
        std::string output_file_native{ "modelA2." + model.native_extension() };
        save_structural_model( model, output_file_native );

        Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
