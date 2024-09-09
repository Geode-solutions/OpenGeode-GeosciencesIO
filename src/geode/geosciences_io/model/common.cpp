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

#include <geode/geosciences_io/model/common.hpp>

#include <gdal_priv.h>

#include <geode/mesh/io/regular_grid_input.hpp>
#include <geode/mesh/io/triangulated_surface_input.hpp>

#include <geode/geosciences/explicit/common.hpp>
#include <geode/geosciences/explicit/representation/io/structural_model_input.hpp>
#include <geode/geosciences/explicit/representation/io/structural_model_output.hpp>
#include <geode/geosciences/implicit/representation/io/horizons_stack_input.hpp>

#include <geode/geosciences_io/model/internal/brep_fem_output.hpp>
#include <geode/geosciences_io/model/internal/horizons_stack_skua_input.hpp>
#include <geode/geosciences_io/model/internal/lso_input.hpp>
#include <geode/geosciences_io/model/internal/lso_output.hpp>
#include <geode/geosciences_io/model/internal/ml_input.hpp>
#include <geode/geosciences_io/model/internal/ml_output_brep.hpp>
#include <geode/geosciences_io/model/internal/ml_output_structural_model.hpp>
#include <geode/geosciences_io/model/internal/shp_input.hpp>

namespace
{
    void register_structural_model_input()
    {
        geode::StructuralModelInputFactory::register_creator<
            geode::internal::MLInput >(
            geode::internal::MLInput::extension().data() );
        geode::StructuralModelInputFactory::register_creator<
            geode::internal::LSOInput >(
            geode::internal::LSOInput::extension().data() );
    }

    void register_structural_model_output()
    {
        geode::StructuralModelOutputFactory::register_creator<
            geode::internal::MLOutputStructuralModel >(
            geode::internal::MLOutputStructuralModel::extension().data() );
        geode::StructuralModelOutputFactory::register_creator<
            geode::internal::LSOOutput >(
            geode::internal::LSOOutput::extension().data() );
    }

    void register_section_input()
    {
        for( const auto& shp_ext : geode::internal::SHPInput::extensions() )
        {
            geode::SectionInputFactory::register_creator<
                geode::internal::SHPInput >( shp_ext );
        }
    }

    void register_brep_output()
    {
        geode::BRepOutputFactory::register_creator<
            geode::internal::MLOutputBRep >(
            geode::internal::MLOutputBRep::extension().data() );
    }

    void register_brep_fem_output()
    {
        geode::BRepOutputFactory::register_creator<
            geode::internal::BRepFemOutput >(
            geode::internal::BRepFemOutput::extension().data() );
    }

    void register_horizons_stack_input()
    {
        geode::HorizonsStackInputFactory< 2 >::register_creator<
            geode::internal::HorizonStackSKUAInput< 2 > >(
            geode::internal::HorizonStackSKUAInput< 2 >::extension().data() );
        geode::HorizonsStackInputFactory< 3 >::register_creator<
            geode::internal::HorizonStackSKUAInput< 3 > >(
            geode::internal::HorizonStackSKUAInput< 3 >::extension().data() );
    }
} // namespace

namespace geode
{
    OPENGEODE_LIBRARY_IMPLEMENTATION( GeosciencesIOModel )
    {
        GeosciencesExplicitLibrary::initialize();
        GeosciencesImplicitLibrary::initialize();
        register_structural_model_input();
        register_structural_model_output();
        register_section_input();
        register_brep_output();
        register_brep_fem_output();
        register_horizons_stack_input();
        GDALAllRegister();
    }
} // namespace geode
