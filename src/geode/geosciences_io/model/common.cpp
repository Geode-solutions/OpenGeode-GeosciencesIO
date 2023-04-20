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

#include <geode/geosciences_io/model/common.h>

#include <gdal_priv.h>

#include <geode/mesh/io/regular_grid_input.h>
#include <geode/mesh/io/triangulated_surface_input.h>

#include <geode/geosciences/explicit/common.h>
#include <geode/geosciences/explicit/representation/io/structural_model_input.h>
#include <geode/geosciences/explicit/representation/io/structural_model_output.h>
#include <geode/geosciences/implicit/representation/io/stratigraphic_units_stack_input.h>

#include <geode/geosciences_io/model/private/lso_input.h>
#include <geode/geosciences_io/model/private/lso_output.h>
#include <geode/geosciences_io/model/private/ml_input.h>
#include <geode/geosciences_io/model/private/ml_output_brep.h>
#include <geode/geosciences_io/model/private/ml_output_structural_model.h>
#include <geode/geosciences_io/model/private/shp_input.h>
#include <geode/geosciences_io/model/private/su_stack_skua_input.h>

namespace
{
    void register_structural_model_input()
    {
        geode::StructuralModelInputFactory::register_creator<
            geode::detail::MLInput >(
            geode::detail::MLInput::extension().data() );
        geode::StructuralModelInputFactory::register_creator<
            geode::detail::LSOInput >(
            geode::detail::LSOInput::extension().data() );
    }

    void register_structural_model_output()
    {
        geode::StructuralModelOutputFactory::register_creator<
            geode::detail::MLOutputStructuralModel >(
            geode::detail::MLOutputStructuralModel::extension().data() );
        geode::StructuralModelOutputFactory::register_creator<
            geode::detail::LSOOutput >(
            geode::detail::LSOOutput::extension().data() );
    }

    void register_section_input()
    {
        for( const auto& shp_ext : geode::detail::SHPInput::extensions() )
        {
            geode::SectionInputFactory::register_creator<
                geode::detail::SHPInput >( shp_ext );
        }
    }

    void register_brep_output()
    {
        geode::BRepOutputFactory::register_creator<
            geode::detail::MLOutputBRep >(
            geode::detail::MLOutputBRep::extension().data() );
    }

    void register_stratigraphic_units_stack_input()
    {
        geode::StratigraphicUnitsStackInputFactory< 2 >::register_creator<
            geode::detail::SUStackSKUAInput< 2 > >(
            geode::detail::SUStackSKUAInput< 2 >::extension().data() );
        geode::StratigraphicUnitsStackInputFactory< 3 >::register_creator<
            geode::detail::SUStackSKUAInput< 3 > >(
            geode::detail::SUStackSKUAInput< 3 >::extension().data() );
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
        register_stratigraphic_units_stack_input();
        GDALAllRegister();
    }
} // namespace geode
