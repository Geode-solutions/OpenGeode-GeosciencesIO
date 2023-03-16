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

#include <geode/io/geosciences/common.h>

#include <gdal_priv.h>

#include <geode/mesh/io/regular_grid_input.h>
#include <geode/mesh/io/triangulated_surface_input.h>

#include <geode/geosciences/common.h>
#include <geode/geosciences/representation/io/structural_model_input.h>
#include <geode/geosciences/representation/io/structural_model_output.h>

#include <geode/io/geosciences/private/lso_input.h>
#include <geode/io/geosciences/private/lso_output.h>
#include <geode/io/geosciences/private/ml_input.h>
#include <geode/io/geosciences/private/ml_output_brep.h>
#include <geode/io/geosciences/private/ml_output_structural_model.h>
#include <geode/io/geosciences/private/pl_output.h>
#include <geode/io/geosciences/private/shp_input.h>
#include <geode/io/geosciences/private/ts_input.h>
#include <geode/io/geosciences/private/vo_input.h>
#include <geode/io/geosciences/private/well_dat_input.h>
#include <geode/io/geosciences/private/well_dev_input.h>
#include <geode/io/geosciences/private/well_txt_input.h>
#include <geode/io/geosciences/private/wl_input.h>

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

    void register_triangulated_surface_input()
    {
        geode::TriangulatedSurfaceInputFactory3D::register_creator<
            geode::detail::TSInput >(
            geode::detail::TSInput::extension().data() );
    }

    void register_edged_curve_output()
    {
        geode::EdgedCurveOutputFactory3D::register_creator<
            geode::detail::PLOutput >(
            geode::detail::PLOutput::extension().data() );
    }

    void register_edged_curve_input()
    {
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::detail::WLInput >(
            geode::detail::WLInput::extension().data() );
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::detail::WellDatInput >(
            geode::detail::WellDatInput::extension().data() );
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::detail::WellTxtInput >(
            geode::detail::WellTxtInput::extension().data() );
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::detail::WellDevInput >(
            geode::detail::WellDevInput::extension().data() );
    }

    void register_regular_grid_input()
    {
        geode::RegularGridInputFactory3D::register_creator<
            geode::detail::VOInput >(
            geode::detail::VOInput::extension().data() );
    }
} // namespace

namespace geode
{
    OPENGEODE_LIBRARY_IMPLEMENTATION( OpenGeodeGeosciencesIOGeosciences )
    {
        OpenGeodeGeosciencesGeosciences::initialize();
        register_structural_model_input();
        register_structural_model_output();
        register_section_input();
        register_brep_output();
        register_triangulated_surface_input();
        register_edged_curve_input();
        register_edged_curve_output();
        register_regular_grid_input();
        GDALAllRegister();
    }
} // namespace geode
