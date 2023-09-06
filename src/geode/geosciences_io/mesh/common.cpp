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

#include <geode/geosciences_io/mesh/common.h>

#include <geode/mesh/io/regular_grid_input.h>
#include <geode/mesh/io/triangulated_surface_input.h>

#include <geode/geosciences_io/mesh/private/pl_output.h>
#include <geode/geosciences_io/mesh/private/ts_input.h>
#include <geode/geosciences_io/mesh/private/ts_output.h>
#include <geode/geosciences_io/mesh/private/vo_input.h>
#include <geode/geosciences_io/mesh/private/well_dat_input.h>
#include <geode/geosciences_io/mesh/private/well_dev_input.h>
#include <geode/geosciences_io/mesh/private/well_txt_input.h>
#include <geode/geosciences_io/mesh/private/wl_input.h>
#include <geode/geosciences_io/mesh/private/grdecl_input.h>

namespace
{
    void register_triangulated_surface_input()
    {
        geode::TriangulatedSurfaceInputFactory3D::register_creator<
            geode::detail::TSInput >(
            geode::detail::TSInput::extension().data() );
    }

    void register_triangulated_surface_output()
    {
        geode::TriangulatedSurfaceOutputFactory3D::register_creator<
            geode::detail::TSOutput >(
            geode::detail::TSOutput::extension().data() );
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
   
    void register_hybrid_solid_input()
    {
        geode::HybridSolidInputFactory3D::register_creator<
        geode::detail::GRDECLInput >(
        geode::detail::GRDECLInput::extension().data() ); 
    }

} // namespace

namespace geode
{
    OPENGEODE_LIBRARY_IMPLEMENTATION( GeosciencesIOMesh )
    {
        OpenGeodeMeshLibrary::initialize();
        register_triangulated_surface_input();
        register_triangulated_surface_output();
        register_edged_curve_input();
        register_edged_curve_output();
        register_regular_grid_input();
        register_hybrid_solid_input();
    }
} // namespace geode
