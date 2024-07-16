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

#include <geode/geosciences_io/mesh/common.hpp>

#include <geode/mesh/io/regular_grid_input.hpp>
#include <geode/mesh/io/triangulated_surface_input.hpp>

#include <geode/geosciences_io/mesh/internal/grdecl_input.hpp>
#include <geode/geosciences_io/mesh/internal/pl_input.hpp>
#include <geode/geosciences_io/mesh/internal/pl_output.hpp>
#include <geode/geosciences_io/mesh/internal/ts_input.hpp>
#include <geode/geosciences_io/mesh/internal/ts_output.hpp>
#include <geode/geosciences_io/mesh/internal/vo_input.hpp>
#include <geode/geosciences_io/mesh/internal/well_dat_input.hpp>
#include <geode/geosciences_io/mesh/internal/well_dev_input.hpp>
#include <geode/geosciences_io/mesh/internal/well_txt_input.hpp>
#include <geode/geosciences_io/mesh/internal/wl_input.hpp>

namespace
{
    void register_triangulated_surface_input()
    {
        geode::TriangulatedSurfaceInputFactory3D::register_creator<
            geode::internal::TSInput >(
            geode::internal::TSInput::extension().data() );
    }

    void register_triangulated_surface_output()
    {
        geode::TriangulatedSurfaceOutputFactory3D::register_creator<
            geode::internal::TSOutput >(
            geode::internal::TSOutput::extension().data() );
    }

    void register_edged_curve_output()
    {
        geode::EdgedCurveOutputFactory3D::register_creator<
            geode::internal::PLOutput >(
            geode::internal::PLOutput::extension().data() );
    }

    void register_edged_curve_input()
    {
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::internal::WLInput >(
            geode::internal::WLInput::extension().data() );
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::internal::WellDatInput >(
            geode::internal::WellDatInput::extension().data() );
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::internal::WellTxtInput >(
            geode::internal::WellTxtInput::extension().data() );
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::internal::WellDevInput >(
            geode::internal::WellDevInput::extension().data() );
        geode::EdgedCurveInputFactory3D::register_creator<
            geode::internal::PLInput >(
            geode::internal::PLInput::extension().data() );
    }

    void register_regular_grid_input()
    {
        geode::RegularGridInputFactory3D::register_creator<
            geode::internal::VOInput >(
            geode::internal::VOInput::extension().data() );
    }

    void register_hybrid_solid_input()
    {
        geode::HybridSolidInputFactory3D::register_creator<
            geode::internal::GRDECLInput >(
            geode::internal::GRDECLInput::extension().data() );
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
