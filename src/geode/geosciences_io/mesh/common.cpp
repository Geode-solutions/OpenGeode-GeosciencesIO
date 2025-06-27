/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

#include <geode/io/image/common.hpp>

#include <geode/geosciences_io/mesh/internal/dem_input.hpp>
#include <geode/geosciences_io/mesh/internal/fem_output.hpp>
#include <geode/geosciences_io/mesh/internal/geotiff_input.hpp>
#include <geode/geosciences_io/mesh/internal/grdecl_input.hpp>
#include <geode/geosciences_io/mesh/internal/pl_input.hpp>
#include <geode/geosciences_io/mesh/internal/pl_output.hpp>
#include <geode/geosciences_io/mesh/internal/ts_input.hpp>
#include <geode/geosciences_io/mesh/internal/ts_output.hpp>
#include <geode/geosciences_io/mesh/internal/vo_input.hpp>
#include <geode/geosciences_io/mesh/internal/vs_input.hpp>
#include <geode/geosciences_io/mesh/internal/vs_output.hpp>
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

    void register_polygonal_surface_input()
    {
        geode::PolygonalSurfaceInputFactory3D::register_creator<
            geode::internal::DEMInput >(
            geode::internal::DEMInput::extension().data() );
    }

    void register_tetrahedral_solid_output()
    {
        geode::TetrahedralSolidOutputFactory3D::register_creator<
            geode::internal::SolidFemOutput >(
            geode::internal::SolidFemOutput::extension().data() );
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

    void register_light_regular_grid_input()
    {
        for( const auto& tif_ext : geode::internal::GEOTIFFInput::extensions() )
        {
            geode::LightRegularGridInputFactory2D::register_creator<
                geode::internal::GEOTIFFInput >( tif_ext );
        }
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

    void register_point_set_input()
    {
        geode::PointSetInputFactory3D::register_creator<
            geode::internal::VSInput >(
            geode::internal::VSInput::extension().data() );
    }

    void register_point_set_output()
    {
        geode::PointSetOutputFactory3D::register_creator<
            geode::internal::VSOutput >(
            geode::internal::VSOutput::extension().data() );
    }
} // namespace

namespace geode
{
    OPENGEODE_LIBRARY_IMPLEMENTATION( GeosciencesIOMesh )
    {
        OpenGeodeMeshLibrary::initialize();
        IOImageLibrary::initialize();
        register_triangulated_surface_input();
        register_triangulated_surface_output();
        register_polygonal_surface_input();
        register_tetrahedral_solid_output();
        register_edged_curve_input();
        register_edged_curve_output();
        register_light_regular_grid_input();
        register_regular_grid_input();
        register_hybrid_solid_input();
        register_point_set_input();
        register_point_set_output();
    }
} // namespace geode
