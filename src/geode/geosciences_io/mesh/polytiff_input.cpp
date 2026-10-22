/*
 * Copyright (c) 2019 - 2026 Geode-solutions
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

#include <geode/geosciences_io/mesh/internal/polytiff_input.hpp>

#include <geode/mesh/core/polygonal_surface.hpp>

#include <gdal_priv.h>

#include <geode/io/image/detail/gdal_file.hpp>

#include <geode/geosciences_io/mesh/internal/dem_input.hpp>

namespace geode::internal
{
    std::unique_ptr< PolygonalSurface3D > PolyTIFFInput::read(
        const MeshImpl& impl )
    {
        auto surface = PolygonalSurface3D::create( impl );
        read_elevation_raster( filename(), *surface );
        return surface;
    }

    auto PolyTIFFInput::additional_files() const -> AdditionalFiles
    {
        detail::GDALFile reader{ filename() };
        return reader.additional_files< AdditionalFiles >();
    }

    Percentage PolyTIFFInput::is_loadable() const
    {
        detail::GDALFile reader{ this->filename() };
        if( !reader.is_coordinate_system_loadable()
            || reader.dataset().GetRasterCount() == 0 )
        {
            return Percentage{ 0 };
        }
        return Percentage{ 1 };
    }
} // namespace geode::internal