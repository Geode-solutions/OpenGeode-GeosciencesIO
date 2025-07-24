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

#include <geode/geosciences_io/mesh/internal/gdal_utils.hpp>

#include <gdal_priv.h>

#include <geode/geometry/coordinate_system.hpp>

namespace geode
{
    namespace internal
    {
        CoordinateSystem2D read_coordinate_system( GDALDataset& dataset )
        {
            std::array< double, 6 > geo_transform;
            const auto status = dataset.GetGeoTransform( geo_transform.data() );
            OPENGEODE_EXCEPTION( status == CE_None,
                "Failed to read geotransform from GDALDataset" );
            Point2D origin;
            Vector2D x_direction;
            Vector2D y_direction;
            for( const auto axis : geode::LRange{ 2 } )
            {
                origin.set_value( axis, geo_transform[3 * axis] );
                x_direction.set_value( axis, geo_transform[3 * axis + 1] );
                y_direction.set_value( axis, geo_transform[3 * axis + 2] );
            }
            return { origin, { x_direction, y_direction } };
        }
    } // namespace internal
} // namespace geode