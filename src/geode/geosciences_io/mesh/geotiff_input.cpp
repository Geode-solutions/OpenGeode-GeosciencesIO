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

#include <geode/geosciences_io/mesh/internal/geotiff_input.hpp>

#include <ogrsf_frmts.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/file.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>

#include <geode/mesh/core/light_regular_grid.hpp>

#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/image/core/raster_image.hpp>
#include <geode/image/core/rgb_color.hpp>
#include <geode/image/io/raster_image_output.hpp>

#include <geode/io/image/internal/tiff_input.hpp>

namespace
{
    class GEOTIFFInputImpl
    {
    public:
        GEOTIFFInputImpl( std::string_view filename )
            : gdal_data_{ GDALDataset::Open(
                geode::to_string( filename ).c_str(), GDAL_OF_READONLY ) }
        {
            OPENGEODE_EXCEPTION(
                gdal_data_, "[GEOTIFFInput] Failed to open file ", filename );
            get_geotransform();
        }

        const geode::Point2D& get_origin() const
        {
            return origin_;
        }

        std::array< geode::Vector2D, 2 > get_directions() const
        {
            return { x_direction_, y_direction_ };
        }

        bool get_geotransform()
        {
            double geoTransform[6];
            if( gdal_data_->GetGeoTransform( geoTransform ) == CE_Failure )
            {
                return false;
            }
            origin_ = geode::Point2D{ { geoTransform[0], geoTransform[3] } };
            x_direction_ =
                geode::Vector2D{ { geoTransform[1], geoTransform[2] } };
            y_direction_ =
                geode::Vector2D{ { geoTransform[4], geoTransform[5] } };
            return true;
        }

        bool get_projection_ref()
        {
            const auto* projection_ref = gdal_data_->GetProjectionRef();
            if( projection_ref == nullptr || strlen( projection_ref ) == 0 )
            {
                return false;
            }
            // apparently there is everthing in this string to load the CRS...
            return true;
        }

    private:
        geode::Point2D origin_{ { 0., 0. } };
        geode::Vector2D x_direction_{ { 1., 0. } };
        geode::Vector2D y_direction_{ { 0., 1. } };

        GDALDatasetUniquePtr gdal_data_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        LightRegularGrid2D GEOTIFFInput::read()
        {
            TIFFInput image_reader( filename() );
            auto image = image_reader.read();
            std::array< index_t, 2 > cells_number{ image.nb_cells_in_direction(
                                                       0 ),
                image.nb_cells_in_direction( 1 ) };
            GEOTIFFInputImpl geo_reader( filename() );
            const auto& origin = geo_reader.get_origin();

            auto cells_directions = geo_reader.get_directions();

            LightRegularGrid2D grid{ origin, cells_number, cells_directions };
            auto color =
                grid.cell_attribute_manager()
                    .find_or_create_attribute< VariableAttribute, RGBColor >(
                        "RGB_data", {} );
            for( const auto cell_id : geode::Range{ image.nb_cells() } )
            {
                color->set_value( cell_id, image.color( cell_id ) );
            }
            return grid;
        }

        bool GEOTIFFInput::is_loadable() const
        {
            GDALDatasetUniquePtr gdal_data{ GDALDataset::Open(
                geode::to_string( filename() ).c_str(), GDAL_OF_VECTOR ) };
            if( !gdal_data )
            {
                return false;
            }
            double geoTransform[6];
            if( gdal_data->GetGeoTransform( geoTransform ) == CE_Failure )
            {
                return false;
            }
            return true;
        }

    } // namespace internal
} // namespace geode