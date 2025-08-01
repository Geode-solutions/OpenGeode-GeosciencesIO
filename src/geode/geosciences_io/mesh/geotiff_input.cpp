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

#include <gdal_priv.h>

#include <geode/basic/file.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/coordinate_system.hpp>
#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/image/core/raster_image.hpp>
#include <geode/image/io/raster_image_input.hpp>

#include <geode/mesh/core/light_regular_grid.hpp>
#include <geode/mesh/helpers/convert_grid.hpp>

#include <geode/io/image/detail/gdal_file.hpp>

namespace
{
    class GEOTIFFInputImpl : public geode::detail::GDALFile
    {
    public:
        GEOTIFFInputImpl( std::string_view filename )
            : geode::detail::GDALFile{ filename }, filename_{ filename }
        {
        }

        geode::LightRegularGrid2D read_file()
        {
            const auto level = geode::Logger::level();
            geode::Logger::set_level( geode::Logger::LEVEL::critical );
            const auto raster = geode::load_raster_image< 2 >( filename_ );
            geode::Logger::set_level( level );
            const auto cooridinate_system = read_coordinate_system();
            return geode::convert_raster_image_into_grid(
                raster, cooridinate_system );
        }

    private:
        std::string_view filename_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        LightRegularGrid2D GEOTIFFInput::read()
        {
            GEOTIFFInputImpl geo_reader( filename() );
            return geo_reader.read_file();
        }

        Percentage GEOTIFFInput::is_loadable() const
        {
            const auto raster_percent =
                is_raster_image_loadable< 2 >( filename() );
            if( raster_percent.value() != 1 )
            {
                return raster_percent;
            }
            detail::GDALFile reader{ this->filename() };
            if( !reader.is_coordinate_system_loadable() )
            {
                return Percentage{ 0 };
            }
            return Percentage{ 1 };
        }

        auto GEOTIFFInput::additional_files() const -> AdditionalFiles
        {
            detail::GDALFile reader{ this->filename() };
            return reader.additional_files< AdditionalFiles >();
        }
    } // namespace internal
} // namespace geode