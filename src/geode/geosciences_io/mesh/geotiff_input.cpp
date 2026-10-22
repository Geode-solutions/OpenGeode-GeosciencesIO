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

#include <geode/geosciences_io/mesh/internal/geotiff_input.hpp>

#include <cmath>
#include <cstdint>
#include <limits>

#include <gdal_priv.h>

#include <absl/strings/str_cat.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/file.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/variable_attribute.hpp>

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
            const auto width = dataset().GetRasterXSize();
            const auto height = dataset().GetRasterYSize();
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                width > 0 && height > 0
                    && static_cast< std::uint64_t >( width )
                               * static_cast< std::uint64_t >( height )
                           < geode::NO_ID,
                nullptr, geode::OpenGeodeException::TYPE::data,
                "[GEOTIFFInput] Raster size ", width, "x", height,
                " exceeds the maximum number of cells" );
            const auto raster = load_raster_quietly();
            auto grid = geode::convert_raster_image_into_grid(
                raster, bottom_up_coordinate_system() );
            for( const auto band_id :
                geode::Range{ 1, dataset().GetRasterCount() + 1 } )
            {
                auto& band = *dataset().GetRasterBand( band_id );
                if( band.GetRasterDataType() == GDT_Byte
                    && band.GetColorInterpretation() != GCI_Undefined )
                {
                    continue;
                }
                read_data_band( grid, band, band_id );
            }
            return grid;
        }

    private:
        /*!
         * Band values that are not colors (e.g. elevations) are stored as a
         * double cell attribute, no data values being NaN
         */
        void read_data_band( geode::LightRegularGrid2D& grid,
            GDALRasterBand& band,
            geode::index_t band_id )
        {
            const auto width = grid.nb_cells_in_direction( 0 );
            const auto height = grid.nb_cells_in_direction( 1 );
            std::vector< double > values(
                static_cast< std::size_t >( width ) * height );
            const auto status = band.RasterIO( GF_Read, 0, 0,
                static_cast< int >( width ), static_cast< int >( height ),
                values.data(), static_cast< int >( width ),
                static_cast< int >( height ), GDT_Float64, 0, 0 );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                status == CE_None, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[GEOTIFFInput] Failed to read band ", band_id );
            int has_no_data{ 0 };
            const auto no_data_value = band.GetNoDataValue( &has_no_data );
            constexpr auto NO_DATA = std::numeric_limits< double >::quiet_NaN();
            auto& manager = grid.cell_attribute_manager();
            std::string name{ band.GetDescription() };
            if( name.empty() || manager.attribute_ids_matching_name( name ) )
            {
                absl::StrAppend(
                    &name, name.empty() ? "" : "_", "band_", band_id );
            }
            geode::AttributeValues< double > default_values;
            default_values.default_value = NO_DATA;
            default_values.no_value = NO_DATA;
            geode::AttributeProperties properties;
            properties.assignable = false;
            properties.interpolable = true;
            properties.transferable = true;
            const auto attribute_id =
                manager.create_attribute< geode::VariableAttribute, double >(
                    name, std::move( default_values ),
                    std::move( properties ) );
            const auto attribute =
                manager.find_attribute< geode::VariableAttribute, double >(
                    attribute_id );
            for( const auto row : geode::Range{ height } )
            {
                for( const auto column : geode::Range{ width } )
                {
                    const auto value = values[row * width + column];
                    if( has_no_data != 0
                        && ( value == no_data_value
                             || ( std::isnan( no_data_value )
                                  && std::isnan( value ) ) ) )
                    {
                        continue;
                    }
                    attribute->set_value(
                        grid.cell_index( { column, height - 1 - row } ),
                        value );
                }
            }
        }

        geode::RasterImage2D load_raster_quietly() const
        {
            class LoggerLevelRestorer
            {
            public:
                LoggerLevelRestorer() : level_{ geode::Logger::level() }
                {
                    geode::Logger::set_level( geode::Logger::LEVEL::critical );
                }
                LoggerLevelRestorer( const LoggerLevelRestorer& ) = delete;
                LoggerLevelRestorer& operator=(
                    const LoggerLevelRestorer& ) = delete;
                LoggerLevelRestorer( LoggerLevelRestorer&& ) = delete;
                LoggerLevelRestorer& operator=(
                    LoggerLevelRestorer&& ) = delete;
                ~LoggerLevelRestorer()
                {
                    geode::Logger::set_level( level_ );
                }

            private:
                geode::Logger::LEVEL level_;
            };
            const LoggerLevelRestorer restorer;
            return geode::load_raster_image< 2 >( filename_ );
        }

        /*!
         * TIFF raster images are read bottom-up (cell 0 is the last image
         * row) while the geotransform origin is the top-left image corner.
         */
        geode::CoordinateSystem2D bottom_up_coordinate_system()
        {
            const auto coordinate_system = read_coordinate_system();
            const auto& row_direction = coordinate_system.direction( 1 );
            const auto height =
                static_cast< double >( dataset().GetRasterYSize() );
            return { { coordinate_system.direction( 0 ), row_direction * -1. },
                coordinate_system.origin() + row_direction * height };
        }

    private:
        std::string_view filename_;
    };
} // namespace

namespace geode::internal
{
    LightRegularGrid2D GEOTIFFInput::read()
    {
        GEOTIFFInputImpl geo_reader( filename() );
        return geo_reader.read_file();
    }

    Percentage GEOTIFFInput::is_loadable() const
    {
        const auto raster_percent = is_raster_image_loadable< 2 >( filename() );
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
} // namespace geode::internal