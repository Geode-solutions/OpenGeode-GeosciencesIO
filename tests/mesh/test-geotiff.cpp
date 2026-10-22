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

#include <geode/tests_config.hpp>

#include <cmath>
#include <fstream>
#include <limits>
#include <vector>

#include <geode/basic/assert.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>

#include <geode/image/core/rgb_color.hpp>

#include <geode/mesh/core/light_regular_grid.hpp>

#include <geode/mesh/io/light_regular_grid_input.hpp>
#include <geode/mesh/io/light_regular_grid_output.hpp>

#include <geode/geosciences_io/mesh/common.hpp>
#include <geode/io/mesh/common.hpp>

std::string write_raster( std::string_view name,
    std::string_view extension,
    geode::index_t width,
    geode::index_t height,
    const std::vector< float >& values,
    std::string_view no_data )
{
    auto filename = absl::StrCat( name, ".", extension );
    std::ofstream vrt{ filename };
    vrt << "<VRTDataset rasterXSize=\"" << width << "\" rasterYSize=\""
        << height << "\">\n"
        << "  <GeoTransform>100, 10, 0, 500, 0, -10</GeoTransform>\n";
    if( values.empty() )
    {
        vrt << "  <VRTRasterBand dataType=\"Float32\" band=\"1\"/>\n"
            << "</VRTDataset>\n";
        return filename;
    }
    const auto raw_filename = absl::StrCat( name, ".raw" );
    {
        std::ofstream raw{ raw_filename, std::ios::binary };
        raw.write( reinterpret_cast< const char* >( values.data() ),
            static_cast< std::streamsize >( values.size() * sizeof( float ) ) );
    }
    vrt << "  <VRTRasterBand dataType=\"Float32\" band=\"1\" "
           "subClass=\"VRTRawRasterBand\">\n";
    if( !no_data.empty() )
    {
        vrt << "    <NoDataValue>" << no_data << "</NoDataValue>\n";
    }
    vrt << "    <SourceFilename relativeToVRT=\"1\">" << raw_filename
        << "</SourceFilename>\n"
        << "    <ImageOffset>0</ImageOffset>\n"
        << "    <PixelOffset>4</PixelOffset>\n"
        << "    <LineOffset>" << 4 * width << "</LineOffset>\n"
        << "    <ByteOrder>LSB</ByteOrder>\n"
        << "  </VRTRasterBand>\n"
        << "</VRTDataset>\n";
    return filename;
}

template < typename T >
std::shared_ptr< geode::ReadOnlyAttribute< T > > find_named_attribute(
    const geode::AttributeManager& manager, std::string_view name )
{
    const auto ids = manager.attribute_ids_matching_name( name );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        ids && ids->size() == 1, "Cannot find attribute ", name );
    return manager.find_read_only_attribute< T >( ids->front() );
}

geode::Grid2D::CellIndices find_pixel_cell(
    const geode::LightRegularGrid2D& grid,
    geode::index_t row,
    geode::index_t column )
{
    const geode::Point2D center{ { 105. + 10. * column, 495. - 10. * row } };
    for( const auto j : geode::Range{ grid.nb_cells_in_direction( 1 ) } )
    {
        for( const auto i : geode::Range{ grid.nb_cells_in_direction( 0 ) } )
        {
            if( grid.cell_barycenter( { i, j } ).inexact_equal( center ) )
            {
                return { i, j };
            }
        }
    }
    geode::OpenGeodeGeosciencesIOMeshException::test(
        false, "Cannot find cell of pixel ", row, " ", column );
    return {};
}

void test_pixel_orientation()
{
    const std::vector< float > values{ 1, 2, 3, 4, 5, 6 };
    const auto filename =
        write_raster( "raster_geotiff", "tif", 3, 2, values, "" );
    const auto grid = geode::load_light_regular_grid< 2 >( filename );
    const auto colors = find_named_attribute< geode::RGBColor >(
        grid.cell_attribute_manager(), "RGB_data" );
    for( const auto row : geode::Range{ 2 } )
    {
        for( const auto column : geode::Range{ 3 } )
        {
            const auto cell =
                grid.cell_index( find_pixel_cell( grid, row, column ) );
            geode::OpenGeodeGeosciencesIOMeshException::test(
                colors->value( cell ).red() == values[row * 3 + column],
                "Wrong color for pixel ", row, " ", column );
        }
    }
}

void test_logger_level_on_error()
{
    const auto filename = "raster_geotiff_missing_source.tif";
    {
        std::ofstream vrt{ filename };
        vrt << "<VRTDataset rasterXSize=\"3\" rasterYSize=\"2\">\n"
            << "  <GeoTransform>100, 10, 0, 500, 0, -10</GeoTransform>\n"
            << "  <VRTRasterBand dataType=\"Byte\" band=\"1\">\n"
            << "    <SimpleSource>\n"
            << "      <SourceFilename relativeToVRT=\"1\">"
               "raster_geotiff_does_not_exist.tif</SourceFilename>\n"
            << "      <SourceBand>1</SourceBand>\n"
            << "      <SourceProperties RasterXSize=\"3\" RasterYSize=\"2\" "
               "DataType=\"Byte\" BlockXSize=\"3\" BlockYSize=\"2\"/>\n"
            << "    </SimpleSource>\n"
            << "  </VRTRasterBand>\n"
            << "</VRTDataset>\n";
    }
    const auto level = geode::Logger::level();
    bool thrown{ false };
    try
    {
        const auto grid = geode::load_light_regular_grid< 2 >( filename );
    }
    catch( const geode::OpenGeodeException& )
    {
        thrown = true;
    }
    geode::OpenGeodeGeosciencesIOMeshException::test(
        thrown, "Loading a raster with a missing source should throw" );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        geode::Logger::level() == level,
        "Logger level should be restored after a loading error" );
}

void test_data_band()
{
    const auto nan = std::numeric_limits< float >::quiet_NaN();
    const std::vector< float > values{ -5.5F, 2, 3, 400, 5, nan };
    const auto filename =
        write_raster( "raster_geotiff_data", "tif", 3, 2, values, "nan" );
    const auto grid = geode::load_light_regular_grid< 2 >( filename );
    const auto band = find_named_attribute< double >(
        grid.cell_attribute_manager(), "band_1" );
    for( const auto row : geode::Range{ 2 } )
    {
        for( const auto column : geode::Range{ 3 } )
        {
            const auto cell =
                grid.cell_index( find_pixel_cell( grid, row, column ) );
            const auto expected = values[row * 3 + column];
            if( std::isnan( expected ) )
            {
                geode::OpenGeodeGeosciencesIOMeshException::test(
                    !band->has_value( cell ), "No data pixel ", row, " ",
                    column, " should have no value" );
                continue;
            }
            geode::OpenGeodeGeosciencesIOMeshException::test(
                band->value( cell ) == expected, "Wrong band value for pixel ",
                row, " ", column );
        }
    }
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        geode::OpenGeodeIOMeshLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::trace );

        auto grid = geode::load_light_regular_grid< 2 >(
            absl::StrCat( geode::DATA_PATH, "cea.tiff" ) );
        geode::save_light_regular_grid( grid, "cea.vti" );

        test_pixel_orientation();
        test_logger_level_on_error();
        test_data_band();

        geode::Logger::info( "[TEST SUCCESS]" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
