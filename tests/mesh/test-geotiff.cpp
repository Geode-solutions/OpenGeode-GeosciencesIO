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

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <optional>
#include <vector>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include <geode/basic/assert.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>

#include <geode/image/core/rgb_color.hpp>

#include <geode/mesh/core/grid_coordinate_reference_system_manager.hpp>
#include <geode/mesh/core/light_regular_grid.hpp>

#include <geode/mesh/io/light_regular_grid_input.hpp>
#include <geode/mesh/io/light_regular_grid_output.hpp>

#include <geode/geosciences/explicit/geometry/geographic_coordinate_system.hpp>

#include <geode/geosciences_io/mesh/common.hpp>
#include <geode/io/mesh/common.hpp>

std::string write_raster( std::string_view name,
    geode::index_t width,
    geode::index_t height,
    const std::vector< float >& values,
    std::optional< double > no_data,
    std::optional< int > epsg_code = std::nullopt )
{
    GDALAllRegister();
    auto filename = absl::StrCat( name, ".tif" );
    auto* driver = GetGDALDriverManager()->GetDriverByName( "GTiff" );
    const std::array< const char*, 2 > options{ "SPARSE_OK=TRUE", nullptr };
    auto* dataset =
        driver->Create( filename.c_str(), static_cast< int >( width ),
            static_cast< int >( height ), 1, GDT_Float32, options.data() );
    std::array< double, 6 > transform{ 100, 10, 0, 500, 0, -10 };
    geode::OpenGeodeGeosciencesIOMeshException::test(
        dataset->SetGeoTransform( transform.data() ) == CE_None,
        "Cannot set raster transform" );
    if( epsg_code )
    {
        OGRSpatialReference reference;
        geode::OpenGeodeGeosciencesIOMeshException::test(
            reference.importFromEPSG( epsg_code.value() ) == OGRERR_NONE
                && dataset->SetSpatialRef( &reference ) == CE_None,
            "Cannot set raster spatial reference" );
    }
    auto* band = dataset->GetRasterBand( 1 );
    if( no_data )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            band->SetNoDataValue( no_data.value() ) == CE_None,
            "Cannot set raster no data value" );
    }
    if( !values.empty() )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            band->RasterIO( GF_Write, 0, 0, static_cast< int >( width ),
                static_cast< int >( height ),
                const_cast< float* >( values.data() ),
                static_cast< int >( width ), static_cast< int >( height ),
                GDT_Float32, 0, 0 )
                == CE_None,
            "Cannot write raster values" );
    }
    GDALClose( dataset );
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
        write_raster( "raster_geotiff", 3, 2, values, std::nullopt );
    const auto grid = geode::load_light_regular_grid< 2 >( filename );
    const auto colors = find_named_attribute< geode::RGBColor >(
        grid.cell_attribute_manager(), "RGB_data" );
    const auto band = find_named_attribute< double >(
        grid.cell_attribute_manager(), "band_1" );
    // Pixel values increase in row major order and colors are rescaled
    // values, so the red channel must increase in the same order
    int previous_red{ -1 };
    for( const auto row : geode::Range{ 2 } )
    {
        for( const auto column : geode::Range{ 3 } )
        {
            const auto cell =
                grid.cell_index( find_pixel_cell( grid, row, column ) );
            geode::OpenGeodeGeosciencesIOMeshException::test(
                band->value( cell ) == values[row * 3 + column],
                "Wrong band value for pixel ", row, " ", column );
            const int red = colors->value( cell ).red();
            geode::OpenGeodeGeosciencesIOMeshException::test(
                red > previous_red, "Wrong color for pixel ", row, " ",
                column );
            previous_red = red;
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
        write_raster( "raster_geotiff_data", 3, 2, values, nan );
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

void test_huge_raster()
{
    const auto filename =
        write_raster( "raster_geotiff_huge", 70000, 70000, {}, std::nullopt );
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
        thrown, "Loading a raster with too many pixels should throw" );
}

void test_coordinate_system()
{
    const std::vector< float > values{ 0, 1, 2, 3, 4, 5 };
    const auto filename =
        write_raster( "raster_geotiff_crs", 3, 2, values, std::nullopt, 32631 );
    const auto grid = geode::load_light_regular_grid< 2 >( filename );
    const auto& crs_manager = grid.grid_coordinate_reference_system_manager();
    geode::OpenGeodeGeosciencesIOMeshException::test(
        crs_manager.active_coordinate_reference_system_name() == "EPSG:32631",
        "Wrong GeoTIFF CRS name" );
    const auto* crs =
        dynamic_cast< const geode::GeographicGridCoordinateSystem2D* >(
            &crs_manager.active_coordinate_reference_system() );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        crs != nullptr && crs->info().authority_code() == "EPSG:32631"
            && crs->info().name == "WGS 84 / UTM zone 31N",
        "Wrong GeoTIFF geographic CRS" );

    const auto no_crs_filename =
        write_raster( "raster_geotiff_no_crs", 3, 2, values, std::nullopt );
    const auto no_crs_grid =
        geode::load_light_regular_grid< 2 >( no_crs_filename );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        no_crs_grid.grid_coordinate_reference_system_manager()
                .active_coordinate_reference_system_name()
            == "points",
        "GeoTIFF without CRS should keep the default CRS" );
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
        test_huge_raster();
        test_coordinate_system();

        geode::Logger::info( "[TEST SUCCESS]" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
