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

#include <fstream>
#include <limits>
#include <tuple>
#include <vector>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/core/polygonal_surface.hpp>

#include <geode/mesh/io/polygonal_surface_input.hpp>
#include <geode/mesh/io/polygonal_surface_output.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

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

void check_normals( const geode::PolygonalSurface3D& surface )
{
    for( const auto p : geode::Range{ surface.nb_polygons() } )
    {
        const auto normal = surface.polygon_normal( p );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            normal && normal->value( 2 ) > 0,
            "Polygon normal should be upward" );
    }
}

void check_polygons( const geode::PolygonalSurface3D& surface )
{
    check_normals( surface );
    for( const auto p : geode::Range{ surface.nb_polygons() } )
    {
        for( const auto e : geode::LRange{ 4 } )
        {
            const auto& point =
                surface.point( surface.polygon_vertex( { p, e } ) );
            const auto& next = surface.point( surface.polygon_vertex(
                { p, static_cast< geode::local_index_t >( ( e + 1 ) % 4 ) } ) );
            const geode::Vector2D edge{ { next.value( 0 ) - point.value( 0 ),
                next.value( 1 ) - point.value( 1 ) } };
            geode::OpenGeodeGeosciencesIOMeshException::test(
                std::fabs( edge.length() - 10 ) < geode::GLOBAL_EPSILON,
                "Polygon edges should join neighbor pixels" );
        }
    }
}

void test_polygons()
{
    const auto filename = write_raster(
        "raster_dem_polygons", "dem", 3, 2, { 1, 2, 3, 4, 5, 6 }, "-9999" );
    const auto surface = geode::load_polygonal_surface< 3 >( filename );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface->nb_vertices() == 6 && surface->nb_polygons() == 2,
        "Wrong number of elements in raster DEM" );
    check_polygons( *surface );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface->point( 0 ).inexact_equal( geode::Point3D{ { 105, 495, 1 } } )
            && surface->point( 5 ).inexact_equal(
                geode::Point3D{ { 125, 485, 6 } } ),
        "DEM vertices should be at pixel centers" );
}

void test_huge_raster()
{
    const auto filename =
        write_raster( "raster_dem_huge", "dem", 70000, 70000, {}, "" );
    bool thrown{ false };
    try
    {
        const auto surface = geode::load_polygonal_surface< 3 >( filename );
    }
    catch( const geode::OpenGeodeException& )
    {
        thrown = true;
    }
    geode::OpenGeodeGeosciencesIOMeshException::test(
        thrown, "Loading a raster with too many pixels should throw" );
}

void test_no_data()
{
    const auto nan = std::numeric_limits< float >::quiet_NaN();
    for( const auto& [name, no_data, value] :
        { std::tuple{ "raster_dem_nan", "nan", nan },
            std::tuple{ "raster_dem_float", "0.1", 0.1F } } )
    {
        const auto filename = write_raster(
            name, "dem", 3, 2, { value, 2, 3, 4, 5, 6 }, no_data );
        const auto surface = geode::load_polygonal_surface< 3 >( filename );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 5 && surface->nb_polygons() == 1,
            "No data pixels should be skipped in ", name );
        check_polygons( *surface );
    }
}

int main()
{
    try
    {
        geode::Logger::set_level( geode::Logger::LEVEL::trace );
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        const auto surface = geode::load_polygonal_surface< 3 >(
            absl::StrCat( geode::DATA_PATH, "bathy_IrishSea_DEM.dem" ) );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 456427,
            "Number of vertices in the loaded Surface is not correct" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_polygons() == 453226,
            "Number of polygons in the loaded Surface is not correct" );
        geode::save_polygonal_surface(
            *surface, "bathy_IrishSea_DEM.og_psf3d" );
        check_normals( *surface );
        test_polygons();
        test_huge_raster();
        test_no_data();

        geode::Logger::info( "[TEST SUCCESS]" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
