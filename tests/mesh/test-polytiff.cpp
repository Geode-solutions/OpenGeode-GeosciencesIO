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
#include <limits>
#include <optional>
#include <tuple>
#include <vector>

#include <gdal_priv.h>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/core/polygonal_surface.hpp>

#include <geode/mesh/io/polygonal_surface_input.hpp>
#include <geode/mesh/io/polygonal_surface_output.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

std::string write_raster( std::string_view name,
    geode::index_t width,
    geode::index_t height,
    const std::vector< float >& values,
    std::optional< double > no_data )
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
        "raster_polytiff_polygons", 3, 2, { 1, 2, 3, 4, 5, 6 }, -9999. );
    const auto surface = geode::load_polygonal_surface< 3 >( filename );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface->nb_vertices() == 6 && surface->nb_polygons() == 2,
        "Wrong number of elements in raster PolyTIFF" );
    check_polygons( *surface );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface->point( 0 ).inexact_equal( geode::Point3D{ { 105, 495, 1 } } )
            && surface->point( 5 ).inexact_equal(
                geode::Point3D{ { 125, 485, 6 } } ),
        "PolyTIFF vertices should be at pixel centers" );
}

void test_huge_raster()
{
    const auto filename =
        write_raster( "raster_polytiff_huge", 70000, 70000, {}, std::nullopt );
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
        { std::tuple{
              "raster_polytiff_nan", static_cast< double >( nan ), nan },
            std::tuple{ "raster_polytiff_float", 0.1, 0.1F } } )
    {
        const auto filename =
            write_raster( name, 3, 2, { value, 2, 3, 4, 5, 6 }, no_data );
        const auto surface = geode::load_polygonal_surface< 3 >( filename );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 5 && surface->nb_polygons() == 1,
            "No data pixels should be skipped in ", name );
        check_polygons( *surface );
    }
}

void test_pixel_values()
{
    const auto nan = std::numeric_limits< float >::quiet_NaN();
    const auto filename =
        write_raster( "raster_polytiff", 3, 2, { 1, 2, 3, 4, nan, 6 }, nan );
    const auto surface = geode::load_polygonal_surface< 3 >( filename );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface->nb_vertices() == 5 && surface->nb_polygons() == 0,
        "Wrong number of elements in raster PolyTIFF" );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface->point( 0 ).inexact_equal( geode::Point3D{ { 105, 495, 1 } } )
            && surface->point( 2 ).inexact_equal(
                geode::Point3D{ { 125, 495, 3 } } )
            && surface->point( 4 ).inexact_equal(
                geode::Point3D{ { 125, 485, 6 } } ),
        "PolyTIFF vertices should be at pixel centers with pixel values" );
    const auto filename2 = write_raster(
        "raster_polytiff_full", 3, 2, { 1, 2, 3, 4, 5, 6 }, std::nullopt );
    const auto surface2 = geode::load_polygonal_surface< 3 >( filename2 );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface2->nb_vertices() == 6 && surface2->nb_polygons() == 2,
        "Wrong number of elements in full raster PolyTIFF" );
}

int main()
{
    try
    {
        geode::Logger::set_level( geode::Logger::LEVEL::trace );
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        const auto surface = geode::load_polygonal_surface< 3 >(
            absl::StrCat( geode::DATA_PATH, "cea.tiff" ) );
        geode::save_polygonal_surface( *surface, "cea.og_psf3d" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 264710,
            "Number of vertices in the loaded Surface is not correct" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_polygons() == 263682,
            "Number of polygons in the loaded Surface is not correct" );

        test_pixel_values();
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
