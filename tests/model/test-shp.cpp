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

#include <string>
#include <vector>

#include <absl/strings/str_cat.h>

#include <ogrsf_frmts.h>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/string.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/corner.hpp>
#include <geode/model/mixin/core/line.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/core/section.hpp>
#include <geode/model/representation/io/section_input.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace
{
    std::string write_shapefile( std::string_view name,
        OGRwkbGeometryType type,
        const std::vector< std::string >& wkts )
    {
        GDALAllRegister();
        const auto filename = absl::StrCat( name, ".shp" );
        auto* driver =
            GetGDALDriverManager()->GetDriverByName( "ESRI Shapefile" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            driver != nullptr, "Shapefile driver not found" );
        driver->Delete( filename.c_str() );
        GDALDatasetUniquePtr dataset{ driver->Create(
            filename.c_str(), 0, 0, 0, GDT_Unknown, nullptr ) };
        geode::OpenGeodeGeosciencesIOModelException::test(
            dataset != nullptr, "Cannot create ", filename );
        auto* layer = dataset->CreateLayer(
            geode::to_string( name ).c_str(), nullptr, type, nullptr );
        for( const auto& wkt : wkts )
        {
            OGRFeature feature{ layer->GetLayerDefn() };
            if( !wkt.empty() )
            {
                OGRGeometry* geometry{ nullptr };
                OGRGeometryFactory::createFromWkt(
                    wkt.c_str(), nullptr, &geometry );
                geode::OpenGeodeGeosciencesIOModelException::test(
                    geometry != nullptr, "Invalid WKT ", wkt );
                feature.SetGeometryDirectly( geometry );
            }
            geode::OpenGeodeGeosciencesIOModelException::test(
                layer->CreateFeature( &feature ) == OGRERR_NONE,
                "Cannot write feature ", wkt );
        }
        return filename;
    }

    void check_section_vertices( const geode::Section& section )
    {
        for( const auto& line : section.lines() )
        {
            const auto& mesh = line.mesh();
            for( const auto v : geode::Range{ mesh.nb_vertices() } )
            {
                geode::OpenGeodeGeosciencesIOModelException::test(
                    section.unique_vertex( { line.component_id(), v } )
                        != geode::NO_ID,
                    "Line vertex without unique vertex" );
            }
            for( const auto e : geode::Range{ mesh.nb_edges() } )
            {
                geode::OpenGeodeGeosciencesIOModelException::test(
                    mesh.edge_vertex( { e, 0 } ) != mesh.edge_vertex( { e, 1 } )
                        && mesh.edge_vertex( { e, 1 } ) < mesh.nb_vertices(),
                    "Degenerated or invalid Line edge" );
            }
        }
        for( const auto& surface : section.surfaces() )
        {
            const auto& mesh = surface.mesh();
            for( const auto v : geode::Range{ mesh.nb_vertices() } )
            {
                geode::OpenGeodeGeosciencesIOModelException::test(
                    section.unique_vertex( { surface.component_id(), v } )
                        != geode::NO_ID,
                    "Surface vertex without unique vertex" );
            }
            for( const auto p : geode::Range{ mesh.nb_polygons() } )
            {
                geode::OpenGeodeGeosciencesIOModelException::test(
                    mesh.nb_polygon_vertices( p ) >= 3,
                    "Degenerated Surface polygon" );
            }
        }
    }

    void test_degenerated_parts()
    {
        const auto lines = geode::load_section(
            write_shapefile( "shp_degenerated_lines", wkbLineString,
                { "LINESTRING (0 0)", "LINESTRING (0 0,0 0)",
                    "LINESTRING (0 0,1 0,2 1)" } ) );
        check_section_vertices( lines );
        geode::index_t nb_edges{ 0 };
        for( const auto& line : lines.lines() )
        {
            nb_edges += line.mesh().nb_edges();
        }
        geode::OpenGeodeGeosciencesIOModelException::test(
            nb_edges == 2, "Wrong number of edges with degenerated lines" );

        const auto polygons = geode::load_section(
            write_shapefile( "shp_degenerated_polygons", wkbPolygon,
                { "POLYGON ((0 0,1 1,0 0))",
                    "POLYGON ((0 0,1 0,1 1,0 1,0 0))" } ) );
        check_section_vertices( polygons );
        geode::index_t nb_polygons{ 0 };
        for( const auto& surface : polygons.surfaces() )
        {
            nb_polygons += surface.mesh().nb_polygons();
        }
        geode::OpenGeodeGeosciencesIOModelException::test( nb_polygons == 1,
            "Wrong number of polygons with degenerated polygons" );
    }
    void test_points()
    {
        const auto section = geode::load_section( write_shapefile( "shp_points",
            wkbPoint, { "POINT (0 0)", "POINT (1 0)", "POINT (2 3)" } ) );
        geode::OpenGeodeGeosciencesIOModelException::test(
            section.nb_corners() == 3, "Wrong number of corners" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            section.nb_unique_vertices() == 3,
            "Wrong number of unique vertices" );
        for( const auto& corner : section.corners() )
        {
            geode::OpenGeodeGeosciencesIOModelException::test(
                corner.mesh().nb_vertices() == 1
                    && section.unique_vertex( { corner.component_id(), 0 } )
                           != geode::NO_ID,
                "Wrong Corner mesh" );
        }
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        const auto section = geode::load_section(
            absl::StrCat( geode::DATA_PATH, "toulouse.shz" ) );
        geode::OpenGeodeGeosciencesIOModelException::test(
            section.nb_corners() == 51718, "Wrong number of corners" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            section.nb_lines() == 10, "Wrong number of lines" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            section.nb_surfaces() == 12, "Wrong number of surfaces" );

        const auto section2 =
            geode::load_section( absl::StrCat( geode::DATA_PATH, "test.shp" ) );
        geode::OpenGeodeGeosciencesIOModelException::test(
            section2.nb_corners() == 0, "Wrong number of corners" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            section2.nb_lines() == 0, "Wrong number of lines" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            section2.nb_surfaces() == 1, "Wrong number of surfaces" );

        test_degenerated_parts();
        test_points();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
