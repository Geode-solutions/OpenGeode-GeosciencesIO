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
#include <sstream>

#include <geode/basic/assert.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/image/core/rgb_color.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/light_regular_grid.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/polygonal_surface.hpp>
#include <geode/mesh/core/triangulated_surface.hpp>
#include <geode/mesh/io/edged_curve_input.hpp>
#include <geode/mesh/io/light_regular_grid_input.hpp>
#include <geode/mesh/io/point_set_input.hpp>
#include <geode/mesh/io/polygonal_surface_input.hpp>
#include <geode/mesh/io/triangulated_surface_input.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

namespace
{
    std::string data_file( std::string_view name )
    {
        return absl::StrCat( geode::DATA_PATH, name, ".omf" );
    }

    template < typename T >
    std::shared_ptr< geode::ReadOnlyAttribute< T > > find_attribute(
        const geode::AttributeManager& manager, std::string_view name )
    {
        const auto ids = manager.attribute_ids_matching_name( name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            ids.has_value() && ids->size() == 1, "Attribute ", name,
            " should exist once" );
        auto attribute = manager.find_read_only_attribute< T >( ids->front() );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            attribute != nullptr, "Attribute ", name, " has wrong type" );
        return attribute;
    }

    void check_no_attribute(
        const geode::AttributeManager& manager, std::string_view name )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !manager.attribute_ids_matching_name( name ).has_value(),
            "Attribute ", name, " should not exist" );
    }

    template < typename T >
    void check_values( const geode::AttributeManager& manager,
        std::string_view name,
        const std::vector< T >& values )
    {
        const auto attribute = find_attribute< T >( manager, name );
        for( const auto element : geode::Indices{ values } )
        {
            geode::OpenGeodeGeosciencesIOMeshException::test(
                attribute->value( element ) == values[element],
                "Wrong value of ", name, " for element ", element );
        }
    }

    void check_doubles( const geode::AttributeManager& manager,
        std::string_view name,
        const std::vector< double >& values )
    {
        const auto attribute = find_attribute< double >( manager, name );
        for( const auto element : geode::Indices{ values } )
        {
            if( std::isnan( values[element] ) )
            {
                geode::OpenGeodeGeosciencesIOMeshException::test(
                    !attribute->has_value( element ), "Element ", element,
                    " of ", name, " should have no value" );
                continue;
            }
            geode::OpenGeodeGeosciencesIOMeshException::test(
                attribute->value( element ) == values[element],
                "Wrong value of ", name, " for element ", element );
        }
    }

    /*
     * Write an OMF v1 file made of a header, binary blocks starting at byte
     * 60 and the JSON index
     */
    void write_omf( std::string_view filename,
        std::string_view json,
        std::string_view binary = {} )
    {
        std::string header( 60, '\0' );
        const std::array< unsigned char, 4 > magic{ 0x84, 0x83, 0x82, 0x81 };
        for( const auto byte : geode::LRange{ 4 } )
        {
            header[byte] = static_cast< char >( magic[byte] );
        }
        header.replace( 4, 10, "OMF-v0.9.0" );
        const auto json_start = static_cast< uint64_t >( 60 + binary.size() );
        for( const auto byte : geode::LRange{ 8 } )
        {
            header[52 + byte] =
                static_cast< char >( ( json_start >> ( 8 * byte ) ) & 0xFF );
        }
        std::ofstream file{ geode::to_string( filename ), std::ios::binary };
        file << header << binary << json;
    }

    void check_point( const geode::Point3D& point,
        const geode::Point3D& expected,
        std::string_view message )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            point.inexact_equal( expected ), message, ": got ", point.string(),
            " instead of ", expected.string() );
    }

    /*
     * Each reader returns the fraction of the project elements it loads
     */
    void check_loadable( const std::string& file,
        double point_set,
        double edged_curve,
        double triangulated_surface,
        double polygonal_surface,
        double grid )
    {
        const auto check = [&file]( geode::Percentage percentage,
                               double expected, std::string_view type ) {
            geode::OpenGeodeGeosciencesIOMeshException::test(
                std::fabs( percentage.value() - expected )
                    <= geode::GLOBAL_EPSILON,
                file, " loadability as ", type, " is ", percentage.value(),
                " instead of ", expected );
        };
        check(
            geode::is_point_set_loadable< 3 >( file ), point_set, "PointSet" );
        check( geode::is_edged_curve_loadable< 3 >( file ), edged_curve,
            "EdgedCurve" );
        check( geode::is_triangulated_surface_loadable< 3 >( file ),
            triangulated_surface, "TriangulatedSurface" );
        check( geode::is_polygonal_surface_loadable< 3 >( file ),
            polygonal_surface, "PolygonalSurface" );
        check( geode::is_light_regular_grid_loadable< 3 >( file ), grid,
            "LightRegularGrid" );
    }

    void test_point_set()
    {
        const auto file = data_file( "omf_points" );
        check_loadable( file, 1, 0, 0, 0, 0 );
        const auto point_set = geode::load_point_set< 3 >( file );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            point_set->nb_vertices() == 4, "Wrong number of points" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            point_set->name() == "collars", "Wrong point set name" );
        // project origin + geometry origin + vertex
        check_point( point_set->point( 0 ), geode::Point3D{ { 110, 200, 300 } },
            "Wrong point 0" );
        check_point( point_set->point( 3 ), geode::Point3D{ { 110, 200, 301 } },
            "Wrong point 3" );
        const auto& manager = point_set->vertex_attribute_manager();
        check_doubles( manager, "grade", { 0.5, NAN, 2.5, 3.5 } );
        check_values< std::string >( manager, "label", { "a", "b", "c", "d" } );
        check_values< geode::RGBColor >( manager, "rgb",
            { { 255, 0, 0 }, { 0, 255, 0 }, { 0, 0, 255 }, { 255, 0, 7 } } );
        check_values< int >( manager, "rock_index", { 1, 0, -1, 1 } );
        check_values< std::string >(
            manager, "rock", { "basalt", "granite", "", "basalt" } );
        check_values< geode::RGBColor >( manager, "rock_color",
            { { 40, 50, 60 }, { 10, 20, 30 }, {}, { 40, 50, 60 } } );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !find_attribute< int >( manager, "rock_index" )->has_value( 2 ),
            "Unmapped rock_index should have no value" );
        check_values< geode::Vector3D >( manager, "dir",
            { geode::Vector3D{ { 1, 0, 0 } }, geode::Vector3D{ { 0, 1, 0 } },
                geode::Vector3D{ { 0, 0, 1 } },
                geode::Vector3D{ { 1, 1, 1 } } } );
        check_no_attribute( manager, "omf_element" );
    }

    void test_edged_curve()
    {
        const auto file = data_file( "omf_lines" );
        check_loadable( file, 0, 1, 0, 0, 0 );
        const auto curve = geode::load_edged_curve< 3 >( file );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 5 && curve->nb_edges() == 3,
            "Wrong number of curve vertices or edges" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->name() == "omf_lines",
            "Merged curve should be named after the project" );
        check_point( curve->point( 3 ), geode::Point3D{ { 5, 0, 0 } },
            "Wrong first point of hole2" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->edge_vertex( { 2, 0 } ) == 3
                && curve->edge_vertex( { 2, 1 } ) == 4,
            "Wrong vertices of hole2 segment" );
        const auto& vertices = curve->vertex_attribute_manager();
        const auto& edges = curve->edge_attribute_manager();
        check_doubles( edges, "cu", { 0.1, 0.2, 0.3 } );
        check_doubles( vertices, "depth", { 0, 1, 3, NAN, NAN } );
        // same name as hole1 depth but other kind of data
        check_values< std::string >(
            vertices, "depth_string", { "", "", "", "top", "bottom" } );
        check_values< geode::index_t >(
            vertices, "omf_element", { 0, 0, 0, 1, 1 } );
        check_values< geode::index_t >( edges, "omf_element", { 0, 0, 1 } );
    }

    void test_corrupted_edged_curve()
    {
        const auto curve =
            geode::load_edged_curve< 3 >( data_file( "omf_lines_corrupted" ) );
        // hole1 second segment references a missing vertex and is skipped
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_vertices() == 5 && curve->nb_edges() == 2,
            "Invalid segment should be skipped" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->edge_vertex( { 1, 0 } ) == 3
                && curve->edge_vertex( { 1, 1 } ) == 4,
            "Wrong vertices of hole2 segment" );
        const auto& edges = curve->edge_attribute_manager();
        // value of the skipped segment is dropped and hole2 cu block is
        // outside the data section
        check_doubles( edges, "cu", { 0.1, NAN } );
        // depth of hole1 has a missing value
        check_doubles(
            curve->vertex_attribute_manager(), "depth", { 0, 1, NAN } );
    }

    void test_triangulated_surface()
    {
        const auto file = data_file( "omf_surfaces" );
        check_loadable( file, 0, 0, 0.5, 0.5, 0 );
        const auto surface = geode::load_triangulated_surface< 3 >( file );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 4 && surface->nb_polygons() == 2,
            "Wrong number of triangulated surface vertices or triangles" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->name() == "tri", "Wrong triangulated surface name" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->polygon_adjacent( { 0, 2 } ) == 1,
            "Triangle adjacencies should be computed" );
        check_doubles(
            surface->polygon_attribute_manager(), "area", { 0.5, 0.25 } );
    }

    void test_polygonal_surface()
    {
        const auto surface =
            geode::load_polygonal_surface< 3 >( data_file( "omf_surfaces" ) );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 6 && surface->nb_polygons() == 2,
            "Wrong number of grid surface vertices or polygons" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->name() == "grid", "Wrong grid surface name" );
        // axis_u = Y, axis_v = -X, so offsets are along Z; u varies fastest
        check_point( surface->point( 0 ), geode::Point3D{ { 10, 20, 30 } },
            "Wrong grid surface node (0, 0)" );
        check_point( surface->point( 2 ), geode::Point3D{ { 10, 23, 32 } },
            "Wrong grid surface node (2, 0)" );
        check_point( surface->point( 4 ), geode::Point3D{ { 7, 21, 34 } },
            "Wrong grid surface node (1, 1)" );
        const std::array< geode::index_t, 4 > second_quad{ 1, 2, 5, 4 };
        for( const auto vertex : geode::LRange{ 4 } )
        {
            geode::OpenGeodeGeosciencesIOMeshException::test(
                surface->polygon_vertex( { 1, vertex } ) == second_quad[vertex],
                "Wrong vertex ", vertex, " of second grid surface quad" );
        }
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->polygon_adjacent( { 0, 1 } ) == 1,
            "Quad adjacencies should be computed" );
        check_doubles( surface->vertex_attribute_manager(), "node_id",
            { 0, 1, 2, 3, 4, 5 } );
        check_doubles(
            surface->polygon_attribute_manager(), "face_id", { 10, 11 } );
    }

    void test_light_regular_grid()
    {
        const auto file = data_file( "omf_volume" );
        check_loadable( file, 0, 0, 0, 0, 1 );
        const auto grid = geode::load_light_regular_grid< 3 >( file );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            grid.nb_cells_in_direction( 0 ) == 3
                && grid.nb_cells_in_direction( 1 ) == 2
                && grid.nb_cells_in_direction( 2 ) == 4,
            "Wrong number of grid cells" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            grid.name() == "block_model", "Wrong grid name" );
        // axis_u = Y with 2m cells, axis_v = -X with 3m cells, 1m along Z
        check_point( grid.point( grid.vertex_index( { 1, 1, 1 } ) ),
            geode::Point3D{ { -2, 4, 4 } }, "Wrong grid vertex (1, 1, 1)" );
        // OMF volume data is in C order over (u, v, w)
        const auto cell_id = find_attribute< double >(
            grid.cell_attribute_manager(), "cell_id" );
        const auto density = find_attribute< double >(
            grid.cell_attribute_manager(), "density" );
        for( const auto u : geode::Range{ 3 } )
        {
            for( const auto v : geode::Range{ 2 } )
            {
                for( const auto w : geode::Range{ 4 } )
                {
                    const auto cell = grid.cell_index( { u, v, w } );
                    const auto value = u * 8 + v * 4 + w;
                    geode::OpenGeodeGeosciencesIOMeshException::test(
                        cell_id->value( cell ) == value,
                        "Wrong cell_id of cell (", u, ", ", v, ", ", w, ")" );
                    geode::OpenGeodeGeosciencesIOMeshException::test(
                        density->value( cell )
                            == ( value % 2 == 0 ? 2.5 : 2.7 ),
                        "Wrong density of cell (", u, ", ", v, ", ", w, ")" );
                }
            }
        }
        const auto vertex_id = find_attribute< double >(
            grid.grid_vertex_attribute_manager(), "vertex_id" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_id->value( grid.vertex_index( { 3, 2, 4 } ) ) == 59
                && vertex_id->value( grid.vertex_index( { 1, 2, 0 } ) ) == 25,
            "Wrong vertex_id" );
    }

    void test_non_uniform_grid()
    {
        check_loadable( data_file( "omf_volume_tensor" ), 0, 0, 0, 0, 0 );
        try
        {
            const auto grid = geode::load_light_regular_grid< 3 >(
                data_file( "omf_volume_tensor" ) );
        }
        catch( const geode::OpenGeodeException& )
        {
            return;
        }
        geode::OpenGeodeGeosciencesIOMeshException::test( false,
            "Non uniform grid should not be loaded as a LightRegularGrid" );
    }

    void test_mixed_loadability()
    {
        // 1 point set, 1 line set, 2 surfaces and 2 volumes: only the first
        // volume is loaded as a LightRegularGrid
        const auto file = data_file( "omf_mixed" );
        check_loadable( file, 1. / 6, 1. / 6, 2. / 6, 0, 1. / 6 );
        const auto grid = geode::load_light_regular_grid< 3 >( file );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            grid.name() == "volume1", "First volume should be loaded" );
    }

    void test_invalid_json_start()
    {
        std::ifstream input{ data_file( "omf_points" ), std::ios::binary };
        std::stringstream buffer;
        buffer << input.rdbuf();
        auto content = buffer.str();
        for( const auto byte : geode::Range{ 8 } )
        {
            content[52 + byte] = static_cast< char >( 0xFF );
        }
        constexpr auto filename = "omf_invalid_json_start.omf";
        std::ofstream{ filename, std::ios::binary } << content;
        check_loadable( filename, 0, 0, 0, 0, 0 );
        try
        {
            const auto point_set = geode::load_point_set< 3 >( filename );
        }
        catch( const geode::OpenGeodeException& )
        {
            return;
        }
        geode::OpenGeodeGeosciencesIOMeshException::test(
            false, "OMF file with invalid JSON start should not be loaded" );
    }
    void test_degenerate_elements()
    {
        constexpr auto filename = "omf_degenerate_elements.omf";
        write_omf( filename, R"({
            "p": { "__class__": "Project", "elements": [ "l", "s" ] },
            "l": { "__class__": "LineSetElement", "name": "lines",
                "geometry": "lg" },
            "lg": { "__class__": "LineSetGeometry", "vertices": "lv",
                "segments": "ls" },
            "lv": { "__class__": "Vector3Array",
                "array": [ 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0 ] },
            "ls": { "__class__": "Int2Array",
                "array": [ [ 0, 1 ], [ 3, 3 ], [ 1, 2 ] ] },
            "s": { "__class__": "SurfaceElement", "name": "surface",
                "geometry": "sg" },
            "sg": { "__class__": "SurfaceGeometry", "vertices": "sv",
                "triangles": "st" },
            "sv": { "__class__": "Vector3Array",
                "array": [ 0, 0, 0, 1, 0, 0, 0, 1, 0 ] },
            "st": { "__class__": "Int3Array",
                "array": [ [ 0, 1, 2 ], [ 0, 0, 1 ] ] }
        })" );
        const auto curve = geode::load_edged_curve< 3 >( filename );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            curve->nb_edges() == 2, "Degenerate segment should be skipped" );
        const auto surface = geode::load_triangulated_surface< 3 >( filename );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_polygons() == 1,
            "Degenerate triangle should be skipped" );
    }
    template < typename Loader >
    void check_load_throws( std::string_view message, Loader&& loader )
    {
        try
        {
            loader();
        }
        catch( const geode::OpenGeodeException& )
        {
            return;
        }
        geode::OpenGeodeGeosciencesIOMeshException::test( false, message );
    }

    void test_collinear_axes()
    {
        constexpr auto surface_file = "omf_collinear_surface_axes.omf";
        write_omf( surface_file, R"({
            "p": { "__class__": "Project", "elements": [ "s" ] },
            "s": { "__class__": "SurfaceElement", "name": "grid",
                "geometry": "sg" },
            "sg": { "__class__": "SurfaceGridGeometry", "tensor_u": [ 1, 1 ],
                "tensor_v": [ 1 ], "axis_u": [ 1, 0, 0 ],
                "axis_v": [ 2, 0, 0 ] }
        })" );
        check_load_throws(
            "Surface grid with collinear axes should throw", [&] {
                return geode::load_polygonal_surface< 3 >( surface_file );
            } );
        constexpr auto volume_file = "omf_coplanar_volume_axes.omf";
        write_omf( volume_file, R"({
            "p": { "__class__": "Project", "elements": [ "v" ] },
            "v": { "__class__": "VolumeElement", "name": "grid",
                "geometry": "vg" },
            "vg": { "__class__": "VolumeGridGeometry", "tensor_u": [ 1 ],
                "tensor_v": [ 1 ], "tensor_w": [ 1 ], "axis_u": [ 1, 0, 0 ],
                "axis_v": [ 0, 1, 0 ], "axis_w": [ 1, 1, 0 ] }
        })" );
        check_load_throws( "Volume grid with coplanar axes should throw", [&] {
            return geode::load_light_regular_grid< 3 >( volume_file );
        } );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::trace );
        test_point_set();
        test_edged_curve();
        test_corrupted_edged_curve();
        test_triangulated_surface();
        test_polygonal_surface();
        test_light_regular_grid();
        test_non_uniform_grid();
        test_invalid_json_start();
        test_mixed_loadability();
        test_degenerate_elements();
        test_collinear_axes();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
