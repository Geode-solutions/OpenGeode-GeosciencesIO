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

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include <pugixml.hpp>

#include <absl/strings/str_cat.h>

#include <geode/basic/logger.hpp>
#include <geode/tests_config.hpp>

#include <geode/geometry/point.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/ogs_vtu.hpp>

namespace
{
    pugi::xml_document load_xml( const std::string& filename )
    {
        pugi::xml_document document;
        const auto result = document.load_file( filename.c_str() );
        geode::OpenGeodeGeosciencesIOModelException::test(
            static_cast< bool >( result ), "[Test] Cannot parse ", filename );
        return document;
    }

    pugi::xml_node data_array( const pugi::xml_document& document,
        std::string_view section,
        std::string_view name )
    {
        const auto xpath = absl::StrCat( "//", section, "/DataArray[@Name='",
            name, "']" );
        return document.select_node( xpath.c_str() ).node();
    }

    std::string text( const pugi::xml_node& node )
    {
        return node.text().as_string();
    }

    void check( bool condition, std::string_view message )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            condition, "[Test] ", message );
    }
} // namespace

void test_ogs_mesh_file()
{
    geode::internal::OgsMeshFile mesh;
    mesh.points = { geode::Point3D{ { 0, 0, 0 } },
        geode::Point3D{ { 1, 0, 0 } }, geode::Point3D{ { 1, 1, 0 } },
        geode::Point3D{ { 0, 1, 0 } } };
    mesh.cells = { { 0, 1, 2 }, { 0, 2, 3 } };
    mesh.bulk_node_ids = { 7, 3, 9, 4 };
    mesh.point_arrays.emplace_back(
        "values", std::vector< double >{ 1.5, 2.5, 3.5, 4.5 } );
    geode::internal::save_ogs_mesh_file( "ogs_mesh_file.vtu", mesh );
    const auto document = load_xml( "ogs_mesh_file.vtu" );
    const auto piece = document.select_node( "//Piece" ).node();
    check( piece.attribute( "NumberOfPoints" ).as_uint() == 4,
        "Wrong number of points" );
    check( piece.attribute( "NumberOfCells" ).as_uint() == 2,
        "Wrong number of cells" );
    const auto ids = data_array( document, "PointData", "bulk_node_ids" );
    check( std::string{ ids.attribute( "type" ).value() } == "UInt64",
        "bulk_node_ids must be UInt64" );
    check( text( ids ) == "7 3 9 4", "Wrong bulk_node_ids" );
    const auto values = data_array( document, "PointData", "values" );
    check( std::string{ values.attribute( "type" ).value() } == "Float64",
        "Point arrays must be Float64" );
    check( text( data_array( document, "Cells", "types" ) ) == "5 5",
        "Triangles must be VTK type 5" );
    check( !data_array( document, "CellData", "bulk_element_ids" ),
        "Empty cell id arrays must not be written" );

    geode::internal::OgsMeshFile vertices;
    vertices.points = { geode::Point3D{ { 0.5, 0.5, 0.5 } } };
    vertices.cells = { { 0 } };
    vertices.bulk_node_ids = { 12 };
    geode::internal::save_ogs_mesh_file( "ogs_vertex_file.vtu", vertices );
    const auto vertex_document = load_xml( "ogs_vertex_file.vtu" );
    check( text( data_array( vertex_document, "Cells", "types" ) ) == "1",
        "Single-vertex cells must be VTK type 1" );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        test_ogs_mesh_file();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
