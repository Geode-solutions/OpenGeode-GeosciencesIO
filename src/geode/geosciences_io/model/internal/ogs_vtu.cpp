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

#include <geode/geosciences_io/model/internal/ogs_vtu.hpp>

#include <string>

#include <pugixml.hpp>

#include <absl/strings/str_cat.h>

#include <geode/basic/range.hpp>

#include <geode/geosciences_io/model/internal/ogs_simulation_rules.hpp>
#include <geode/geosciences_io/model/internal/simulation_xml.hpp>

namespace
{
    // VTK cell types
    constexpr std::uint8_t VTK_VERTEX{ 1 };
    constexpr std::uint8_t VTK_TRIANGLE{ 5 };
    constexpr std::uint8_t VTK_POLYGON{ 7 };
    constexpr std::uint8_t VTK_QUAD{ 9 };
    constexpr std::uint8_t VTK_TETRA{ 10 };
    constexpr std::uint8_t VTK_HEXAHEDRON{ 12 };
    constexpr std::uint8_t VTK_WEDGE{ 13 };
    constexpr std::uint8_t VTK_PYRAMID{ 14 };
    constexpr geode::index_t NB_TRIANGLE_VERTICES{ 3 };
    constexpr geode::index_t NB_QUAD_VERTICES{ 4 };

    std::uint8_t vtk_polyhedron_type( geode::index_t nb_vertices )
    {
        switch( nb_vertices )
        {
            case geode::internal::NB_TETRAHEDRON_VERTICES:
                return VTK_TETRA;
            case geode::internal::NB_PYRAMID_VERTICES:
                return VTK_PYRAMID;
            case geode::internal::NB_PRISM_VERTICES:
                return VTK_WEDGE;
            default:
                return VTK_HEXAHEDRON;
        }
    }

    std::uint8_t vtk_cell_type( geode::index_t nb_vertices, bool volumetric )
    {
        if( volumetric )
        {
            return vtk_polyhedron_type( nb_vertices );
        }
        if( nb_vertices == 1 )
        {
            return VTK_VERTEX;
        }
        if( nb_vertices == NB_TRIANGLE_VERTICES )
        {
            return VTK_TRIANGLE;
        }
        if( nb_vertices == NB_QUAD_VERTICES )
        {
            return VTK_QUAD;
        }
        return VTK_POLYGON;
    }

    pugi::xml_node add_data_array(
        pugi::xml_node& parent, std::string_view type, std::string_view name )
    {
        auto array = parent.append_child( "DataArray" );
        geode::internal::set_attribute( array, "type", type );
        if( !name.empty() )
        {
            geode::internal::set_attribute( array, "Name", name );
        }
        geode::internal::set_attribute( array, "format", "ascii" );
        return array;
    }

    template < typename Container, typename Format >
    void set_values(
        pugi::xml_node& array, const Container& values, const Format& format )
    {
        std::string text;
        for( const auto& value : values )
        {
            if( !text.empty() )
            {
                text += ' ';
            }
            absl::StrAppend( &text, format( value ) );
        }
        array.text().set( text.c_str() );
    }

    void add_ids( pugi::xml_node& parent,
        std::string_view name,
        const std::vector< std::uint64_t >& ids )
    {
        if( ids.empty() )
        {
            return;
        }
        auto array = add_data_array( parent, "UInt64", name );
        set_values( array, ids, []( std::uint64_t value ) {
            return value;
        } );
    }

    void add_points(
        pugi::xml_node& piece, const geode::internal::OgsMeshFile& mesh )
    {
        auto points_node = piece.append_child( "Points" );
        auto array = add_data_array( points_node, "Float64", "" );
        geode::internal::set_attribute( array, "NumberOfComponents", "3" );
        std::vector< double > coordinates;
        coordinates.reserve( 3 * mesh.points.size() );
        for( const auto& point : mesh.points )
        {
            for( const auto axis : geode::LRange{ 3 } )
            {
                coordinates.push_back( point.value( axis ) );
            }
        }
        set_values( array, coordinates, geode::internal::xml_real );
    }

    void add_cells(
        pugi::xml_node& piece, const geode::internal::OgsMeshFile& mesh )
    {
        auto cells_node = piece.append_child( "Cells" );
        std::vector< geode::index_t > connectivity;
        std::vector< geode::index_t > offsets;
        std::vector< geode::index_t > types;
        for( const auto& cell : mesh.cells )
        {
            connectivity.insert( connectivity.end(), cell.begin(), cell.end() );
            offsets.push_back( connectivity.size() );
            types.push_back( vtk_cell_type( cell.size(), mesh.volumetric ) );
        }
        const auto as_is = []( geode::index_t value ) {
            return value;
        };
        auto connectivity_array =
            add_data_array( cells_node, "Int64", "connectivity" );
        set_values( connectivity_array, connectivity, as_is );
        auto offsets_array = add_data_array( cells_node, "Int64", "offsets" );
        set_values( offsets_array, offsets, as_is );
        auto types_array = add_data_array( cells_node, "UInt8", "types" );
        set_values( types_array, types, as_is );
    }
} // namespace

namespace geode::internal
{
    void save_ogs_mesh_file(
        std::string_view filename, const OgsMeshFile& mesh )
    {
        pugi::xml_document document;
        auto root = document.append_child( "VTKFile" );
        set_attribute( root, "type", "UnstructuredGrid" );
        set_attribute( root, "version", "1.0" );
        set_attribute( root, "byte_order", "LittleEndian" );
        set_attribute( root, "header_type", "UInt64" );
        auto grid = root.append_child( "UnstructuredGrid" );
        auto piece = grid.append_child( "Piece" );
        set_attribute(
            piece, "NumberOfPoints", absl::StrCat( mesh.points.size() ) );
        set_attribute(
            piece, "NumberOfCells", absl::StrCat( mesh.cells.size() ) );
        auto point_data = piece.append_child( "PointData" );
        add_ids( point_data, "bulk_node_ids", mesh.bulk_node_ids );
        for( const auto& [name, values] : mesh.point_arrays )
        {
            auto array = add_data_array( point_data, "Float64", name );
            set_values( array, values, xml_real );
        }
        auto cell_data = piece.append_child( "CellData" );
        add_ids( cell_data, "bulk_element_ids", mesh.bulk_element_ids );
        add_ids( cell_data, "bulk_face_ids", mesh.bulk_face_ids );
        for( const auto& [name, values] : mesh.cell_arrays )
        {
            auto array = add_data_array( cell_data, "Float64", name );
            set_values( array, values, xml_real );
        }
        add_points( piece, mesh );
        add_cells( piece, mesh );
        const auto file = std::string{ filename };
        OpenGeodeGeosciencesIOModelException::check_exception(
            document.save_file( file.c_str(), "  " ), nullptr,
            OpenGeodeException::TYPE::data, "[OgsExporter] Cannot write file ",
            file );
    }
} // namespace geode::internal
