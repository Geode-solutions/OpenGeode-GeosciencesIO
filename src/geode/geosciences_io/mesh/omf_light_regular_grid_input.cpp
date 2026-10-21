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

#include <geode/geosciences_io/mesh/internal/omf_light_regular_grid_input.hpp>

#include <cmath>

#include <absl/algorithm/container.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/identifier_builder.hpp>
#include <geode/basic/logger.hpp>

#include <geode/mesh/core/light_regular_grid.hpp>

#include <geode/geosciences_io/mesh/internal/omf_common.hpp>

namespace
{
    constexpr auto OMF_LIGHT_REGULAR_GRID_ELEMENT_CLASS = "VolumeElement";
    constexpr auto OMF_LIGHT_REGULAR_GRID_GEOMETRY_CLASS = "VolumeGridGeometry";
    constexpr std::array< std::string_view, 3 > OMF_TENSOR_KEYS{ "tensor_u",
        "tensor_v", "tensor_w" };
    constexpr std::array< std::string_view, 3 > OMF_AXIS_KEYS{ "axis_u",
        "axis_v", "axis_w" };

    /*
     * LightRegularGrid only supports a constant cell size along each axis
     */
    bool omf_is_uniform( const std::vector< double >& tensor )
    {
        const auto length = tensor.front();
        return absl::c_all_of( tensor, [length]( double cell_length ) {
            return std::fabs( cell_length - length ) <= geode::GLOBAL_EPSILON;
        } );
    }

    double omf_uniform_cell_length( const nlohmann::json& element,
        const std::vector< double >& tensor,
        std::string_view key )
    {
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            omf_is_uniform( tensor ), nullptr,
            geode::OpenGeodeException::TYPE::data, "[OMFInput] Volume \"",
            geode::internal::OMFFile::name( element ), "\" has a non uniform ",
            key, ", which cannot be represented by a LightRegularGrid" );
        return tensor.front();
    }

    /*
     * Data index of OMF volumes: C order over (u, v, w), w varying fastest
     */
    std::array< geode::index_t, 3 > omf_grid_indices(
        geode::index_t value, const std::array< geode::index_t, 3 >& sizes )
    {
        const auto w = value % sizes[2];
        const auto v = ( value / sizes[2] ) % sizes[1];
        const auto u = value / ( sizes[2] * sizes[1] );
        return { u, v, w };
    }
} // namespace

namespace geode::internal
{
    LightRegularGrid3D OMFLightRegularGridInput::read()
    {
        const OMFFile file{ filename() };
        const auto elements =
            file.elements( OMF_LIGHT_REGULAR_GRID_ELEMENT_CLASS,
                OMF_LIGHT_REGULAR_GRID_GEOMETRY_CLASS );
        OpenGeodeGeosciencesIOMeshException::check_exception( !elements.empty(),
            nullptr, OpenGeodeException::TYPE::data, "[OMFInput] ", filename(),
            " has no grid volume" );
        const auto& element = *elements.front();
        if( elements.size() > 1 )
        {
            Logger::warning( "[OMFInput] ", filename(), " has ",
                elements.size(), " grid volumes, only the first one (\"",
                OMFFile::name( element ), "\") is loaded" );
        }
        const auto& geometry = file.geometry( element );
        std::array< index_t, 3 > cells_number;
        std::array< Vector3D, 3 > directions;
        uint64_t nb_vertices{ 1 };
        for( const auto axis : LRange{ 3 } )
        {
            const auto tensor =
                OMFFile::read_tensor( geometry, OMF_TENSOR_KEYS[axis] );
            nb_vertices *= tensor.size() + 1;
            OpenGeodeGeosciencesIOMeshException::check_exception(
                nb_vertices < NO_ID, nullptr, OpenGeodeException::TYPE::data,
                "[OMFInput] Volume \"", OMFFile::name( element ),
                "\" has too many cells" );
            cells_number[axis] = static_cast< index_t >( tensor.size() );
            directions[axis] =
                OMFFile::axis( geometry, OMF_AXIS_KEYS[axis], axis )
                * omf_uniform_cell_length(
                    element, tensor, OMF_TENSOR_KEYS[axis] );
        }
        LightRegularGrid3D grid{ Point3D{ file.origin( geometry ) },
            cells_number, directions };
        IdentifierBuilder{ grid }.set_name( OMFFile::name( element ) );
        file.read_data( element, "cells", grid.nb_cells(),
            grid.cell_attribute_manager(),
            [&grid, &cells_number]( index_t value ) {
                return grid.cell_index(
                    omf_grid_indices( value, cells_number ) );
            } );
        const std::array< index_t, 3 > vertices_number{ cells_number[0] + 1,
            cells_number[1] + 1, cells_number[2] + 1 };
        file.read_data( element, "vertices", grid.nb_grid_vertices(),
            grid.grid_vertex_attribute_manager(),
            [&grid, &vertices_number]( index_t value ) {
                return grid.vertex_index(
                    omf_grid_indices( value, vertices_number ) );
            } );
        return grid;
    }

    Percentage OMFLightRegularGridInput::is_loadable() const
    {
        // Only the first volume is loaded, if its cells are uniform
        return OMFFile::is_loadable( filename(), []( const OMFFile& file ) {
            const auto elements =
                file.elements( OMF_LIGHT_REGULAR_GRID_ELEMENT_CLASS,
                    OMF_LIGHT_REGULAR_GRID_GEOMETRY_CLASS );
            if( elements.empty() )
            {
                return std::size_t{ 0 };
            }
            const auto& element = *elements.front();
            const auto& geometry = file.geometry( element );
            for( const auto& key : OMF_TENSOR_KEYS )
            {
                if( !omf_is_uniform( OMFFile::read_tensor( geometry, key ) ) )
                {
                    return std::size_t{ 0 };
                }
            }
            return std::size_t{ 1 };
        } );
    }
} // namespace geode::internal
