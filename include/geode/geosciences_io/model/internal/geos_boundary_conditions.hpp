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

#pragma once

#include <array>
#include <optional>
#include <string_view>
#include <vector>

#include <absl/container/linked_hash_map.h>
#include <absl/types/span.h>

#include <geode/geometry/point.hpp>

#include <geode/model/mixin/core/physical_properties.hpp>

#include <geode/geosciences/simulation/core/fluid.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    class PhysicalProblem;
} // namespace geode

namespace geode::internal
{
    struct GeosBoundaryField
    {
        std::string_view field_name;
        std::string_view name_prefix;
    };

    /// GEOS face field of the pressure and temperature boundary conditions
    [[nodiscard]] const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME,
        GeosBoundaryField >&
        geos_boundary_condition_fields();

    struct GeosFluidFraction
    {
        std::string_view fluid;
        PHASE phase;
    };

    /// Phase of each boundary fluid fraction property
    [[nodiscard]] const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME,
        GeosFluidFraction >&
        geos_fluid_fractions();

    /// GEOS component of a boundary phase fraction: the component carried by
    /// the phase in the problem fluid
    [[nodiscard]] std::optional< local_index_t >
        geos_boundary_fraction_component(
            const PhysicalProblem& problem, const GeosFluidFraction& fraction );

    /*!
     * Boundary condition properties the GEOS solver of the problem uses:
     * pressure, temperature if thermal or multiphase, and the fraction of
     * every fluid phase if multiphase.
     */
    [[nodiscard]] std::vector< PHYSICAL_PROPERTY_NAME >
        geos_boundary_properties( const PhysicalProblem& problem );

    using TableCoordinates = std::array< std::vector< double >, 3 >;

    /*!
     * Table coordinates looking up every polygon center on its own node with
     * GEOS nearest interpolation. No value if two centers cannot be
     * separated.
     */
    [[nodiscard]] std::optional< TableCoordinates >
        separating_table_coordinates( absl::Span< const Point3D > centers );

    /// Node of each polygon center, in Fortran order
    [[nodiscard]] std::vector< index_t > table_polygon_nodes(
        const TableCoordinates& coordinates,
        absl::Span< const Point3D > centers );

    /*!
     * Polygon whose value is stored at each table node, in Fortran order.
     * Nodes looked up by no polygon center take the value of the closest one.
     */
    [[nodiscard]] std::vector< index_t > table_node_polygons(
        const TableCoordinates& coordinates,
        absl::Span< const index_t > polygon_nodes,
        absl::Span< const Point3D > centers );
} // namespace geode::internal
