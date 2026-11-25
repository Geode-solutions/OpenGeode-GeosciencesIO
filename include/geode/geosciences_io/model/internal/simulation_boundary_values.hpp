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

#include <absl/container/flat_hash_map.h>
#include <absl/container/linked_hash_map.h>
#include <absl/types/span.h>

#include <geode/basic/uuid.hpp>

#include <geode/geosciences/simulation/core/fluid.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    class BRep;
    FORWARD_DECLARATION_DIMENSION_CLASS( Surface );
    ALIAS_3D( Surface );
} // namespace geode

// Boundary condition values of a model, shared by the simulator exporters
namespace geode::internal
{
    inline constexpr std::array< std::string_view, 3 > AXIS_NAMES{ "x", "y",
        "z" };

    // Value of the attribute on each Surface polygon, reduced to a single value
    // for a ConstantAttribute. context prefixes the exception messages.
    template < typename T >
    [[nodiscard]] std::optional< std::vector< T > > surface_values(
        const Surface3D& surface,
        const uuid& attribute_id,
        std::string_view context );

    template < typename T >
    struct VariableSurface
    {
        const geode::Surface3D* surface;
        geode::index_t region_id;
        std::vector< T > values;
    };

    // Boundary condition values of the Surfaces: regions of the uniform
    // Surfaces grouped by value, and the other Surfaces with their values
    template < typename T >
    struct SurfaceValueGroups
    {
        absl::linked_hash_map< T, std::vector< geode::index_t > >
            uniform_regions;
        std::vector< VariableSurface< T > > variable_surfaces;
    };

    template < typename T >
    [[nodiscard]] SurfaceValueGroups< T > group_surface_values(
        const BRep& model,
        const absl::flat_hash_map< uuid, index_t >& surface_regions,
        const uuid& attribute_id,
        std::string_view context );

    [[nodiscard]] std::vector< double > axis_values(
        absl::Span< const std::array< double, 3 > > values,
        local_index_t axis );

    /// Boundary condition PhysicalProperties of the model
    struct ModelBoundaryConditions
    {
        /// Phases of the boundary fluid fraction PhysicalProperties
        std::vector< PHASE > fraction_phases;
        bool has_pressure{ false };
        bool has_temperature{ false };
    };

    [[nodiscard]] ModelBoundaryConditions model_boundary_conditions(
        const BRep& model );
} // namespace geode::internal
