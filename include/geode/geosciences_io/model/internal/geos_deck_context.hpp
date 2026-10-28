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
#include <string>
#include <string_view>
#include <vector>

#include <absl/container/flat_hash_map.h>
#include <absl/types/span.h>

#include <geode/basic/uuid.hpp>

#include <geode/geosciences/simulation/core/fluid.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    class BRep;
} // namespace geode

namespace geode::internal
{
    /// GEOS cell region made of the model cells whose region attribute is id
    struct GeosRegion
    {
        index_t id{ NO_ID };
        std::string name;
        std::string component_name;
        uuid component_id;
        /// Blocks merged into this region
        std::vector< uuid > blocks;
    };

    /// Cell attribute of the exported mesh
    struct GeosCellArray
    {
        uuid attribute_id;
        std::string name;
    };

    /// Model data needed to translate a PhysicalProblem into a GEOS deck
    struct GeosDeckContext
    {
        std::vector< GeosRegion > solid_regions;
        /// GEOS region id of each Block
        absl::flat_hash_map< uuid, index_t > block_regions;
        /// Regions guessed from a BRep, not given by the exporter
        bool regions_are_inferred{ false };
        /// Cell array and mean value of the porosity PhysicalProperty
        std::optional< GeosCellArray > porosity_array;
        std::optional< double > mean_porosity;
        /// Cell array and mean value of the permeability PhysicalProperty
        std::optional< GeosCellArray > permeability_array;
        std::optional< std::array< double, 3 > > mean_permeability;
        /// Phases of the boundary fluid fractions carried by Surfaces
        std::vector< PHASE > boundary_fraction_phases;
        /// Boundary pressure carried by at least one Surface
        bool has_boundary_pressure{ false };
        /// Boundary temperature carried by at least one Surface
        bool has_boundary_temperature{ false };
    };

    /// Model array imported in GEOS, possibly copied under a new VTU name
    struct GeosImportedField
    {
        [[nodiscard]] bool is_copy() const
        {
            return vtu_name != source.name;
        }

        GeosCellArray source;
        std::string vtu_name;
        std::string geos_name;
    };

    /// Item (material or initial condition) assigned to each solid region
    struct GeosRegionAssignment
    {
        std::vector< std::optional< index_t > > region_items;
        std::vector< std::string > issues;
    };

    /*!
     * Assigns items to regions: items with regions first (the first one
     * wins), then the default item (no region) to the remaining regions. An
     * item region must cover whole GEOS regions.
     */
    [[nodiscard]] GeosRegionAssignment assign_geos_regions( const BRep& model,
        const GeosDeckContext& context,
        absl::Span< const absl::Span< const uuid > > item_regions,
        absl::Span< const std::string > item_names,
        std::string_view item_kind );
} // namespace geode::internal
