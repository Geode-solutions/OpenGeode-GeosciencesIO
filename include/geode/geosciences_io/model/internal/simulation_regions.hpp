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

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <absl/container/flat_hash_map.h>

#include <geode/basic/uuid.hpp>

#include <geode/geosciences/simulation/helpers/simulation_report.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    class BRep;
    class PhysicalProblem;
    class StructuralModel;
} // namespace geode

namespace geode::internal
{
    /// Simulator cell region made of the model cells whose region attribute is
    /// id
    struct SimulationRegion
    {
        index_t id{ NO_ID };
        std::string name;
        std::string component_name;
        uuid component_id;
    };

    /// Simulator cell regions of the exported Blocks
    struct SimulationRegions
    {
        /// Regions indexed by id
        std::vector< SimulationRegion > solid_regions;
        /// Region id of each Block
        absl::flat_hash_map< uuid, index_t > block_regions;
        /// Block, StratigraphicUnit or collection uuid -> region ids
        absl::flat_hash_map< uuid, std::vector< index_t > > regions_by_uuid;
        /// Blocks merged with other Blocks into one region (e.g. the
        /// Blocks of a StratigraphicUnit) -> that region id
        absl::flat_hash_map< uuid, index_t > merged_blocks;
    };

    /// Item (material or initial condition) assigned to each solid region
    struct RegionAssignment
    {
        std::vector< std::optional< index_t > > region_items;
        /// Assignments the simulator cannot write, in the regions stage
        std::vector< SimulationIssue > issues;
    };

    /// One region per Block
    [[nodiscard]] SimulationRegions brep_block_regions( const BRep& model );

    /// One region per StratigraphicUnit (merging its Blocks), then one per
    /// Block of no StratigraphicUnit. context prefixes the warning about
    /// Blocks in several StratigraphicUnits.
    [[nodiscard]] SimulationRegions structural_model_unit_regions(
        const StructuralModel& model, std::string_view context );

    /// Maps the uuid of each BlockCollection to the regions of its Blocks.
    /// A collection with a Block merged into another region is skipped.
    /// regions must come from brep_block_regions or
    /// structural_model_unit_regions.
    void add_block_collection_regions(
        const BRep& model, SimulationRegions& regions );

    /*!
     * Material of each solid region, with the materials the simulator cannot
     * assign. Blocks covered by several materials or by none are not issues
     * here: SimulationChecker reports them. simulator_name is used in the
     * issue texts.
     */
    [[nodiscard]] RegionAssignment assign_region_materials(
        const PhysicalProblem& problem,
        const BRep& model,
        const SimulationRegions& regions,
        std::string_view simulator_name );

    /// Initial condition of each solid region, assigned as materials are
    [[nodiscard]] RegionAssignment assign_region_initial_conditions(
        const PhysicalProblem& problem,
        const BRep& model,
        const SimulationRegions& regions,
        std::string_view simulator_name );
} // namespace geode::internal
