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
#include <string_view>

#include <absl/types/span.h>

#include <geode/geosciences/simulation/helpers/physical_problem_simulation.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSet );
    ALIAS_3D( PointSet );
    class BRep;
    class PhysicalProblem;
} // namespace geode

namespace geode::internal
{
    inline constexpr auto GEOS_SIMULATOR_KEY = "geos";

    /// Name of the solid cell attribute holding the GEOS region ids
    inline constexpr auto GEOS_REGION_ATTRIBUTE_NAME = "attribute";

    struct GeosWell
    {
        std::string_view name;
        const PointSet3D& perforations;
    };

    struct GeosExportCheck
    {
        SimulationReport report;
        /// Regions and model data of the deck, mean values excepted
        GeosDeckContext context;
        GeosRegionAssignment materials;
        GeosRegionAssignment initial_conditions;
    };

    /*!
     * Every check of a GEOS export, done on the model without converting it:
     * model meshes and properties, then the problem. Without wells, source
     * term wells are not checked.
     * @param[in] regions Context holding the GEOS regions of the exporter
     */
    [[nodiscard]] GeosExportCheck check_geos_export( const BRep& model,
        GeosDeckContext regions,
        const PhysicalProblem& problem,
        std::optional< absl::Span< const GeosWell > > wells );

    /*!
     * Simulation checker registered as GEOS_SIMULATOR_KEY. StratigraphicUnits
     * of the model are merged into one GEOS region, as by the StructuralModel
     * exporter: an item covering part of a unit is only a warning.
     */
    class GeosSimulationChecker final : public SimulationChecker
    {
    public:
        [[nodiscard]] SimulationReport check(
            const PhysicalProblem& problem, const BRep& model ) const final;
    };
} // namespace geode::internal
