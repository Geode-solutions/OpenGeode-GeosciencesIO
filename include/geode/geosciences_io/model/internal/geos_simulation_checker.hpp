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
#include <vector>

#include <geode/geosciences/simulation/helpers/simulation_checker.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>

namespace geode::internal
{
    inline constexpr auto GEOS_SIMULATOR_KEY = "geos";

    /// What only the GEOS exporter knows about the export
    struct GeosExportData
    {
        /// GEOS regions of the exporter
        GeosRegions regions;
        /// Names of the well perforation sets given to the exporter
        std::vector< std::string > well_names;
    };

    /*!
     * Checks that a GEOS deck can be written for a problem on a model.
     * Registered as GEOS_SIMULATOR_KEY, it does not know the exporter wells
     * and regions: source term wells and the assignment of materials to GEOS
     * regions merging several Blocks (e.g. StratigraphicUnits) are only
     * checked with GeosExportData.
     */
    class GeosSimulationChecker final : public SimulationChecker
    {
    public:
        GeosSimulationChecker(
            const PhysicalProblem& problem, const BRep& model );
        GeosSimulationChecker( const PhysicalProblem& problem,
            const BRep& model,
            GeosExportData export_data );

    private:
        void check_simulator() final;

        void check_rock_properties();

        void check_boundary_properties();

    private:
        std::optional< GeosExportData > export_data_;
    };
} // namespace geode::internal
