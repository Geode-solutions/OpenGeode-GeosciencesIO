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

namespace geode::internal
{
    /// What only the OGS exporter knows about the export
    struct OgsExportData
    {
        /// Names of the well perforation sets given to the exporter
        std::vector< std::string > well_names;
    };

    /*!
     * Checks that an OGS project can be written for a problem on a model.
     * Registered as OGS_SIMULATOR_KEY, it does not know the exporter wells:
     * source term wells are only checked with OgsExportData.
     */
    class OgsSimulationChecker final : public SimulationChecker
    {
    public:
        OgsSimulationChecker(
            const PhysicalProblem& problem, const BRep& model );
        OgsSimulationChecker( const PhysicalProblem& problem,
            const BRep& model,
            OgsExportData export_data );

    private:
        void check_simulator() final;

        void check_cells();

        void check_properties();

        /// Warns once per pair of Surfaces giving different values to a
        /// shared vertex
        void check_conflicting_boundaries( PHYSICAL_PROPERTY_NAME property );

    private:
        std::optional< OgsExportData > export_data_;
    };
} // namespace geode::internal
