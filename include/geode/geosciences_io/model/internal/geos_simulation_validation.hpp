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
#include <absl/types/span.h>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    class PhysicalProblem;
    class SimulationReport;
    struct SourceTerm;
} // namespace geode

namespace geode::internal
{
    struct GeosDeckContext;
    struct GeosRegionAssignment;
} // namespace geode::internal

namespace geode::internal
{
    /// XPath of the generated nodes that "geos" solver options can modify
    [[nodiscard]] const absl::flat_hash_map< std::string, std::string >&
        geos_solver_option_roles();

    /*!
     * Temperature of an isothermal GEOS simulation: the one shared by every
     * initial condition, if any
     */
    [[nodiscard]] std::optional< double > geos_isothermal_temperature(
        const PhysicalProblem& problem );

    /*!
     * GEOS compositional solvers need a temperature on every boundary face
     * set. Returns the temperature to impose when the model has none.
     */
    [[nodiscard]] std::optional< double > geos_default_boundary_temperature(
        const PhysicalProblem& problem, const GeosDeckContext& context );

    /// Name of the GEOS SourceFlux of a source term
    [[nodiscard]] std::string geos_source_term_name(
        const SourceTerm& source_term );

    /// Material of each solid region, with the assignment issues
    [[nodiscard]] GeosRegionAssignment assign_geos_materials(
        const PhysicalProblem& problem, const GeosDeckContext& context );

    /// Initial condition of each solid region, with the assignment issues
    [[nodiscard]] GeosRegionAssignment assign_geos_initial_conditions(
        const PhysicalProblem& problem, const GeosDeckContext& context );

    /*!
     * GEOS-specific checks of a problem that need no model data, in the
     * simulator stage. Source term wells are checked against the names of
     * the well perforations given to the exporter, if known.
     */
    void check_geos_problem( const PhysicalProblem& problem,
        std::optional< absl::Span< const std::string > > well_names,
        SimulationReport& report );

    /*!
     * GEOS rules on the boundary conditions of the context, in the
     * properties stage: boundary pressures, temperatures and fluid fractions
     * GEOS needs together
     */
    void check_geos_boundary_conditions( const PhysicalProblem& problem,
        const GeosDeckContext& context,
        SimulationReport& report );
} // namespace geode::internal
