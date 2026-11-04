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

#include <absl/container/flat_hash_map.h>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    class BRep;
    class PhysicalProblem;
    struct SourceTerm;
} // namespace geode

namespace geode::internal
{
    struct GeosBoundaryConditions;
    struct GeosRegionAssignment;
    struct GeosRegions;
} // namespace geode::internal

// GEOS rules shared by the deck writer and GeosSimulationChecker
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
        const PhysicalProblem& problem,
        const GeosBoundaryConditions& boundary_conditions );

    /// Name of the GEOS SourceFlux of a source term
    [[nodiscard]] std::string geos_source_term_name(
        const SourceTerm& source_term );

    /*!
     * Material of each solid region, with the materials GEOS cannot assign.
     * Blocks covered by several materials or by none are not issues here:
     * SimulationChecker reports them.
     */
    [[nodiscard]] GeosRegionAssignment assign_geos_materials(
        const PhysicalProblem& problem,
        const BRep& model,
        const GeosRegions& regions );

    /// Initial condition of each solid region, assigned as materials are
    [[nodiscard]] GeosRegionAssignment assign_geos_initial_conditions(
        const PhysicalProblem& problem,
        const BRep& model,
        const GeosRegions& regions );
} // namespace geode::internal
