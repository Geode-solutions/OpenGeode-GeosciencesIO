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

#include <geode/geosciences_io/model/internal/geos_simulation_rules.hpp>

#include <optional>
#include <string>

#include <absl/strings/str_cat.h>

#include <geode/basic/range.hpp>

#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>
#include <geode/geosciences/simulation/helpers/internal/simulation_checker_helpers.hpp>

#include <geode/geosciences_io/model/internal/simulation_boundary_values.hpp>

namespace geode::internal
{
    const SolverOptionRoles& geos_solver_option_roles()
    {
        static const SolverOptionRoles roles{
            { "solvers", "/Problem/Solvers" },
            { "flow", "/Problem/Solvers/*[@name='flow']" },
            { "mechanics", "/Problem/Solvers/*[@name='mechanics']" },
            { "poromechanics", "/Problem/Solvers/*[@name='poromechanics']" },
            { "poromechanicsNonlinear", "/Problem/Solvers/"
                                        "*[@name='poromechanics']/"
                                        "NonlinearSolverParameters" },
            { "poromechanicsLinear", "/Problem/Solvers/"
                                     "*[@name='poromechanics']/"
                                     "LinearSolverParameters" },
            { "nonlinear", "/Problem/Solvers/*/NonlinearSolverParameters" },
            { "linear", "/Problem/Solvers/*/LinearSolverParameters" },
            { "mesh", "/Problem/Mesh/VTKMesh" },
            { "tpfa", "/Problem/NumericalMethods/FiniteVolume/"
                      "TwoPointFluxApproximation" },
            { "fluid", "/Problem/Constitutive/*[@name='fluid']" },
            { "events", "/Problem/Events" },
            { "solverEvent",
                "/Problem/Events/PeriodicEvent[@name='solverApplications']" },
            { "vtk", "/Problem/Outputs/VTK" },
        };
        return roles;
    }

    SimulatorSolverOptions geos_solver_options()
    {
        return { "geos", "GEOS", "[GeosExporter]", geos_solver_option_roles() };
    }

    std::optional< double > geos_isothermal_temperature(
        const PhysicalProblem& problem )
    {
        std::optional< double > temperature;
        for( const auto condition_id :
            Range{ problem.nb_initial_conditions() } )
        {
            const auto* condition =
                dynamic_cast< const UniformInitialCondition* >(
                    &problem.initial_condition( condition_id ) );
            if( condition == nullptr || !condition->temperature()
                || ( temperature
                     && temperature.value()
                            != condition->temperature().value() ) )
            {
                return std::nullopt;
            }
            temperature = condition->temperature();
        }
        return temperature;
    }

    std::optional< double > geos_default_boundary_temperature(
        const PhysicalProblem& problem,
        const ModelBoundaryConditions& boundary_conditions )
    {
        if( !boundary_conditions.has_pressure
            || boundary_conditions.has_temperature || problem.thermal()
            || !is_multiphase( problem ) )
        {
            return std::nullopt;
        }
        return geos_isothermal_temperature( problem );
    }

    std::string geos_source_term_name( const SourceTerm& source_term )
    {
        return absl::StrCat(
            "sourceFlux_", source_term.well, "_", source_term.component );
    }
} // namespace geode::internal
