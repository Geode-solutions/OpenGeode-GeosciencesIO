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

#include <geode/geosciences_io/model/internal/geos_simulation_checker.hpp>

#include <array>
#include <utility>

#include <absl/strings/str_cat.h>

#include <geode/mesh/core/regular_grid_solid.hpp>
#include <geode/mesh/core/solid_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/helpers/internal/simulation_checker_helpers.hpp>

#include <geode/geosciences_io/model/internal/geos_simulation_validation.hpp>

namespace
{
    using geode::PHYSICAL_PROPERTY_NAME;
    using geode::SIMULATION_CHECK_STAGE;

    constexpr std::array< std::pair< PHYSICAL_PROPERTY_NAME, geode::PHASE >, 3 >
        BOUNDARY_FRACTIONS{ {
            { PHYSICAL_PROPERTY_NAME::boundary_oil_fraction,
                geode::PHASE::oil },
            { PHYSICAL_PROPERTY_NAME::boundary_gas_fraction,
                geode::PHASE::gas },
            { PHYSICAL_PROPERTY_NAME::boundary_water_fraction,
                geode::PHASE::water },
        } };

    void check_regular_grids(
        const geode::BRep& model, geode::SimulationReport& report )
    {
        if( model.nb_blocks() == 0 )
        {
            return;
        }
        for( const auto& block : model.blocks() )
        {
            if( block.mesh().type_name()
                != geode::RegularGrid3D::type_name_static() )
            {
                return;
            }
        }
        geode::SimulationIssue issue{ geode::SIMULATION_ISSUE_SEVERITY::error,
            SIMULATION_CHECK_STAGE::mesh, "Blocks",
            "all Blocks are RegularGrids, GEOS needs an unstructured mesh" };
        issue.fix = "convert the Block meshes into HybridSolids";
        report.add_issue( std::move( issue ) );
    }

    // Boundary conditions written in the deck
    geode::internal::GeosDeckContext boundary_context(
        const geode::BRep& model )
    {
        geode::internal::GeosDeckContext context;
        context.has_boundary_pressure = model.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_pressure );
        context.has_boundary_temperature = model.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_temperature );
        for( const auto& [property, phase] : BOUNDARY_FRACTIONS )
        {
            if( model.has_physical_property( property ) )
            {
                context.boundary_fraction_phases.push_back( phase );
            }
        }
        return context;
    }

    void report_assignment_issues(
        const geode::internal::GeosRegionAssignment& assignment,
        std::string_view subject,
        geode::SimulationReport& report )
    {
        for( const auto& issue : assignment.issues )
        {
            report.add_error( SIMULATION_CHECK_STAGE::regions,
                std::string{ subject }, issue );
        }
    }
} // namespace

namespace geode::internal
{
    GeosSimulationChecker::GeosSimulationChecker() = default;

    GeosSimulationChecker::GeosSimulationChecker( GeosExportData export_data )
        : export_data_( std::move( export_data ) )
    {
    }

    void GeosSimulationChecker::check_simulator( const PhysicalProblem& problem,
        const BRep& model,
        SimulationReport& report ) const
    {
        report.set_simulator( "GEOS" );
        std::optional< absl::Span< const std::string > > well_names;
        if( export_data_ )
        {
            well_names = export_data_->well_names;
        }
        check_geos_problem( problem, well_names, report );
        check_blocks( model, report );
        check_unique_vertices_in_blocks( model, report );
        check_regular_grids( model, report );
        if( !problem.has_flow() || !problem.flow().has_fluid() )
        {
            return;
        }
        check_rock_properties( problem, model, report );
        check_boundary_properties( problem, model, report );
        check_geos_boundary_conditions(
            problem, boundary_context( model ), report );
        if( export_data_ )
        {
            report_assignment_issues(
                assign_geos_materials( problem, export_data_->regions ),
                "Materials", report );
            report_assignment_issues( assign_geos_initial_conditions(
                                          problem, export_data_->regions ),
                "Initial conditions", report );
        }
    }

    // GEOS imports a model porosity or permeability on every cell, merged
    // from the Block meshes. Without such a property, the deck uses the
    // value of each material.
    void GeosSimulationChecker::check_rock_properties(
        const PhysicalProblem& problem,
        const BRep& model,
        SimulationReport& report )
    {
        for( const auto property : { PHYSICAL_PROPERTY_NAME::porosity,
                 PHYSICAL_PROPERTY_NAME::permeability } )
        {
            PhysicalPropertyRequirement requirement{ property };
            requirement.merged = true;
            const auto coverage =
                check_physical_property( model, requirement, report );
            check_rock_property_availability(
                problem, model, coverage, report );
        }
    }

    void GeosSimulationChecker::check_boundary_properties(
        const PhysicalProblem& problem,
        const BRep& model,
        SimulationReport& report )
    {
        const auto pressure = check_physical_property( model,
            PhysicalPropertyRequirement{
                PHYSICAL_PROPERTY_NAME::boundary_pressure },
            report );
        check_physical_property_presence( model, pressure,
            "nothing drives the flow in the GEOS deck", report );
        static_cast< void >( check_physical_property( model,
            PhysicalPropertyRequirement{
                PHYSICAL_PROPERTY_NAME::boundary_temperature },
            report ) );
        if( !is_multiphase( problem ) )
        {
            return;
        }
        for( const auto& fraction : BOUNDARY_FRACTIONS )
        {
            static_cast< void >( check_physical_property( model,
                PhysicalPropertyRequirement{ fraction.first }, report ) );
        }
    }
} // namespace geode::internal
