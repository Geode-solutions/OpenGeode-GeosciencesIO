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

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <absl/strings/str_cat.h>
#include <absl/types/span.h>

#include <geode/basic/range.hpp>
#include <geode/basic/uuid.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/solver_options.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>
#include <geode/geosciences/simulation/helpers/internal/physical_problem_blocks.hpp>
#include <geode/geosciences/simulation/helpers/internal/simulation_checker_helpers.hpp>
#include <geode/geosciences/simulation/helpers/simulation_report.hpp>

#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>

namespace
{
    using geode::PHYSICAL_PROPERTY_NAME;

    constexpr std::array< std::pair< PHYSICAL_PROPERTY_NAME, geode::PHASE >, 3 >
        BOUNDARY_FRACTION_PHASES{ {
            { PHYSICAL_PROPERTY_NAME::boundary_oil_fraction,
                geode::PHASE::oil },
            { PHYSICAL_PROPERTY_NAME::boundary_gas_fraction,
                geode::PHASE::gas },
            { PHYSICAL_PROPERTY_NAME::boundary_water_fraction,
                geode::PHASE::water },
        } };

    // Assigns items (materials or initial conditions) to the GEOS regions:
    // items with regions first, then the default item (no region) to the
    // remaining regions
    class GeosRegionAssigner
    {
    public:
        GeosRegionAssigner( const geode::BRep& model,
            const geode::internal::GeosRegions& regions )
            : model_( model ), regions_( regions )
        {
            assignment_.region_items.resize( regions_.solid_regions.size() );
        }

        void assign( geode::index_t item,
            const std::string& subject,
            absl::Span< const geode::uuid > components )
        {
            for( const auto& component : components )
            {
                const auto regions = regions_.regions_by_uuid.find( component );
                if( regions == regions_.regions_by_uuid.end() )
                {
                    report_unknown_region( subject, component );
                    continue;
                }
                for( const auto region : regions->second )
                {
                    assignment_.region_items[region] = item;
                }
            }
        }

        geode::internal::GeosRegionAssignment finish(
            std::optional< geode::index_t > default_item )
        {
            for( auto& item : assignment_.region_items )
            {
                if( !item )
                {
                    item = default_item;
                }
            }
            return std::move( assignment_ );
        }

    private:
        void add_issue( geode::SimulationIssue issue )
        {
            assignment_.issues.push_back( std::move( issue ) );
        }

        // Regions which are not Blocks nor collections of Blocks are
        // reported by SimulationChecker::check_physical_problem_regions
        void report_unknown_region(
            const std::string& subject, const geode::uuid& component )
        {
            if( const auto merged = regions_.merged_blocks.find( component );
                merged != regions_.merged_blocks.end() )
            {
                report_merged_block( subject, component, merged->second );
                return;
            }
            const auto blocks = geode::internal::physical_problem_region_blocks(
                model_, component );
            if( !blocks )
            {
                return;
            }
            geode::SimulationIssue issue{
                geode::SIMULATION_ISSUE_SEVERITY::error,
                geode::SIMULATION_CHECK_STAGE::regions, subject,
                absl::StrCat( "region ", component.string(),
                    " contains Blocks merged into other GEOS regions" )
            };
            issue.fix = "assign the StratigraphicUnits, or Blocks of no "
                        "StratigraphicUnit";
            for( const auto& block : blocks.value() )
            {
                issue.add_location( model_.block( block ) );
            }
            add_issue( std::move( issue ) );
        }

        void report_merged_block( const std::string& subject,
            const geode::uuid& block,
            geode::index_t region )
        {
            const auto& component_name =
                regions_.solid_regions[region].component_name;
            geode::SimulationIssue issue{
                geode::SIMULATION_ISSUE_SEVERITY::error,
                geode::SIMULATION_CHECK_STAGE::regions, subject,
                absl::StrCat( "the Block is part of the GEOS region of ",
                    component_name, ", which merges several Blocks" )
            };
            issue.fix = absl::StrCat(
                "assign ", component_name, " instead of its Blocks" );
            issue.add_location( model_.block( block ) );
            add_issue( std::move( issue ) );
        }

    private:
        const geode::BRep& model_;
        const geode::internal::GeosRegions& regions_;
        geode::internal::GeosRegionAssignment assignment_;
    };

    geode::internal::GeosRegionAssignment assign_geos_regions(
        const geode::BRep& model,
        const geode::internal::GeosRegions& regions,
        absl::Span< const absl::Span< const geode::uuid > > item_regions,
        absl::Span< const std::string > item_subjects )
    {
        GeosRegionAssigner assigner{ model, regions };
        std::optional< geode::index_t > default_item;
        for( const auto item : geode::Indices{ item_regions } )
        {
            if( item_regions[item].empty() )
            {
                default_item = item;
                continue;
            }
            assigner.assign( item, item_subjects[item], item_regions[item] );
        }
        return assigner.finish( default_item );
    }
} // namespace

namespace geode::internal
{
    const absl::flat_hash_map< std::string, std::string >&
        geos_solver_option_roles()
    {
        static const absl::flat_hash_map< std::string, std::string > roles{
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

    GeosBoundaryConditions geos_boundary_conditions( const BRep& model )
    {
        GeosBoundaryConditions conditions;
        conditions.has_pressure = model.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_pressure );
        conditions.has_temperature = model.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_temperature );
        for( const auto& [property, phase] : BOUNDARY_FRACTION_PHASES )
        {
            if( model.has_physical_property( property ) )
            {
                conditions.fraction_phases.push_back( phase );
            }
        }
        return conditions;
    }

    GeosRegionAssignment assign_geos_materials( const PhysicalProblem& problem,
        const BRep& model,
        const GeosRegions& regions )
    {
        std::vector< absl::Span< const uuid > > material_regions;
        std::vector< std::string > subjects;
        for( const auto& material : problem.materials() )
        {
            material_regions.emplace_back( material.regions );
            subjects.push_back( material_subject( material.name ) );
        }
        return assign_geos_regions(
            model, regions, material_regions, subjects );
    }

    GeosRegionAssignment assign_geos_initial_conditions(
        const PhysicalProblem& problem,
        const BRep& model,
        const GeosRegions& regions )
    {
        std::vector< absl::Span< const uuid > > condition_regions;
        std::vector< std::string > subjects;
        for( const auto condition_id :
            Range{ problem.nb_initial_conditions() } )
        {
            condition_regions.emplace_back(
                problem.initial_condition( condition_id ).regions() );
            subjects.push_back( initial_condition_subject( condition_id ) );
        }
        return assign_geos_regions(
            model, regions, condition_regions, subjects );
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
        const GeosBoundaryConditions& boundary_conditions )
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
