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

#include <geode/geosciences_io/model/internal/geos_simulation_validation.hpp>

#include <algorithm>
#include <limits>
#include <optional>

#include <absl/algorithm/container.h>
#include <absl/container/btree_set.h>
#include <absl/container/flat_hash_set.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>

#include <geode/basic/range.hpp>

#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/schedule.hpp>
#include <geode/geosciences/simulation/core/solver_options.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>
#include <geode/geosciences/simulation/helpers/physical_problem_blocks.hpp>
#include <geode/geosciences/simulation/helpers/physical_problem_simulation.hpp>

#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_writer.hpp>
#include <geode/geosciences_io/model/internal/geos_xml.hpp>

namespace
{
    template < typename... Types, typename Base >
    bool is_one_of( const Base& object )
    {
        return (
            ... || ( dynamic_cast< const Types* >( &object ) != nullptr ) );
    }

    bool is_geos_multiphase( const geode::PhysicalProblem& problem )
    {
        return dynamic_cast< const geode::MultiphaseFlow* >( &problem.flow() )
               != nullptr;
    }

    class GeosValidator
    {
    public:
        GeosValidator( const geode::PhysicalProblem& problem,
            geode::SimulationReport& report )
            : problem_( problem ), report_( report )
        {
        }

        void check_problem(
            std::optional< absl::Span< const std::string > > well_names )
        {
            check_outputs();
            check_solver_options();
            check_source_term_wells( well_names );
            if( !has_fluid() )
            {
                return;
            }
            check_supported_flow();
            check_supported_materials();
            check_supported_initial_conditions();
            check_thermal();
            check_source_terms();
            check_ignored_fields();
            check_validity_ranges();
        }

        void check_deck( const geode::internal::GeosDeckContext& context )
        {
            if( !has_fluid() )
            {
                return;
            }
            check_boundary_conditions( context );
        }

    private:
        template < typename... Args >
        void add_issue( const Args&... args )
        {
            report_.add_error( args... );
        }

        template < typename... Args >
        void add_warning( const Args&... args )
        {
            report_.add_warning( args... );
        }

        template < typename Object >
        void add_unsupported( const Object& object )
        {
            add_issue(
                object.type_name(), " is not supported by the GEOS writer" );
        }

        [[nodiscard]] bool has_fluid() const
        {
            return problem_.has_flow() && problem_.flow().has_fluid();
        }

        [[nodiscard]] const geode::FluidModel& fluid() const
        {
            return problem_.fluid();
        }

        [[nodiscard]] bool is_multiphase() const
        {
            return is_geos_multiphase( problem_ );
        }

        void check_outputs()
        {
            const auto& output = problem_.output();
            if( !geode::internal::geos_writes_vtk_output( output )
                && !geode::internal::geos_writes_restart_output( output ) )
            {
                add_issue( "The problem has no output: set an output "
                           "frequency, output times or restarts, otherwise "
                           "GEOS writes no result" );
            }
        }

        void check_supported_flow()
        {
            if( !is_one_of< geode::SinglePhaseFlow, geode::MultiphaseFlow >(
                    problem_.flow() ) )
            {
                add_unsupported( problem_.flow() );
            }
            if( !is_one_of< geode::SinglePhaseFluid, geode::DeadOilFluid,
                    geode::Co2BrineFluid >( fluid() ) )
            {
                add_unsupported( fluid() );
            }
            if( !is_multiphase() )
            {
                return;
            }
            for( const auto& component : fluid().components() )
            {
                if( !geode::internal::is_geos_name( component.name ) )
                {
                    add_issue( "Fluid component name \"", component.name,
                        "\" must only contain letters, digits, '_' and '-'" );
                }
            }
        }

        void check_supported_materials()
        {
            for( const auto& material : problem_.materials() )
            {
                if( !geode::internal::is_geos_name( material.name ) )
                {
                    add_issue( "Material name \"", material.name,
                        "\" must only contain letters, digits, '_' and '-'" );
                }
                if( material.relative_permeability )
                {
                    check_supported_relperm( *material.relative_permeability );
                }
                if( material.capillary_pressure )
                {
                    check_supported_capillary( *material.capillary_pressure );
                }
            }
        }

        void check_supported_relperm(
            const geode::RelativePermeability& relperm )
        {
            if( const auto* table =
                    dynamic_cast< const geode::TableRelativePermeability* >(
                        &relperm ) )
            {
                if( table->tables().size() != 1 )
                {
                    add_issue( "The GEOS writer only supports two-phase "
                               "relative permeability tables" );
                }
                return;
            }
            if( !is_one_of< geode::BrooksCoreyRelativePermeability >(
                    relperm ) )
            {
                add_unsupported( relperm );
            }
        }

        void check_supported_capillary(
            const geode::CapillaryPressure& capillary )
        {
            if( const auto* table =
                    dynamic_cast< const geode::TableCapillaryPressure* >(
                        &capillary ) )
            {
                if( table->tables().size() != 1 )
                {
                    add_issue( "The GEOS writer only supports two-phase "
                               "capillary pressure tables" );
                }
                return;
            }
            if( !is_one_of< geode::BrooksCoreyCapillaryPressure >( capillary ) )
            {
                add_unsupported( capillary );
            }
        }

        void check_supported_initial_conditions()
        {
            for( const auto condition_id :
                geode::Range{ problem_.nb_initial_conditions() } )
            {
                const auto& condition =
                    problem_.initial_condition( condition_id );
                if( !is_one_of< geode::UniformInitialCondition >( condition ) )
                {
                    add_unsupported( condition );
                }
            }
        }

        void check_thermal()
        {
            if( problem_.thermal() )
            {
                if( is_one_of< geode::DeadOilFluid >( fluid() ) )
                {
                    add_issue( "The GEOS writer does not support thermal "
                               "simulations with DeadOilFluid" );
                }
                return;
            }
            if( is_multiphase()
                && !geode::internal::geos_isothermal_temperature( problem_ ) )
            {
                add_issue( "GEOS isothermal multiphase simulations need the "
                           "same temperature in every initial condition" );
            }
        }

        void check_source_term_wells(
            std::optional< absl::Span< const std::string > > well_names )
        {
            absl::btree_set< std::string > unchecked_wells;
            for( const auto& source : problem_.source_terms() )
            {
                if( !well_names )
                {
                    unchecked_wells.insert( source.well );
                }
                else if( absl::c_find( well_names.value(), source.well )
                         == well_names->end() )
                {
                    add_issue( "Source term well ", source.well,
                        " has no perforations given to the exporter" );
                }
            }
            for( const auto& well : unchecked_wells )
            {
                add_warning( "Source term well ", well,
                    " is not checked: its perforations are given to the GEOS "
                    "exporter" );
            }
        }

        void check_source_terms()
        {
            for( const auto& source : problem_.source_terms() )
            {
                if( !geode::internal::is_geos_name(
                        geode::internal::geos_source_term_name( source ) ) )
                {
                    add_issue( "Source term well ", source.well,
                        " and component ", source.component,
                        " names must only contain letters, digits, '_' "
                        "and '-'" );
                }
            }
            if( problem_.thermal() && !problem_.source_terms().empty() )
            {
                add_issue( "GEOS source fluxes carry no energy: the GEOS "
                           "writer does not support source terms in thermal "
                           "simulations" );
            }
        }

        void check_ignored_fields()
        {
            check_ignored_material_fields();
            const auto* single_phase =
                dynamic_cast< const geode::SinglePhaseFluid* >( &fluid() );
            if( !problem_.thermal() && single_phase != nullptr
                && single_phase->thermal() )
            {
                add_warning( "SinglePhaseFluid thermal properties are ignored "
                             "by isothermal GEOS simulations" );
            }
            if( const auto* co2_brine =
                    dynamic_cast< const geode::Co2BrineFluid* >( &fluid() ) )
            {
                check_ezrokhi_coefficients( *co2_brine );
            }
        }

        void check_ignored_material_fields()
        {
            for( const auto& material : problem_.materials() )
            {
                if( !is_multiphase()
                    && ( material.relative_permeability
                         || material.capillary_pressure ) )
                {
                    add_warning( "Material ", material.name,
                        " saturation functions are ignored by single phase "
                        "GEOS simulations" );
                }
                if( !problem_.thermal() && material.thermal )
                {
                    add_warning( "Material ", material.name,
                        " thermal properties are ignored by isothermal GEOS "
                        "simulations" );
                }
            }
        }

        void check_ezrokhi_coefficients( const geode::Co2BrineFluid& co2_brine )
        {
            const auto has_density =
                !co2_brine.ezrokhi_density_coefficients().empty();
            const auto has_viscosity =
                !co2_brine.ezrokhi_viscosity_coefficients().empty();
            if( co2_brine.correlation()
                    == geode::CO2_BRINE_CORRELATION::phillips
                && ( has_density || has_viscosity ) )
            {
                add_warning( "Co2BrineFluid Ezrokhi coefficients are ignored "
                             "by the Phillips correlation" );
            }
            if( co2_brine.correlation() == geode::CO2_BRINE_CORRELATION::ezrokhi
                && !( has_density && has_viscosity ) )
            {
                add_warning( "Co2BrineFluid Ezrokhi correlation without "
                             "density or viscosity coefficients: they are "
                             "omitted from the GEOS deck" );
            }
        }

        void check_validity_ranges()
        {
            const auto pressure_range = valid_pressures();
            const auto temperature_range = valid_temperatures();
            for( const auto condition_id :
                geode::Range{ problem_.nb_initial_conditions() } )
            {
                const auto* condition =
                    dynamic_cast< const geode::UniformInitialCondition* >(
                        &problem_.initial_condition( condition_id ) );
                if( condition == nullptr )
                {
                    continue;
                }
                if( pressure_range
                    && !is_inside( condition->pressure(), *pressure_range ) )
                {
                    add_warning( "Initial condition #", condition_id,
                        " pressure is outside the fluid validity range" );
                }
                if( temperature_range && condition->temperature()
                    && !is_inside(
                        condition->temperature().value(), *temperature_range ) )
                {
                    add_warning( "Initial condition #", condition_id,
                        " temperature is outside the fluid validity range" );
                }
            }
        }

        static bool is_inside(
            double value, const std::array< double, 2 >& range )
        {
            return value >= range[0] && value <= range[1];
        }

        [[nodiscard]] std::optional< std::array< double, 2 > >
            valid_pressures() const
        {
            return geode::internal::geos_fluid_pressure_range( fluid() );
        }

        [[nodiscard]] std::optional< std::array< double, 2 > >
            valid_temperatures() const
        {
            return geode::internal::geos_fluid_temperature_range( fluid() );
        }

        void check_solver_options()
        {
            const auto& roles = geode::internal::geos_solver_option_roles();
            const auto& options = problem_.solver_options().options( "geos" );
            for( const auto& [key, value] : options )
            {
                const std::vector< std::string_view > parts =
                    absl::StrSplit( key, absl::MaxSplits( '/', 1 ) );
                if( parts.size() != 2 || parts[1].empty()
                    || !roles.contains( parts[0] ) )
                {
                    add_issue( "GEOS solver option \"", key,
                        R"(" must be "<role>/<attribute>" with a known role)" );
                    continue;
                }
                if( !geode::internal::is_xml_attribute_name( parts[1] ) )
                {
                    add_issue( "GEOS solver option \"", key, "\": \"", parts[1],
                        "\" is not an XML attribute name" );
                }
                if( parts[0] == "vtk"
                    && !geode::internal::geos_writes_vtk_output(
                        problem_.output() ) )
                {
                    add_warning( "GEOS solver option \"", key,
                        "\" is ignored: the problem has no VTK output" );
                }
            }
        }

        void check_boundary_conditions(
            const geode::internal::GeosDeckContext& context )
        {
            if( !is_multiphase() )
            {
                if( problem_.thermal()
                    && context.has_boundary_pressure
                           != context.has_boundary_temperature )
                {
                    add_issue( "GEOS thermal single phase boundary "
                               "conditions need both a boundary pressure and "
                               "a boundary temperature" );
                }
                return;
            }
            const auto has_fractions =
                !context.boundary_fraction_phases.empty();
            if( !context.has_boundary_pressure )
            {
                if( has_fractions || context.has_boundary_temperature )
                {
                    add_issue( "GEOS boundary fluid fractions and "
                               "temperatures need a boundary pressure" );
                }
                return;
            }
            check_boundary_components( context );
            if( problem_.thermal() && !context.has_boundary_temperature )
            {
                add_issue( "GEOS thermal multiphase boundary pressures need a "
                           "boundary temperature" );
            }
        }

        void check_boundary_components(
            const geode::internal::GeosDeckContext& context )
        {
            std::vector< bool > has_fraction(
                fluid().components().size(), false );
            for( const auto phase : context.boundary_fraction_phases )
            {
                const auto component = geode::phase_component( fluid(), phase );
                if( !component )
                {
                    add_issue( "Boundary ", geode::phase_name( phase ),
                        " fraction: the fluid has no ",
                        geode::phase_name( phase ),
                        " phase with its own component" );
                    continue;
                }
                has_fraction[component.value()] = true;
            }
            for( const auto component : geode::Indices{ has_fraction } )
            {
                if( !has_fraction[component] )
                {
                    add_issue( "GEOS boundary conditions need a fluid "
                               "fraction for every component: ",
                        fluid().components()[component].name, " is missing" );
                }
            }
        }

    private:
        const geode::PhysicalProblem& problem_;
        geode::SimulationReport& report_;
    };

    class GeosRegionAssigner
    {
    public:
        GeosRegionAssigner( const geode::BRep& model,
            const geode::internal::GeosDeckContext& context,
            std::string_view item_kind )
            : model_( model ), context_( context ), item_kind_( item_kind )
        {
            assignment_.region_items.resize( context_.solid_regions.size() );
        }

        void assign( geode::index_t item,
            std::string_view item_name,
            absl::Span< const geode::uuid > components )
        {
            for( const auto& component : components )
            {
                const auto blocks =
                    geode::physical_problem_region_blocks( model_, component );
                if( blocks )
                {
                    assign_blocks( item, item_name, component, *blocks );
                }
            }
        }

        // Regions without item get the default one, if any
        geode::internal::GeosRegionAssignment finish(
            std::optional< geode::index_t > default_item )
        {
            for( auto& assigned : assignment_.region_items )
            {
                if( !assigned )
                {
                    assigned = default_item;
                }
            }
            return std::move( assignment_ );
        }

    private:
        void assign_blocks( geode::index_t item,
            std::string_view item_name,
            const geode::uuid& component,
            absl::Span< const geode::uuid > blocks )
        {
            const absl::flat_hash_set< geode::uuid > block_set{ blocks.begin(),
                blocks.end() };
            absl::btree_set< geode::index_t > regions;
            for( const auto& block : blocks )
            {
                const auto region = context_.block_regions.find( block );
                if( region != context_.block_regions.end() )
                {
                    regions.insert( region->second );
                }
            }
            for( const auto region_id : regions )
            {
                const auto& region = context_.solid_regions[region_id];
                if( !absl::c_all_of(
                        region.blocks, [&block_set]( const auto& block ) {
                            return block_set.contains( block );
                        } ) )
                {
                    assignment_.issues.emplace_back( absl::StrCat( "Region ",
                        component.string(), " of ", item_kind_, " ", item_name,
                        " covers only part of the GEOS region of ",
                        region.component_name,
                        ": its Blocks are merged into one region, assign it "
                        "as a whole" ) );
                    continue;
                }
                auto& assigned = assignment_.region_items[region_id];
                if( !assigned )
                {
                    assigned = item;
                }
            }
        }

    private:
        const geode::BRep& model_;
        const geode::internal::GeosDeckContext& context_;
        std::string_view item_kind_;
        geode::internal::GeosRegionAssignment assignment_;
    };
} // namespace

namespace geode::internal
{
    const absl::flat_hash_map< std::string, std::string >&
        geos_solver_option_roles()
    {
        static const absl::flat_hash_map< std::string, std::string > roles{
            { "solvers", "/Problem/Solvers" },
            { "flow", "/Problem/Solvers/*[@name='flow']" },
            { "nonlinear",
                "/Problem/Solvers/*[@name='flow']/NonlinearSolverParameters" },
            { "linear",
                "/Problem/Solvers/*[@name='flow']/LinearSolverParameters" },
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

    GeosRegionAssignment assign_geos_regions( const BRep& model,
        const GeosDeckContext& context,
        absl::Span< const absl::Span< const uuid > > item_regions,
        absl::Span< const std::string > item_names,
        std::string_view item_kind )
    {
        GeosRegionAssigner assigner{ model, context, item_kind };
        std::optional< index_t > default_item;
        for( const auto item : Indices{ item_regions } )
        {
            if( item_regions[item].empty() )
            {
                if( !default_item )
                {
                    default_item = item;
                }
                continue;
            }
            assigner.assign( item, item_names[item], item_regions[item] );
        }
        return assigner.finish( default_item );
    }

    GeosRegionAssignment assign_geos_materials( const PhysicalProblem& problem,
        const BRep& model,
        const GeosDeckContext& context )
    {
        std::vector< absl::Span< const uuid > > material_regions;
        std::vector< std::string > material_names;
        for( const auto& material : problem.materials() )
        {
            material_regions.emplace_back( material.regions );
            material_names.push_back( material.name );
        }
        return assign_geos_regions(
            model, context, material_regions, material_names, "material" );
    }

    GeosRegionAssignment assign_geos_initial_conditions(
        const PhysicalProblem& problem,
        const BRep& model,
        const GeosDeckContext& context )
    {
        std::vector< absl::Span< const uuid > > condition_regions;
        std::vector< std::string > condition_names;
        for( const auto condition_id :
            Range{ problem.nb_initial_conditions() } )
        {
            condition_regions.emplace_back(
                problem.initial_condition( condition_id ).regions() );
            condition_names.push_back( absl::StrCat( "#", condition_id ) );
        }
        return assign_geos_regions( model, context, condition_regions,
            condition_names, "initial condition" );
    }

    std::optional< std::array< double, 2 > > geos_fluid_pressure_range(
        const FluidModel& fluid )
    {
        if( const auto* co2_brine =
                dynamic_cast< const Co2BrineFluid* >( &fluid ) )
        {
            return std::array< double, 2 >{ co2_brine->pressure_range().min,
                co2_brine->pressure_range().max };
        }
        const auto* dead_oil = dynamic_cast< const DeadOilFluid* >( &fluid );
        if( dead_oil == nullptr || dead_oil->hydrocarbon_pvt().empty() )
        {
            return std::nullopt;
        }
        std::array< double, 2 > range{ std::numeric_limits< double >::lowest(),
            std::numeric_limits< double >::max() };
        for( const auto& table : dead_oil->hydrocarbon_pvt() )
        {
            if( table.pressures.empty() )
            {
                return std::nullopt;
            }
            range[0] = std::max( range[0], table.pressures.front() );
            range[1] = std::min( range[1], table.pressures.back() );
        }
        return range;
    }

    std::optional< std::array< double, 2 > > geos_fluid_temperature_range(
        const FluidModel& fluid )
    {
        if( const auto* co2_brine =
                dynamic_cast< const Co2BrineFluid* >( &fluid ) )
        {
            return std::array< double, 2 >{ co2_brine->temperature_range().min,
                co2_brine->temperature_range().max };
        }
        return std::nullopt;
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
        const PhysicalProblem& problem, const GeosDeckContext& context )
    {
        if( !context.has_boundary_pressure || context.has_boundary_temperature
            || problem.thermal() || !is_geos_multiphase( problem ) )
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

    void validate_geos_problem( const PhysicalProblem& problem,
        std::optional< absl::Span< const std::string > > well_names,
        SimulationReport& report )
    {
        GeosValidator{ problem, report }.check_problem( well_names );
    }

    void validate_geos_deck( const PhysicalProblem& problem,
        const GeosDeckContext& context,
        SimulationReport& report )
    {
        GeosValidator{ problem, report }.check_deck( context );
    }
} // namespace geode::internal
