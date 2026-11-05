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
#include <cctype>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <absl/algorithm/container.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>
#include <absl/types/span.h>

#include <geode/basic/range.hpp>

#include <geode/mesh/core/regular_grid_solid.hpp>
#include <geode/mesh/core/solid_mesh.hpp>

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
#include <geode/geosciences/simulation/helpers/internal/simulation_checker_helpers.hpp>

#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_rules.hpp>

namespace
{
    using geode::PHYSICAL_PROPERTY_NAME;
    using geode::SIMULATION_CHECK_STAGE;

    constexpr std::array< PHYSICAL_PROPERTY_NAME, 3 > BOUNDARY_FRACTIONS{
        PHYSICAL_PROPERTY_NAME::boundary_oil_fraction,
        PHYSICAL_PROPERTY_NAME::boundary_gas_fraction,
        PHYSICAL_PROPERTY_NAME::boundary_water_fraction
    };

    template < typename... Types, typename Base >
    bool is_one_of( const Base& object )
    {
        return (
            ... || ( dynamic_cast< const Types* >( &object ) != nullptr ) );
    }

    std::string fraction_property_name( geode::PHASE phase )
    {
        return absl::StrCat(
            "boundary_", geode::phase_name( phase ), "_fraction" );
    }

    bool is_geos_name( std::string_view name )
    {
        return !name.empty()
               && std::all_of( name.begin(), name.end(), []( char character ) {
                      return std::isalnum(
                                 static_cast< unsigned char >( character ) )
                                 != 0
                             || character == '_' || character == '-';
                  } );
    }

    class GeosProblemChecker
    {
    public:
        GeosProblemChecker( const geode::PhysicalProblem& problem,
            geode::SimulationReport& report )
            : problem_( problem ), report_( report )
        {
        }

        void check_problem(
            std::optional< absl::Span< const std::string > > well_names )
        {
            check_solver_options();
            if( !has_fluid() )
            {
                return;
            }
            check_supported_flow();
            check_supported_materials();
            check_supported_initial_conditions();
            check_thermal();
            check_source_terms( well_names );
        }

        void check_boundary_conditions(
            const geode::internal::GeosBoundaryConditions& boundary_conditions )
        {
            if( !has_fluid() )
            {
                return;
            }
            if( !is_multiphase() )
            {
                check_single_phase_boundary_conditions( boundary_conditions );
                return;
            }
            const auto has_fractions =
                !boundary_conditions.fraction_phases.empty();
            if( !boundary_conditions.has_pressure )
            {
                if( has_fractions || boundary_conditions.has_temperature )
                {
                    add_boundary_issue( "boundary_pressure",
                        "GEOS boundary fluid fractions and temperatures need "
                        "a boundary pressure" );
                }
                return;
            }
            check_boundary_components( boundary_conditions );
            if( problem_.thermal() && !boundary_conditions.has_temperature )
            {
                add_boundary_issue( "boundary_pressure",
                    "GEOS thermal multiphase boundary pressures need a "
                    "boundary temperature" );
            }
        }

    private:
        void check_single_phase_boundary_conditions(
            const geode::internal::GeosBoundaryConditions& boundary_conditions )
        {
            if( problem_.thermal() && boundary_conditions.has_pressure
                && !boundary_conditions.has_temperature )
            {
                add_boundary_issue( "boundary_pressure",
                    "GEOS thermal boundary pressures need a boundary "
                    "temperature" );
            }
        }

        template < typename... Args >
        void add_issue( std::string subject, const Args&... message )
        {
            report_.add_error( geode::SIMULATION_CHECK_STAGE::simulator,
                std::move( subject ), message... );
        }

        template < typename... Args >
        void add_boundary_issue( std::string subject, const Args&... message )
        {
            report_.add_error( geode::SIMULATION_CHECK_STAGE::properties,
                std::move( subject ), message... );
        }

        template < typename Object >
        void add_unsupported( std::string subject, const Object& object )
        {
            add_issue( std::move( subject ), object.type_name(),
                " is not supported by the GEOS writer" );
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
            return geode::internal::is_multiphase( problem_ );
        }

        void check_supported_flow()
        {
            if( !is_one_of< geode::SinglePhaseFlow, geode::MultiphaseFlow >(
                    problem_.flow() ) )
            {
                add_unsupported( "Flow model", problem_.flow() );
            }
            if( !is_one_of< geode::SinglePhaseFluid, geode::DeadOilFluid,
                    geode::Co2BrineFluid >( fluid() ) )
            {
                add_unsupported( "Fluid", fluid() );
            }
        }

        void check_supported_materials()
        {
            for( const auto& material : problem_.materials() )
            {
                const auto subject =
                    geode::internal::material_subject( material.name );
                if( !is_geos_name( material.name ) )
                {
                    add_issue( subject, "name must only contain letters, "
                                        "digits, '_' and '-'" );
                }
                if( material.relative_permeability )
                {
                    check_supported_relperm(
                        *material.relative_permeability, subject );
                }
                if( material.capillary_pressure )
                {
                    check_supported_capillary(
                        *material.capillary_pressure, subject );
                }
            }
        }

        void check_supported_relperm(
            const geode::RelativePermeability& relperm,
            const std::string& subject )
        {
            if( const auto* table =
                    dynamic_cast< const geode::TableRelativePermeability* >(
                        &relperm ) )
            {
                if( table->tables().size() != 1 )
                {
                    add_issue( subject, "the GEOS writer only supports "
                                        "two-phase relative permeability "
                                        "tables" );
                }
                return;
            }
            if( !is_one_of< geode::BrooksCoreyRelativePermeability >(
                    relperm ) )
            {
                add_unsupported( subject, relperm );
            }
        }

        void check_supported_capillary(
            const geode::CapillaryPressure& capillary,
            const std::string& subject )
        {
            if( const auto* table =
                    dynamic_cast< const geode::TableCapillaryPressure* >(
                        &capillary ) )
            {
                if( table->tables().size() != 1 )
                {
                    add_issue( subject, "the GEOS writer only supports "
                                        "two-phase capillary pressure tables" );
                }
                return;
            }
            if( !is_one_of< geode::BrooksCoreyCapillaryPressure >( capillary ) )
            {
                add_unsupported( subject, capillary );
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
                    add_unsupported( geode::internal::initial_condition_subject(
                                         condition_id ),
                        condition );
                }
            }
        }

        void check_thermal()
        {
            if( problem_.thermal() )
            {
                if( is_one_of< geode::DeadOilFluid >( fluid() ) )
                {
                    add_issue( "Fluid", "the GEOS writer does not support "
                                        "thermal simulations with "
                                        "DeadOilFluid" );
                }
                return;
            }
            if( is_multiphase()
                && !geode::internal::geos_isothermal_temperature( problem_ ) )
            {
                add_issue( "Initial conditions",
                    "GEOS isothermal multiphase simulations need the same "
                    "temperature in every initial condition" );
            }
        }

        void check_source_terms(
            std::optional< absl::Span< const std::string > > well_names )
        {
            for( const auto& source : problem_.source_terms() )
            {
                const auto subject =
                    absl::StrCat( "Source term on well \"", source.well, "\"" );
                if( well_names
                    && absl::c_find( well_names.value(), source.well )
                           == well_names->end() )
                {
                    geode::SimulationIssue issue{
                        geode::SIMULATION_ISSUE_SEVERITY::error,
                        geode::SIMULATION_CHECK_STAGE::simulator, subject,
                        "the well has no perforations given to the exporter"
                    };
                    issue.fix = "add the well perforations to the exporter";
                    report_.add_issue( std::move( issue ) );
                }
                if( !is_geos_name(
                        geode::internal::geos_source_term_name( source ) ) )
                {
                    add_issue( subject, "well and component ", source.component,
                        " names must only contain letters, digits, '_' and "
                        "'-'" );
                }
            }
            if( problem_.thermal() && !problem_.source_terms().empty() )
            {
                add_issue( "Source terms",
                    "GEOS source fluxes carry no energy: the GEOS writer does "
                    "not support source terms in thermal simulations" );
            }
        }

        void check_solver_options()
        {
            const auto& roles = geode::internal::geos_solver_option_roles();
            for( const auto& [key, value] :
                problem_.solver_options().options( "geos" ) )
            {
                const std::vector< std::string_view > parts =
                    absl::StrSplit( key, absl::MaxSplits( '/', 1 ) );
                if( parts.size() != 2 || parts[1].empty()
                    || !roles.contains( parts[0] ) )
                {
                    add_issue( "Solver options", "GEOS solver option \"", key,
                        R"(" must be "<role>/<attribute>" with a known role)" );
                }
            }
        }

        void check_boundary_components(
            const geode::internal::GeosBoundaryConditions& boundary_conditions )
        {
            std::vector< bool > has_fraction(
                fluid().components().size(), false );
            for( const auto phase : boundary_conditions.fraction_phases )
            {
                const auto component = geode::phase_component( fluid(), phase );
                if( !component )
                {
                    add_boundary_issue( fraction_property_name( phase ),
                        "the fluid has no ", geode::phase_name( phase ),
                        " phase with its own component" );
                    continue;
                }
                has_fraction[component.value()] = true;
            }
            for( const auto component : geode::Indices{ has_fraction } )
            {
                if( !has_fraction[component] )
                {
                    report_missing_fraction( component );
                }
            }
        }

        // The component of each phase has the phase index, when every phase
        // has its own component
        void report_missing_fraction( geode::index_t component )
        {
            const auto& name = fluid().components()[component].name;
            const auto phases = fluid().phases();
            if( phases.size() != fluid().components().size() )
            {
                add_boundary_issue( "Boundary fluid fractions",
                    "GEOS boundary conditions need a fluid fraction for every "
                    "component: ",
                    name, " has none" );
                return;
            }
            const auto property = fraction_property_name( phases[component] );
            geode::SimulationIssue issue{
                geode::SIMULATION_ISSUE_SEVERITY::error,
                geode::SIMULATION_CHECK_STAGE::properties, property,
                absl::StrCat( "missing: GEOS boundary conditions need a fluid "
                              "fraction for every component, ",
                    name, " has none" )
            };
            issue.fix = absl::StrCat( "add a ", property,
                " physical property on the boundary Surfaces" );
            report_.add_issue( std::move( issue ) );
        }

    private:
        const geode::PhysicalProblem& problem_;
        geode::SimulationReport& report_;
    };

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

    void add_issues( std::vector< geode::SimulationIssue >&& issues,
        geode::SimulationReport& report )
    {
        for( auto& issue : issues )
        {
            report.add_issue( std::move( issue ) );
        }
    }
} // namespace

namespace geode::internal
{
    GeosSimulationChecker::GeosSimulationChecker(
        const PhysicalProblem& problem, const BRep& model )
        : SimulationChecker( problem, model )
    {
    }

    GeosSimulationChecker::GeosSimulationChecker(
        const PhysicalProblem& problem,
        const BRep& model,
        GeosExportData export_data )
        : SimulationChecker( problem, model ),
          export_data_( std::move( export_data ) )
    {
    }

    void GeosSimulationChecker::check_simulator()
    {
        const auto& problem = this->problem();
        const auto& model = this->model();
        report().set_simulator( "GEOS" );
        std::optional< absl::Span< const std::string > > well_names;
        if( export_data_ )
        {
            well_names = export_data_->well_names;
        }
        GeosProblemChecker problem_checker{ problem, report() };
        problem_checker.check_problem( well_names );
        check_blocks();
        check_unique_vertices_in_blocks();
        check_regular_grids( model, report() );
        if( !problem.has_flow() || !problem.flow().has_fluid() )
        {
            return;
        }
        check_rock_properties();
        check_boundary_properties();
        problem_checker.check_boundary_conditions(
            geos_boundary_conditions( model ) );
        if( export_data_ )
        {
            add_issues(
                assign_geos_materials( problem, model, export_data_->regions )
                    .issues,
                report() );
            add_issues( assign_geos_initial_conditions(
                            problem, model, export_data_->regions )
                            .issues,
                report() );
        }
    }

    // GEOS imports a model porosity or permeability on every cell, merged
    // from the Block meshes. Without such a property, the deck uses the
    // value of each material.
    void GeosSimulationChecker::check_rock_properties()
    {
        std::vector< PhysicalPropertyCoverage > coverages;
        for( const auto property : { PHYSICAL_PROPERTY_NAME::porosity,
                 PHYSICAL_PROPERTY_NAME::permeability } )
        {
            PhysicalPropertyRequirement requirement{ property };
            requirement.merged = true;
            coverages.push_back( check_physical_property( requirement ) );
        }
        check_rock_property_availability( coverages );
    }

    void GeosSimulationChecker::check_boundary_properties()
    {
        const auto pressure =
            check_physical_property( PhysicalPropertyRequirement{
                PHYSICAL_PROPERTY_NAME::boundary_pressure } );
        check_physical_property_presence(
            pressure, "nothing drives the flow in the GEOS deck" );
        check_physical_property( PhysicalPropertyRequirement{
            PHYSICAL_PROPERTY_NAME::boundary_temperature } );
        if( !is_multiphase( problem() ) )
        {
            return;
        }
        for( const auto fraction : BOUNDARY_FRACTIONS )
        {
            check_physical_property( PhysicalPropertyRequirement{ fraction } );
        }
    }
} // namespace geode::internal
