#include <geode/geosciences_io/model/internal/geos_simulation_validation.hpp>

#include <algorithm>
#include <optional>

#include <absl/algorithm/container.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>

#include <geode/basic/range.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/solver_options.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>

#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>

namespace
{
    template < typename... Types, typename Base >
    bool is_one_of( const Base& object )
    {
        return (
            ... || ( dynamic_cast< const Types* >( &object ) != nullptr ) );
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

    bool is_geos_multiphase( const geode::PhysicalProblem& problem )
    {
        return dynamic_cast< const geode::MultiphaseFlow* >( &problem.flow() )
               != nullptr;
    }

    std::vector< std::vector< geode::uuid > > material_regions(
        const geode::PhysicalProblem& problem )
    {
        std::vector< std::vector< geode::uuid > > regions;
        for( const auto& material : problem.materials() )
        {
            regions.push_back( material.regions );
        }
        return regions;
    }

    std::vector< std::vector< geode::uuid > > initial_condition_regions(
        const geode::PhysicalProblem& problem )
    {
        std::vector< std::vector< geode::uuid > > regions;
        for( const auto condition_id :
            geode::Range{ problem.nb_initial_conditions() } )
        {
            regions.push_back(
                problem.initial_condition( condition_id ).regions );
        }
        return regions;
    }

    class GeosValidator
    {
    public:
        GeosValidator( const geode::PhysicalProblem& problem,
            const geode::internal::GeosDeckContext& context )
            : problem_( problem ), context_( context )
        {
        }

        std::vector< std::string > run()
        {
            if( !problem_.has_flow() || !problem_.flow().fluid )
            {
                return std::move( issues_ );
            }
            check_supported_flow();
            check_supported_materials();
            check_supported_initial_conditions();
            check_thermal();
            check_source_terms();
            check_region_assignments();
            check_rock_properties();
            check_solver_options();
            check_boundary_conditions();
            return std::move( issues_ );
        }

    private:
        template < typename... Args >
        void add_issue( const Args&... args )
        {
            issues_.emplace_back( absl::StrCat( args... ) );
        }

        template < typename Object >
        void add_unsupported( const Object& object )
        {
            add_issue(
                object.type_name(), " is not supported by the GEOS writer" );
        }

        [[nodiscard]] const geode::FluidModel& fluid() const
        {
            return *problem_.flow().fluid;
        }

        [[nodiscard]] bool is_multiphase() const
        {
            return is_geos_multiphase( problem_ );
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
        }

        void check_supported_materials()
        {
            for( const auto& material : problem_.materials() )
            {
                if( !is_geos_name( material.name ) )
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
                if( table->tables.size() != 1 )
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
                if( table->tables.size() != 1 )
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

        void check_source_terms()
        {
            for( const auto& source : problem_.source_terms() )
            {
                if( absl::c_find( context_.well_names, source.well )
                    == context_.well_names.end() )
                {
                    add_issue( "Source term well ", source.well,
                        " has no perforations given to the exporter" );
                }
                if( !is_geos_name(
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

        void check_region_assignments()
        {
            check_assignment( material_regions( problem_ ), "material" );
            check_assignment(
                initial_condition_regions( problem_ ), "initial condition" );
        }

        void check_assignment(
            const std::vector< std::vector< geode::uuid > >& regions,
            std::string_view kind )
        {
            std::vector< absl::Span< const geode::uuid > > item_regions{
                regions.begin(), regions.end()
            };
            const auto assignment = geode::internal::assign_geos_regions(
                context_, item_regions, kind );
            issues_.insert( issues_.end(), assignment.issues.begin(),
                assignment.issues.end() );
        }

        void check_rock_properties()
        {
            for( const auto& material : problem_.materials() )
            {
                if( !context_.porosity_attribute && !material.porosity )
                {
                    add_issue( "Material ", material.name,
                        " needs a porosity: the model has no porosity "
                        "property" );
                }
                if( !context_.permeability_attribute && !material.permeability )
                {
                    add_issue( "Material ", material.name,
                        " needs a permeability: the model has no "
                        "permeability property" );
                }
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
                    add_issue( "GEOS solver option \"", key,
                        R"(" must be "<role>/<attribute>" with a known role)" );
                }
            }
        }

        void check_boundary_conditions()
        {
            const auto has_fractions =
                !context_.boundary_fraction_phases.empty();
            if( !is_multiphase() )
            {
                if( has_fractions )
                {
                    add_issue(
                        "Boundary fluid fractions need a multiphase flow" );
                }
                return;
            }
            if( !context_.has_boundary_pressure )
            {
                if( has_fractions || context_.has_boundary_temperature )
                {
                    add_issue( "GEOS boundary fluid fractions and "
                               "temperatures need a boundary pressure" );
                }
                return;
            }
            check_boundary_components();
            if( problem_.thermal() && !context_.has_boundary_temperature )
            {
                add_issue( "GEOS thermal multiphase boundary pressures need a "
                           "boundary temperature" );
            }
        }

        void check_boundary_components()
        {
            std::vector< bool > has_fraction(
                fluid().components.size(), false );
            for( const auto phase : context_.boundary_fraction_phases )
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
                        fluid().components[component].name, " is missing" );
                }
            }
        }

    private:
        const geode::PhysicalProblem& problem_;
        const geode::internal::GeosDeckContext& context_;
        std::vector< std::string > issues_;
    };
    class GeosRegionAssigner
    {
    public:
        GeosRegionAssigner( const geode::internal::GeosDeckContext& context,
            std::string_view item_kind )
            : context_( context ), item_kind_( item_kind )
        {
            assignment_.region_items.resize( context_.solid_regions.size() );
        }

        void assign(
            geode::index_t item, absl::Span< const geode::uuid > components )
        {
            for( const auto& component : components )
            {
                const auto regions = context_.regions_by_uuid.find( component );
                if( regions == context_.regions_by_uuid.end() )
                {
                    add_issue( "Region ", component.string(), " of ",
                        item_kind_, " ", item, " is not a model region" );
                    continue;
                }
                assign_regions( item, regions->second );
            }
        }

        // Regions without item get the default one, if any
        geode::internal::GeosRegionAssignment finish(
            std::optional< geode::index_t > default_item )
        {
            for( const auto region :
                geode::Indices{ assignment_.region_items } )
            {
                auto& assigned = assignment_.region_items[region];
                if( !assigned )
                {
                    assigned = default_item;
                }
                if( !assigned )
                {
                    add_issue( "Region ",
                        context_.solid_regions[region].component_name,
                        " has no ", item_kind_ );
                }
            }
            return std::move( assignment_ );
        }

    private:
        template < typename... Args >
        void add_issue( const Args&... args )
        {
            assignment_.issues.emplace_back( absl::StrCat( args... ) );
        }

        void assign_regions(
            geode::index_t item, absl::Span< const geode::index_t > regions )
        {
            for( const auto region : regions )
            {
                auto& assigned = assignment_.region_items[region];
                if( assigned && assigned.value() != item )
                {
                    add_issue( "Region ",
                        context_.solid_regions[region].component_name,
                        " has several ", item_kind_, "s" );
                }
                assigned = item;
            }
        }

    private:
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

    GeosRegionAssignment assign_geos_regions( const GeosDeckContext& context,
        absl::Span< const absl::Span< const uuid > > item_regions,
        std::string_view item_kind )
    {
        GeosRegionAssigner assigner{ context, item_kind };
        std::optional< index_t > default_item;
        for( const auto item : Indices{ item_regions } )
        {
            if( item_regions[item].empty() )
            {
                default_item = item;
                continue;
            }
            assigner.assign( item, item_regions[item] );
        }
        return assigner.finish( default_item );
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
            if( condition == nullptr || !condition->temperature
                || ( temperature
                     && temperature.value()
                            != condition->temperature.value() ) )
            {
                return std::nullopt;
            }
            temperature = condition->temperature;
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

    std::vector< std::string > validate_geos_simulation(
        const PhysicalProblem& problem, const GeosDeckContext& context )
    {
        return GeosValidator{ problem, context }.run();
    }
} // namespace geode::internal
