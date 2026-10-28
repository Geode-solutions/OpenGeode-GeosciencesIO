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

#include <geode/geosciences_io/model/internal/geos_simulation_writer.hpp>

#include <array>
#include <optional>
#include <string>
#include <vector>

#include <pugixml.hpp>

#include <absl/algorithm/container.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>

#include <geode/basic/logger.hpp>
#include <geode/basic/pimpl_impl.hpp>
#include <geode/basic/range.hpp>

#include <geode/geometry/vector.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/schedule.hpp>
#include <geode/geosciences/simulation/core/solver_options.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>

#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_validation.hpp>
#include <geode/geosciences_io/model/internal/geos_xml.hpp>

namespace
{
    constexpr auto FLUID_MODEL_NAME = "fluid";
    constexpr auto FLOW_SOLVER_NAME = "flow";
    constexpr auto TPFA_NAME = "tpfa";
    constexpr auto NULL_SOLID_NAME = "nullSolid";
    constexpr auto VTK_OUTPUT_NAME = "vtkOutput";
    constexpr auto RESTART_OUTPUT_NAME = "restartOutput";

    std::vector< std::string > geos_phase_names(
        absl::Span< const geode::PHASE > phases )
    {
        std::vector< std::string > names;
        names.reserve( phases.size() );
        for( const auto phase : phases )
        {
            names.emplace_back( geode::phase_name( phase ) );
        }
        return names;
    }

    // Per-phase values given in the law phase order, reordered in the fluid
    // phase order expected by GEOS
    std::vector< double > in_fluid_phase_order(
        absl::Span< const geode::PHASE > law_phases,
        absl::Span< const double > values,
        absl::Span< const geode::PHASE > fluid_phases )
    {
        std::vector< double > ordered;
        ordered.reserve( fluid_phases.size() );
        for( const auto phase : fluid_phases )
        {
            const auto* const law_phase = absl::c_find( law_phases, phase );
            ordered.push_back( values[static_cast< size_t >(
                law_phase - law_phases.begin() )] );
        }
        return ordered;
    }

    std::vector< std::string > component_names( const geode::FluidModel& fluid )
    {
        std::vector< std::string > names;
        for( const auto& component : fluid.components() )
        {
            names.push_back( component.name );
        }
        return names;
    }

    std::vector< double > component_molar_weights(
        const geode::FluidModel& fluid )
    {
        std::vector< double > weights;
        for( const auto& component : fluid.components() )
        {
            weights.push_back( component.molar_weight );
        }
        return weights;
    }

    std::string pvt_table_name( geode::PHASE phase, std::string_view quantity )
    {
        return absl::StrCat(
            "pvt_", geode::phase_name( phase ), "_", quantity );
    }

    pugi::xml_node add_named_child( pugi::xml_node& parent,
        std::string_view element,
        std::string_view name )
    {
        auto node = parent.append_child( std::string{ element }.c_str() );
        geode::internal::set_attribute( node, "name", name );
        return node;
    }

    pugi::xml_node add_field_specification( pugi::xml_node& parent,
        std::string_view name,
        std::string_view region_name,
        std::string_view field_name )
    {
        auto node = add_named_child( parent, "FieldSpecification", name );
        geode::internal::set_attribute( node, "initialCondition", "1" );
        geode::internal::set_attribute( node, "setNames", "{ all }" );
        geode::internal::set_attribute( node, "objectPath",
            absl::StrCat( "ElementRegions/", region_name ) );
        geode::internal::set_attribute( node, "fieldName", field_name );
        return node;
    }
} // namespace

namespace geode::internal
{
    struct GeosTableFunction
    {
        std::string name;
        std::vector< double > coordinates;
        std::vector< double > values;
    };

    class GeosSimulationWriter::Impl
    {
    public:
        Impl( const PhysicalProblem& problem,
            const GeosDeckContext& context,
            const GeosRegionAssignment& materials,
            const GeosRegionAssignment& initial_conditions )
            : problem_( problem ),
              context_( context ),
              material_of_region_( assigned_items( materials ) ),
              condition_of_region_( assigned_items( initial_conditions ) )
        {
            build_porosity_groups();
            build_permeability_groups();
            build_imported_fields();
            build_tables();
        }

        [[nodiscard]] absl::Span< const GeosImportedField >
            imported_fields() const
        {
            return imported_fields_;
        }

        void write_physics( pugi::xml_node& problem_node ) const
        {
            auto solvers = problem_node.append_child( "Solvers" );
            const auto& gravity = problem_.gravity();
            set_attribute( solvers, "gravityVector",
                geos_array( { gravity.value( 0 ), gravity.value( 1 ),
                    gravity.value( 2 ) } ) );
            write_flow_solver( solvers );
            auto finite_volume = problem_node.append_child( "NumericalMethods" )
                                     .append_child( "FiniteVolume" );
            add_named_child(
                finite_volume, "TwoPointFluxApproximation", TPFA_NAME );
            write_element_regions( problem_node );
            write_constitutive( problem_node );
        }

        void write_initial_conditions(
            pugi::xml_node& field_specifications ) const
        {
            for( const auto& region : context_.solid_regions )
            {
                const auto* condition =
                    dynamic_cast< const UniformInitialCondition* >(
                        &problem_.initial_condition(
                            condition_of_region_[region.id] ) );
                write_uniform_condition(
                    field_specifications, *condition, region.name );
            }
        }

        void write_source_terms( pugi::xml_node& field_specifications ) const
        {
            for( const auto& source : problem_.source_terms() )
            {
                auto node = add_named_child( field_specifications, "SourceFlux",
                    geos_source_term_name( source ) );
                set_attribute( node, "objectPath", "ElementRegions" );
                set_attribute(
                    node, "setNames", absl::StrCat( "{ ", source.well, " }" ) );
                if( is_multiphase() )
                {
                    set_attribute( node, "component",
                        absl::StrCat( component_index( source.component ) ) );
                }
                // GEOS source fluxes are positive when producing
                set_real_attribute( node, "scale", -source.rate );
            }
        }

        [[nodiscard]] bool has_functions() const
        {
            return !tables_.empty();
        }

        void write_functions( pugi::xml_node& functions ) const
        {
            for( const auto& table : tables_ )
            {
                auto node =
                    add_named_child( functions, "TableFunction", table.name );
                set_attribute(
                    node, "coordinates", geos_array( table.coordinates ) );
                set_attribute( node, "values", geos_array( table.values ) );
            }
        }

        void write_events_and_outputs( pugi::xml_node& problem_node ) const
        {
            const auto& schedule = problem_.schedule();
            auto events = problem_node.append_child( "Events" );
            set_real_attribute( events, "minTime", schedule.start_time );
            set_real_attribute( events, "maxTime", schedule.end_time );
            // GEOS runs the events of a cycle in order: outputs must come
            // before the solver to write the state at their own time
            write_output_events( events );
            auto solver_event = add_named_child(
                events, "PeriodicEvent", "solverApplications" );
            set_attribute( solver_event, "target",
                absl::StrCat( "/Solvers/", FLOW_SOLVER_NAME ) );
            if( schedule.max_time_step )
            {
                set_real_attribute( solver_event, "maxEventDt",
                    schedule.max_time_step.value() );
            }
            write_outputs( problem_node );
        }

        void apply_solver_options( pugi::xml_node& problem_node ) const
        {
            const auto& roles = geos_solver_option_roles();
            for( const auto& [key, value] :
                problem_.solver_options().options( "geos" ) )
            {
                const std::vector< std::string > parts =
                    absl::StrSplit( key, absl::MaxSplits( '/', 1 ) );
                const auto& xpath = roles.at( parts[0] );
                // Options without a deck node (vtk without VTK output) are
                // reported by the GEOS simulation check
                for( const auto& selected :
                    problem_node.select_nodes( xpath.c_str() ) )
                {
                    auto node = selected.node();
                    set_attribute( node, parts[1], value );
                }
            }
        }

    private:
        [[nodiscard]] const FluidModel& fluid() const
        {
            return problem_.fluid();
        }

        [[nodiscard]] index_t component_index(
            std::string_view component_name ) const
        {
            const auto& components = fluid().components();
            return static_cast< index_t >( std::distance( components.begin(),
                absl::c_find_if( components,
                    [component_name]( const FluidComponent& component ) {
                        return component.name == component_name;
                    } ) ) );
        }

        [[nodiscard]] bool is_multiphase() const
        {
            return dynamic_cast< const MultiphaseFlow* >( &problem_.flow() )
                   != nullptr;
        }

        [[nodiscard]] bool is_thermal() const
        {
            return problem_.thermal().has_value();
        }

        [[nodiscard]] static std::vector< index_t > assigned_items(
            const GeosRegionAssignment& assignment )
        {
            std::vector< index_t > items;
            items.reserve( assignment.region_items.size() );
            for( const auto& item : assignment.region_items )
            {
                items.push_back( item.value() );
            }
            return items;
        }

        void build_porosity_groups()
        {
            auto groups = geos_porosity_groups( problem_, context_ );
            porosity_groups_ = std::move( groups.groups );
            material_porosity_group_ = std::move( groups.material_groups );
        }

        void build_permeability_groups()
        {
            for( const auto& material : problem_.materials() )
            {
                const auto components =
                    context_.permeability_array
                        ? material.permeability.value_or(
                            context_.mean_permeability.value_or(
                                std::array< double, 3 >{ 0, 0, 0 } ) )
                        : material.permeability.value();
                material_permeability_group_.push_back(
                    permeability_group( components ) );
            }
        }

        index_t permeability_group( const std::array< double, 3 >& components )
        {
            if( context_.permeability_array && !permeability_groups_.empty() )
            {
                return 0;
            }
            for( const auto group_id : Indices{ permeability_groups_ } )
            {
                if( permeability_groups_[group_id] == components )
                {
                    return group_id;
                }
            }
            permeability_groups_.push_back( components );
            return static_cast< index_t >( permeability_groups_.size() - 1 );
        }

        void build_imported_fields()
        {
            if( context_.permeability_array )
            {
                const auto& permeability = context_.permeability_array.value();
                imported_fields_.push_back( { permeability, permeability.name,
                    "permeability_0_permeability" } );
            }
            if( !context_.porosity_array )
            {
                return;
            }
            // Each porosity model imports its own copy of the porosity array
            const auto& porosity = context_.porosity_array.value();
            for( const auto group_id : Indices{ porosity_groups_ } )
            {
                imported_fields_.push_back( { porosity,
                    geos_porosity_vtu_name( porosity.name, group_id,
                        static_cast< index_t >( porosity_groups_.size() ) ),
                    absl::StrCat(
                        "porosity_", group_id, "_referencePorosity" ) } );
            }
        }

        void build_tables()
        {
            if( const auto* dead_oil =
                    dynamic_cast< const DeadOilFluid* >( &fluid() ) )
            {
                build_dead_oil_tables( *dead_oil );
            }
            if( !is_multiphase() )
            {
                return;
            }
            for( const auto& material : problem_.materials() )
            {
                build_saturation_tables( material );
            }
        }

        void build_dead_oil_tables( const DeadOilFluid& fluid_model )
        {
            for( const auto phase : fluid_model.phases() )
            {
                for( const auto& table : fluid_model.hydrocarbon_pvt() )
                {
                    if( table.phase != phase )
                    {
                        continue;
                    }
                    tables_.push_back( { pvt_table_name( phase, "fvf" ),
                        table.pressures, table.formation_volume_factors } );
                    tables_.push_back( { pvt_table_name( phase, "visc" ),
                        table.pressures, table.viscosities } );
                }
            }
        }

        void build_saturation_tables( const Material& material )
        {
            if( const auto* relperm =
                    dynamic_cast< const TableRelativePermeability* >(
                        material.relative_permeability.get() ) )
            {
                const auto tables = relperm->tables();
                const auto& table = tables.front();
                tables_.push_back( { absl::StrCat( "relperm_", material.name,
                                         "_wetting" ),
                    table.wetting_kr.saturations, table.wetting_kr.values } );
                tables_.push_back(
                    { absl::StrCat( "relperm_", material.name, "_non_wetting" ),
                        table.non_wetting_kr.saturations,
                        table.non_wetting_kr.values } );
            }
            if( !material.capillary_pressure )
            {
                return;
            }
            if( const auto* capillary =
                    dynamic_cast< const TableCapillaryPressure* >(
                        material.capillary_pressure.get() ) )
            {
                const auto tables = capillary->tables();
                const auto& table = tables.front();
                tables_.push_back(
                    { absl::StrCat( "cappres_", material.name, "_pc" ),
                        table.pressure.saturations, table.pressure.values } );
            }
        }

        [[nodiscard]] std::vector< std::string > region_names() const
        {
            std::vector< std::string > names;
            for( const auto& region : context_.solid_regions )
            {
                names.push_back( region.name );
            }
            return names;
        }

        void write_flow_solver( pugi::xml_node& solvers ) const
        {
            auto solver = add_named_child( solvers,
                is_multiphase() ? "CompositionalMultiphaseFVM"
                                : "SinglePhaseFVM",
                FLOW_SOLVER_NAME );
            set_attribute( solver, "logLevel", "1" );
            set_attribute( solver, "discretization", TPFA_NAME );
            set_attribute(
                solver, "targetRegions", geos_names( region_names() ) );
            const auto temperature =
                is_thermal() ? problem_.thermal()->reference_temperature
                             : geos_isothermal_temperature( problem_ );
            if( temperature )
            {
                set_real_attribute(
                    solver, "temperature", temperature.value() );
            }
            if( const auto* multiphase =
                    dynamic_cast< const MultiphaseFlow* >( &problem_.flow() ) )
            {
                set_attribute( solver, "useMass",
                    multiphase->basis() == COMPOSITION_BASIS::mass ? "1"
                                                                   : "0" );
            }
            if( is_thermal() )
            {
                set_attribute( solver, "isThermal", "1" );
            }
            set_real_attribute(
                solver, "initialDt", problem_.schedule().initial_time_step );
            write_solver_parameters( solver );
        }

        void write_solver_parameters( pugi::xml_node& solver ) const
        {
            auto nonlinear = solver.append_child( "NonlinearSolverParameters" );
            set_attribute( nonlinear, "newtonTol", "1e-6" );
            set_attribute( nonlinear, "newtonMaxIter", "15" );
            set_attribute( nonlinear, "maxTimeStepCuts", "10" );
            auto linear = solver.append_child( "LinearSolverParameters" );
            if( is_multiphase() )
            {
                set_attribute( linear, "solverType", "fgmres" );
                set_attribute( linear, "preconditionerType", "mgr" );
                set_attribute( linear, "krylovTol", "1e-5" );
                return;
            }
            set_attribute( linear, "solverType", "gmres" );
            set_attribute( linear, "preconditionerType", "amg" );
        }

        [[nodiscard]] std::vector< std::string > material_list(
            const Material& material ) const
        {
            std::vector< std::string > materials{ FLUID_MODEL_NAME,
                absl::StrCat( "rock_", material.name ) };
            if( is_multiphase() )
            {
                materials.push_back(
                    absl::StrCat( "relperm_", material.name ) );
                if( material.capillary_pressure )
                {
                    materials.push_back(
                        absl::StrCat( "cappres_", material.name ) );
                }
            }
            if( is_thermal() )
            {
                materials.push_back(
                    absl::StrCat( "thermalCond_", material.name ) );
            }
            return materials;
        }

        void write_element_regions( pugi::xml_node& problem_node ) const
        {
            auto regions = problem_node.append_child( "ElementRegions" );
            for( const auto& region : context_.solid_regions )
            {
                add_comment( regions,
                    absl::StrCat( region.name, ": ", region.component_name,
                        " (", region.component_id.string(), ")" ) );
                auto node = add_named_child(
                    regions, "CellElementRegion", region.name );
                set_attribute(
                    node, "cellBlocks", absl::StrCat( "{ ", region.id, " }" ) );
                const auto materials = problem_.materials();
                const auto& material =
                    materials[material_of_region_[region.id]];
                set_attribute( node, "materialList",
                    geos_names( material_list( material ) ) );
            }
        }

        void write_constitutive( pugi::xml_node& problem_node ) const
        {
            auto constitutive = problem_node.append_child( "Constitutive" );
            write_fluid( constitutive );
            add_named_child( constitutive, "NullModel", NULL_SOLID_NAME );
            write_porosity_models( constitutive );
            write_permeability_models( constitutive );
            for( const auto material_id : Indices{ problem_.materials() } )
            {
                write_material( constitutive, material_id );
            }
        }

        void write_fluid( pugi::xml_node& constitutive ) const
        {
            if( const auto* single_phase =
                    dynamic_cast< const SinglePhaseFluid* >( &fluid() ) )
            {
                write_single_phase_fluid( constitutive, *single_phase );
            }
            else if( const auto* dead_oil =
                         dynamic_cast< const DeadOilFluid* >( &fluid() ) )
            {
                write_dead_oil_fluid( constitutive, *dead_oil );
            }
            else if( const auto* co2_brine =
                         dynamic_cast< const Co2BrineFluid* >( &fluid() ) )
            {
                write_co2_brine_fluid( constitutive, *co2_brine );
            }
        }

        void write_single_phase_fluid( pugi::xml_node& constitutive,
            const SinglePhaseFluid& fluid_model ) const
        {
            auto node = add_named_child( constitutive,
                is_thermal() ? "ThermalCompressibleSinglePhaseFluid"
                             : "CompressibleSinglePhaseFluid",
                FLUID_MODEL_NAME );
            set_real_attribute(
                node, "defaultDensity", fluid_model.reference_density() );
            set_real_attribute(
                node, "defaultViscosity", fluid_model.reference_viscosity() );
            set_real_attribute(
                node, "referenceDensity", fluid_model.reference_density() );
            set_real_attribute(
                node, "referenceViscosity", fluid_model.reference_viscosity() );
            set_real_attribute(
                node, "referencePressure", fluid_model.reference_pressure() );
            set_real_attribute(
                node, "compressibility", fluid_model.compressibility() );
            set_real_attribute(
                node, "viscosibility", fluid_model.viscosibility() );
            if( is_thermal() && fluid_model.thermal() )
            {
                const auto& thermal = fluid_model.thermal().value();
                set_real_attribute( node, "referenceTemperature",
                    problem_.thermal()->reference_temperature );
                set_real_attribute( node, "specificHeatCapacity",
                    thermal.specific_heat_capacity );
                set_real_attribute(
                    node, "thermalExpansionCoeff", thermal.thermal_expansion );
                set_real_attribute( node, "temperatureViscosityCoefficient",
                    thermal.viscosity_temperature_coefficient );
            }
        }

        void write_dead_oil_fluid( pugi::xml_node& constitutive,
            const DeadOilFluid& fluid_model ) const
        {
            auto node = add_named_child(
                constitutive, "DeadOilFluid", FLUID_MODEL_NAME );
            set_attribute( node, "phaseNames",
                geos_names( geos_phase_names( fluid_model.phases() ) ) );
            set_attribute( node, "componentNames",
                geos_names( component_names( fluid_model ) ) );
            set_attribute( node, "componentMolarWeight",
                geos_array( component_molar_weights( fluid_model ) ) );
            set_attribute( node, "surfaceDensities",
                geos_array( fluid_model.surface_densities() ) );
            std::vector< std::string > fvf_tables;
            std::vector< std::string > viscosity_tables;
            for( const auto phase : fluid_model.phases() )
            {
                if( phase == PHASE::water )
                {
                    continue;
                }
                fvf_tables.push_back( pvt_table_name( phase, "fvf" ) );
                viscosity_tables.push_back( pvt_table_name( phase, "visc" ) );
            }
            set_attribute( node, "hydrocarbonFormationVolFactorTableNames",
                geos_names( fvf_tables ) );
            set_attribute( node, "hydrocarbonViscosityTableNames",
                geos_names( viscosity_tables ) );
            if( fluid_model.water_pvt() )
            {
                write_water_pvt( node, fluid_model.water_pvt().value() );
            }
        }

        static void write_water_pvt(
            pugi::xml_node& node, const WaterPvt& water )
        {
            set_real_attribute(
                node, "waterReferencePressure", water.reference_pressure );
            set_real_attribute( node, "waterFormationVolumeFactor",
                water.formation_volume_factor );
            set_real_attribute(
                node, "waterCompressibility", water.compressibility );
            set_real_attribute( node, "waterViscosity", water.viscosity );
        }

        void write_co2_brine_fluid( pugi::xml_node& constitutive,
            const Co2BrineFluid& fluid_model ) const
        {
            const auto is_ezrokhi =
                fluid_model.correlation() == CO2_BRINE_CORRELATION::ezrokhi;
            auto node = add_named_child( constitutive,
                absl::StrCat( "CO2Brine", is_ezrokhi ? "Ezrokhi" : "Phillips",
                    is_thermal() ? "Thermal" : "", "Fluid" ),
                FLUID_MODEL_NAME );
            set_attribute( node, "phaseNames",
                geos_names( geos_phase_names( fluid_model.phases() ) ) );
            set_attribute( node, "componentNames",
                geos_names( component_names( fluid_model ) ) );
            set_attribute( node, "componentMolarWeight",
                geos_array( component_molar_weights( fluid_model ) ) );
            const auto& pressure = fluid_model.pressure_range();
            set_attribute( node, "pressureCoordinates",
                geos_array( { pressure.min, pressure.max } ) );
            set_real_attribute( node, "pressureInterval", pressure.step );
            const auto& temperature = fluid_model.temperature_range();
            set_attribute( node, "temperatureCoordinates",
                geos_array( { temperature.min, temperature.max } ) );
            set_real_attribute( node, "temperatureInterval", temperature.step );
            set_real_attribute( node, "salinity", fluid_model.salinity() );
            if( is_ezrokhi
                && !fluid_model.ezrokhi_density_coefficients().empty() )
            {
                set_attribute( node, "ezrokhiDensityCoefficients",
                    geos_array( fluid_model.ezrokhi_density_coefficients() ) );
            }
            if( is_ezrokhi
                && !fluid_model.ezrokhi_viscosity_coefficients().empty() )
            {
                set_attribute( node, "ezrokhiViscosityCoefficients",
                    geos_array(
                        fluid_model.ezrokhi_viscosity_coefficients() ) );
            }
        }

        void write_porosity_models( pugi::xml_node& constitutive ) const
        {
            for( const auto group_id : Indices{ porosity_groups_ } )
            {
                const auto& group = porosity_groups_[group_id];
                auto node = add_named_child( constitutive, "PressurePorosity",
                    absl::StrCat( "porosity_", group_id ) );
                set_real_attribute(
                    node, "defaultReferencePorosity", group.default_porosity );
                set_real_attribute(
                    node, "referencePressure", group.reference_pressure );
                set_real_attribute(
                    node, "compressibility", group.compressibility );
            }
        }

        void write_permeability_models( pugi::xml_node& constitutive ) const
        {
            for( const auto group_id : Indices{ permeability_groups_ } )
            {
                auto node =
                    add_named_child( constitutive, "ConstantPermeability",
                        absl::StrCat( "permeability_", group_id ) );
                set_attribute( node, "permeabilityComponents",
                    geos_array( permeability_groups_[group_id] ) );
            }
        }

        void write_material(
            pugi::xml_node& constitutive, index_t material_id ) const
        {
            const auto materials = problem_.materials();
            const auto& material = materials[material_id];
            auto rock = add_named_child( constitutive,
                "CompressibleSolidConstantPermeability",
                absl::StrCat( "rock_", material.name ) );
            set_attribute( rock, "solidModelName", NULL_SOLID_NAME );
            set_attribute( rock, "porosityModelName",
                absl::StrCat(
                    "porosity_", material_porosity_group_[material_id] ) );
            set_attribute( rock, "permeabilityModelName",
                absl::StrCat( "permeability_",
                    material_permeability_group_[material_id] ) );
            if( is_thermal() )
            {
                set_attribute( rock, "solidInternalEnergyModelName",
                    absl::StrCat( "rockEnergy_", material.name ) );
                write_rock_thermal( constitutive, material );
            }
            if( is_multiphase() )
            {
                write_relative_permeability( constitutive, material );
                write_capillary_pressure( constitutive, material );
            }
        }

        void write_rock_thermal(
            pugi::xml_node& constitutive, const Material& material ) const
        {
            const auto& thermal = material.thermal.value();
            auto energy = add_named_child( constitutive, "SolidInternalEnergy",
                absl::StrCat( "rockEnergy_", material.name ) );
            set_real_attribute( energy, "referenceVolumetricHeatCapacity",
                thermal.volumetric_heat_capacity );
            set_real_attribute( energy, "referenceTemperature",
                problem_.thermal()->reference_temperature );
            set_attribute( energy, "referenceInternalEnergy", "0" );
            const auto conductivity_name =
                absl::StrCat( "thermalCond_", material.name );
            if( is_multiphase() )
            {
                auto conductivity = add_named_child( constitutive,
                    "MultiPhaseConstantThermalConductivity",
                    conductivity_name );
                set_attribute( conductivity, "phaseNames",
                    geos_names( geos_phase_names( fluid().phases() ) ) );
                set_attribute( conductivity, "thermalConductivityComponents",
                    geos_array( thermal.thermal_conductivity ) );
                return;
            }
            auto conductivity = add_named_child( constitutive,
                "SinglePhaseThermalConductivity", conductivity_name );
            set_attribute( conductivity, "defaultThermalConductivityComponents",
                geos_array( thermal.thermal_conductivity ) );
        }

        void write_relative_permeability(
            pugi::xml_node& constitutive, const Material& material ) const
        {
            const auto name = absl::StrCat( "relperm_", material.name );
            const auto& fluid_phases = fluid().phases();
            if( const auto* brooks_corey =
                    dynamic_cast< const BrooksCoreyRelativePermeability* >(
                        material.relative_permeability.get() ) )
            {
                const auto& phases = brooks_corey->phases();
                auto node = add_named_child(
                    constitutive, "BrooksCoreyRelativePermeability", name );
                set_attribute( node, "phaseNames",
                    geos_names( geos_phase_names( fluid_phases ) ) );
                set_attribute( node, "phaseMinVolumeFraction",
                    geos_array( in_fluid_phase_order( phases,
                        brooks_corey->residual_saturations(),
                        fluid_phases ) ) );
                set_attribute( node, "phaseRelPermExponent",
                    geos_array( in_fluid_phase_order(
                        phases, brooks_corey->exponents(), fluid_phases ) ) );
                set_attribute( node, "phaseRelPermMaxValue",
                    geos_array( in_fluid_phase_order(
                        phases, brooks_corey->max_values(), fluid_phases ) ) );
                return;
            }
            auto node = add_named_child(
                constitutive, "TableRelativePermeability", name );
            set_attribute( node, "phaseNames",
                geos_names( geos_phase_names( fluid_phases ) ) );
            set_attribute( node, "wettingNonWettingRelPermTableNames",
                geos_names( { absl::StrCat( name, "_wetting" ),
                    absl::StrCat( name, "_non_wetting" ) } ) );
        }

        void write_capillary_pressure(
            pugi::xml_node& constitutive, const Material& material ) const
        {
            if( !material.capillary_pressure )
            {
                return;
            }
            const auto name = absl::StrCat( "cappres_", material.name );
            const auto& fluid_phases = fluid().phases();
            if( const auto* brooks_corey =
                    dynamic_cast< const BrooksCoreyCapillaryPressure* >(
                        material.capillary_pressure.get() ) )
            {
                const auto& phases = brooks_corey->phases();
                auto node = add_named_child(
                    constitutive, "BrooksCoreyCapillaryPressure", name );
                set_attribute( node, "phaseNames",
                    geos_names( geos_phase_names( fluid_phases ) ) );
                set_attribute( node, "phaseMinVolumeFraction",
                    geos_array( in_fluid_phase_order( phases,
                        brooks_corey->residual_saturations(),
                        fluid_phases ) ) );
                set_attribute( node, "phaseEntryPressure",
                    geos_array( in_fluid_phase_order( phases,
                        brooks_corey->entry_pressures(), fluid_phases ) ) );
                set_attribute( node, "phaseCapPressureExponentInv",
                    geos_array( in_fluid_phase_order( phases,
                        brooks_corey->pore_size_distribution_indices(),
                        fluid_phases ) ) );
                return;
            }
            auto node =
                add_named_child( constitutive, "TableCapillaryPressure", name );
            set_attribute( node, "phaseNames",
                geos_names( geos_phase_names( fluid_phases ) ) );
            set_attribute( node, "wettingNonWettingCapPressureTableName",
                absl::StrCat( name, "_pc" ) );
        }

        void write_uniform_condition( pugi::xml_node& field_specifications,
            const UniformInitialCondition& condition,
            std::string_view region_name ) const
        {
            auto pressure = add_field_specification( field_specifications,
                absl::StrCat( "initialPressure_", region_name ), region_name,
                "pressure" );
            set_real_attribute( pressure, "scale", condition.pressure() );
            if( condition.temperature() )
            {
                auto temperature =
                    add_field_specification( field_specifications,
                        absl::StrCat( "initialTemperature_", region_name ),
                        region_name, "temperature" );
                set_real_attribute(
                    temperature, "scale", condition.temperature().value() );
            }
            if( !is_multiphase() )
            {
                return;
            }
            for( const auto component_id : Indices{ condition.composition() } )
            {
                auto composition =
                    add_field_specification( field_specifications,
                        absl::StrCat( "initialComposition_",
                            fluid().components()[component_id].name, "_",
                            region_name ),
                        region_name, "globalCompFraction" );
                set_attribute(
                    composition, "component", absl::StrCat( component_id ) );
                set_real_attribute( composition, "scale",
                    condition.composition()[component_id] );
            }
        }

        void write_output_events( pugi::xml_node& events ) const
        {
            const auto& output = problem_.output();
            const auto vtk_target =
                absl::StrCat( "/Outputs/", VTK_OUTPUT_NAME );
            if( output.frequency )
            {
                auto event =
                    add_named_child( events, "PeriodicEvent", "vtkOutputs" );
                set_real_attribute(
                    event, "timeFrequency", output.frequency.value() );
                set_attribute( event, "target", vtk_target );
            }
            for( const auto time_id : Indices{ output.times } )
            {
                auto event = add_named_child( events, "SoloEvent",
                    absl::StrCat( "vtkOutput_", time_id ) );
                set_real_attribute(
                    event, "targetTime", output.times[time_id] );
                set_attribute( event, "target", vtk_target );
            }
            const auto restart_target =
                absl::StrCat( "/Outputs/", RESTART_OUTPUT_NAME );
            if( output.restart_frequency )
            {
                auto event =
                    add_named_child( events, "PeriodicEvent", "restarts" );
                set_real_attribute(
                    event, "timeFrequency", output.restart_frequency.value() );
                // Restarts must not shorten the time steps of the run
                set_attribute( event, "targetExactTimestep", "0" );
                set_attribute( event, "target", restart_target );
            }
            if( output.restart )
            {
                auto event = add_named_child( events, "SoloEvent", "restart" );
                set_real_attribute(
                    event, "targetTime", problem_.schedule().end_time );
                set_attribute( event, "target", restart_target );
            }
        }

        void write_outputs( pugi::xml_node& problem_node ) const
        {
            const auto& output = problem_.output();
            const auto has_vtk = geos_writes_vtk_output( output );
            const auto has_restart = geos_writes_restart_output( output );
            if( !has_vtk && !has_restart )
            {
                return;
            }
            auto outputs = problem_node.append_child( "Outputs" );
            if( has_vtk )
            {
                auto vtk = add_named_child( outputs, "VTK", VTK_OUTPUT_NAME );
                set_attribute(
                    vtk, "format", output.binary ? "binary" : "ascii" );
            }
            if( has_restart )
            {
                add_named_child( outputs, "Restart", RESTART_OUTPUT_NAME );
            }
        }

    private:
        const PhysicalProblem& problem_;
        const GeosDeckContext& context_;
        std::vector< index_t > material_of_region_;
        std::vector< index_t > condition_of_region_;
        std::vector< GeosPorosityGroup > porosity_groups_;
        std::vector< index_t > material_porosity_group_;
        std::vector< std::array< double, 3 > > permeability_groups_;
        std::vector< index_t > material_permeability_group_;
        std::vector< GeosImportedField > imported_fields_;
        std::vector< GeosTableFunction > tables_;
    };

    GeosPorosityGroups geos_porosity_groups(
        const PhysicalProblem& problem, const GeosDeckContext& context )
    {
        GeosPorosityGroups result;
        const auto imported = context.porosity_array.has_value();
        for( const auto& material : problem.materials() )
        {
            const GeosPorosityGroup group{ material.pore_compressibility,
                material.reference_pressure,
                material.porosity.value_or(
                    context.mean_porosity.value_or( 0 ) ) };
            const auto existing = absl::c_find_if(
                result.groups, [&group, imported]( const auto& other ) {
                    return other.compressibility == group.compressibility
                           && other.reference_pressure
                                  == group.reference_pressure
                           && ( imported
                                || other.default_porosity
                                       == group.default_porosity );
                } );
            result.material_groups.push_back(
                static_cast< index_t >( existing - result.groups.begin() ) );
            if( existing == result.groups.end() )
            {
                result.groups.push_back( group );
            }
        }
        return result;
    }

    bool geos_writes_vtk_output( const OutputSettings& output )
    {
        return output.frequency || !output.times.empty();
    }

    bool geos_writes_restart_output( const OutputSettings& output )
    {
        return output.restart || output.restart_frequency.has_value();
    }

    std::string geos_porosity_vtu_name(
        std::string_view name, index_t group, index_t nb_groups )
    {
        if( nb_groups == 1 )
        {
            return to_string( name );
        }
        return absl::StrCat( name, "_", group );
    }

    GeosSimulationWriter::GeosSimulationWriter( const PhysicalProblem& problem,
        const GeosDeckContext& context,
        const GeosRegionAssignment& materials,
        const GeosRegionAssignment& initial_conditions )
        : impl_( problem, context, materials, initial_conditions )
    {
    }

    GeosSimulationWriter::~GeosSimulationWriter() = default;

    absl::Span< const GeosImportedField >
        GeosSimulationWriter::imported_fields() const
    {
        return impl_->imported_fields();
    }

    void GeosSimulationWriter::write_physics(
        pugi::xml_node& problem_node ) const
    {
        impl_->write_physics( problem_node );
    }

    void GeosSimulationWriter::write_initial_conditions(
        pugi::xml_node& field_specifications ) const
    {
        impl_->write_initial_conditions( field_specifications );
    }

    void GeosSimulationWriter::write_source_terms(
        pugi::xml_node& field_specifications ) const
    {
        impl_->write_source_terms( field_specifications );
    }

    bool GeosSimulationWriter::has_functions() const
    {
        return impl_->has_functions();
    }

    void GeosSimulationWriter::write_functions(
        pugi::xml_node& functions ) const
    {
        impl_->write_functions( functions );
    }

    void GeosSimulationWriter::write_events_and_outputs(
        pugi::xml_node& problem_node ) const
    {
        impl_->write_events_and_outputs( problem_node );
    }

    void GeosSimulationWriter::apply_solver_options(
        pugi::xml_node& problem_node ) const
    {
        impl_->apply_solver_options( problem_node );
    }
} // namespace geode::internal
