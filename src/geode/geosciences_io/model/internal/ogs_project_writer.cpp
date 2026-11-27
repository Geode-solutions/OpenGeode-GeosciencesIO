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

#include <geode/geosciences_io/model/internal/ogs_project_writer.hpp>

#include <array>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#include <pugixml.hpp>

#include <absl/algorithm/container.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_split.h>

#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/heat_source.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/schedule.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>

#include <geode/geosciences_io/model/internal/simulation_xml.hpp>

namespace
{
    // OGS solver defaults, overridable through the "ogs" solver options
    constexpr double LIQUID_FLOW_RELTOL{ 1e-10 };
    constexpr geode::index_t PICARD_MAX_ITERATIONS{ 20 };
    constexpr geode::index_t NEWTON_MAX_ITERATIONS{ 60 };
    // McWhorter benchmark values, capillary then gas pressure
    constexpr std::string_view TWO_PHASE_ABSTOLS{ "1e-05 1e-05" };
    constexpr std::string_view TWO_PHASE_RELTOLS{ "1e-08 1e-06" };
    constexpr std::string_view JACOBIAN_EPSILONS{ "0.001 0.01" };
    // HeatTransportInStationaryFlow benchmark, temperature then pressure
    constexpr std::string_view HT_RELTOLS{ "1e-10 1e-10" };
    constexpr geode::index_t HT_MAX_ITERATIONS{ 50 };
    // Newton of the mechanics processes: one reltol per component, and a
    // tiny abstol stopping on components that stay exactly 0 (a relative
    // criterion never holds on them)
    constexpr double MECHANICS_RELTOL{ 1e-10 };
    constexpr double MECHANICS_ABSTOL{ 1e-15 };
    constexpr geode::index_t MECHANICS_MAX_ITERATIONS{ 50 };
    // IterationNumberBasedTimeStepping values of the OGS benchmarks
    constexpr std::string_view NUMBER_ITERATIONS{ "1 4 10" };
    constexpr std::string_view MULTIPLIERS{ "2 1 0.5" };
    // Newton of the mechanics processes needs 3 to 5 iterations to reach
    // MECHANICS_RELTOL on nonlinear (thermal) problems: the steps grow up to
    // 5 iterations
    constexpr std::string_view MECHANICS_NUMBER_ITERATIONS{ "5 8 15" };
    // Output only at the fixed output times and at the last step
    constexpr geode::index_t OUTPUT_EACH_STEPS{ 1000000000 };

    constexpr std::string_view PROCESS_NAME{ "flow" };
    constexpr std::string_view NONLINEAR_SOLVER{ "nonlinear_solver" };
    constexpr std::string_view LINEAR_SOLVER{ "linear_solver" };

    pugi::xml_node add_text(
        pugi::xml_node parent, std::string_view name, std::string_view value )
    {
        auto node = parent.append_child( std::string{ name }.c_str() );
        node.text().set( std::string{ value }.c_str() );
        return node;
    }

    pugi::xml_node add_real(
        pugi::xml_node parent, std::string_view name, double value )
    {
        return add_text( parent, name, geode::internal::xml_real( value ) );
    }

    std::string reals( absl::Span< const double > values )
    {
        return absl::StrJoin(
            values, " ", []( std::string* out, double value ) {
                absl::StrAppend( out, geode::internal::xml_real( value ) );
            } );
    }

    // Real term of an OGS (exprtk) expression, negative ones in parentheses
    std::string expression_real( double value )
    {
        const auto real = geode::internal::xml_real( value );
        return value < 0. ? absl::StrCat( "(", real, ")" ) : real;
    }

    // Value of an exponential law in the OGS variables
    std::string law_expression( const geode::internal::OgsExponentialLaw& law )
    {
        return absl::StrCat( expression_real( law.reference_value ), "*exp(",
            expression_real( law.pressure_factor ), "*(liquid_phase_pressure-",
            expression_real( law.reference_pressure ), ")+",
            expression_real( law.temperature_factor ), "*(temperature-",
            expression_real( law.reference_temperature ), "))" );
    }

    // Function property of a law of pressure and temperature: its value and
    // both derivatives (OGS stops when HT asks a missing one)
    void write_function_law( pugi::xml_node& property,
        const geode::internal::OgsExponentialLaw& law )
    {
        add_text( property, "type", "Function" );
        const auto value = law_expression( law );
        add_text( property.append_child( "value" ), "expression", value );
        for( const auto& [variable, factor] :
            std::array< std::pair< std::string_view, double >, 2 >{
                { { "liquid_phase_pressure", law.pressure_factor },
                    { "temperature", law.temperature_factor } } } )
        {
            auto dvalue = property.append_child( "dvalue" );
            add_text( dvalue, "variable_name", variable );
            add_text( dvalue, "expression",
                absl::StrCat( expression_real( factor ), "*", value ) );
        }
    }

    // Sets the text of the descendant at path ("a/b"), created if missing
    void apply_element_text(
        pugi::xml_node& node, std::string_view path, std::string_view value )
    {
        auto current = node;
        for( const auto name : absl::StrSplit( path, '/' ) )
        {
            const std::string child_name{ name };
            auto child = current.child( child_name.c_str() );
            if( !child )
            {
                child = current.append_child( child_name.c_str() );
            }
            current = child;
        }
        current.text().set( std::string{ value }.c_str() );
    }

    // Parameter read from a mesh array
    struct OgsMeshParameter
    {
        std::string_view name;
        // MeshElement or MeshNode
        std::string_view type;
        std::string_view mesh;
        std::string_view field;
    };

    // Medium or phase property interpolated in a curve
    struct OgsCurveProperty
    {
        std::string_view name;
        std::string_view variable;
        std::string_view curve;
    };

    class OgsProjectBuilder
    {
    public:
        OgsProjectBuilder( const geode::PhysicalProblem& problem,
            const geode::internal::OgsProjectContext& context,
            pugi::xml_node& root )
            : problem_( problem ), context_( context ), root_( root )
        {
            parameters_ = root_.append_child( "parameters" );
        }

        void write()
        {
            write_meshes();
            write_process();
            write_media();
            write_time_loop();
            write_process_variables();
            write_solvers();
            if( is_two_phase() || has_heat_curves() )
            {
                auto curves = root_.append_child( "curves" );
                if( is_two_phase() )
                {
                    write_two_phase_curves( curves );
                }
                write_heat_curves( curves );
            }
        }

    private:
        [[nodiscard]] bool is_two_phase() const
        {
            return context_.process
                   == geode::internal::OGS_PROCESS::two_phase_flow_pp;
        }

        [[nodiscard]] bool is_ht() const
        {
            return context_.process == geode::internal::OGS_PROCESS::ht;
        }

        [[nodiscard]] bool is_mechanics() const
        {
            return geode::internal::has_displacement( context_.process );
        }

        [[nodiscard]] bool is_small_deformation() const
        {
            return context_.process
                   == geode::internal::OGS_PROCESS::small_deformation;
        }

        // Primary variable components: temperature, pressure, displacement
        [[nodiscard]] geode::index_t nb_components() const
        {
            geode::index_t count{ 3 };
            count += is_small_deformation() ? 0 : 1;
            count +=
                geode::internal::has_temperature( context_.process ) ? 1 : 0;
            return count;
        }

        // Uniform value as a Constant, a bulk cell array as a MeshElement
        void add_cell_parameter(
            std::string_view name, const geode::internal::OgsCellValue& value )
        {
            if( value.value )
            {
                add_constant_parameter( name, value.value.value() );
                return;
            }
            add_mesh_parameter(
                { name, "MeshElement", context_.prefix, value.field_name } );
        }

        // Constant property of a uniform value, Parameter property of the
        // bulk cell parameter otherwise
        void write_cell_property( pugi::xml_node& properties,
            std::string_view name,
            std::string_view parameter,
            const geode::internal::OgsCellValue& value ) const
        {
            if( value.value )
            {
                write_constant( properties, name, value.value.value() );
                return;
            }
            write_parameter_property( properties, name, parameter );
        }

        void write_mechanics_process( pugi::xml_node& process )
        {
            add_text( process, "type", mechanics_process_type() );
            add_text( process, "integration_order", "2" );
            if( geode::internal::has_temperature( context_.process ) )
            {
                add_text( process, "equation_balance_type", "mass" );
            }
            auto relation = process.append_child( "constitutive_relation" );
            add_text( relation, "type", "LinearElasticIsotropic" );
            add_text( relation, "youngs_modulus", "youngs_modulus" );
            add_text( relation, "poissons_ratio", "poissons_ratio" );
            add_cell_parameter( "youngs_modulus", context_.youngs_modulus );
            add_cell_parameter( "poissons_ratio", context_.poissons_ratio );
            write_mechanics_process_variables(
                process.append_child( "process_variables" ) );
            auto secondary = process.append_child( "secondary_variables" );
            for( const auto* name : { "sigma", "epsilon" } )
            {
                secondary.append_child( "secondary_variable" )
                    .append_attribute( "name" )
                    .set_value( name );
            }
            add_text( process, "specific_body_force", body_force() );
            write_initial_stress( process );
        }

        [[nodiscard]] std::string_view mechanics_process_type() const
        {
            if( is_small_deformation() )
            {
                return "SMALL_DEFORMATION";
            }
            return geode::internal::has_temperature( context_.process )
                       ? "THERMO_HYDRO_MECHANICS"
                       : "HYDRO_MECHANICS";
        }

        void write_mechanics_process_variables( pugi::xml_node variables ) const
        {
            if( is_small_deformation() )
            {
                add_text( variables, "process_variable", "displacement" );
                return;
            }
            if( geode::internal::has_temperature( context_.process ) )
            {
                add_text( variables, "temperature", "temperature" );
            }
            add_text( variables, "pressure", "pressure" );
            add_text( variables, "displacement", "displacement" );
        }

        // SMALL_DEFORMATION reads an effective stress, the poromechanics
        // processes a total one
        void write_initial_stress( pugi::xml_node& process )
        {
            if( context_.initial_stress_field.empty() )
            {
                return;
            }
            auto stress =
                add_text( process, "initial_stress", "initial_stress" );
            if( !is_small_deformation() )
            {
                stress.append_attribute( "type" ).set_value( "total" );
            }
            add_mesh_parameter( { "initial_stress", "MeshElement",
                context_.prefix, context_.initial_stress_field } );
        }

        // Bulk cell parameters shared by every medium
        void add_media_cell_parameters()
        {
            for( const auto& [name, value] :
                std::array< std::pair< std::string_view,
                                const geode::internal::OgsCellValue* >,
                    4 >{ { { "solid_density", &context_.solid_density },
                    { "biot_coefficient", &context_.biot_coefficient },
                    { "thermal_expansivity", &context_.thermal_expansivity },
                    { "solid_heat_capacity",
                        &context_.solid_heat_capacity } } } )
            {
                if( !value->value && !value->field_name.empty() )
                {
                    add_cell_parameter( name, *value );
                }
            }
        }

        // Without flow, a single medium with the solid density
        void write_mechanics_media( pugi::xml_node& media )
        {
            add_media_cell_parameters();
            if( is_small_deformation() )
            {
                auto medium = media.append_child( "medium" );
                medium.append_attribute( "id" ).set_value( 0 );
                write_mechanics_solid_phase( medium.append_child( "phases" ) );
                medium.append_child( "properties" );
                return;
            }
            for( const auto material_id :
                geode::Indices{ problem_.materials() } )
            {
                const auto& material = problem_.materials()[material_id];
                auto medium = media.append_child( "medium" );
                medium.append_attribute( "id" ).set_value( material_id );
                auto phases = medium.append_child( "phases" );
                write_liquid_phase( phases );
                write_mechanics_solid_phase( phases );
                auto properties = medium.append_child( "properties" );
                write_permeability( properties, material, material_id );
                write_porosity( properties, material );
                write_cell_property( properties, "biot_coefficient",
                    "biot_coefficient", context_.biot_coefficient );
                if( geode::internal::has_temperature( context_.process ) )
                {
                    write_constant_values( properties, "thermal_conductivity",
                        material.thermal.value().thermal_conductivity );
                    continue;
                }
                write_constant( properties, "reference_temperature",
                    context_.reference_temperature );
            }
        }

        void write_mechanics_solid_phase( pugi::xml_node phases ) const
        {
            auto solid = phases.append_child( "phase" );
            add_text( solid, "type", "Solid" );
            auto properties = solid.append_child( "properties" );
            write_cell_property( properties, "density", "solid_density",
                context_.solid_density );
            if( !geode::internal::has_temperature( context_.process ) )
            {
                return;
            }
            // The rock heat capacity is rho_s c_s (1 - phi) with the grain
            // density, also the density of the gravity
            write_cell_property( properties, "specific_heat_capacity",
                "solid_heat_capacity", context_.solid_heat_capacity );
            write_cell_property( properties, "thermal_expansivity",
                "thermal_expansivity", context_.thermal_expansivity );
        }

        void write_mechanics_convergence( pugi::xml_node& convergence ) const
        {
            add_text( convergence, "type", "PerComponentDeltaX" );
            add_text( convergence, "norm_type", "NORM2" );
            const std::vector< double > abstols(
                nb_components(), MECHANICS_ABSTOL );
            const std::vector< double > reltols(
                nb_components(), MECHANICS_RELTOL );
            add_text( convergence, "abstols", reals( abstols ) );
            add_text( convergence, "reltols", reals( reltols ) );
        }

        void write_mechanics_output_variables( pugi::xml_node& variables ) const
        {
            if( geode::internal::has_temperature( context_.process ) )
            {
                add_text( variables, "variable", "temperature" );
            }
            if( !is_small_deformation() )
            {
                add_text( variables, "variable", "pressure" );
            }
            for( const auto* name : { "displacement", "sigma", "epsilon" } )
            {
                add_text( variables, "variable", name );
            }
        }

        void write_mechanics_variables( pugi::xml_node& variables )
        {
            if( geode::internal::has_temperature( context_.process ) )
            {
                write_temperature_variable( variables );
            }
            if( !is_small_deformation() )
            {
                write_pressure_variable( variables );
            }
            write_displacement_variable( variables );
        }

        void write_displacement_variable( pugi::xml_node& variables )
        {
            auto initial = parameters_.append_child( "parameter" );
            add_text( initial, "name", "initial_displacement" );
            add_text( initial, "type", "Constant" );
            add_text( initial, "values", "0 0 0" );
            auto variable = variables.append_child( "process_variable" );
            add_text( variable, "name", "displacement" );
            add_text( variable, "components", "3" );
            add_text( variable, "order", "1" );
            add_text( variable, "initial_condition", "initial_displacement" );
            auto conditions = variable.append_child( "boundary_conditions" );
            for( const auto& boundary : context_.boundaries )
            {
                if( const auto type = mechanics_condition_type( boundary ) )
                {
                    add_mechanics_condition( conditions, boundary, *type );
                }
            }
        }

        static std::optional< std::string_view > mechanics_condition_type(
            const geode::internal::OgsBoundaryMesh& boundary )
        {
            switch( boundary.property )
            {
                case geode::PHYSICAL_PROPERTY_NAME::boundary_displacement:
                    return "Dirichlet";
                case geode::PHYSICAL_PROPERTY_NAME::boundary_traction:
                    return "Neumann";
                case geode::PHYSICAL_PROPERTY_NAME::boundary_normal_traction:
                    return "NormalTraction";
                default:
                    return std::nullopt;
            }
        }

        // Uniform value as a Constant, per polygon values as a MeshElement
        // parameter of the boundary mesh
        void add_mechanics_condition( pugi::xml_node& conditions,
            const geode::internal::OgsBoundaryMesh& boundary,
            std::string_view type )
        {
            if( boundary.value )
            {
                add_constant_parameter( boundary.stem, boundary.value.value() );
            }
            else
            {
                add_mesh_parameter(
                    { boundary.stem, "MeshElement", boundary.stem, "values" } );
            }
            auto condition = conditions.append_child( "boundary_condition" );
            add_text( condition, "mesh", boundary.stem );
            add_text( condition, "type", type );
            if( boundary.component )
            {
                add_text( condition, "component",
                    absl::StrCat( boundary.component.value() ) );
            }
            add_text( condition, "parameter", boundary.stem );
        }

        void add_constant_parameter( std::string_view name, double value )
        {
            auto node = parameters_.append_child( "parameter" );
            add_text( node, "name", name );
            add_text( node, "type", "Constant" );
            add_real( node, "value", value );
        }

        void add_mesh_parameter( const OgsMeshParameter& parameter )
        {
            auto node = parameters_.append_child( "parameter" );
            add_text( node, "name", parameter.name );
            add_text( node, "type", parameter.type );
            add_text( node, "mesh", parameter.mesh );
            add_text( node, "field_name", parameter.field );
        }

        void add_initial_parameter( std::string_view name,
            const geode::internal::OgsInitialValue& initial )
        {
            if( initial.value )
            {
                add_constant_parameter( name, initial.value.value() );
                return;
            }
            add_mesh_parameter(
                { name, "MeshNode", context_.prefix, initial.field_name } );
        }

        void write_meshes()
        {
            auto meshes = root_.append_child( "meshes" );
            add_text( meshes, "mesh", absl::StrCat( context_.prefix, ".vtu" ) );
            for( const auto& boundary : context_.boundaries )
            {
                add_text(
                    meshes, "mesh", absl::StrCat( boundary.stem, ".vtu" ) );
            }
            for( const auto& well : context_.wells )
            {
                add_text( meshes, "mesh", absl::StrCat( well.stem, ".vtu" ) );
            }
            for( const auto& heat_source : context_.heat_sources )
            {
                add_text(
                    meshes, "mesh", absl::StrCat( heat_source.stem, ".vtu" ) );
            }
        }

        void write_process()
        {
            auto process =
                root_.append_child( "processes" ).append_child( "process" );
            add_text( process, "name", PROCESS_NAME );
            switch( context_.process )
            {
                case geode::internal::OGS_PROCESS::two_phase_flow_pp:
                    write_two_phase_process( process );
                    return;
                case geode::internal::OGS_PROCESS::ht:
                    write_ht_process( process );
                    return;
                case geode::internal::OGS_PROCESS::small_deformation:
                case geode::internal::OGS_PROCESS::hydro_mechanics:
                case geode::internal::OGS_PROCESS::thermo_hydro_mechanics:
                    write_mechanics_process( process );
                    return;
                case geode::internal::OGS_PROCESS::liquid_flow:
                    write_liquid_flow_process( process );
                    return;
            }
        }

        [[nodiscard]] std::string body_force() const
        {
            const auto& gravity = problem_.gravity();
            return reals( std::array< double, 3 >{
                gravity.value( 0 ), gravity.value( 1 ), gravity.value( 2 ) } );
        }

        void write_two_phase_process( pugi::xml_node& process )
        {
            add_text( process, "type", "TWOPHASE_FLOW_PP" );
            add_text( process, "integration_order", "2" );
            auto jacobian = process.append_child( "jacobian_assembler" );
            add_text( jacobian, "type", "CentralDifferences" );
            add_text( jacobian, "epsilons", JACOBIAN_EPSILONS );
            auto variables = process.append_child( "process_variables" );
            add_text( variables, "capillary_pressure", "capillary_pressure" );
            add_text( variables, "gas_pressure", "gas_pressure" );
            process.append_child( "secondary_variables" )
                .append_child( "secondary_variable" )
                .append_attribute( "name" )
                .set_value( "saturation" );
            add_text( process, "specific_body_force", body_force() );
            add_text( process, "mass_lumping", "true" );
            add_text( process, "temperature", "temperature" );
            add_constant_parameter(
                "temperature", context_.reference_temperature );
        }

        void write_ht_process( pugi::xml_node& process ) const
        {
            add_text( process, "type", "HT" );
            add_text( process, "integration_order", "2" );
            auto variables = process.append_child( "process_variables" );
            add_text( variables, "temperature", "temperature" );
            add_text( variables, "pressure", "pressure" );
            add_text( process, "specific_body_force", body_force() );
            add_text( process, "equation_balance_type", "mass" );
        }

        void write_liquid_flow_process( pugi::xml_node& process ) const
        {
            add_text( process, "type", "LIQUID_FLOW" );
            add_text( process, "integration_order", "2" );
            add_text( process.append_child( "process_variables" ),
                "process_variable", "pressure" );
            add_text( process, "specific_body_force", body_force() );
            add_text( process, "equation_balance_type", "mass" );
        }

        // Variable of the exponent of an Exponential property
        struct OgsExponent
        {
            std::string_view variable;
            double reference_condition{ 0 };
            double factor{ 0 };
        };

        // Constant without factor, Exponential with one, Function with both
        void write_fluid_law( pugi::xml_node& properties,
            std::string_view name,
            const geode::internal::OgsExponentialLaw& law ) const
        {
            const auto on_pressure = law.pressure_factor != 0.;
            const auto on_temperature = law.temperature_factor != 0.;
            if( !on_pressure && !on_temperature )
            {
                write_constant( properties, name, law.reference_value );
                return;
            }
            auto property = properties.append_child( "property" );
            add_text( property, "name", name );
            if( on_pressure && on_temperature )
            {
                write_function_law( property, law );
                return;
            }
            add_text( property, "type", "Exponential" );
            add_real( property, "reference_value", law.reference_value );
            add_text( property, "offset", "0" );
            const auto exponent =
                on_pressure
                    ? OgsExponent{ "liquid_phase_pressure",
                          law.reference_pressure, law.pressure_factor }
                    : OgsExponent{ "temperature", law.reference_temperature,
                          law.temperature_factor };
            auto node = property.append_child( "exponent" );
            add_text( node, "variable_name", exponent.variable );
            add_real(
                node, "reference_condition", exponent.reference_condition );
            add_real( node, "factor", exponent.factor );
        }

        void write_constant( pugi::xml_node& properties,
            std::string_view name,
            double value ) const
        {
            auto property = properties.append_child( "property" );
            add_text( property, "name", name );
            add_text( property, "type", "Constant" );
            add_real( property, "value", value );
        }

        void write_parameter_property( pugi::xml_node& properties,
            std::string_view name,
            std::string_view parameter ) const
        {
            auto property = properties.append_child( "property" );
            add_text( property, "name", name );
            add_text( property, "type", "Parameter" );
            add_text( property, "parameter_name", parameter );
        }

        void write_curve_property(
            pugi::xml_node& properties, const OgsCurveProperty& curve ) const
        {
            auto property = properties.append_child( "property" );
            add_text( property, "name", curve.name );
            add_text( property, "type", "Curve" );
            add_text( property, "independent_variable", curve.variable );
            add_text( property, "curve", curve.curve );
        }

        void write_liquid_phase( pugi::xml_node& phases ) const
        {
            auto phase = phases.append_child( "phase" );
            add_text( phase, "type", "AqueousLiquid" );
            auto properties = phase.append_child( "properties" );
            if( is_two_phase() )
            {
                const auto& fluid = dynamic_cast< const geode::DeadOilFluid& >(
                    problem_.fluid() );
                const auto water = water_phase_index( fluid );
                write_constant( properties, "density",
                    fluid.surface_densities()[water]
                        / fluid.water_pvt()->formation_volume_factor );
                write_constant(
                    properties, "viscosity", fluid.water_pvt()->viscosity );
                return;
            }
            write_fluid_law( properties, "density",
                geode::internal::ogs_density_law( problem_ ) );
            write_fluid_law( properties, "viscosity",
                geode::internal::ogs_viscosity_law( problem_ ) );
            if( geode::internal::has_temperature( context_.process ) )
            {
                const auto& fluid =
                    dynamic_cast< const geode::SinglePhaseFluid& >(
                        problem_.fluid() );
                write_constant( properties, "specific_heat_capacity",
                    fluid.thermal().value().specific_heat_capacity );
            }
            if( is_ht() )
            {
                // Unread: the medium has its own conductivity
                write_constant( properties, "thermal_conductivity", 0. );
            }
        }

        // HT reads the rock heat capacity as rho_s * c_s * (1 - phi): a unit
        // density and the volumetric heat capacity give it per solid volume
        void write_solid_phase(
            pugi::xml_node& phases, const geode::Material& material ) const
        {
            auto phase = phases.append_child( "phase" );
            add_text( phase, "type", "Solid" );
            auto properties = phase.append_child( "properties" );
            write_storage( properties, material );
            write_constant( properties, "density", 1. );
            write_constant( properties, "specific_heat_capacity",
                material.thermal.value().volumetric_heat_capacity );
            // Unread: the medium has its own conductivity
            write_constant( properties, "thermal_conductivity", 0. );
        }

        void write_storage(
            pugi::xml_node& properties, const geode::Material& material ) const
        {
            if( context_.has_porosity_field )
            {
                write_parameter_property( properties, "storage", "storage" );
                return;
            }
            write_constant( properties, "storage",
                geode::internal::ogs_storage( material.porosity.value_or( 0. ),
                    material.pore_compressibility ) );
        }

        static geode::index_t water_phase_index(
            const geode::DeadOilFluid& fluid )
        {
            for( const auto phase : geode::Indices{ fluid.phases() } )
            {
                if( fluid.phases()[phase] == geode::PHASE::water )
                {
                    return phase;
                }
            }
            return 0;
        }

        static geode::index_t oil_component( const geode::DeadOilFluid& fluid )
        {
            return geode::phase_component( fluid, geode::PHASE::oil ).value();
        }

        void write_oil_phase( pugi::xml_node& phases ) const
        {
            const auto& fluid =
                dynamic_cast< const geode::DeadOilFluid& >( problem_.fluid() );
            auto phase = phases.append_child( "phase" );
            add_text( phase, "type", "Gas" );
            auto properties = phase.append_child( "properties" );
            write_curve_property( properties,
                { "density", "gas_phase_pressure", "oil_density" } );
            write_curve_property( properties,
                { "viscosity", "gas_phase_pressure", "oil_viscosity" } );
            write_constant( properties, "molar_mass",
                fluid.components()[oil_component( fluid )].molar_weight );
        }

        void write_media()
        {
            auto media = root_.append_child( "media" );
            if( is_mechanics() )
            {
                write_mechanics_media( media );
                return;
            }
            for( const auto material_id :
                geode::Indices{ problem_.materials() } )
            {
                const auto& material = problem_.materials()[material_id];
                auto medium = media.append_child( "medium" );
                medium.append_attribute( "id" ).set_value( material_id );
                auto phases = medium.append_child( "phases" );
                write_liquid_phase( phases );
                if( is_two_phase() )
                {
                    write_oil_phase( phases );
                }
                if( is_ht() )
                {
                    write_solid_phase( phases, material );
                }
                auto properties = medium.append_child( "properties" );
                write_rock_properties( properties, material, material_id );
            }
        }

        void write_rock_properties( pugi::xml_node& properties,
            const geode::Material& material,
            geode::index_t material_id )
        {
            write_permeability( properties, material, material_id );
            write_porosity( properties, material );
            switch( context_.process )
            {
                case geode::internal::OGS_PROCESS::two_phase_flow_pp:
                    write_saturation_properties( properties, material_id );
                    return;
                case geode::internal::OGS_PROCESS::ht:
                    write_thermal_medium( properties, material );
                    return;
                case geode::internal::OGS_PROCESS::small_deformation:
                case geode::internal::OGS_PROCESS::hydro_mechanics:
                case geode::internal::OGS_PROCESS::thermo_hydro_mechanics:
                    // Written by write_mechanics_media
                    return;
                case geode::internal::OGS_PROCESS::liquid_flow:
                    write_storage( properties, material );
                    write_constant( properties, "reference_temperature",
                        context_.reference_temperature );
                    return;
            }
        }

        void write_porosity(
            pugi::xml_node& properties, const geode::Material& material ) const
        {
            if( context_.has_porosity_field )
            {
                write_parameter_property( properties, "porosity", "porosity" );
                return;
            }
            write_constant(
                properties, "porosity", material.porosity.value_or( 0. ) );
        }

        void write_constant_values( pugi::xml_node& properties,
            std::string_view name,
            absl::Span< const double > values ) const
        {
            auto property = properties.append_child( "property" );
            add_text( property, "name", name );
            add_text( property, "type", "Constant" );
            add_text( property, "value", reals( values ) );
        }

        // Diagonal bulk conductivity, no dispersion
        void write_thermal_medium(
            pugi::xml_node& properties, const geode::Material& material ) const
        {
            write_constant_values( properties, "thermal_conductivity",
                material.thermal.value().thermal_conductivity );
            write_constant(
                properties, "thermal_longitudinal_dispersivity", 0. );
            write_constant(
                properties, "thermal_transversal_dispersivity", 0. );
        }

        void write_saturation_properties(
            pugi::xml_node& properties, geode::index_t material_id ) const
        {
            const auto curves = absl::StrCat( "_", material_id );
            const auto saturation = absl::StrCat( "saturation", curves );
            const auto wetting =
                absl::StrCat( "relative_permeability", curves );
            const auto non_wetting = absl::StrCat(
                "relative_permeability_nonwetting_phase", curves );
            write_curve_property( properties,
                { "saturation", "capillary_pressure", saturation } );
            write_curve_property( properties,
                { "relative_permeability", "liquid_saturation", wetting } );
            write_curve_property(
                properties, { "relative_permeability_nonwetting_phase",
                                "liquid_saturation", non_wetting } );
        }

        void write_permeability( pugi::xml_node& properties,
            const geode::Material& material,
            geode::index_t material_id )
        {
            if( context_.has_permeability_field )
            {
                write_parameter_property(
                    properties, "permeability", "permeability" );
            }
            else
            {
                const auto name = absl::StrCat( "permeability_", material_id );
                const auto permeability = material.permeability.value_or(
                    std::array< double, 3 >{ 0., 0., 0. } );
                auto parameter = parameters_.append_child( "parameter" );
                add_text( parameter, "name", name );
                add_text( parameter, "type", "Constant" );
                // TWOPHASE_FLOW_PP reads a scalar permeability
                add_text( parameter, "values",
                    is_two_phase()
                        ? geode::internal::xml_real( permeability[0] )
                        : reals( permeability ) );
                write_parameter_property( properties, "permeability", name );
            }
        }

        void write_time_loop()
        {
            auto time_loop = root_.append_child( "time_loop" );
            auto process =
                time_loop.append_child( "processes" ).append_child( "process" );
            process.append_attribute( "ref" ).set_value(
                std::string{ PROCESS_NAME }.c_str() );
            add_text( process, "nonlinear_solver", NONLINEAR_SOLVER );
            write_convergence(
                process.append_child( "convergence_criterion" ) );
            add_text( process.append_child( "time_discretization" ), "type",
                "BackwardEuler" );
            write_time_stepping( process );
            write_output( time_loop );
        }

        void write_convergence( pugi::xml_node convergence ) const
        {
            switch( context_.process )
            {
                case geode::internal::OGS_PROCESS::two_phase_flow_pp:
                    add_text( convergence, "type", "PerComponentDeltaX" );
                    add_text( convergence, "norm_type", "NORM2" );
                    add_text( convergence, "abstols", TWO_PHASE_ABSTOLS );
                    add_text( convergence, "reltols", TWO_PHASE_RELTOLS );
                    return;
                case geode::internal::OGS_PROCESS::ht:
                    add_text( convergence, "type", "PerComponentDeltaX" );
                    add_text( convergence, "norm_type", "NORM2" );
                    add_text( convergence, "reltols", HT_RELTOLS );
                    return;
                case geode::internal::OGS_PROCESS::small_deformation:
                case geode::internal::OGS_PROCESS::hydro_mechanics:
                case geode::internal::OGS_PROCESS::thermo_hydro_mechanics:
                    write_mechanics_convergence( convergence );
                    return;
                case geode::internal::OGS_PROCESS::liquid_flow:
                    add_text( convergence, "type", "DeltaX" );
                    add_text( convergence, "norm_type", "NORM2" );
                    add_real( convergence, "reltol", LIQUID_FLOW_RELTOL );
                    return;
            }
        }

        void write_time_stepping( pugi::xml_node& process ) const
        {
            const auto& schedule = problem_.schedule();
            auto stepping = process.append_child( "time_stepping" );
            if( schedule.max_time_step )
            {
                add_text(
                    stepping, "type", "IterationNumberBasedTimeStepping" );
                add_real( stepping, "t_initial", schedule.start_time );
                add_real( stepping, "t_end", schedule.end_time );
                add_real( stepping, "initial_dt", schedule.initial_time_step );
                add_real( stepping, "minimum_dt", schedule.initial_time_step );
                add_real(
                    stepping, "maximum_dt", schedule.max_time_step.value() );
                add_text( stepping, "number_iterations",
                    is_mechanics() ? MECHANICS_NUMBER_ITERATIONS
                                   : NUMBER_ITERATIONS );
                add_text( stepping, "multiplier", MULTIPLIERS );
                return;
            }
            add_text( stepping, "type", "FixedTimeStepping" );
            add_real( stepping, "t_initial", schedule.start_time );
            add_real( stepping, "t_end", schedule.end_time );
            const auto nb_steps = static_cast< geode::index_t >(
                std::ceil( ( schedule.end_time - schedule.start_time )
                               / schedule.initial_time_step
                           - geode::GLOBAL_EPSILON ) );
            auto pair =
                stepping.append_child( "timesteps" ).append_child( "pair" );
            add_text( pair, "repeat", absl::StrCat( nb_steps ) );
            add_real( pair, "delta_t", schedule.initial_time_step );
        }

        void write_output( pugi::xml_node& time_loop ) const
        {
            auto output = time_loop.append_child( "output" );
            add_text( output, "type", "VTK" );
            add_text( output, "prefix", context_.prefix );
            auto pair =
                output.append_child( "timesteps" ).append_child( "pair" );
            add_text( pair, "repeat", "1" );
            add_text( pair, "each_steps", absl::StrCat( OUTPUT_EACH_STEPS ) );
            add_text( output, "fixed_output_times",
                reals( geode::internal::ogs_output_times(
                    problem_.schedule(), problem_.output() ) ) );
            write_output_variables( output.append_child( "variables" ) );
            add_text( output, "suffix", "_ts_{:timestep}_t_{:time}" );
        }

        void write_output_variables( pugi::xml_node variables ) const
        {
            switch( context_.process )
            {
                case geode::internal::OGS_PROCESS::two_phase_flow_pp:
                    add_text( variables, "variable", "capillary_pressure" );
                    add_text( variables, "variable", "gas_pressure" );
                    add_text( variables, "variable", "saturation" );
                    return;
                case geode::internal::OGS_PROCESS::ht:
                    add_text( variables, "variable", "temperature" );
                    add_text( variables, "variable", "pressure" );
                    return;
                case geode::internal::OGS_PROCESS::small_deformation:
                case geode::internal::OGS_PROCESS::hydro_mechanics:
                case geode::internal::OGS_PROCESS::thermo_hydro_mechanics:
                    write_mechanics_output_variables( variables );
                    return;
                case geode::internal::OGS_PROCESS::liquid_flow:
                    add_text( variables, "variable", "pressure" );
                    return;
            }
        }

        pugi::xml_node add_process_variable( pugi::xml_node& variables,
            std::string_view name,
            std::string_view initial_parameter )
        {
            auto variable = variables.append_child( "process_variable" );
            add_text( variable, "name", name );
            add_text( variable, "components", "1" );
            add_text( variable, "order", "1" );
            add_text( variable, "initial_condition", initial_parameter );
            return variable;
        }

        void add_dirichlet_conditions(
            pugi::xml_node& variable, geode::PHYSICAL_PROPERTY_NAME property )
        {
            auto conditions = variable.append_child( "boundary_conditions" );
            for( const auto& boundary : context_.boundaries )
            {
                if( boundary.property != property )
                {
                    continue;
                }
                if( boundary.value )
                {
                    add_constant_parameter(
                        boundary.stem, boundary.value.value() );
                }
                else
                {
                    add_mesh_parameter( { boundary.stem, "MeshNode",
                        boundary.stem, "values" } );
                }
                auto condition =
                    conditions.append_child( "boundary_condition" );
                add_text( condition, "mesh", boundary.stem );
                add_text( condition, "type", "Dirichlet" );
                add_text( condition, "parameter", boundary.stem );
            }
        }

        void add_source_terms( pugi::xml_node& variable )
        {
            auto sources = variable.append_child( "source_terms" );
            for( const auto source_id :
                geode::Indices{ problem_.source_terms() } )
            {
                const auto& source = problem_.source_terms()[source_id];
                const geode::internal::OgsWellMesh* well{ nullptr };
                for( const auto& candidate : context_.wells )
                {
                    if( candidate.name == source.well )
                    {
                        well = &candidate;
                    }
                }
                if( well == nullptr )
                {
                    geode::Logger::warning( "[OgsExporter] Source term on "
                                            "well ",
                        source.well, " is ignored: the well has no vertex." );
                    continue;
                }
                const auto name = absl::StrCat( "source_", source_id );
                add_constant_parameter( name,
                    source.rate
                        / static_cast< double >( well->mesh.points.size() ) );
                auto term = sources.append_child( "source_term" );
                add_text( term, "mesh", well->stem );
                add_text( term, "type", "Nodal" );
                add_text( term, "parameter", name );
            }
        }

        void write_process_variables()
        {
            auto variables = root_.append_child( "process_variables" );
            if( context_.has_permeability_field )
            {
                add_mesh_parameter( { "permeability", "MeshElement",
                    context_.prefix, "permeability" } );
            }
            if( context_.has_porosity_field )
            {
                add_mesh_parameter( { "porosity", "MeshElement",
                    context_.prefix, "porosity" } );
                // Only LIQUID_FLOW and HT read a storage
                if( !is_two_phase() && !is_mechanics() )
                {
                    add_mesh_parameter( { "storage", "MeshElement",
                        context_.prefix, "storage" } );
                }
            }
            switch( context_.process )
            {
                case geode::internal::OGS_PROCESS::two_phase_flow_pp:
                    write_two_phase_variables( variables );
                    return;
                case geode::internal::OGS_PROCESS::ht:
                    write_temperature_variable( variables );
                    write_pressure_variable( variables );
                    return;
                case geode::internal::OGS_PROCESS::small_deformation:
                case geode::internal::OGS_PROCESS::hydro_mechanics:
                case geode::internal::OGS_PROCESS::thermo_hydro_mechanics:
                    write_mechanics_variables( variables );
                    return;
                case geode::internal::OGS_PROCESS::liquid_flow:
                    write_pressure_variable( variables );
                    return;
            }
        }

        void write_two_phase_variables( pugi::xml_node& variables )
        {
            add_initial_parameter( "initial_capillary_pressure",
                context_.initial_capillary_pressure );
            auto capillary = add_process_variable(
                variables, "capillary_pressure", "initial_capillary_pressure" );
            add_dirichlet_conditions( capillary,
                geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction );
            add_initial_parameter(
                "initial_gas_pressure", context_.initial_pressure );
            auto gas = add_process_variable(
                variables, "gas_pressure", "initial_gas_pressure" );
            add_dirichlet_conditions(
                gas, geode::PHYSICAL_PROPERTY_NAME::boundary_pressure );
        }

        void write_temperature_variable( pugi::xml_node& variables )
        {
            add_initial_parameter(
                "initial_temperature", context_.initial_temperature );
            auto temperature = add_process_variable(
                variables, "temperature", "initial_temperature" );
            add_dirichlet_conditions( temperature,
                geode::PHYSICAL_PROPERTY_NAME::boundary_temperature );
            add_heat_source_terms( temperature );
        }

        // Nodal power (W) shared by the vertices of a point heat source,
        // Volumetric power density (W/m3) of a Block heat source. A power
        // table is a curve scaled by 1 / (number of vertices or volume).
        void add_heat_source_terms( pugi::xml_node& variable )
        {
            if( context_.heat_sources.empty() )
            {
                return;
            }
            auto terms = variable.append_child( "source_terms" );
            for( const auto& heat_mesh : context_.heat_sources )
            {
                const auto& source =
                    problem_.heat_sources()[heat_mesh.heat_source];
                const auto divisor =
                    source.block
                        ? heat_mesh.volume
                        : static_cast< double >( heat_mesh.mesh.points.size() );
                const auto name = absl::StrCat( "heat_", source.name );
                if( source.times.empty() )
                {
                    add_constant_parameter(
                        name, source.powers.front() / divisor );
                }
                else
                {
                    add_curve_scaled_parameter( name, 1. / divisor );
                }
                auto term = terms.append_child( "source_term" );
                add_text( term, "mesh", heat_mesh.stem );
                add_text( term, "type", source.block ? "Volumetric" : "Nodal" );
                add_text( term, "parameter", name );
            }
        }

        // The curve of the name times a Constant <name>_scale
        void add_curve_scaled_parameter( std::string_view name, double scale )
        {
            const auto scale_name = absl::StrCat( name, "_scale" );
            add_constant_parameter( scale_name, scale );
            auto node = parameters_.append_child( "parameter" );
            add_text( node, "name", name );
            add_text( node, "type", "CurveScaled" );
            add_text( node, "curve", name );
            add_text( node, "parameter", scale_name );
        }

        // Source terms on pressure only: in HT, injected fluid enters at the
        // local temperature
        void write_pressure_variable( pugi::xml_node& variables )
        {
            add_initial_parameter(
                "initial_pressure", context_.initial_pressure );
            auto pressure = add_process_variable(
                variables, "pressure", "initial_pressure" );
            add_dirichlet_conditions(
                pressure, geode::PHYSICAL_PROPERTY_NAME::boundary_pressure );
            add_source_terms( pressure );
        }

        void write_solvers()
        {
            auto nonlinear = root_.append_child( "nonlinear_solvers" )
                                 .append_child( "nonlinear_solver" );
            add_text( nonlinear, "name", NONLINEAR_SOLVER );
            write_nonlinear_type( nonlinear );
            add_text( nonlinear, "linear_solver", LINEAR_SOLVER );
            auto linear = root_.append_child( "linear_solvers" )
                              .append_child( "linear_solver" );
            add_text( linear, "name", LINEAR_SOLVER );
            auto eigen = linear.append_child( "eigen" );
            add_text( eigen, "solver_type", "SparseLU" );
            add_text( eigen, "scaling", "true" );
        }

        void write_nonlinear_type( pugi::xml_node& nonlinear ) const
        {
            switch( context_.process )
            {
                case geode::internal::OGS_PROCESS::two_phase_flow_pp:
                    add_text( nonlinear, "type", "Newton" );
                    add_text( nonlinear, "max_iter",
                        absl::StrCat( NEWTON_MAX_ITERATIONS ) );
                    return;
                case geode::internal::OGS_PROCESS::ht:
                    add_text( nonlinear, "type", "Picard" );
                    add_text( nonlinear, "max_iter",
                        absl::StrCat( HT_MAX_ITERATIONS ) );
                    return;
                case geode::internal::OGS_PROCESS::small_deformation:
                case geode::internal::OGS_PROCESS::hydro_mechanics:
                case geode::internal::OGS_PROCESS::thermo_hydro_mechanics:
                    add_text( nonlinear, "type", "Newton" );
                    add_text( nonlinear, "max_iter",
                        absl::StrCat( MECHANICS_MAX_ITERATIONS ) );
                    return;
                case geode::internal::OGS_PROCESS::liquid_flow:
                    add_text( nonlinear, "type", "Picard" );
                    add_text( nonlinear, "max_iter",
                        absl::StrCat( PICARD_MAX_ITERATIONS ) );
                    return;
            }
        }

        void add_curve( pugi::xml_node& curves,
            std::string_view name,
            const geode::internal::OgsCurve& curve ) const
        {
            auto node = curves.append_child( "curve" );
            add_text( node, "name", name );
            add_text( node, "coords", reals( curve.coords ) );
            add_text( node, "values", reals( curve.values ) );
        }

        [[nodiscard]] bool has_heat_curves() const
        {
            return absl::c_any_of(
                problem_.heat_sources(), []( const geode::HeatSource& source ) {
                    return !source.times.empty();
                } );
        }

        // Power (W) of each time (s)
        void write_heat_curves( pugi::xml_node& curves ) const
        {
            for( const auto& source : problem_.heat_sources() )
            {
                if( !source.times.empty() )
                {
                    add_curve( curves, absl::StrCat( "heat_", source.name ),
                        { source.times, source.powers } );
                }
            }
        }

        void write_two_phase_curves( pugi::xml_node& curves ) const
        {
            const auto& fluid =
                dynamic_cast< const geode::DeadOilFluid& >( problem_.fluid() );
            const auto oil = geode::internal::oil_curves( fluid );
            add_curve( curves, "oil_density", oil[0] );
            add_curve( curves, "oil_viscosity", oil[1] );
            for( const auto material_id :
                geode::Indices{ problem_.materials() } )
            {
                const auto& material = problem_.materials()[material_id];
                const auto suffix = absl::StrCat( "_", material_id );
                if( material.capillary_pressure )
                {
                    add_curve( curves, absl::StrCat( "saturation", suffix ),
                        geode::internal::saturation_curve(
                            *material.capillary_pressure ) );
                }
                if( material.relative_permeability )
                {
                    const auto relperm =
                        geode::internal::relative_permeability_curves(
                            *material.relative_permeability );
                    add_curve( curves,
                        absl::StrCat( "relative_permeability", suffix ),
                        relperm[0] );
                    add_curve( curves,
                        absl::StrCat(
                            "relative_permeability_nonwetting_phase", suffix ),
                        relperm[1] );
                }
            }
        }

    private:
        const geode::PhysicalProblem& problem_;
        const geode::internal::OgsProjectContext& context_;
        pugi::xml_node& root_;
        pugi::xml_node parameters_;
    };
} // namespace

namespace geode::internal
{
    OgsProjectWriter::OgsProjectWriter(
        const PhysicalProblem& problem, const OgsProjectContext& context )
        : problem_( problem ), context_( context )
    {
    }

    void OgsProjectWriter::write( pugi::xml_document& document ) const
    {
        auto root = document.append_child( "OpenGeoSysProject" );
        OgsProjectBuilder builder{ problem_, context_, root };
        builder.write();
        apply_solver_options( root, problem_, ogs_solver_options(),
            []( pugi::xml_node& node, std::string_view name,
                std::string_view value ) {
                apply_element_text( node, name, value );
            } );
    }
} // namespace geode::internal
