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

#include <geode/geosciences_io/model/internal/ogs_simulation_rules.hpp>

#include <algorithm>
#include <cmath>

#include <absl/algorithm/container.h>
#include <absl/strings/ascii.h>

#include <geode/basic/range.hpp>

#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/schedule.hpp>

namespace
{
    constexpr double DEFAULT_REFERENCE_TEMPERATURE{ 293.15 };

    bool is_valid_ogs_name( std::string_view name )
    {
        return !name.empty()
               && absl::c_all_of( name, []( char character ) {
                      return absl::ascii_isalnum(
                                 static_cast< unsigned char >( character ) )
                             || character == '_' || character == '-';
                  } );
    }

    geode::index_t phase_index(
        absl::Span< const geode::PHASE > phases, geode::PHASE phase )
    {
        const auto position = absl::c_find( phases, phase );
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            position != phases.end(), nullptr,
            geode::OpenGeodeException::TYPE::data, "[OgsExporter] Phase ",
            geode::phase_name( phase ), " is missing" );
        return static_cast< geode::index_t >( position - phases.begin() );
    }

    bool is_water_oil( absl::Span< const geode::PHASE > phases )
    {
        return phases.size() == 2
               && absl::c_linear_search( phases, geode::PHASE::water )
               && absl::c_linear_search( phases, geode::PHASE::oil );
    }

    const geode::HydrocarbonPvtTable& oil_table(
        const geode::DeadOilFluid& fluid )
    {
        for( const auto& table : fluid.hydrocarbon_pvt() )
        {
            if( table.phase == geode::PHASE::oil )
            {
                return table;
            }
        }
        throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[OgsExporter] The fluid has no oil PVT table" };
    }

    double clamp_unit( double value )
    {
        return std::clamp( value, 0., 1. );
    }

    // Samples in [0, 1], first 0, last 1
    double sample( geode::index_t index )
    {
        return static_cast< double >( index )
               / static_cast< double >( geode::internal::CURVE_SAMPLES - 1 );
    }

    geode::internal::OgsCurve brooks_corey_saturation(
        const geode::BrooksCoreyCapillaryPressure& law )
    {
        const auto water = phase_index( law.phases(), geode::PHASE::water );
        const auto oil = phase_index( law.phases(), geode::PHASE::oil );
        const auto water_residual = law.residual_saturations()[water];
        const auto max_saturation = 1. - law.residual_saturations()[oil];
        const auto mobile = max_saturation - water_residual;
        const auto entry_pressure = law.entry_pressures()[water];
        const auto inverse_lambda =
            1. / law.pore_size_distribution_indices()[water];
        geode::internal::OgsCurve curve;
        curve.coords.push_back( 0. );
        curve.values.push_back( max_saturation );
        // Effective saturation from 1 down to its first non-zero sample
        for( auto index = geode::internal::CURVE_SAMPLES - 1; index > 0;
            index-- )
        {
            const auto effective = sample( index );
            curve.coords.push_back(
                entry_pressure * std::pow( effective, -inverse_lambda ) );
            curve.values.push_back( water_residual + mobile * effective );
        }
        return curve;
    }

    geode::internal::OgsCurve table_saturation(
        const geode::TableCapillaryPressure& law )
    {
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            law.tables().size() == 1, nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[OgsExporter] A two-phase capillary pressure needs one table" );
        const auto& table = law.tables()[0].pressure;
        geode::internal::OgsCurve curve;
        for( auto index = table.values.size(); index > 0; index-- )
        {
            const auto pressure = table.values[index - 1];
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                curve.coords.empty() || pressure > curve.coords.back(),
                nullptr, geode::OpenGeodeException::TYPE::data,
                "[OgsExporter] The capillary pressure table must strictly "
                "decrease with the water saturation" );
            curve.coords.push_back( pressure );
            curve.values.push_back( table.saturations[index - 1] );
        }
        if( !curve.coords.empty() && curve.coords.front() > 0. )
        {
            curve.coords.insert( curve.coords.begin(), 0. );
            curve.values.insert( curve.values.begin(), curve.values.front() );
        }
        return curve;
    }

    std::array< geode::internal::OgsCurve, 2 > brooks_corey_relperm(
        const geode::BrooksCoreyRelativePermeability& law )
    {
        const auto water = phase_index( law.phases(), geode::PHASE::water );
        const auto oil = phase_index( law.phases(), geode::PHASE::oil );
        const auto water_residual = law.residual_saturations()[water];
        const auto oil_residual = law.residual_saturations()[oil];
        const auto mobile = 1. - water_residual - oil_residual;
        std::array< geode::internal::OgsCurve, 2 > curves;
        for( const auto index : geode::Range{ geode::internal::CURVE_SAMPLES } )
        {
            const auto saturation = sample( index );
            const auto water_effective =
                clamp_unit( ( saturation - water_residual ) / mobile );
            const auto oil_effective =
                clamp_unit( ( 1. - saturation - oil_residual ) / mobile );
            for( auto& curve : curves )
            {
                curve.coords.push_back( saturation );
            }
            curves[0].values.push_back( law.max_values()[water]
                                        * std::pow( water_effective,
                                            law.exponents()[water] ) );
            curves[1].values.push_back(
                law.max_values()[oil]
                * std::pow( oil_effective, law.exponents()[oil] ) );
        }
        return curves;
    }

    std::array< geode::internal::OgsCurve, 2 > table_relperm(
        const geode::TableRelativePermeability& law )
    {
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            law.tables().size() == 1, nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[OgsExporter] A two-phase relative permeability needs one "
            "table" );
        const auto& table = law.tables()[0];
        std::array< geode::internal::OgsCurve, 2 > curves;
        curves[0].coords = table.wetting_kr.saturations;
        curves[0].values = table.wetting_kr.values;
        // The non-wetting table is versus the oil saturation
        const auto& oil = table.non_wetting_kr;
        for( auto index = oil.saturations.size(); index > 0; index-- )
        {
            curves[1].coords.push_back( 1. - oil.saturations[index - 1] );
            curves[1].values.push_back( oil.values[index - 1] );
        }
        return curves;
    }
} // namespace

namespace geode::internal
{
    SimulatorSolverOptions ogs_solver_options()
    {
        static const SolverOptionRoles roles{
            { "process", "/OpenGeoSysProject/processes/process" },
            { "nonlinear",
                "/OpenGeoSysProject/nonlinear_solvers/nonlinear_solver" },
            { "linear",
                "/OpenGeoSysProject/linear_solvers/linear_solver/eigen" },
            { "time_stepping", "/OpenGeoSysProject/time_loop/processes/"
                               "process/time_stepping" },
            { "convergence", "/OpenGeoSysProject/time_loop/processes/"
                             "process/convergence_criterion" },
            { "output", "/OpenGeoSysProject/time_loop/output" },
        };
        return { OGS_SIMULATOR_KEY, "OGS", "[OgsExporter]", roles };
    }

    const SimulationNaming& ogs_naming()
    {
        static const SimulationNaming naming{ "[OgsExporter]", "OGS",
            &is_valid_ogs_name, " /\\:*?\"<>|&'.,;=()[]{}" };
        return naming;
    }

    std::optional< OGS_PROCESS > ogs_process( const PhysicalProblem& problem )
    {
        if( !problem.has_flow() || problem.thermal()
            || problem.solid_mechanics() || !problem.flow().has_fluid() )
        {
            return std::nullopt;
        }
        const auto& flow = problem.flow();
        const auto& fluid = flow.fluid();
        if( dynamic_cast< const SinglePhaseFlow* >( &flow )
            && dynamic_cast< const SinglePhaseFluid* >( &fluid ) )
        {
            return OGS_PROCESS::liquid_flow;
        }
        if( dynamic_cast< const MultiphaseFlow* >( &flow )
            && dynamic_cast< const DeadOilFluid* >( &fluid )
            && is_water_oil( fluid.phases() ) )
        {
            return OGS_PROCESS::two_phase_flow_pp;
        }
        return std::nullopt;
    }

    double ogs_reference_temperature( const PhysicalProblem& problem )
    {
        std::optional< double > temperature;
        for( const auto condition_id :
            Range{ problem.nb_initial_conditions() } )
        {
            const auto* condition =
                dynamic_cast< const UniformInitialCondition* >(
                    &problem.initial_condition( condition_id ) );
            if( condition == nullptr || !condition->temperature() )
            {
                continue;
            }
            const auto value = condition->temperature().value();
            if( temperature && temperature.value() != value )
            {
                return DEFAULT_REFERENCE_TEMPERATURE;
            }
            temperature = value;
        }
        return temperature.value_or( DEFAULT_REFERENCE_TEMPERATURE );
    }

    double ogs_storage( double porosity, double pore_compressibility )
    {
        return porosity * pore_compressibility;
    }

    std::vector< double > ogs_output_times(
        const Schedule& schedule, const OutputSettings& output )
    {
        std::vector< double > times = output.times;
        if( output.frequency && output.frequency.value() > 0. )
        {
            const auto frequency = output.frequency.value();
            const auto nb_outputs = static_cast< index_t >(
                std::floor( ( schedule.end_time - schedule.start_time )
                                / frequency
                            + GLOBAL_EPSILON ) );
            for( const auto index : Range{ 1, nb_outputs + 1 } )
            {
                times.push_back( schedule.start_time + index * frequency );
            }
            if( times.empty() || times.back() < schedule.end_time )
            {
                times.push_back( schedule.end_time );
            }
        }
        if( times.empty() )
        {
            times.push_back( schedule.end_time );
        }
        absl::c_sort( times );
        times.erase( std::unique( times.begin(), times.end() ), times.end() );
        return times;
    }

    double water_saturation( const DeadOilFluid& fluid,
        COMPOSITION_BASIS basis,
        absl::Span< const double > composition,
        double pressure )
    {
        const auto water = phase_index( fluid.phases(), PHASE::water );
        const auto water_component =
            phase_component( fluid, PHASE::water ).value();
        const auto oil_component = phase_component( fluid, PHASE::oil ).value();
        auto water_mass = composition[water_component];
        auto oil_mass = composition[oil_component];
        if( basis == COMPOSITION_BASIS::mole )
        {
            water_mass *= fluid.components()[water_component].molar_weight;
            oil_mass *= fluid.components()[oil_component].molar_weight;
        }
        const auto water_density =
            fluid.surface_densities()[water]
            / fluid.water_pvt().value().formation_volume_factor;
        const auto oil_density =
            curve_value( oil_curves( fluid )[0], pressure );
        const auto water_volume = water_mass / water_density;
        const auto oil_volume = oil_mass / oil_density;
        return water_volume / ( water_volume + oil_volume );
    }

    std::array< OgsCurve, 2 > oil_curves( const DeadOilFluid& fluid )
    {
        const auto oil = phase_index( fluid.phases(), PHASE::oil );
        const auto& table = oil_table( fluid );
        std::array< OgsCurve, 2 > curves;
        for( const auto index : Indices{ table.pressures } )
        {
            for( auto& curve : curves )
            {
                curve.coords.push_back( table.pressures[index] );
            }
            curves[0].values.push_back( fluid.surface_densities()[oil]
                                        / table.formation_volume_factors[index] );
            curves[1].values.push_back( table.viscosities[index] );
        }
        return curves;
    }

    double curve_value( const OgsCurve& curve, double coord )
    {
        if( coord <= curve.coords.front() )
        {
            return curve.values.front();
        }
        for( const auto index : Range{ 1, curve.coords.size() } )
        {
            if( coord <= curve.coords[index] )
            {
                const auto ratio =
                    ( coord - curve.coords[index - 1] )
                    / ( curve.coords[index] - curve.coords[index - 1] );
                return curve.values[index - 1]
                       + ratio
                             * ( curve.values[index]
                                 - curve.values[index - 1] );
            }
        }
        return curve.values.back();
    }

    OgsCurve saturation_curve( const CapillaryPressure& capillary_pressure )
    {
        if( const auto* brooks_corey =
                dynamic_cast< const BrooksCoreyCapillaryPressure* >(
                    &capillary_pressure ) )
        {
            return brooks_corey_saturation( *brooks_corey );
        }
        if( const auto* table = dynamic_cast< const TableCapillaryPressure* >(
                &capillary_pressure ) )
        {
            return table_saturation( *table );
        }
        throw OpenGeodeGeosciencesIOModelException{ nullptr,
            OpenGeodeException::TYPE::data, "[OgsExporter] ",
            capillary_pressure.type_name(), " is not supported by OGS" };
    }

    double capillary_pressure(
        const CapillaryPressure& capillary_pressure, double water_saturation )
    {
        const auto curve = saturation_curve( capillary_pressure );
        if( water_saturation >= curve.values.front() )
        {
            return 0.;
        }
        // Saturations do not increase along the curve
        for( const auto index : Range{ 1, curve.values.size() } )
        {
            if( water_saturation >= curve.values[index] )
            {
                const auto delta =
                    curve.values[index - 1] - curve.values[index];
                if( delta <= 0. )
                {
                    return curve.coords[index];
                }
                const auto ratio =
                    ( curve.values[index - 1] - water_saturation ) / delta;
                return curve.coords[index - 1]
                       + ratio
                             * ( curve.coords[index]
                                 - curve.coords[index - 1] );
            }
        }
        return curve.coords.back();
    }

    std::array< OgsCurve, 2 > relative_permeability_curves(
        const RelativePermeability& relative_permeability )
    {
        if( const auto* brooks_corey =
                dynamic_cast< const BrooksCoreyRelativePermeability* >(
                    &relative_permeability ) )
        {
            return brooks_corey_relperm( *brooks_corey );
        }
        if( const auto* table =
                dynamic_cast< const TableRelativePermeability* >(
                    &relative_permeability ) )
        {
            return table_relperm( *table );
        }
        throw OpenGeodeGeosciencesIOModelException{ nullptr,
            OpenGeodeException::TYPE::data, "[OgsExporter] ",
            relative_permeability.type_name(), " is not supported by OGS" };
    }
} // namespace geode::internal
