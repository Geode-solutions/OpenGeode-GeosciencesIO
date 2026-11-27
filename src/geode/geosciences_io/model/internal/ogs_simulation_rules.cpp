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

#include <geode/mesh/core/solid_mesh.hpp>

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
        return !name.empty() && absl::c_all_of( name, []( char character ) {
            return absl::ascii_isalnum(
                       static_cast< unsigned char >( character ) )
                   || character == '_' || character == '-';
        } );
    }

    geode::index_t phase_index(
        absl::Span< const geode::PHASE > phases, geode::PHASE phase )
    {
        const auto* const position = absl::c_find( phases, phase );
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

    using Vector3 = std::array< double, 3 >;

    // Gauss points of 2 points per axis on [-1, 1], weight 1
    constexpr double GAUSS_POINT{ 0.57735026918962576 };

    double determinant( const std::array< Vector3, 3 >& columns )
    {
        const auto& [a, b, c] = columns;
        return a[0] * ( b[1] * c[2] - b[2] * c[1] )
               - b[0] * ( a[1] * c[2] - a[2] * c[1] )
               + c[0] * ( a[1] * b[2] - a[2] * b[1] );
    }

    Vector3 edge( const geode::Point3D& start, const geode::Point3D& end )
    {
        return { end.value( 0 ) - start.value( 0 ),
            end.value( 1 ) - start.value( 1 ),
            end.value( 2 ) - start.value( 2 ) };
    }

    // Jacobian determinant of the map sum_i N_i x_i, given the derivatives
    // of the shape functions along the 3 reference axes
    double jacobian( absl::Span< const geode::Point3D > points,
        absl::Span< const Vector3 > derivatives )
    {
        std::array< Vector3, 3 > columns{};
        for( const auto vertex : geode::Indices{ points } )
        {
            for( const auto axis : geode::LRange{ 3 } )
            {
                for( const auto coordinate : geode::LRange{ 3 } )
                {
                    columns[axis][coordinate] +=
                        derivatives[vertex][axis]
                        * points[vertex].value( coordinate );
                }
            }
        }
        return determinant( columns );
    }

    // Reference corners of the VTK hexahedron
    constexpr std::array< Vector3, 8 > HEXAHEDRON_CORNERS{ { { -1, -1, -1 },
        { 1, -1, -1 }, { 1, 1, -1 }, { -1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 },
        { 1, 1, 1 }, { -1, 1, 1 } } };

    // Derivatives of the trilinear shape function of the corner at the
    // reference point
    Vector3 trilinear_derivatives(
        const Vector3& corner, const Vector3& reference )
    {
        const Vector3 factors{ 1. + reference[0] * corner[0],
            1. + reference[1] * corner[1], 1. + reference[2] * corner[2] };
        return { corner[0] * factors[1] * factors[2] / 8.,
            corner[1] * factors[0] * factors[2] / 8.,
            corner[2] * factors[0] * factors[1] / 8. };
    }

    // Trilinear hexahedron, VTK order, 2 x 2 x 2 Gauss points (exact); a
    // pyramid is the hexahedron whose top vertices are its apex
    double hexahedron_volume( absl::Span< const geode::Point3D > points )
    {
        double volume{ 0 };
        for( const auto& gauss : HEXAHEDRON_CORNERS )
        {
            const Vector3 reference{ gauss[0] * GAUSS_POINT,
                gauss[1] * GAUSS_POINT, gauss[2] * GAUSS_POINT };
            std::array< Vector3, 8 > derivatives{};
            for( const auto vertex : geode::Indices{ HEXAHEDRON_CORNERS } )
            {
                derivatives[vertex] = trilinear_derivatives(
                    HEXAHEDRON_CORNERS[vertex], reference );
            }
            volume += jacobian( points, derivatives );
        }
        return std::fabs( volume );
    }

    // Triangle (0 1 2) times segment (bottom, top), VTK order: exact with 3
    // triangle points (weight 1/6) times 2 Gauss points
    double prism_volume( absl::Span< const geode::Point3D > points )
    {
        static constexpr std::array< std::array< double, 2 >, 3 > TRIANGLE{
            { { 1. / 6., 1. / 6. }, { 2. / 3., 1. / 6. }, { 1. / 6., 2. / 3. } }
        };
        double volume{ 0 };
        for( const auto& [xi, eta] : TRIANGLE )
        {
            for( const auto zeta : { -GAUSS_POINT, GAUSS_POINT } )
            {
                const std::array< double, 3 > area{ 1. - xi - eta, xi, eta };
                const Vector3 d_xi{ -1., 1., 0. };
                const Vector3 d_eta{ -1., 0., 1. };
                std::array< Vector3, 6 > derivatives{};
                for( const auto vertex : geode::LRange{ 3 } )
                {
                    for( const auto& [level, sign] :
                        { std::pair< geode::local_index_t, double >{ 0, -1. },
                            { 3, 1. } } )
                    {
                        const auto height = ( 1. + sign * zeta ) / 2.;
                        derivatives[vertex + level] = { d_xi[vertex] * height,
                            d_eta[vertex] * height, sign * area[vertex] / 2. };
                    }
                }
                volume += jacobian( points, derivatives ) / 6.;
            }
        }
        return std::fabs( volume );
    }

    bool is_single_phase( const geode::FlowModel& flow )
    {
        return dynamic_cast< const geode::SinglePhaseFlow* >( &flow ) != nullptr
               && dynamic_cast< const geode::SinglePhaseFluid* >(
                      &flow.fluid() )
                      != nullptr;
    }

    bool is_water_oil_flow( const geode::FlowModel& flow )
    {
        return dynamic_cast< const geode::MultiphaseFlow* >( &flow ) != nullptr
               && dynamic_cast< const geode::DeadOilFluid* >( &flow.fluid() )
                      != nullptr
               && is_water_oil( flow.fluid().phases() );
    }

    // Process of a SinglePhaseFlow of a SinglePhaseFluid
    geode::internal::OGS_PROCESS single_phase_process(
        const geode::PhysicalProblem& problem )
    {
        using geode::internal::OGS_PROCESS;
        const auto is_thermal = problem.thermal().has_value();
        if( problem.solid_mechanics() )
        {
            return is_thermal ? OGS_PROCESS::thermo_hydro_mechanics
                              : OGS_PROCESS::hydro_mechanics;
        }
        return is_thermal ? OGS_PROCESS::ht : OGS_PROCESS::liquid_flow;
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

    const geode::SinglePhaseFluid& single_phase_fluid(
        const geode::PhysicalProblem& problem )
    {
        const auto* fluid =
            problem.has_flow() && problem.flow().has_fluid()
                ? dynamic_cast< const geode::SinglePhaseFluid* >(
                      &problem.flow().fluid() )
                : nullptr;
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            fluid != nullptr, nullptr,
            geode::OpenGeodeException::TYPE::internal,
            "[OgsExporter] Fluid laws need a SinglePhaseFluid" );
        return *fluid;
    }

    const geode::FluidThermalProperties& fluid_thermal(
        const geode::SinglePhaseFluid& fluid )
    {
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            fluid.thermal().has_value(), nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[OgsExporter] The fluid needs thermal properties for a thermal "
            "problem" );
        return fluid.thermal().value();
    }

    // Fluid property following an exponential law of pressure and
    // temperature
    struct OgsFluidProperty
    {
        double reference_value{ 0 };
        double pressure_factor{ 0 };
        // Thermal coefficient, minus the temperature factor of the law
        double geode::FluidThermalProperties::* temperature_coefficient{
            nullptr
        };
    };

    // Law of the property, at the fluid reference pressure and the reference
    // temperature of the problem, with no temperature factor when the
    // problem is not thermal
    geode::internal::OgsExponentialLaw ogs_fluid_law(
        const geode::PhysicalProblem& problem,
        const geode::SinglePhaseFluid& fluid,
        const OgsFluidProperty& property )
    {
        geode::internal::OgsExponentialLaw law;
        law.reference_value = property.reference_value;
        law.pressure_factor = property.pressure_factor;
        law.reference_pressure = fluid.reference_pressure();
        law.reference_temperature =
            geode::internal::ogs_reference_temperature( problem );
        if( problem.thermal() )
        {
            law.temperature_factor =
                -( fluid_thermal( fluid ).*property.temperature_coefficient );
        }
        return law;
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
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            entry_pressure > 0., nullptr, geode::OpenGeodeException::TYPE::data,
            "[OgsExporter] The water entry pressure of a Brooks-Corey "
            "capillary pressure must be positive: OGS needs strictly "
            "increasing capillary pressures" );
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
            curve.values.push_back(
                max_saturation - mobile * ( 1. - effective ) );
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
                curve.coords.empty() || pressure > curve.coords.back(), nullptr,
                geode::OpenGeodeException::TYPE::data,
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
            curves[0].values.push_back(
                law.max_values()[water]
                * std::pow( water_effective, law.exponents()[water] ) );
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
        const auto has_mechanics = problem.solid_mechanics().has_value();
        if( has_mechanics && !problem.has_flow() )
        {
            return OGS_PROCESS::small_deformation;
        }
        if( !problem.has_flow() || !problem.flow().has_fluid() )
        {
            return std::nullopt;
        }
        const auto& flow = problem.flow();
        if( is_single_phase( flow ) )
        {
            return single_phase_process( problem );
        }
        if( problem.thermal() || has_mechanics || !is_water_oil_flow( flow ) )
        {
            return std::nullopt;
        }
        return OGS_PROCESS::two_phase_flow_pp;
    }

    bool has_temperature( OGS_PROCESS process )
    {
        return process == OGS_PROCESS::ht
               || process == OGS_PROCESS::thermo_hydro_mechanics;
    }

    bool has_displacement( OGS_PROCESS process )
    {
        return process == OGS_PROCESS::small_deformation
               || process == OGS_PROCESS::hydro_mechanics
               || process == OGS_PROCESS::thermo_hydro_mechanics;
    }

    OgsElasticity ogs_elasticity( const OgsModuli& moduli )
    {
        const auto denominator = 3. * moduli.bulk + moduli.shear;
        return { 9. * moduli.bulk * moduli.shear / denominator,
            ( 3. * moduli.bulk - 2. * moduli.shear ) / ( 2. * denominator ) };
    }

    OgsStressTensor ogs_stress_tensor(
        const StressTensor& tensor, double factor )
    {
        return { factor * tensor[0], factor * tensor[1], factor * tensor[2],
            factor * tensor[5], factor * tensor[3], factor * tensor[4] };
    }

    double ogs_reference_temperature( const PhysicalProblem& problem )
    {
        if( problem.thermal() )
        {
            return problem.thermal()->reference_temperature;
        }
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

    OgsExponentialLaw ogs_density_law( const PhysicalProblem& problem )
    {
        const auto& fluid = single_phase_fluid( problem );
        return ogs_fluid_law( problem, fluid,
            { fluid.reference_density(), fluid.compressibility(),
                &FluidThermalProperties::thermal_expansion } );
    }

    OgsExponentialLaw ogs_viscosity_law( const PhysicalProblem& problem )
    {
        const auto& fluid = single_phase_fluid( problem );
        return ogs_fluid_law( problem, fluid,
            { fluid.reference_viscosity(), fluid.viscosibility(),
                &FluidThermalProperties::viscosity_temperature_coefficient } );
    }

    double ogs_law_value(
        const OgsExponentialLaw& law, const OgsFluidState& state )
    {
        return law.reference_value
               * std::exp(
                   law.pressure_factor
                       * ( state.pressure - law.reference_pressure )
                   + law.temperature_factor
                         * ( state.temperature - law.reference_temperature ) );
    }

    bool is_ogs_polyhedron( const SolidMesh3D& mesh, index_t polyhedron )
    {
        constexpr local_index_t NB_PYRAMID_FACETS{ 5 };
        constexpr local_index_t NB_PRISM_FACETS{ 5 };
        constexpr local_index_t NB_HEXAHEDRON_FACETS{ 6 };
        const auto nb_vertices = mesh.nb_polyhedron_vertices( polyhedron );
        const auto nb_facets = mesh.nb_polyhedron_facets( polyhedron );
        return ( nb_vertices == NB_TETRAHEDRON_VERTICES && nb_facets == 4 )
               || ( nb_vertices == NB_PYRAMID_VERTICES
                    && nb_facets == NB_PYRAMID_FACETS )
               || ( nb_vertices == NB_PRISM_VERTICES
                    && nb_facets == NB_PRISM_FACETS )
               || ( nb_vertices == NB_HEXAHEDRON_VERTICES
                    && nb_facets == NB_HEXAHEDRON_FACETS );
    }

    double ogs_cell_volume( const SolidMesh3D& mesh, index_t polyhedron )
    {
        std::vector< Point3D > points;
        for( const auto vertex : mesh.polyhedron_vertices( polyhedron ) )
        {
            points.push_back( mesh.point( vertex ) );
        }
        switch( points.size() )
        {
            case NB_TETRAHEDRON_VERTICES:
                return std::fabs( determinant( { edge( points[0], points[1] ),
                           edge( points[0], points[2] ),
                           edge( points[0], points[3] ) } ) )
                       / 6.;
            case NB_PYRAMID_VERTICES:
                points.insert( points.end(), 3, points.back() );
                return hexahedron_volume( points );
            case NB_PRISM_VERTICES:
                return prism_volume( points );
            default:
                return hexahedron_volume( points );
        }
    }

    bool is_isotropic( const std::array< double, 3 >& permeability )
    {
        return permeability[1] == permeability[0]
               && permeability[2] == permeability[0];
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
            const auto nb_outputs = static_cast< index_t >( std::floor(
                ( schedule.end_time - schedule.start_time ) / frequency
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
            curves[0].values.push_back(
                fluid.surface_densities()[oil]
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
            // End of the maximum saturation plateau: the entry pressure
            index_t last{ 0 };
            while( last + 1 < curve.values.size()
                   && curve.values[last + 1] >= curve.values.front() )
            {
                last++;
            }
            return curve.coords[last];
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
