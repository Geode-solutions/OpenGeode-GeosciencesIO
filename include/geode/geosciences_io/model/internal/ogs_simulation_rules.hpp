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

#pragma once

#include <array>
#include <optional>
#include <string_view>
#include <vector>

#include <absl/types/span.h>

#include <geode/basic/types.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/simulation_cell_properties.hpp>
#include <geode/geosciences_io/model/internal/simulation_xml.hpp>

namespace geode
{
    class CapillaryPressure;
    class DeadOilFluid;
    class PhysicalProblem;
    class RelativePermeability;
    struct OutputSettings;
    struct Schedule;
} // namespace geode

namespace geode::internal
{
    inline constexpr std::string_view OGS_SIMULATOR_KEY{ "ogs" };

    /// Number of samples of a curve written from an analytical law
    inline constexpr index_t CURVE_SAMPLES{ 1001 };

    /*!
     * "ogs" solver options: "<role>/<path>" sets the text of the element
     * <path> (created if missing) below each node of the role. Roles:
     * process, nonlinear, linear, time_stepping, convergence, output.
     */
    [[nodiscard]] SimulatorSolverOptions opengeode_geosciencesio_model_api
        ogs_solver_options();

    /// OGS names (materials, wells, fields): non-empty, [A-Za-z0-9_-] only
    [[nodiscard]] const SimulationNaming& opengeode_geosciencesio_model_api
        ogs_naming();

    enum struct OGS_PROCESS
    {
        liquid_flow,
        two_phase_flow_pp
    };

    /*!
     * LIQUID_FLOW for a SinglePhaseFlow of a SinglePhaseFluid,
     * TWOPHASE_FLOW_PP for a MultiphaseFlow of a water-oil DeadOilFluid,
     * nothing otherwise (thermal and solid mechanics included)
     */
    [[nodiscard]] std::optional< OGS_PROCESS >
        opengeode_geosciencesio_model_api ogs_process(
            const PhysicalProblem& problem );

    /*!
     * Temperature (K) shared by every initial condition giving one, else
     * 293.15 K. OGS needs it as LIQUID_FLOW reference_temperature and as
     * TWOPHASE_FLOW_PP temperature, both unused by isothermal flow.
     */
    [[nodiscard]] double opengeode_geosciencesio_model_api
        ogs_reference_temperature( const PhysicalProblem& problem );

    /// Same permeability along x, y and z: TWOPHASE_FLOW_PP reads a scalar
    [[nodiscard]] bool opengeode_geosciencesio_model_api is_isotropic(
        const std::array< double, 3 >& permeability );

    /// OGS storage (1/Pa) of a porosity with a relative pore compressibility
    [[nodiscard]] double opengeode_geosciencesio_model_api ogs_storage(
        double porosity, double pore_compressibility );

    /*!
     * Sorted output times: the explicit times, and start + k * frequency up
     * to the end time (the end time included). The end time alone when
     * neither is given.
     */
    [[nodiscard]] std::vector< double > opengeode_geosciencesio_model_api
        ogs_output_times(
            const Schedule& schedule, const OutputSettings& output );

    /*!
     * Water saturation of a water-oil composition (global fractions in the
     * fluid component order) at a pressure, with the water density of its
     * PVT reference and the oil density surface density / Bo(pressure)
     */
    [[nodiscard]] double opengeode_geosciencesio_model_api water_saturation(
        const DeadOilFluid& fluid,
        COMPOSITION_BASIS basis,
        absl::Span< const double > composition,
        double pressure );

    /// OGS curve: values linearly interpolated between strictly increasing
    /// coords
    struct OgsCurve
    {
        std::vector< double > coords;
        std::vector< double > values;
    };

    /// Oil density (kg/m3) and viscosity (Pa.s) versus oil pressure, at the
    /// PVT table pressures
    [[nodiscard]] std::array< OgsCurve, 2 > opengeode_geosciencesio_model_api
        oil_curves( const DeadOilFluid& fluid );

    /// Interpolates linearly, constant outside the coords
    [[nodiscard]] double opengeode_geosciencesio_model_api curve_value(
        const OgsCurve& curve, double coord );

    /*!
     * Water saturation versus capillary pressure (Pa), from pc = 0. Throws if
     * a table capillary pressure is not strictly decreasing with the water
     * saturation.
     */
    [[nodiscard]] OgsCurve opengeode_geosciencesio_model_api saturation_curve(
        const CapillaryPressure& capillary_pressure );

    /// Capillary pressure (Pa) of a water saturation, through the law
    [[nodiscard]] double opengeode_geosciencesio_model_api capillary_pressure(
        const CapillaryPressure& capillary_pressure, double water_saturation );

    /// Water then oil relative permeability versus the water saturation
    [[nodiscard]] std::array< OgsCurve, 2 > opengeode_geosciencesio_model_api
        relative_permeability_curves(
            const RelativePermeability& relative_permeability );
} // namespace geode::internal
