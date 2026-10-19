#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <absl/container/flat_hash_map.h>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    class PhysicalProblem;
    struct SourceTerm;
} // namespace geode

namespace geode::internal
{
    struct GeosDeckContext;
} // namespace geode::internal

namespace geode::internal
{
    /// XPath of the generated nodes that "geos" solver options can modify
    [[nodiscard]] const absl::flat_hash_map< std::string, std::string >&
        geos_solver_option_roles();

    /*!
     * Temperature of an isothermal GEOS simulation: the one shared by every
     * initial condition, if any
     */
    [[nodiscard]] std::optional< double > geos_isothermal_temperature(
        const PhysicalProblem& problem );

    /*!
     * GEOS compositional solvers need a temperature on every boundary face
     * set. Returns the temperature to impose when the model has none.
     */
    [[nodiscard]] std::optional< double > geos_default_boundary_temperature(
        const PhysicalProblem& problem, const GeosDeckContext& context );

    /// Name of the GEOS SourceFlux of a source term
    [[nodiscard]] std::string geos_source_term_name(
        const SourceTerm& source_term );

    /*!
     * GEOS-specific checks of a problem on a model (simulator independent
     * checks are done by validate_physical_problem).
     * Returns one message per issue, empty if the deck can be written.
     */
    [[nodiscard]] std::vector< std::string > validate_geos_simulation(
        const PhysicalProblem& problem, const GeosDeckContext& context );
} // namespace geode::internal
