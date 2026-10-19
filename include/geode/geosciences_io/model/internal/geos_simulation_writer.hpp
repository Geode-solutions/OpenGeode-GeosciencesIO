#pragma once

#include <absl/types/span.h>

#include <geode/basic/pimpl.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace pugi
{
    class xml_node;
} // namespace pugi

namespace geode
{
    class PhysicalProblem;
} // namespace geode

namespace geode::internal
{
    struct GeosDeckContext;
    struct GeosImportedField;
} // namespace geode::internal

namespace geode::internal
{
    /*!
     * Translates a validated PhysicalProblem into the GEOS deck sections that
     * do not depend on the mesh export.
     */
    class GeosSimulationWriter
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( GeosSimulationWriter );

    public:
        GeosSimulationWriter(
            const PhysicalProblem& problem, const GeosDeckContext& context );
        ~GeosSimulationWriter();

        /// Model arrays to write in the VTU file and import in GEOS
        [[nodiscard]] absl::Span< const GeosImportedField >
            imported_fields() const;

        /// Solvers, NumericalMethods, ElementRegions and Constitutive
        void write_physics( pugi::xml_node& problem_node ) const;

        void write_initial_conditions(
            pugi::xml_node& field_specifications ) const;

        /// SourceFlux on the well perforation sets
        void write_source_terms( pugi::xml_node& field_specifications ) const;

        [[nodiscard]] bool has_functions() const;

        void write_functions( pugi::xml_node& functions ) const;

        void write_events_and_outputs( pugi::xml_node& problem_node ) const;

        /// Applies the "geos" solver options to the generated nodes
        void apply_solver_options( pugi::xml_node& problem_node ) const;

    private:
        IMPLEMENTATION_MEMBER( impl_ );
    };
} // namespace geode::internal
