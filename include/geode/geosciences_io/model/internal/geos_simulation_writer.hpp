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

#include <string>
#include <string_view>
#include <vector>

#include <absl/types/span.h>

#include <geode/basic/pimpl.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace pugi
{
    class xml_node;
} // namespace pugi

namespace geode
{
    struct OutputSettings;
    class PhysicalProblem;
} // namespace geode

namespace geode::internal
{
    struct GeosDeckContext;
    struct GeosImportedField;
    struct GeosRegionAssignment;
} // namespace geode::internal

namespace geode::internal
{
    struct GeosPorosityGroup
    {
        double compressibility{ 0 };
        double reference_pressure{ 0 };
        double default_porosity{ 0 };
    };

    struct GeosPorosityGroups
    {
        /// Distinct GEOS porosity models
        std::vector< GeosPorosityGroup > groups;
        /// Porosity model of each material
        std::vector< index_t > material_groups;
    };

    /// Materials sharing compressibility and reference pressure share a model
    [[nodiscard]] GeosPorosityGroups geos_porosity_groups(
        const PhysicalProblem& problem, const GeosDeckContext& context );

    /// The deck has a VTK output node
    [[nodiscard]] bool geos_writes_vtk_output( const OutputSettings& output );

    /// The deck has a restart output node
    [[nodiscard]] bool geos_writes_restart_output(
        const OutputSettings& output );

    /// VTU name of the porosity array imported by a porosity model
    [[nodiscard]] std::string geos_porosity_vtu_name(
        std::string_view name, index_t group, index_t nb_groups );

    /*!
     * Translates a validated PhysicalProblem into the GEOS deck sections that
     * do not depend on the mesh export. Every solid region must have a
     * material and an initial condition in the given assignments.
     */
    class GeosSimulationWriter
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( GeosSimulationWriter );

    public:
        GeosSimulationWriter( const PhysicalProblem& problem,
            const GeosDeckContext& context,
            const GeosRegionAssignment& materials,
            const GeosRegionAssignment& initial_conditions );
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
