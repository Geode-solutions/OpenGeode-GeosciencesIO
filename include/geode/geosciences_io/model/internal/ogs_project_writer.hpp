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

#include <optional>
#include <string>
#include <string_view>

#include <absl/types/span.h>

#include <geode/geosciences_io/model/internal/ogs_simulation_rules.hpp>
#include <geode/geosciences_io/model/internal/ogs_vtu.hpp>

namespace pugi
{
    class xml_document;
} // namespace pugi

namespace geode
{
    class PhysicalProblem;
} // namespace geode

namespace geode::internal
{
    /// Initial value of a primary variable: uniform, or a bulk point array
    struct OgsInitialValue
    {
        std::optional< double > value;
        /// Name of the bulk mesh point array used without uniform value
        std::string field_name;
    };

    /// What the exporter wrote that the project file refers to
    struct OgsProjectContext
    {
        OGS_PROCESS process{ OGS_PROCESS::liquid_flow };
        /// Bulk mesh file stem, also the output prefix
        std::string prefix;
        /// Bulk mesh cell arrays of the model properties, if written
        bool has_porosity_field{ false };
        bool has_permeability_field{ false };
        /// Boundary meshes; boundary_water_fraction ones carry capillary
        /// pressures
        absl::Span< const OgsBoundaryMesh > boundaries;
        absl::Span< const OgsWellMesh > wells;
        OgsInitialValue initial_pressure;
        /// Two-phase only
        OgsInitialValue initial_capillary_pressure;
        double reference_temperature{ 0 };
    };

    /// Writes the OGS project file of a problem
    class OgsProjectWriter
    {
    public:
        OgsProjectWriter(
            const PhysicalProblem& problem, const OgsProjectContext& context );

        void write( pugi::xml_document& document ) const;

    private:
        const PhysicalProblem& problem_;
        const OgsProjectContext& context_;
    };
} // namespace geode::internal
