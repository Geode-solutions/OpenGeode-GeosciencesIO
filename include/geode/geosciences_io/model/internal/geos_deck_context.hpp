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
#include <string>
#include <string_view>
#include <vector>

#include <absl/container/flat_hash_map.h>
#include <absl/types/span.h>

#include <geode/basic/uuid.hpp>

#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/helpers/simulation_report.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/simulation_boundary_values.hpp>
#include <geode/geosciences_io/model/internal/simulation_regions.hpp>

namespace geode
{
    class BRep;
} // namespace geode

namespace geode::internal
{
    /// Cell attribute of the exported mesh
    struct GeosCellArray
    {
        uuid attribute_id;
        std::string name;
    };

    /// Cell array of a PhysicalProperty with its mean value
    struct GeosCellProperty
    {
        GeosCellArray array;
        double mean{ 0 };
    };

    /// Model data needed to translate a PhysicalProblem into a GEOS deck
    struct GeosDeckContext
    {
        SimulationRegions regions;
        /// Cell array and mean value of the porosity PhysicalProperty
        std::optional< GeosCellArray > porosity_array;
        std::optional< double > mean_porosity;
        /// Cell array and mean value of the permeability PhysicalProperty
        std::optional< GeosCellArray > permeability_array;
        std::optional< std::array< double, 3 > > mean_permeability;
        ModelBoundaryConditions boundary_conditions;
        std::optional< GeosCellProperty > bulk_modulus;
        std::optional< GeosCellProperty > shear_modulus;
        /// Bulk density, or grain density with flow
        std::optional< GeosCellProperty > density;
        /// Poromechanics only (flow and solid mechanics)
        std::optional< GeosCellProperty > biot_coefficient;
        std::optional< GeosCellProperty > grain_bulk_modulus;
        /// Drained linear thermal expansion coefficient, thermal
        /// poromechanics only
        std::optional< GeosCellProperty > thermal_expansion;
        /// Initial stress tensor (xx, yy, zz, yz, xz, xy) in the GEOS sign
        /// convention: negative values are compressions. With flow, it is
        /// the stress of the solid: total stress + Biot coefficient x pore
        /// pressure
        std::optional< GeosCellArray > initial_stress;
        /// Names of the well perforation sets given to the exporter
        std::vector< std::string > well_names;
    };

    /// Model array imported in GEOS, possibly copied under a new VTU name
    struct GeosImportedField
    {
        [[nodiscard]] bool is_copy() const
        {
            return vtu_name != source.name;
        }

        GeosCellArray source;
        std::string vtu_name;
        std::string geos_name;
    };
} // namespace geode::internal
