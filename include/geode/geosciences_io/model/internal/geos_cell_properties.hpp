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

#include <absl/types/span.h>

#include <geode/basic/variable_attribute.hpp>

#include <geode/model/mixin/core/physical_properties.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    class AttributeManager;
    class BRep;
    class PhysicalProblem;
    FORWARD_DECLARATION_DIMENSION_CLASS( SolidMesh );
    ALIAS_3D( SolidMesh );
    struct ModelToMeshMappings;
} // namespace geode

// Block physical properties copied on the GEOS solid cells
namespace geode::internal
{
    using StressTensor = std::array< double, 6 >;

    // Factor converting model stresses and normal tractions to GEOS, which
    // counts tensions positively
    [[nodiscard]] double geos_stress_factor( const PhysicalProblem& problem );

    // Block values are read one by one: importing one ConstantAttribute per
    // Block onto the solid would keep only the last Block value.
    void transfer_stress_components( const BRep& model,
        const ModelToMeshMappings& model2solid,
        const PhysicalProperties::Info& property_info,
        std::string_view name,
        local_index_t first_component,
        double factor,
        VariableAttribute< StressTensor >& tensor );

    // Name of a new solid attribute, valid in GEOS lists and different from
    // the already exported ones
    [[nodiscard]] std::string unique_attribute_name(
        const AttributeManager& manager,
        std::string_view attribute_name,
        std::string_view geos_name );

    // Copies the property attribute of the Blocks on the solid cells: name
    // of the solid attribute, none if a Block lacks the attribute
    [[nodiscard]] std::optional< std::string > transfer_block_attribute(
        const BRep& model,
        SolidMesh3D& solid,
        absl::Span< const index_t > block_polyhedra_offsets,
        const PhysicalProperties::Info& property_info,
        std::string_view geos_name );
} // namespace geode::internal
