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

// Block physical properties copied on the simulator solid cells
namespace geode::internal
{
    using StressTensor = std::array< double, 6 >;

    // Simulator conventions of the messages and attribute names
    struct SimulationNaming
    {
        // Prefix of the exception and warning messages
        std::string_view context;
        // Simulator name used in the warnings
        std::string_view simulator_name;
        // Whether a solid attribute name can be used by the simulator
        bool ( *is_valid_name )( std::string_view );
        // Characters of an invalid name replaced by '_'
        std::string_view invalid_characters;
    };

    // Factor converting model stresses and normal tractions to the simulator
    // convention (tensions positive or compressions positive)
    [[nodiscard]] double stress_factor(
        const PhysicalProblem& problem, bool simulator_tension_positive );

    // Block values are read one by one: importing one ConstantAttribute per
    // Block onto the solid would keep only the last Block value.
    void transfer_stress_components( const SimulationNaming& naming,
        const BRep& model,
        const ModelToMeshMappings& model2solid,
        const PhysicalProperties::Info& property_info,
        std::string_view name,
        local_index_t first_component,
        double factor,
        VariableAttribute< StressTensor >& tensor );

    // Name of a new solid attribute, valid for the simulator and different
    // from the already exported ones. simulator_field_name is the name of the
    // simulator field reading the attribute: fallback for an empty name and
    // suffix making the name unique.
    [[nodiscard]] std::string unique_attribute_name(
        const SimulationNaming& naming,
        const AttributeManager& manager,
        std::string_view attribute_name,
        std::string_view simulator_field_name );

    // Copies the property attribute of the Blocks on the solid cells: name
    // of the solid attribute, none if a Block lacks the attribute
    [[nodiscard]] std::optional< std::string > transfer_block_attribute(
        const SimulationNaming& naming,
        const BRep& model,
        SolidMesh3D& solid,
        absl::Span< const index_t > block_polyhedra_offsets,
        const PhysicalProperties::Info& property_info,
        std::string_view simulator_field_name );
} // namespace geode::internal
