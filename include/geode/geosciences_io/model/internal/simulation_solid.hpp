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

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <absl/types/span.h>

#include <geode/basic/uuid.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/model/helpers/convert_to_mesh.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/simulation_cell_properties.hpp>

namespace geode
{
    class BRep;
    FORWARD_DECLARATION_DIMENSION_CLASS( SolidMesh );
    ALIAS_3D( SolidMesh );
} // namespace geode

namespace geode::internal
{
    struct SimulationRegions;

    /// Solid mesh of the Blocks of a model, with the simulator region of each
    /// polyhedron stored in a polyhedron attribute
    class opengeode_geosciencesio_model_api SimulationSolid
    {
    public:
        /// Converts the model Blocks into a solid and fills the region
        /// attribute from the regions of the Blocks.
        /// @param[in] context Prefix of the exception messages
        SimulationSolid( const BRep& model,
            const SimulationRegions& regions,
            std::string_view context );
        ~SimulationSolid();
        OPENGEODE_DISABLE_COPY_AND_MOVE( SimulationSolid );

        [[nodiscard]] SolidMesh3D& mesh();
        [[nodiscard]] const SolidMesh3D& mesh() const;
        [[nodiscard]] const ModelToMeshMappings& mappings() const;

        /// First solid polyhedron of each Block, in the model Block order,
        /// followed by the number of solid polyhedra
        [[nodiscard]] absl::Span< const index_t >
            block_polyhedra_offsets() const;

        /// Id of the region attribute in the solid polyhedron attributes
        [[nodiscard]] const uuid& region_attribute_id() const;

        /// Region of the polyhedron (index of its Block region)
        [[nodiscard]] index_t region( index_t polyhedron ) const;

        /// Copies the property attribute of the Blocks on the solid cells:
        /// name of the solid attribute, none if a Block lacks the attribute.
        /// simulator_field_name is the name of the simulator field reading
        /// the attribute.
        [[nodiscard]] std::optional< std::string > transfer_block_property(
            const SimulationNaming& naming,
            const BRep& model,
            const PhysicalProperties::Info& info,
            std::string_view simulator_field_name );

    private:
        std::unique_ptr< SolidMesh3D > mesh_;
        ModelToMeshMappings mappings_;
        std::vector< index_t > block_polyhedra_offsets_;
        uuid region_attribute_id_;
        std::shared_ptr< VariableAttribute< index_t > > region_attribute_;
    };
} // namespace geode::internal
