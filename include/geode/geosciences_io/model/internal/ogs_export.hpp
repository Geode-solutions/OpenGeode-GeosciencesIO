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
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <absl/container/flat_hash_map.h>

#include <geode/basic/uuid.hpp>

#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/ogs_vtu.hpp>
#include <geode/geosciences_io/model/internal/simulation_export.hpp>
#include <geode/geosciences_io/model/internal/simulation_regions.hpp>
#include <geode/geosciences_io/model/internal/simulation_solid.hpp>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( SolidMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( Surface );
    ALIAS_3D( SolidMesh );
    ALIAS_3D( Surface );
} // namespace geode

namespace geode::internal
{
    /*!
     * Writes an OpenGeoSys project: the bulk mesh with its MaterialIDs, one
     * boundary mesh per group of Surfaces sharing a boundary condition, one
     * point mesh per well and, with a physical problem, the project file.
     */
    template < typename Model >
    class opengeode_geosciencesio_model_api OgsExporterImpl
        : public SimulationExporterImpl< Model >
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( OgsExporterImpl );

    public:
        OgsExporterImpl() = delete;
        /// @param[in] context Prefix of the exception and warning messages
        OgsExporterImpl( std::string_view files_directory,
            const Model& model,
            std::string_view context );
        ~OgsExporterImpl() override;

        void prepare_export();

        void write_files() const;

    protected:
        [[nodiscard]] virtual SimulationRegions create_region_definitions(
            const Model& model ) const = 0;

    private:
        void check_well_name( std::string_view name ) const final;

        void clear_export();

        void check_problem() const;

        void build_bulk();

        void compute_material_ids();

        void transfer_cell_properties();

        void build_boundary_meshes();

        void build_boundary_meshes( PHYSICAL_PROPERTY_NAME property );

        void build_well_meshes();

        [[nodiscard]] const SolidMesh3D& bulk() const;

        // Model whose Blocks are converted into the bulk mesh
        [[nodiscard]] const BRep& solid_model() const;

        [[nodiscard]] std::string bulk_file() const;

        void write_bulk() const;

    private:
        std::string context_;
        SimulationRegions regions_;
        // Copy of a model whose Blocks are all RegularGrids, with HybridSolid
        // Blocks: the model solid conversion needs unstructured Blocks
        std::optional< BRep > grid_model_;
        std::optional< SimulationSolid > solid_;
        // Bulk mesh converted into a HybridSolid, when the merged Blocks are
        // not a TetrahedralSolid nor a HybridSolid
        std::unique_ptr< SolidMesh3D > hybrid_;
        std::vector< std::int32_t > material_ids_;
        std::optional< std::vector< double > > porosity_;
        std::optional< std::vector< std::array< double, 3 > > > permeability_;
        std::vector< const Surface3D* > surfaces_;
        absl::flat_hash_map< uuid, index_t > surface_indices_;
        std::vector< OgsBoundaryMesh > boundary_meshes_;
        std::vector< OgsWellMesh > well_meshes_;
    };
} // namespace geode::internal
