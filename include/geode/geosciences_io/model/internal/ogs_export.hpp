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
#include <utility>
#include <vector>

#include <absl/container/flat_hash_map.h>
#include <absl/types/span.h>

#include <geode/basic/uuid.hpp>

#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/ogs_project_writer.hpp>
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
    /// Water fraction in a bulk cell
    struct OgsCellWaterFraction
    {
        index_t cell;
        double water_fraction;
    };

    /// How the warnings name a point of a set: "<point> <id> of <owner>"
    struct OgsPointDescription
    {
        std::string_view point;
        std::string owner;
    };

    /// Mechanical parameters of the project, see OgsProjectContext
    struct OgsMechanicsValues
    {
        OgsCellValue youngs_modulus;
        OgsCellValue poissons_ratio;
        OgsCellValue solid_density;
        OgsCellValue biot_coefficient;
        OgsCellValue thermal_expansivity;
        OgsCellValue solid_heat_capacity;
    };

    /*!
     * Writes an OpenGeoSys project: the bulk mesh with its MaterialIDs, one
     * boundary mesh per group of Surfaces sharing a boundary condition, one
     * point mesh per well and, with a physical problem, one mesh per heat
     * source and the project file.
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

        // Gives every bulk cell the vertex order OGS reads with a positive
        // Jacobian
        void orient_bulk();

        void compute_material_ids();

        void transfer_cell_properties();

        // Name of the solid cell attribute of the property, none if a Block
        // lacks it
        [[nodiscard]] std::optional< std::string > transfer_property(
            PHYSICAL_PROPERTY_NAME property );

        // Process of the problem, none without problem
        [[nodiscard]] std::optional< OGS_PROCESS > process() const;

        // Value of a double Block property on each bulk cell, zero if the
        // model lacks it
        [[nodiscard]] std::vector< double > property_cell_values(
            PHYSICAL_PROPERTY_NAME property );

        // Uniform value, or a new bulk cell array named field
        [[nodiscard]] OgsCellValue cell_value(
            std::string_view field, std::vector< double > values );

        void transfer_mechanics_properties();

        // C_v of the cell material / grain density of the cell
        void compute_solid_heat_capacity(
            absl::Span< const double > grain_densities );

        // OGS initial stress of every cell: written when the model has one,
        // and always (zeros) with flow since OGS starts from it
        void transfer_initial_stress();

        void build_boundary_meshes();

        void build_boundary_meshes( PHYSICAL_PROPERTY_NAME property );

        void add_boundary_mesh( PHYSICAL_PROPERTY_NAME property,
            index_t nb_skipped_polygons,
            OgsBoundaryMesh&& boundary );

        void build_mechanics_boundary_meshes();

        // One Dirichlet mesh per (axis, value) group of uniform Surfaces
        void build_displacement_meshes( const uuid& attribute_id );

        // One mesh per traction component, faces ordered outward
        void build_traction_meshes( const uuid& attribute_id );

        // Per polygon tractions: one mesh per non-zero component
        void build_variable_traction_meshes( const Surface3D& surface,
            absl::Span< const std::array< double, 3 > > tractions );

        // Faces ordered outward, values in the OGS convention
        void build_normal_traction_meshes( const uuid& attribute_id );

        void build_well_meshes();

        [[nodiscard]] OgsWellMesh build_well_mesh(
            index_t well_id, const AABBTree3D& aabb ) const;

        // Bulk vertex nearest to each point in its containing cell, without
        // duplicates; warns about the ignored points outside the model
        [[nodiscard]] OgsMeshFile point_vertices( const PointSet3D& points,
            const AABBTree3D& aabb,
            const OgsPointDescription& description ) const;

        void build_heat_source_meshes();

        [[nodiscard]] OgsHeatSourceMesh build_point_heat_source_mesh(
            index_t heat_source_id ) const;

        [[nodiscard]] OgsHeatSourceMesh build_block_heat_source_mesh(
            index_t heat_source_id ) const;

        void compute_initial_values();

        // Uniform value, or a bulk point array of the mean of the values of
        // the cells around each vertex
        [[nodiscard]] OgsInitialValue initial_value(
            std::string_view field, const std::vector< double >& cell_values );

        [[nodiscard]] SolidMesh3D& bulk_mesh();

        [[nodiscard]] bool is_two_phase() const;

        // Capillary pressure of a water fraction in a bulk cell
        [[nodiscard]] double cell_capillary_pressure(
            const OgsCellWaterFraction& cell_fraction ) const;

        void compute_two_phase_values(
            absl::Span< const index_t > cell_conditions );

        // Replaces the water fractions of the boundary meshes by capillary
        // pressures
        void convert_boundary_water_fractions();

        void convert_water_fractions( OgsBoundaryMesh& boundary ) const;

        // Mean capillary pressure of the polygons around each vertex
        [[nodiscard]] std::vector< double > vertex_capillary_pressures(
            const OgsMeshFile& mesh,
            absl::Span< const double > fractions ) const;

        void write_project() const;

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
        std::optional< std::vector< double > > storage_;
        std::vector< double > cell_pressures_;
        OgsInitialValue initial_pressure_;
        OgsInitialValue initial_capillary_pressure_;
        // HT and THERMO_HYDRO_MECHANICS only
        OgsInitialValue initial_temperature_;
        // Mechanics only
        OgsMechanicsValues mechanics_;
        std::optional< std::vector< OgsStressTensor > > initial_stress_;
        // Bulk cell arrays of the cell values
        std::vector< std::pair< std::string, std::vector< double > > >
            cell_arrays_;
        std::vector< const Surface3D* > surfaces_;
        absl::flat_hash_map< uuid, index_t > surface_indices_;
        std::vector< OgsBoundaryMesh > boundary_meshes_;
        std::vector< OgsWellMesh > well_meshes_;
        std::vector< OgsHeatSourceMesh > heat_source_meshes_;
    };
} // namespace geode::internal
