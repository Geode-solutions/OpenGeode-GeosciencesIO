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
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <absl/container/flat_hash_map.h>
#include <absl/types/span.h>

#include <geode/basic/uuid.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/helpers/convert_to_mesh.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>

namespace pugi
{
    class xml_node;
} // namespace pugi

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSet );
    FORWARD_DECLARATION_DIMENSION_CLASS( SolidMesh );
    ALIAS_3D( PointSet );
    ALIAS_3D( SolidMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( Surface );
    ALIAS_3D( Surface );
    struct ModelToMeshMappings;
    class AttributeManager;
    class PhysicalProblem;
} // namespace geode

namespace geode::internal
{
    class GeosSimulationWriter;
} // namespace geode::internal

namespace geode::internal
{
    template < typename Model >
    class GeosExporterImpl
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( GeosExporterImpl );

    public:
        GeosExporterImpl() = delete;
        GeosExporterImpl(
            std::string_view files_directory, const Model& model );
        virtual ~GeosExporterImpl();

        void prepare_export();
        void write_files() const;

        void add_well_perforations(
            const PointSet3D& perforations, std::string_view name );

        /// Writes a full simulation deck instead of the mesh and boundary
        /// conditions only. The problem must outlive the export.
        void set_physical_problem( const PhysicalProblem& problem );

    protected:
        std::string_view files_directory() const;
        std::string_view prefix() const;

        void check_export();
        void convert_model();
        index_t initialize_solid_region_attribute();
        void initialize_surface_cells( index_t first_surface_region_id );

        void write_well_perforations_boxes( pugi::xml_node& root ) const;
        void write_mesh_files( pugi::xml_node& root ) const;
        void write_boundary_conditions( pugi::xml_node& root ) const;
        void write_boundary_functions( pugi::xml_node& root ) const;

        void transfer_physical_properties();
        void transfer_boundary_conditions();
        void transfer_boundary_condition( const uuid& attribute_id,
            std::string_view field_name,
            std::string_view name_prefix,
            std::optional< local_index_t > component );

        std::string write_solid_file() const;
        void duplicate_imported_fields(
            AttributeManager& cell_attributes ) const;
        void prepare_simulation();

        void add_isothermal_boundary_temperature( double temperature );
        void compute_mean_rock_properties();
        void write_simulation_sections( pugi::xml_node& problem_node ) const;
        void write_well_perforation_file() const;
        void write_boundary_table_files() const;

    private:
        struct BoundaryCondition
        {
            std::string name;
            std::string field_name;
            std::optional< local_index_t > component;
            std::vector< index_t > region_ids;
            double scale{ 0. };
            std::optional< std::string > function_name;
        };

        // GEOS TableFunction grid of a Surface, built from its polygon centers
        // so that each polygon is looked up on its own node
        struct BoundaryGrid
        {
            std::array< std::vector< double >, 3 > coordinates;
            std::array< std::string, 3 > files;
            // Surface polygon whose value is stored at each node, in Fortran
            // order
            std::vector< index_t > node_polygons;
        };

        struct BoundaryTable
        {
            std::string function_name;
            uuid surface_id;
            std::string voxel_file;
            std::vector< double > values;
        };

        const BoundaryGrid& boundary_grid( const Surface3D& surface );
        std::string add_boundary_table( const Surface3D& surface,
            absl::Span< const double > values,
            std::string_view name );

    private:
        const Model& model_;
        std::unique_ptr< SolidMesh3D > model_solid_{};
        ModelToMeshMappings model2solid_;

        uuid region_attribute_id_;
        std::shared_ptr< VariableAttribute< index_t > > region_attribute_{};

        std::vector< PolygonVertices > surface_cells_{};
        std::vector< index_t > surface_cells_region_{};
        absl::flat_hash_map< uuid, index_t > surface_regions_{};

        std::string files_directory_;
        std::string prefix_;

        std::vector< std::pair< std::string, std::string > > imported_fields_{};
        std::vector< BoundaryCondition > boundary_conditions_{};
        absl::flat_hash_map< uuid, BoundaryGrid > boundary_grids_{};
        std::vector< BoundaryTable > boundary_tables_{};

        std::vector< std::unique_ptr< PointSet3D > > well_perforations_{};
        std::vector< std::string > well_names_{};

        const PhysicalProblem* problem_{ nullptr };
        GeosDeckContext context_{};
        std::optional< GeosRegionAssignment > materials_{};
        std::optional< GeosRegionAssignment > initial_conditions_{};
        std::unique_ptr< GeosSimulationWriter > writer_{};
    };
} // namespace geode::internal