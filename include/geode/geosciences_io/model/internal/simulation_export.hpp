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
#include <string>
#include <string_view>
#include <vector>

#include <absl/types/span.h>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSet );
    FORWARD_DECLARATION_DIMENSION_CLASS( SolidMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( AABBTree );
    FORWARD_DECLARATION_DIMENSION_CLASS( Point );
    ALIAS_3D( PointSet );
    ALIAS_3D( SolidMesh );
    ALIAS_3D( AABBTree );
    ALIAS_3D( Point );
    class PhysicalProblem;
} // namespace geode

namespace geode::internal
{
    /// Index of the solid polyhedron containing the point, the nearest to the
    /// point among the polyhedra containing it. NO_ID if the point is outside
    /// the solid.
    [[nodiscard]] index_t containing_cell( const SolidMesh3D& solid,
        const AABBTree3D& aabb,
        const Point3D& point );

    /// Simulator-neutral part of a model exporter: output directory and file
    /// prefix, physical problem, well perforations and heat source points.
    template < typename Model >
    class SimulationExporterImpl
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( SimulationExporterImpl );

    public:
        SimulationExporterImpl() = delete;
        /// The directory is created. Its name is used as file prefix.
        /// @param[in] context Prefix of the exception messages
        SimulationExporterImpl( std::string_view files_directory,
            const Model& model,
            std::string_view context );
        virtual ~SimulationExporterImpl();

        /// Writes a full simulation deck instead of the mesh and boundary
        /// conditions only. The problem must outlive the export.
        void set_physical_problem( const PhysicalProblem& problem );

        void add_well_perforations(
            const PointSet3D& perforations, std::string_view name );

        /// Throws if points of this heat source name were already added
        void add_heat_source_points(
            const PointSet3D& points, std::string_view name );

    protected:
        [[nodiscard]] const Model& model() const;
        [[nodiscard]] const PhysicalProblem* problem() const;
        [[nodiscard]] std::string_view files_directory() const;
        [[nodiscard]] std::string_view prefix() const;
        [[nodiscard]] absl::Span< const std::unique_ptr< PointSet3D > >
            well_perforations() const;
        [[nodiscard]] absl::Span< const std::string > well_names() const;
        [[nodiscard]] absl::Span< const std::unique_ptr< PointSet3D > >
            heat_source_points() const;
        [[nodiscard]] absl::Span< const std::string > heat_source_names() const;

    private:
        /// Throws if the simulator cannot use the name for a well
        virtual void check_well_name( std::string_view name ) const = 0;

    private:
        const Model& model_;
        std::string context_;
        std::string files_directory_;
        std::string prefix_;
        const PhysicalProblem* problem_{ nullptr };
        std::vector< std::unique_ptr< PointSet3D > > well_perforations_{};
        std::vector< std::string > well_names_{};
        std::vector< std::unique_ptr< PointSet3D > > heat_source_points_{};
        std::vector< std::string > heat_source_names_{};
    };
} // namespace geode::internal
