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

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <geode/basic/types.hpp>

#include <geode/geometry/point.hpp>

#include <geode/model/mixin/core/physical_properties.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace geode::internal
{
    /// Boundary or point mesh read by OGS next to the bulk mesh
    struct OgsMeshFile
    {
        std::vector< Point3D > points;
        /// Polygon (triangle, quad, polygon) or single-vertex cells
        std::vector< std::vector< index_t > > cells;
        /// Bulk mesh vertex of each point
        std::vector< std::uint64_t > bulk_node_ids;
        /// Bulk mesh cell owning each face cell, empty if not needed
        std::vector< std::uint64_t > bulk_element_ids;
        /// Facet of the bulk cell matching each face cell, empty if not
        /// needed
        std::vector< std::uint64_t > bulk_face_ids;
        /// Named point values, one per point
        std::vector< std::pair< std::string, std::vector< double > > >
            point_arrays;
    };

    /// Boundary mesh of a group of Surfaces sharing a boundary property
    struct OgsBoundaryMesh
    {
        PHYSICAL_PROPERTY_NAME property;
        /// File name without extension, used by the project to name it
        std::string stem;
        /// Value of a uniform group, none if the mesh has a "values" array
        std::optional< double > value;
        OgsMeshFile mesh;
    };

    /// Point mesh of the vertices receiving a well source term
    struct OgsWellMesh
    {
        std::string name;
        /// File name without extension, used by the project to name it
        std::string stem;
        OgsMeshFile mesh;
    };

    /*!
     * Writes an ASCII VTU file: bulk_node_ids as UInt64 point data, the
     * non-empty bulk_element_ids and bulk_face_ids as UInt64 cell data and
     * the point arrays as Float64 point data
     */
    void opengeode_geosciencesio_model_api save_ogs_mesh_file(
        std::string_view filename, const OgsMeshFile& mesh );
} // namespace geode::internal
