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

#include <fstream>
#include <memory>

#include <absl/functional/function_ref.h>
#include <absl/types/span.h>

#include <nlohmann/json.hpp>

#include <geode/basic/percentage.hpp>

#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

namespace geode
{
    class AttributeManager;
} // namespace geode

namespace geode::internal
{
    /*!
     * Reader of Open Mining Format v1 files (binary layout written by the
     * omf Python package 1.x, version string "OMF-v0.9.0").
     * The file starts with a 60 bytes header (magic number, version, project
     * uuid, JSON start), followed by zlib compressed arrays and ends with a
     * JSON index of all the objects keyed by uuid.
     */
    class OMFFile
    {
    public:
        static constexpr auto EXTENSION = "omf";
        static constexpr auto OMF_ELEMENT_ATTRIBUTE_NAME = "omf_element";
        /*!
         * Default bound of the number of values of an array: meshes are
         * limited to NO_ID points or elements of at most 3 values
         */
        static constexpr uint64_t MAX_ARRAY_VALUES = 3 * uint64_t{ NO_ID };

        explicit OMFFile( std::string_view filename );

        /*!
         * Fraction of the project elements of an OMF v1 file having the given
         * class (and geometry class if not empty), 0 if the file is not an
         * OMF v1 file
         */
        [[nodiscard]] static Percentage is_loadable( std::string_view filename,
            std::string_view element_class,
            std::string_view geometry_class );

        /*!
         * Fraction of the project elements of an OMF v1 file that a reader
         * loads, 0 if the file is not an OMF v1 file
         * @param[in] nb_loadable_elements Number of elements the reader loads
         */
        [[nodiscard]] static Percentage is_loadable( std::string_view filename,
            absl::FunctionRef< std::size_t( const OMFFile& ) >
                nb_loadable_elements );

        /*!
         * Number of project elements present in the file
         */
        [[nodiscard]] std::size_t nb_elements() const;

        [[nodiscard]] std::string_view project_name() const;

        /*!
         * Elements of the project with the given class (and geometry class if
         * not empty), in the project order
         */
        [[nodiscard]] std::vector< const nlohmann::json* > elements(
            std::string_view element_class,
            std::string_view geometry_class ) const;

        [[nodiscard]] const nlohmann::json& object(
            const nlohmann::json& uid ) const;

        [[nodiscard]] const nlohmann::json& geometry(
            const nlohmann::json& element ) const;

        [[nodiscard]] static std::string_view name(
            const nlohmann::json& object );

        /*!
         * Name of the element if there is only one, project name otherwise
         */
        [[nodiscard]] std::string_view mesh_name(
            absl::Span< const nlohmann::json* const > elements ) const;

        /*!
         * Read an inline grid tensor (cell sizes along one axis)
         */
        [[nodiscard]] static std::vector< double > read_tensor(
            const nlohmann::json& geometry, std::string_view key );

        /*!
         * Read a grid axis, normalized
         */
        [[nodiscard]] static Vector3D axis( const nlohmann::json& geometry,
            std::string_view key,
            local_index_t default_axis );

        /*!
         * Project origin plus geometry origin
         */
        [[nodiscard]] Vector3D origin( const nlohmann::json& geometry ) const;

        /*!
         * Read the Vector3Array referenced by geometry[key] and translate
         * the points by the geometry origin
         */
        [[nodiscard]] std::vector< Point3D > read_points(
            const nlohmann::json& geometry, std::string_view key ) const;

        /*!
         * Read the Int2Array/Int3Array referenced by geometry[key].
         * Elements referencing a vertex outside [0, nb_vertices) are
         * filled with NO_ID and reported with a warning.
         */
        template < index_t nb_vertices_per_element >
        [[nodiscard]] std::vector<
            std::array< index_t, nb_vertices_per_element > >
            read_elements( const nlohmann::json& geometry,
                std::string_view key,
                index_t nb_vertices ) const;

        /*!
         * Read the ScalarArray referenced by geometry[key], or an empty
         * vector if there is no such array
         */
        [[nodiscard]] std::vector< double > read_optional_scalars(
            const nlohmann::json& geometry,
            std::string_view key,
            uint64_t max_values = MAX_ARRAY_VALUES ) const;

        /*!
         * Create attributes in the manager for all element data defined at
         * the given location.
         * @param[in] nb_values Number of values expected for this location.
         * @param[in] attribute_element Index in the manager of the i-th
         * value, NO_ID to skip this value.
         */
        void read_data( const nlohmann::json& element,
            std::string_view location,
            index_t nb_values,
            AttributeManager& manager,
            absl::FunctionRef< index_t( index_t ) > attribute_element ) const;

        /*!
         * Read a scalar array. Values of compressed arrays beyond max_values
         * are not decompressed and are ignored with a warning.
         */
        [[nodiscard]] std::vector< double > read_scalars(
            const nlohmann::json& array,
            uint64_t max_values = MAX_ARRAY_VALUES ) const;

        [[nodiscard]] std::vector< int64_t > read_integers(
            const nlohmann::json& array,
            uint64_t max_values = MAX_ARRAY_VALUES ) const;

    private:
        /*!
         * Decompress an array block of 8 bytes values directly into T values
         */
        template < typename T >
        [[nodiscard]] std::vector< T > read_binary(
            const nlohmann::json& binary, uint64_t max_values ) const;

    private:
        /*!
         * Parsed file used by is_loadable, cached per thread while the file
         * is unchanged since all OMF readers check the same file
         */
        [[nodiscard]] static std::shared_ptr< const OMFFile > loadability_file(
            std::string_view filename );

    private:
        std::string filename_;
        mutable std::ifstream file_;
        uint64_t json_start_{ 0 };
        nlohmann::json index_;
        const nlohmann::json* project_{ nullptr };
        Vector3D project_origin_;
        bool warn_missing_elements_{ true };
    };

    /*!
     * Set the mesh name and, when several OMF elements are merged into the
     * same mesh, the omf_element attribute giving the element index
     */
    void set_omf_element_attribute( AttributeManager& manager,
        index_t element_id,
        index_t begin,
        index_t end );
} // namespace geode::internal
