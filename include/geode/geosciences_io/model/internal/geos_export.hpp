/*
 * Copyright (c) 2019 - 2024 Geode-solutions
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

#include <geode/basic/attribute.hpp>

#include <geode/model/helpers/convert_to_mesh.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace pugi
{
    class xml_node;
} // namespace pugi

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSet );
    FORWARD_DECLARATION_DIMENSION_CLASS( EdgedCurve );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( SolidMesh );
    ALIAS_3D( PointSet );
    ALIAS_3D( EdgedCurve );
    ALIAS_3D( SurfaceMesh );
    ALIAS_3D( SolidMesh );
    struct uuid;
    struct ModelToMeshMappings;
} // namespace geode

namespace geode
{
    namespace internal
    {
        template < typename Model >
        class GeosExporterImpl
        {
            OPENGEODE_DISABLE_COPY_AND_MOVE( GeosExporterImpl );

        public:
            GeosExporterImpl() = delete;
            GeosExporterImpl(
                std::string_view files_directory, const Model& model );
            virtual ~GeosExporterImpl() = default;

            void prepare_export();
            void write_files() const;

            void add_well_perforations( const PointSet3D& perforations );
            void add_cell_property1d( std::string_view property_name );
            void add_cell_property2d( std::string_view property_name );
            void add_cell_property3d( std::string_view property_name );

        protected:
            std::string_view files_directory() const;
            std::string_view prefix() const;

            index_t initialize_solid_region_attribute();
            virtual absl::flat_hash_map< uuid, index_t >
                create_region_attribute_map( const Model& model ) const = 0;

            void write_well_perforations_boxes( pugi::xml_node& root ) const;
            void write_mesh_files( pugi::xml_node& root ) const;

            bool check_property_name( std::string_view property_name ) const;
            void transfer_cell_properties();
            void delete_mapping_attributes();

            std::string write_solid_file() const;
            void write_well_perforation_file() const;

        private:
            const Model& model_;
            std::unique_ptr< EdgedCurve3D > model_curve_{};
            std::unique_ptr< SurfaceMesh3D > model_surface_{};
            std::unique_ptr< SolidMesh3D > model_solid_{};
            ModelToMeshMappings model2solid_;

            std::shared_ptr< VariableAttribute< index_t > > region_attribute_{};

            DEBUG_CONST std::string files_directory_;
            DEBUG_CONST std::string prefix_;

            std::vector< std::string > cell_1Dproperty_names_{};
            std::vector< std::string > cell_2Dproperty_names_{};
            std::vector< std::string > cell_3Dproperty_names_{};

            std::vector< std::unique_ptr< PointSet3D > > well_perforations_{};
        };
    } // namespace internal
} // namespace geode