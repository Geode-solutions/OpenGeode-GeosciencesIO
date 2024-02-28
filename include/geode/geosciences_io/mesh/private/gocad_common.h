/*
 * Copyright (c) 2019 Geode-solutions
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

#include <geode/basic/attribute_manager.h>

#include <geode/geometry/point.h>

#include <geode/geosciences_io/mesh/common.h>

namespace geode
{
    namespace detail
    {
        struct HeaderData
        {
            std::string name{ "unknown" };
        };
        HeaderData opengeode_geosciencesio_mesh_api read_header(
            std::ifstream& file );

        std::string opengeode_geosciencesio_mesh_api read_name(
            absl::Span< const absl::string_view > tokens );

        void opengeode_geosciencesio_mesh_api write_header(
            std::ofstream& file, const HeaderData& data );

        struct CRSData
        {
            std::string name{ "Default" };
            std::array< std::string, 3 > axis_names{ { "X", "Y", "Z" } };
            std::array< std::string, 3 > axis_units{ { "m", "m", "m" } };
            int z_sign{ -1 };
        };
        CRSData opengeode_geosciencesio_mesh_api read_CRS(
            std::ifstream& file );

        void opengeode_geosciencesio_mesh_api write_CRS(
            std::ofstream& file, const CRSData& data );

        struct PropHeaderData
        {
            std::vector< std::string > names;
            std::vector< std::pair< std::string, std::string > >
                prop_legal_ranges;
            std::vector< double > no_data_values;
            std::vector< std::string > property_classes;
            std::vector< std::string > kinds;
            std::vector< std::pair< std::string, std::string > >
                property_subclass;
            std::vector< index_t > esizes;
            std::vector< std::string > units;

            bool empty() const
            {
                return names.empty();
            }
        };
        PropHeaderData opengeode_geosciencesio_mesh_api read_prop_header(
            std::ifstream& file, absl::string_view prefix );

        void opengeode_geosciencesio_mesh_api read_properties(
            const PropHeaderData& properties_header,
            std::vector< std::vector< double > >& attribute_values,
            absl::Span< const absl::string_view > tokens,
            geode::index_t line_properties_position );

        void opengeode_geosciencesio_mesh_api create_attributes(
            const PropHeaderData& attributes_header,
            absl::Span< const std::vector< double > > attributes_values,
            geode::AttributeManager& attribute_manager,
            geode::index_t nb_vertices,
            absl::Span< const geode::index_t > inverse_vertex_mapping );

        void opengeode_geosciencesio_mesh_api write_prop_header(
            std::ofstream& file, const PropHeaderData& data );

        struct PropClassHeaderData
        {
            std::string name{ "Default" };
            std::string kind{ "Real Number" };
            std::string unit{ "unitless" };
            bool is_z{ false };
        };
        void opengeode_geosciencesio_mesh_api write_property_class_header(
            std::ofstream& file, const PropClassHeaderData& data );

        struct TSurfBorderData
        {
            TSurfBorderData( index_t corner_id_in, index_t next_id_in )
                : corner_id( corner_id_in ), next_id( next_id_in )
            {
            }
            index_t corner_id;
            index_t next_id;
        };

        struct TSurfData
        {
            index_t tface_id( index_t vertex_id ) const
            {
                for( const auto i : Range{ 1, tface_vertices_offset.size() } )
                {
                    if( vertex_id < tface_vertices_offset[i] )
                    {
                        return i - 1;
                    }
                }
                return tface_vertices_offset.size() - 1;
            }

            index_t OFFSET_START{ 1 };
            HeaderData header;
            CRSData crs;
            PropHeaderData vertices_properties_header;
            std::deque< Point3D > points;
            std::deque< std::array< index_t, 3 > > triangles;
            std::deque< index_t > tface_triangles_offset{ 0 };
            std::deque< index_t > tface_vertices_offset{ 0 };
            std::deque< index_t > bstones;
            std::deque< TSurfBorderData > borders;
            std::vector< std::vector< double > > vertices_attribute_values;
        };
        absl::optional< TSurfData > opengeode_geosciencesio_mesh_api read_tsurf(
            std::ifstream& file );
    } // namespace detail
} // namespace geode
