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

#include <absl/container/flat_hash_map.h>

#include <geode/geometry/point.h>

#include <geode/geosciences/private/utils.h>

namespace geode
{
    class BRep;
} // namespace geode

namespace geode
{
    namespace detail
    {
        struct HeaderData
        {
            std::string name{ "unknown" };
        };
        HeaderData read_header( std::ifstream& file );

        std::string read_name( absl::Span< const absl::string_view > tokens );

        void write_header( std::ofstream& file, const HeaderData& data );

        struct CRSData
        {
            std::string name{ "Default" };
            std::array< std::string, 3 > axis_names{ { "X", "Y", "Z" } };
            std::array< std::string, 3 > axis_units{ { "m", "m", "m" } };
            int z_sign{ 1 };
        };
        CRSData read_CRS( std::ifstream& file );

        void write_CRS( std::ofstream& file, const CRSData& data );

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
        PropHeaderData read_prop_header( std::ifstream& file );

        void write_prop_header(
            std::ofstream& file, const PropHeaderData& data );

        struct PropClassHeaderData
        {
            std::string name{ "Default" };
            std::string kind{ "Real Number" };
            std::string unit{ "unitless" };
            bool is_z{ false };
        };
        void write_property_class_header(
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
            std::deque< Point3D > points;
            std::deque< std::array< index_t, 3 > > triangles;
            std::deque< index_t > tface_triangles_offset{ 0 };
            std::deque< index_t > tface_vertices_offset{ 0 };
            std::deque< index_t > bstones;
            std::deque< TSurfBorderData > borders;
        };
        absl::optional< TSurfData > read_tsurf( std::ifstream& file );

        struct RegionSurfaceSide
        {
            absl::flat_hash_map< uuid, bool > universe_surface_sides;
            absl::flat_hash_map< std::pair< uuid, uuid >, bool >
                regions_surface_sides;
        };
        RegionSurfaceSide determine_surface_to_regions_sides(
            const BRep& brep );
    } // namespace detail
} // namespace geode
