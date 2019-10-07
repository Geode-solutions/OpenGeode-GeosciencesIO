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

#include <array>
#include <deque>
#include <string>

#include <geode/basic/point.h>

#include <geode/geosciences/detail/common.h>

namespace geode
{
    bool opengeode_geosciencesio_geosciences_api string_starts_with(
        const std::string& string, const std::string& check );

    void opengeode_geosciencesio_geosciences_api check_keyword(
        std::ifstream& file, const std::string& keyword );

    std::string opengeode_geosciencesio_geosciences_api read_name(
        std::istringstream& iss );

    struct opengeode_geosciencesio_geosciences_api HeaderData
    {
        std::string name = std::string( "unknown" );
    };
    HeaderData opengeode_geosciencesio_geosciences_api read_header(
        std::ifstream& file );

    struct opengeode_geosciencesio_geosciences_api CRSData
    {
        int z_sign{ 1 };
    };
    CRSData opengeode_geosciencesio_geosciences_api read_CRS(
        std::ifstream& file );

    void opengeode_geosciencesio_geosciences_api goto_keyword(
        std::ifstream& file, const std::string& word );

    struct opengeode_geosciencesio_geosciences_api TSurfBorderData
    {
        TSurfBorderData( index_t corner_id_in, index_t next_id_in )
            : corner_id( corner_id_in ), next_id( next_id_in )
        {
        }
        index_t corner_id;
        index_t next_id;
    };

    struct opengeode_geosciencesio_geosciences_api TSurfData
    {
        geode::index_t tface_id( geode::index_t vertex_id ) const
        {
            for( auto i : geode::Range{ 1, tface_vertices_offset.size() } )
            {
                if( vertex_id < tface_vertices_offset[i] )
                {
                    return i - 1;
                }
            }
            return tface_vertices_offset.size() - 1;
        }

        static constexpr index_t OFFSET_START{ 1 };
        HeaderData header;
        CRSData crs;
        std::deque< Point3D > points;
        std::deque< std::array< index_t, 3 > > triangles;
        std::deque< index_t > tface_triangles_offset{ 0 };
        std::deque< index_t > tface_vertices_offset{ 0 };
        std::deque< index_t > bstones;
        std::deque< TSurfBorderData > borders;
    };
    TSurfData opengeode_geosciencesio_geosciences_api read_tsurf(
        std::ifstream& file );
} // namespace geode