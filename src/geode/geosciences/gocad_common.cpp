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

#include <geode/geosciences/detail/gocad_common.h>

#include <fstream>

namespace
{
void read_tfaces( std::ifstream& file, geode::TSurfData& tsurf )
        {
            geode::goto_keyword( file, "TFACE" );
            std::string line;
            while( std::getline( file, line ) )
            {
                std::istringstream iss{ line };
                std::string keyword;
                iss >> keyword;

                if( keyword == "VRTX" || keyword == "PVRTX" )
                {
                    geode::index_t dummy;
                    double x, y, z;
                    iss >> dummy >> x >> y >> z;
                    tsurf.points.push_back( 
                        geode::Point3D{ { x, y, tsurf.crs.z_sign * z } } );
                }
                else if( keyword == "ATOM" || keyword == "PATOM" )
                {
                    geode::index_t dummy, atom_id;
                    iss >> dummy >> atom_id;
                    tsurf.points.push_back( tsurf.points.at( atom_id - tsurf.OFFSET_START ) );
                }
                else if( keyword == "TRGL" )
                {
                    geode::index_t v0, v1, v2;
                    iss >> v0 >> v1 >> v2;
                    tsurf.triangles.push_back( { v0 - tsurf.OFFSET_START,
                        v1 - tsurf.OFFSET_START, v2 - tsurf.OFFSET_START } );
                }
                else if( keyword == "BSTONE" )
                {
                    geode::index_t bstone_id;
                    iss >> bstone_id;
                    tsurf.bstones.push_back( bstone_id - tsurf.OFFSET_START );
                }
                else if( keyword == "BORDER" )
                {
                    geode::index_t dummy, corner, next;
                    iss >> dummy >> corner >> next;
                    tsurf.borders.emplace_back( corner - tsurf.OFFSET_START, next - tsurf.OFFSET_START );
                }
                else if( keyword == "TFACE" )
                {
                    tsurf.tface_triangles_offset.push_back( tsurf.triangles.size() );
                    tsurf.tface_vertices_offset.push_back( tsurf.points.size() );
                }
                else if( keyword == "END" )
                {
                    tsurf.tface_triangles_offset.push_back( tsurf.triangles.size() );
                    tsurf.tface_vertices_offset.push_back( tsurf.points.size() );
                    return;
                }
            }
        throw geode::OpenGeodeException(
            "Cannot find the end of TSurf section" );
        }
} // namespace

namespace geode
{
    bool string_starts_with(
        const std::string& string, const std::string& check )
    {
        return !string.compare( 0, check.length(), check );
    }

    void check_keyword( std::ifstream& file, const std::string& keyword )
    {
        std::string line;
        std::getline( file, line );
        OPENGEODE_EXCEPTION( string_starts_with( line, keyword ),
            "Line should starts with \"" + keyword + "\"" );
    }

    HeaderData read_header( std::ifstream& file )
    {
        check_keyword( file, "HEADER" );
        HeaderData header;
        std::string line;
        while( std::getline( file, line ) )
        {
            if( string_starts_with( line, "}" ) )
            {
                return header;
            }
            std::istringstream iss{ line };
            std::string keyword;
            iss >> keyword;
            if( keyword == "name:" )
            {
                header.name = read_name( iss );
            }
        }
        throw geode::OpenGeodeException(
            "Cannot find the end of \"HEADER\" section" );
    }

    CRSData read_CRS( std::ifstream& file )
    {
        check_keyword( file, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" );
        CRSData crs;
        std::string line;
        while( std::getline( file, line ) )
        {
            if( string_starts_with(
                    line, "END_ORIGINAL_COORDINATE_SYSTEM" ) )
            {
                return crs;
            }
            std::istringstream iss{ line };
            std::string keyword;
            iss >> keyword;
            if( keyword == "ZPOSITIVE" )
            {
                std::string sign;
                iss >> sign;
                crs.z_sign = sign == "Elevation" ? 1 : -1;
            }
        }
        throw geode::OpenGeodeException{
            "Cannot find the end of CRS section"
        };
    }
    
    void goto_keyword( std::ifstream& file, const std::string& word )
    {
        std::string line;
        while( std::getline( file, line ) )
        {
            if( string_starts_with( line, word ) )
            {
                return;
            }
        }
        throw geode::OpenGeodeException(
            "Cannot find the requested word: ", word );
    }


        std::string read_name( std::istringstream& iss )
        {
            std::string name;
            iss >> name;
            std::string token;
            while( iss >> token )
            {
                name += "_" + token;
            }
            return name;
        }

    TSurfData read_tsurf( std::ifstream& file )
    {
        check_keyword( file, "GOCAD TSurf" );
        TSurfData tsurf;
        tsurf.header = read_header( file );
        tsurf.crs = read_CRS( file );
        read_tfaces( file, tsurf );
        return tsurf;
    }
} // namespace geode
