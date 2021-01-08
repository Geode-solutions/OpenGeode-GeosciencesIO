/*
 * Copyright (c) 2019 - 2020 Geode-solutions
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

#include <geode/geosciences/private/gocad_common.h>

#include <fstream>

#include <absl/strings/match.h>
#include <absl/strings/str_replace.h>

#include <geode/basic/logger.h>

namespace
{
    static constexpr char EOL{ '\n' };
    static constexpr char SPACE{ ' ' };

    void read_tfaces( std::ifstream& file, geode::detail::TSurfData& tsurf )
    {
        geode::detail::goto_keyword( file, "TFACE" );
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
                tsurf.points.push_back(
                    tsurf.points.at( atom_id - tsurf.OFFSET_START ) );
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
                tsurf.borders.emplace_back(
                    corner - tsurf.OFFSET_START, next - tsurf.OFFSET_START );
            }
            else if( keyword == "TFACE" )
            {
                tsurf.tface_triangles_offset.push_back(
                    tsurf.triangles.size() );
                tsurf.tface_vertices_offset.push_back( tsurf.points.size() );
            }
            else if( keyword == "END" )
            {
                tsurf.tface_triangles_offset.push_back(
                    tsurf.triangles.size() );
                tsurf.tface_vertices_offset.push_back( tsurf.points.size() );
                return;
            }
        }
        throw geode::OpenGeodeException{
            "[read_tfaces] Cannot find the end of TSurf section"
        };
    }
} // namespace

namespace geode
{
    namespace detail
    {
        bool string_starts_with(
            absl::string_view string, absl::string_view check )
        {
            return absl::StartsWith( string, check );
        }

        bool line_starts_with( std::ifstream& file, absl::string_view check )
        {
            std::string line;
            std::getline( file, line );
            return string_starts_with( line, check );
        }

        void check_keyword( std::ifstream& file, absl::string_view keyword )
        {
            OPENGEODE_EXCEPTION( line_starts_with( file, keyword ),
                absl::StrCat( "Line should starts with \"", keyword, "\"" ) );
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
            throw geode::OpenGeodeException{
                "[read_header] Cannot find the end of \"HEADER\" section"
            };
        }

        void write_header( std::ofstream& file, const HeaderData& data )
        {
            file << "HEADER {" << EOL;
            file << "name:" << data.name << EOL;
            file << "}" << EOL;
        }

        CRSData read_CRS( std::ifstream& file )
        {
            CRSData crs;
            if( !line_starts_with( file, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) )
            {
                return crs;
            }
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

        void write_CRS( std::ofstream& file, const CRSData& data )
        {
            file << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << EOL;
            file << "NAME " << data.name << EOL;
            file << "AXIS_NAME " << data.axis_names[0] << SPACE
                 << data.axis_names[1] << SPACE << data.axis_names[2] << EOL;
            file << "AXIS_UNIT " << data.axis_units[0] << SPACE
                 << data.axis_units[1] << SPACE << data.axis_units[2] << EOL;
            file << "ZPOSITIVE " << ( data.z_sign == 1 ? "Elevation" : "Depth" )
                 << EOL;
            file << "END_ORIGINAL_COORDINATE_SYSTEM" << EOL;
        }
        /*   PropHeaderData read_prop_header( std::ifstream& file )
           {
               check_keyword( file, "PROPERTY_CLASS_HEADER" );
               PropHeaderData header;
               // header.name = next to  "PROPERTY_CLASS_HEADER"
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
                   if( keyword == "kind:" )
                   {
                       header.kind = read_name( iss );
                   }
                   if( keyword == "unit:" )
                   {
                       header.unit = read_name( iss );
                   }
                   if( keyword == "is_z:" )
                   {
                       header.is_z = true;
                   }
               }
               throw geode::OpenGeodeException{
                   "[read_header] Cannot find the end of "
                   "\"PROPERTY_CLASS_HEADER\" section"
               };
           }*/

        void write_prop_header(
            std::ofstream& file, const PropHeaderData& data )
        {
            file << "PROPERTIES";
            for( const auto& name : data.names )
            {
                file << SPACE << name;
            }
            file << EOL;
            file << "PROP_LEGAL_RANGES";
            for( const auto& prop_range : data.prop_legal_ranges )
            {
                file << SPACE << prop_range.first << SPACE << prop_range.second;
            }
            file << EOL;
            file << "NO_DATA_VALUES";
            for( const auto prop_ndv : data.no_data_values )
            {
                file << SPACE << prop_ndv;
            }
            file << EOL;
            file << "PROPERTY_CLASSES";
            for( const auto prop_class : data.property_classes )
            {
                file << SPACE << prop_class;
            }
            file << EOL;
            file << "PROPERTY_KINDS";
            for( const auto kind : data.kinds )
            {
                file << SPACE << kind;
            }
            file << EOL;
            file << "PROPERTY_SUBCLASSES";
            for( const auto prop_subclasse : data.property_subclass )
            {
                file << SPACE << prop_subclasse.first << SPACE
                     << prop_subclasse.second;
            }
            file << EOL;
            file << "ESIZES";
            for( const auto esize : data.esizes )
            {
                file << SPACE << esize;
            }
            file << EOL;
            file << "UNITS";
            for( const auto unit : data.units )
            {
                file << SPACE << unit;
            }
            file << EOL;
        }
        void write_property_class_header(
            std::ofstream& file, const PropClassHeaderData& data )
        {
            file << "PROPERTY_CLASS_HEADER" << SPACE << data.name << SPACE
                 << "{" << EOL;
            file << "kind:" << data.kind << EOL;
            file << "unit:" << data.unit << EOL;
            if( data.is_z )
            {
                file << "is_Z: on" << EOL;
            }
            file << "}" << EOL;
        }
        std::string goto_keyword( std::ifstream& file, absl::string_view word )
        {
            std::string line;
            while( std::getline( file, line ) )
            {
                if( string_starts_with( line, word ) )
                {
                    return line;
                }
            }
            throw geode::OpenGeodeException{
                "[goto_keyword] Cannot find the requested keyword: ", word
            };
            return "";
        }

        std::string goto_keywords(
            std::ifstream& file, absl::Span< const absl::string_view > words )
        {
            std::string line;
            while( std::getline( file, line ) )
            {
                for( const auto word : words )
                {
                    if( string_starts_with( line, word ) )
                    {
                        return line;
                    }
                }
            }
            throw geode::OpenGeodeException{
                "[goto_keywords] Cannot find one of the requested keywords"
            };
            return "";
        }

        std::string read_name( std::istringstream& iss )
        {
            std::string name;
            iss >> name;
            std::string token;
            while( iss >> token )
            {
                absl::StrAppend( &name, " ", token );
            }
            return absl::StrReplaceAll( name, { { "\"", "" } } );
        }

        absl::optional< TSurfData > read_tsurf( std::ifstream& file )
        {
            if( !line_starts_with( file, "GOCAD TSurf" ) )
            {
                return absl::nullopt;
            }
            TSurfData tsurf;
            tsurf.header = read_header( file );
            tsurf.crs = read_CRS( file );
            read_tfaces( file, tsurf );
            return tsurf;
        }
    } // namespace detail
} // namespace geode
