/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

#include <fstream>
#include <optional>
#include <string>

#include <absl/strings/match.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_replace.h>

#include <geode/basic/file.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/string.hpp>
#include <geode/basic/variable_attribute.hpp>

namespace
{
    std::string get_string_between_quote(
        const std::vector< std::string_view > tokens,
        geode::index_t& from_to_id )
    {
        std::string string_between_quote;
        while(
            from_to_id++ < tokens.size() && tokens[from_to_id].back() != '\"' )
        {
            absl::StrAppend( &string_between_quote, tokens[from_to_id], " " );
        }
        if( from_to_id == tokens.size() )
        {
            throw geode::OpenGeodeException{
                "[Reading Inputs From Skua-Gocad] missing a closing "
                "quote character."
            };
        }
        absl::StrAppend( &string_between_quote, tokens[from_to_id] );
        string_between_quote.erase( string_between_quote.size() - 1, 1 );
        return string_between_quote;
    }

    std::vector< std::string > split_string_considering_quotes(
        std::string_view string_to_split )
    {
        std::vector< std::string > merged;
        const auto tokens = geode::string_split( string_to_split );
        geode::index_t token_id{ 0 };
        while( token_id < tokens.size() )
        {
            if( tokens[token_id] == "\"" )
            {
                merged.emplace_back(
                    get_string_between_quote( tokens, token_id ) );
            }
            else
            {
                merged.emplace_back( tokens[token_id] );
            }
            token_id++;
        }
        return merged;
    }
    std::string write_string_with_quotes( std::string_view string )
    {
        const auto tokens = geode::string_split( string );
        if( tokens.size() > 1 )
        {
            return absl::StrCat(
                "\" ", geode::internal::read_name( tokens ), "\"" );
        }
        return geode::internal::read_name( tokens );
    }

    void read_ilines( std::ifstream& file, geode::internal::ECurveData& ecurve )
    {
        geode::goto_keyword( file, "ILINE" );
        std::string line;
        while( std::getline( file, line ) )
        {
            const auto tokens = geode::string_split( line );
            const auto& keyword = tokens.front();
            if( keyword == "VRTX" || keyword == "PVRTX" )
            {
                if( ecurve.points.empty() )
                {
                    ecurve.OFFSET_START = geode::string_to_index( tokens[1] );
                }
                ecurve.points.emplace_back( std::array< double, 3 >{
                    geode::string_to_double( tokens[2] ),
                    geode::string_to_double( tokens[3] ),
                    geode::string_to_double( tokens[4] )
                        * ( ecurve.crs.z_sign_positive ? 1. : -1. ) } );
            }
            else if( keyword == "SEG" )
            {
                ecurve.edges.emplace_back( std::array< geode::index_t, 2 >{
                    geode::string_to_index( tokens[1] ) - ecurve.OFFSET_START,
                    geode::string_to_index( tokens[2] )
                        - ecurve.OFFSET_START } );
            }
            else if( keyword == "END" )
            {
                return;
            }
        }
    }

    void read_tfaces( std::ifstream& file, geode::internal::TSurfData& tsurf )
    {
        geode::goto_keyword( file, "TFACE" );
        std::string line;
        while( std::getline( file, line ) )
        {
            const auto tokens = geode::string_split( line );
            const auto& keyword = tokens.front();
            if( keyword == "VRTX" || keyword == "PVRTX" )
            {
                if( tsurf.points.empty() )
                {
                    tsurf.OFFSET_START = geode::string_to_index( tokens[1] );
                }
                tsurf.points.emplace_back( std::array< double, 3 >{
                    geode::string_to_double( tokens[2] ),
                    geode::string_to_double( tokens[3] ),
                    geode::string_to_double( tokens[4] )
                        * ( tsurf.crs.z_sign_positive ? 1. : -1. ) } );
                geode::internal::read_properties(
                    tsurf.vertices_properties_header,
                    tsurf.vertices_attribute_values, tokens, 5 );
            }
            else if( keyword == "ATOM" || keyword == "PATOM" )
            {
                tsurf.points.push_back(
                    tsurf.points.at( geode::string_to_index( tokens[2] )
                                     - tsurf.OFFSET_START ) );
                geode::internal::read_properties(
                    tsurf.vertices_properties_header,
                    tsurf.vertices_attribute_values, tokens, 3 );
            }
            else if( keyword == "TRGL" )
            {
                tsurf.triangles.emplace_back( std::array< geode::index_t, 3 >{
                    geode::string_to_index( tokens[1] ) - tsurf.OFFSET_START,
                    geode::string_to_index( tokens[2] ) - tsurf.OFFSET_START,
                    geode::string_to_index( tokens[3] )
                        - tsurf.OFFSET_START } );
            }
            else if( keyword == "BSTONE" )
            {
                tsurf.bstones.push_back(
                    geode::string_to_index( tokens[1] ) - tsurf.OFFSET_START );
            }
            else if( keyword == "BORDER" )
            {
                tsurf.borders.emplace_back(
                    geode::string_to_index( tokens[2] ) - tsurf.OFFSET_START,
                    geode::string_to_index( tokens[3] ) - tsurf.OFFSET_START );
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

    void read_VSet_vertices(
        std::ifstream& file, geode::internal::VSetData& vertex_set )
    {
        auto line = geode::goto_keywords(
            file, std::array< std::string_view, 2 >{ "VRTX", "PVRTX" } );
        do
        {
            const auto tokens = geode::string_split( line );
            const auto& keyword = tokens.front();
            if( keyword == "END" )
            {
                return;
            }
            if( keyword != "VRTX" && keyword != "PVRTX" )
            {
                continue;
            }
            if( vertex_set.points.empty() )
            {
                vertex_set.OFFSET_START = geode::string_to_index( tokens[1] );
            }
            vertex_set.points.emplace_back(
                std::array< double, 3 >{ geode::string_to_double( tokens[2] ),
                    geode::string_to_double( tokens[3] ),
                    geode::string_to_double( tokens[4] )
                        * ( vertex_set.crs.z_sign_positive ? 1. : -1. ) } );
            geode::internal::read_properties(
                vertex_set.vertices_properties_header,
                vertex_set.vertices_attribute_values, tokens, 5 );
        } while( std::getline( file, line ) );
        throw geode::OpenGeodeException{
            "[read_tfaces] Cannot find the end of VSet section"
        };
    }

    void read_property_keyword_with_one_string( std::ifstream& file,
        std::string_view keyword,
        std::vector< std::string >& keyword_data,
        geode::index_t nb_attributes )
    {
        auto line = geode::goto_keyword( file, keyword );
        const auto split_line = split_string_considering_quotes( line );
        keyword_data.resize( nb_attributes );
        for( const auto attr_id : geode::Range{ nb_attributes } )
        {
            keyword_data[attr_id] = geode::to_string( split_line[attr_id + 1] );
        }
    }

    void read_property_keyword_with_two_strings( std::ifstream& file,
        std::string_view keyword,
        std::vector< std::pair< std::string, std::string > >& keyword_data,
        geode::index_t nb_attributes )
    {
        auto line = geode::goto_keyword( file, keyword );
        const auto split_line = geode::string_split( line );
        keyword_data.resize( nb_attributes );
        geode::index_t counter{ 0 };
        for( const auto attr_id : geode::Range{ nb_attributes } )
        {
            keyword_data[attr_id] = {
                geode::to_string( split_line[2 * attr_id + 1 + counter] ),
                geode::to_string( split_line[2 * attr_id + 2 + counter] )
            };
            if( keyword_data[attr_id].first == "LINEARFUNCTION" )
            {
                counter += 2;
            }
        }
    }

    void read_property_keyword_with_one_double( std::ifstream& file,
        std::string_view keyword,
        std::vector< double >& keyword_data,
        geode::index_t nb_attributes )
    {
        auto line = geode::goto_keyword( file, keyword );
        const auto split_line = geode::string_split( line );
        keyword_data.resize( nb_attributes );
        for( const auto attr_id : geode::Range{ nb_attributes } )
        {
            keyword_data[attr_id] =
                geode::string_to_double( split_line[attr_id + 1] );
        }
    }

    void read_property_keyword_with_one_index_t( std::ifstream& file,
        std::string_view keyword,
        std::vector< geode::index_t >& keyword_data,
        geode::index_t nb_attributes )
    {
        auto line = geode::goto_keyword( file, keyword );
        const auto split_line = geode::string_split( line );
        keyword_data.resize( nb_attributes );
        for( const auto attr_id : geode::Range{ nb_attributes } )
        {
            keyword_data[attr_id] =
                geode::string_to_index( split_line[attr_id + 1] );
        }
    }

    template < typename Container >
    void add_vertices_container_attribute( std::string_view attribute_name,
        absl::Span< const double > attribute_values,
        geode::AttributeManager& attribute_manager,
        geode::index_t nb_vertices,
        absl::Span< const geode::index_t > inverse_mapping,
        Container value_array )
    {
        auto attribute = attribute_manager.template find_or_create_attribute<
            geode::VariableAttribute, Container >(
            attribute_name, value_array );
        for( const auto pt_id : geode::Range{ nb_vertices } )
        {
            for( const auto item_id : geode::LRange{ value_array.size() } )
            {
                value_array[item_id] =
                    attribute_values[inverse_mapping[pt_id] * value_array.size()
                                     + item_id];
            }
            attribute->set_value( pt_id, value_array );
        }
    }
} // namespace

namespace geode
{
    namespace internal
    {
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
                constexpr std::string_view NAME_PREFFIX = "name:";
                const auto name_it = line.find( NAME_PREFFIX );
                if( name_it != std::string::npos )
                {
                    std::string_view name_line{ line };
                    name_line.remove_prefix( name_it + NAME_PREFFIX.size() );
                    header.name = read_name( geode::string_split( name_line ) );
                }
            }
            throw geode::OpenGeodeException{
                "[read_header] Cannot find the end of \"HEADER\" section"
            };
        }

        void write_header( std::ofstream& file, const HeaderData& data )
        {
            file << "HEADER {" << EOL;
            if( data.name )
            {
                file << "name:" << data.name.value() << EOL;
            }
            file << "}" << EOL;
        }

        CRSData read_CRS( std::ifstream& file )
        {
            CRSData crs;
            if( !next_keyword_if_it_exists(
                    file, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) )
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
                const auto tokens = split_string_considering_quotes( line );
                if( tokens[0] == "ZPOSITIVE" )
                {
                    crs.z_sign_positive = ( tokens[1] == "Elevation" );
                }
                else if( tokens[0] == "PROJECTION" )
                {
                    crs.projection = tokens[1];
                }
                else if( tokens[0] == "DATUM" )
                {
                    crs.datum = tokens[1];
                }
                else if( tokens[0] == "NAME" )
                {
                    crs.name = tokens[1];
                }
            }
            throw geode::OpenGeodeException{
                "Cannot find the end of CRS section"
            };
        }

        void write_CRS( std::ofstream& file, const CRSData& data )
        {
            file << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << EOL;
            file << "NAME " << write_string_with_quotes( data.name ) << EOL;
            file << "PROJECTION " << data.projection << EOL;
            file << "DATUM " << data.datum << EOL;
            file << "AXIS_NAME " << data.axis_names[0] << SPACE
                 << data.axis_names[1] << SPACE << data.axis_names[2] << EOL;
            file << "AXIS_UNIT " << data.axis_units[0] << SPACE
                 << data.axis_units[1] << SPACE << data.axis_units[2] << EOL;
            file << "ZPOSITIVE "
                 << ( data.z_sign_positive ? "Elevation" : "Depth" ) << EOL;
            file << "END_ORIGINAL_COORDINATE_SYSTEM" << EOL;
        }

        PropHeaderData read_prop_header(
            std::ifstream& file, std::string_view prefix )
        {
            PropHeaderData header;
            const auto opt_line = geode::next_keyword_if_it_exists(
                file, absl::StrCat( prefix, "PROPERTIES" ) );
            if( !opt_line )
            {
                return header;
            }
            const auto split_line =
                split_string_considering_quotes( opt_line.value() );
            const auto nb_attributes = split_line.size() - 1;
            if( nb_attributes == 0 )
            {
                return header;
            }
            header.names.resize( nb_attributes );
            for( const auto attr_id : geode::Range{ nb_attributes } )
            {
                header.names[attr_id] =
                    geode::to_string( split_line[attr_id + 1] );
            }
            read_property_keyword_with_two_strings( file,
                absl::StrCat( prefix, "PROP_LEGAL_RANGES" ),
                header.prop_legal_ranges, nb_attributes );
            read_property_keyword_with_one_double( file,
                absl::StrCat( prefix, "NO_DATA_VALUES" ), header.no_data_values,
                nb_attributes );
            read_property_keyword_with_one_string( file,
                absl::StrCat( prefix, "PROPERTY_CLASSES" ),
                header.property_classes, nb_attributes );
            read_property_keyword_with_one_string( file,
                absl::StrCat( prefix, "PROPERTY_KINDS" ), header.kinds,
                nb_attributes );
            read_property_keyword_with_two_strings( file,
                absl::StrCat( prefix, "PROPERTY_SUBCLASSES" ),
                header.property_subclass, nb_attributes );
            read_property_keyword_with_one_index_t( file,
                absl::StrCat( prefix, "ESIZES" ), header.esizes,
                nb_attributes );
            read_property_keyword_with_one_string( file,
                absl::StrCat( prefix, "UNITS" ), header.units, nb_attributes );
            return header;
        }

        void read_properties( const PropHeaderData& properties_header,
            std::vector< std::vector< double > >& attribute_values,
            absl::Span< const std::string_view > tokens,
            geode::index_t line_properties_position )
        {
            for( const auto attr_id :
                geode::Indices{ properties_header.names } )
            {
                for( const auto item :
                    geode::LRange{ properties_header.esizes[attr_id] } )
                {
                    geode_unused( item );
                    OPENGEODE_ASSERT( line_properties_position < tokens.size(),
                        "[GocadInput::read_point_properties] Cannot read "
                        "properties: number of property items is higher than "
                        "number of tokens." );
                    attribute_values[attr_id].push_back(
                        geode::string_to_double(
                            tokens[line_properties_position] ) );
                    line_properties_position++;
                }
            }
        }

        void create_attributes( const PropHeaderData& attributes_header,
            absl::Span< const std::vector< double > > attributes_values,
            geode::AttributeManager& attribute_manager,
            geode::index_t nb_vertices,
            absl::Span< const geode::index_t > inverse_vertex_mapping )
        {
            for( const auto attr_id :
                geode::Indices{ attributes_header.names } )
            {
                const auto nb_attribute_items =
                    attributes_header.esizes[attr_id];
                if( nb_attribute_items == 1 )
                {
                    auto attribute = attribute_manager.find_or_create_attribute<
                        geode::VariableAttribute, double >(
                        attributes_header.names[attr_id],
                        attributes_header.no_data_values[attr_id] );
                    for( const auto pt_id : geode::Range{ nb_vertices } )
                    {
                        attribute->set_value( pt_id,
                            attributes_values[attr_id]
                                             [inverse_vertex_mapping[pt_id]] );
                    }
                }
                else if( nb_attribute_items == 2 )
                {
                    std::array< double, 2 > container;
                    container.fill( attributes_header.no_data_values[attr_id] );
                    add_vertices_container_attribute(
                        attributes_header.names[attr_id],
                        attributes_values[attr_id], attribute_manager,
                        nb_vertices, inverse_vertex_mapping, container );
                }
                else if( nb_attribute_items == 3 )
                {
                    std::array< double, 3 > container;
                    container.fill( attributes_header.no_data_values[attr_id] );
                    add_vertices_container_attribute(
                        attributes_header.names[attr_id],
                        attributes_values[attr_id], attribute_manager,
                        nb_vertices, inverse_vertex_mapping, container );
                }
                else
                {
                    std::vector< double > container( nb_attribute_items,
                        attributes_header.no_data_values[attr_id] );
                    add_vertices_container_attribute<>(
                        attributes_header.names[attr_id],
                        attributes_values[attr_id], attribute_manager,
                        nb_vertices, inverse_vertex_mapping, container );
                }
            }
        }

        void write_prop_header(
            std::ofstream& file, const PropHeaderData& data )
        {
            file << "PROPERTIES";
            for( const auto& name : data.names )
            {
                file << SPACE << write_string_with_quotes( name );
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
            for( const auto& prop_class : data.property_classes )
            {
                file << SPACE << write_string_with_quotes( prop_class );
            }
            file << EOL;
            file << "PROPERTY_KINDS";
            for( const auto& kind : data.kinds )
            {
                file << SPACE << write_string_with_quotes( kind );
            }
            file << EOL;
            file << "PROPERTY_SUBCLASSES";
            for( const auto& prop_subclasse : data.property_subclass )
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
            for( const auto& unit : data.units )
            {
                file << SPACE << write_string_with_quotes( unit );
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
            file << "name:" << data.name << EOL;
            if( data.is_z )
            {
                file << "is_Z: on" << EOL;
            }
            file << "}" << EOL;
        }

        std::string read_name( absl::Span< const std::string_view > tokens )
        {
            return absl::StrReplaceAll(
                absl::StrJoin( tokens.begin(), tokens.end(), " " ),
                { { "\"", "" } } );
        }

        std::optional< TSurfData > read_tsurf( std::ifstream& file )
        {
            if( !goto_keyword_if_it_exists( file, "GOCAD TSurf" ) )
            {
                return std::nullopt;
            }
            TSurfData tsurf;
            tsurf.header = read_header( file );
            tsurf.crs = read_CRS( file );
            tsurf.vertices_properties_header =
                geode::internal::read_prop_header( file, "" );
            tsurf.vertices_attribute_values.resize(
                tsurf.vertices_properties_header.names.size() );
            read_tfaces( file, tsurf );
            return tsurf;
        }

        std::optional< ECurveData > read_ecurve( std::ifstream& file )
        {
            if( !goto_keyword_if_it_exists( file, "GOCAD PLine" ) )
            {
                return std::nullopt;
            }
            ECurveData ecurve;
            ecurve.header = read_header( file );
            ecurve.crs = read_CRS( file );
            read_ilines( file, ecurve );
            return ecurve;
        }

        std::optional< VSetData > read_vs_points( std::ifstream& file )
        {
            if( !goto_keyword_if_it_exists( file, "GOCAD VSet" ) )
            {
                return std::nullopt;
            }
            VSetData vertex_set;
            vertex_set.header = read_header( file );
            vertex_set.crs = read_CRS( file );
            vertex_set.vertices_properties_header =
                geode::internal::read_prop_header( file, "" );
            vertex_set.vertices_attribute_values.resize(
                vertex_set.vertices_properties_header.names.size() );
            read_VSet_vertices( file, vertex_set );
            return vertex_set;
        }
    } // namespace internal
} // namespace geode
