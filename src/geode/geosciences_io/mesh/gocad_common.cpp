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

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

#include <fstream>
#include <optional>
#include <string>

#include <absl/algorithm/container.h>
#include <absl/strings/match.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_replace.h>
#include <absl/strings/strip.h>

#include <geode/basic/file.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/string.hpp>
#include <geode/basic/variable_attribute.hpp>

namespace
{
    std::string get_string_between_quote(
        absl::Span< const std::string_view > tokens,
        geode::index_t& from_to_id )
    {
        auto token = tokens[from_to_id];
        token.remove_prefix( 1 );
        if( !token.empty() && token.back() == '\"' )
        {
            token.remove_suffix( 1 );
            return geode::to_string( token );
        }
        std::string string_between_quote{ token };
        while( ++from_to_id < tokens.size() )
        {
            token = tokens[from_to_id];
            const auto is_closing = token.back() == '\"';
            if( is_closing )
            {
                token.remove_suffix( 1 );
            }
            if( !string_between_quote.empty() && !token.empty() )
            {
                string_between_quote.push_back( ' ' );
            }
            absl::StrAppend( &string_between_quote, token );
            if( is_closing )
            {
                return string_between_quote;
            }
        }
        geode::Logger::warning( "[Reading Inputs From Skua-Gocad] Missing a "
                                "closing quote character, the end of the line ",
            "is used: \"", string_between_quote, "\"" );
        return string_between_quote;
    }

    std::vector< std::string > split_string_considering_quotes(
        std::string_view string_to_split )
    {
        std::vector< std::string > merged;
        const auto tokens = geode::string_split( string_to_split );
        merged.reserve( tokens.size() );
        for( geode::index_t token_id = 0; token_id < tokens.size(); token_id++ )
        {
            if( tokens[token_id].front() == '\"' )
            {
                merged.emplace_back(
                    get_string_between_quote( tokens, token_id ) );
            }
            else
            {
                merged.emplace_back( tokens[token_id] );
            }
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

    bool is_empty_or_comment( absl::Span< const std::string_view > tokens )
    {
        return tokens.empty() || tokens.front().front() == '#';
    }

    bool has_enough_tokens( absl::Span< const std::string_view > tokens,
        geode::index_t nb_tokens,
        std::string_view line )
    {
        if( tokens.size() >= nb_tokens )
        {
            return true;
        }
        geode::Logger::warning( "[GocadInput] Line \"", line,
            "\" should have at least ", nb_tokens, " tokens, it is skipped" );
        return false;
    }

    void check_enough_tokens( absl::Span< const std::string_view > tokens,
        geode::index_t nb_tokens,
        std::string_view line )
    {
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            tokens.size() >= nb_tokens, nullptr,
            geode::OpenGeodeException::TYPE::data, "[GocadInput] Line \"", line,
            "\" should have at least ", nb_tokens, " tokens" );
    }

    std::optional< geode::index_t > vertex_index( geode::index_t file_id,
        geode::index_t offset_start,
        geode::index_t nb_points )
    {
        if( file_id < offset_start || file_id - offset_start >= nb_points )
        {
            return std::nullopt;
        }
        return file_id - offset_start;
    }

    void resolve_ecurve_indices( geode::internal::ECurveData& ecurve )
    {
        geode::index_t nb_valid{ 0 };
        for( auto edge : ecurve.edges )
        {
            const auto v0 = vertex_index(
                edge[0], ecurve.OFFSET_START, ecurve.points.size() );
            const auto v1 = vertex_index(
                edge[1], ecurve.OFFSET_START, ecurve.points.size() );
            if( !v0 || !v1 || v0 == v1 )
            {
                geode::Logger::warning( "[PLInput] SEG ", edge[0], " ", edge[1],
                    " references an unknown vertex or is degenerated, "
                    "it is skipped" );
                continue;
            }
            ecurve.edges[nb_valid++] = { v0.value(), v1.value() };
        }
        ecurve.edges.resize( nb_valid );
    }

    void resolve_tsurf_triangles( geode::internal::TSurfData& tsurf )
    {
        geode::index_t nb_valid{ 0 };
        geode::index_t nb_invalid{ 0 };
        for( const auto tface :
            geode::Range{ tsurf.tface_triangles_offset.size() - 1 } )
        {
            const auto begin = tsurf.tface_triangles_offset[tface];
            const auto end = tsurf.tface_triangles_offset[tface + 1];
            const auto first_vertex = tsurf.tface_vertices_offset[tface];
            const auto nb_tface_vertices =
                tsurf.tface_vertices_offset[tface + 1] - first_vertex;
            tsurf.tface_triangles_offset[tface] = nb_valid;
            for( const auto triangle_id : geode::Range{ begin, end } )
            {
                auto triangle = tsurf.triangles[triangle_id];
                auto is_valid = true;
                for( auto& vertex : triangle )
                {
                    const auto index = vertex_index( vertex,
                        tsurf.OFFSET_START + first_vertex, nb_tface_vertices );
                    if( !index )
                    {
                        is_valid = false;
                        break;
                    }
                    vertex = index.value() + first_vertex;
                }
                if( !is_valid )
                {
                    nb_invalid++;
                    continue;
                }
                tsurf.triangles[nb_valid++] = triangle;
            }
        }
        tsurf.tface_triangles_offset.back() = nb_valid;
        tsurf.triangles.resize( nb_valid );
        if( nb_invalid > 0 )
        {
            geode::Logger::warning( "[TSInput] ", nb_invalid,
                " TRGL referencing vertices outside of their TFACE are "
                "skipped" );
        }
    }

    void resolve_tsurf_indices( geode::internal::TSurfData& tsurf )
    {
        resolve_tsurf_triangles( tsurf );
        const auto nb_points = tsurf.points.size();
        geode::index_t nb_valid{ 0 };
        for( const auto bstone : tsurf.bstones )
        {
            const auto index =
                vertex_index( bstone, tsurf.OFFSET_START, nb_points );
            if( !index )
            {
                geode::Logger::warning( "[TSInput] BSTONE ", bstone,
                    " references an unknown vertex, it is skipped" );
                continue;
            }
            tsurf.bstones[nb_valid++] = index.value();
        }
        tsurf.bstones.resize( nb_valid );
        nb_valid = 0;
        for( const auto& border : tsurf.borders )
        {
            const auto corner =
                vertex_index( border.corner_id, tsurf.OFFSET_START, nb_points );
            const auto next =
                vertex_index( border.next_id, tsurf.OFFSET_START, nb_points );
            if( !corner || !next
                || tsurf.tface_id( corner.value() )
                       != tsurf.tface_id( next.value() ) )
            {
                geode::Logger::warning( "[TSInput] BORDER ", border.corner_id,
                    " ", border.next_id,
                    " references unknown vertices or vertices from different "
                    "TFACE, it is skipped" );
                continue;
            }
            tsurf.borders[nb_valid++] = { corner.value(), next.value() };
        }
        tsurf.borders.resize( nb_valid, { 0, 0 } );
    }

    bool is_object_end( std::ifstream& file,
        std::string_view line,
        std::string_view keyword,
        std::string_view object )
    {
        if( keyword == "END" )
        {
            return true;
        }
        if( keyword == "GOCAD" )
        {
            geode::Logger::warning(
                "[GocadInput] Missing END keyword in ", object, " section" );
            file.seekg( file.tellg()
                        - static_cast< std::streamoff >( line.size() + 1 ) );
            return true;
        }
        return false;
    }

    void warn_missing_end( std::ifstream& file, std::string_view object )
    {
        geode::Logger::warning(
            "[GocadInput] Missing END keyword in ", object, " section" );
        file.clear();
    }

    void read_ilines( std::ifstream& file, geode::internal::ECurveData& ecurve )
    {
        std::string line;
        while( std::getline( file, line ) )
        {
            const auto tokens = geode::string_split( line );
            if( is_empty_or_comment( tokens ) )
            {
                continue;
            }
            const auto& keyword = tokens.front();
            if( is_object_end( file, line, keyword, "PLine" ) )
            {
                resolve_ecurve_indices( ecurve );
                return;
            }
            if( keyword == "VRTX" || keyword == "PVRTX" )
            {
                check_enough_tokens( tokens, 5, line );
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
                if( !has_enough_tokens( tokens, 3, line ) )
                {
                    continue;
                }
                ecurve.edges.emplace_back( std::array< geode::index_t, 2 >{
                    geode::string_to_index( tokens[1] ),
                    geode::string_to_index( tokens[2] ) } );
            }
        }
        warn_missing_end( file, "PLine" );
        resolve_ecurve_indices( ecurve );
    }

    void end_tsurf( geode::internal::TSurfData& tsurf )
    {
        tsurf.tface_triangles_offset.push_back( tsurf.triangles.size() );
        tsurf.tface_vertices_offset.push_back( tsurf.points.size() );
        resolve_tsurf_indices( tsurf );
    }

    void read_tfaces( std::ifstream& file, geode::internal::TSurfData& tsurf )
    {
        std::string line;
        auto tface_started = false;
        while( std::getline( file, line ) )
        {
            const auto tokens = geode::string_split( line );
            if( is_empty_or_comment( tokens ) )
            {
                continue;
            }
            const auto& keyword = tokens.front();
            if( is_object_end( file, line, keyword, "TSurf" ) )
            {
                end_tsurf( tsurf );
                return;
            }
            if( keyword == "VRTX" || keyword == "PVRTX" )
            {
                check_enough_tokens( tokens, 5, line );
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
                check_enough_tokens( tokens, 3, line );
                tsurf.points.push_back(
                    tsurf.points.at( geode::string_to_index( tokens[2] )
                                     - tsurf.OFFSET_START ) );
                geode::internal::read_properties(
                    tsurf.vertices_properties_header,
                    tsurf.vertices_attribute_values, tokens, 3 );
            }
            else if( keyword == "TRGL" )
            {
                if( !has_enough_tokens( tokens, 4, line ) )
                {
                    continue;
                }
                tsurf.triangles.emplace_back( std::array< geode::index_t, 3 >{
                    geode::string_to_index( tokens[1] ),
                    geode::string_to_index( tokens[2] ),
                    geode::string_to_index( tokens[3] ) } );
            }
            else if( keyword == "BSTONE" )
            {
                if( !has_enough_tokens( tokens, 2, line ) )
                {
                    continue;
                }
                tsurf.bstones.push_back( geode::string_to_index( tokens[1] ) );
            }
            else if( keyword == "BORDER" )
            {
                if( !has_enough_tokens( tokens, 4, line ) )
                {
                    continue;
                }
                tsurf.borders.emplace_back( geode::string_to_index( tokens[2] ),
                    geode::string_to_index( tokens[3] ) );
            }
            else if( keyword == "TFACE" )
            {
                if( tface_started || !tsurf.points.empty()
                    || !tsurf.triangles.empty() )
                {
                    tsurf.tface_triangles_offset.push_back(
                        tsurf.triangles.size() );
                    tsurf.tface_vertices_offset.push_back(
                        tsurf.points.size() );
                }
                tface_started = true;
            }
        }
        warn_missing_end( file, "TSurf" );
        end_tsurf( tsurf );
    }

    void read_VSet_vertices(
        std::ifstream& file, geode::internal::VSetData& vertex_set )
    {
        std::string line;
        while( std::getline( file, line ) )
        {
            const auto tokens = geode::string_split( line );
            if( is_empty_or_comment( tokens ) )
            {
                continue;
            }
            const auto& keyword = tokens.front();
            if( is_object_end( file, line, keyword, "VSet" ) )
            {
                return;
            }
            if( keyword != "VRTX" && keyword != "PVRTX" )
            {
                continue;
            }
            check_enough_tokens( tokens, 5, line );
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
        }
        warn_missing_end( file, "VSet" );
    }

    std::vector< std::string > read_property_strings( std::string_view line )
    {
        auto tokens = split_string_considering_quotes( line );
        tokens.erase( tokens.begin() );
        return tokens;
    }

    std::vector< std::pair< std::string, std::string > >
        read_property_string_pairs( std::string_view line )
    {
        const auto tokens = geode::string_split( line );
        std::vector< std::pair< std::string, std::string > > values;
        for( geode::index_t token_id = 1; token_id + 1 < tokens.size();
            token_id += 2 )
        {
            values.emplace_back( geode::to_string( tokens[token_id] ),
                geode::to_string( tokens[token_id + 1] ) );
            if( values.back().first == "LINEARFUNCTION" )
            {
                token_id += 2;
            }
        }
        return values;
    }

    template < typename T, typename Converter >
    std::vector< T > read_property_numbers(
        std::string_view line, Converter converter )
    {
        const auto tokens = geode::string_split( line );
        std::vector< T > values;
        values.reserve( tokens.size() - 1 );
        for( const auto token_id : geode::Range{ 1, tokens.size() } )
        {
            values.push_back( converter( tokens[token_id] ) );
        }
        return values;
    }

    template < typename T >
    void complete_property_values( std::vector< T >& values,
        geode::index_t nb_properties,
        std::string_view keyword,
        const T& default_value )
    {
        if( values.size() == nb_properties )
        {
            return;
        }
        geode::Logger::warning( "[GocadInput] ", keyword, " has ",
            values.size(), " values for ", nb_properties,
            " properties, missing values are set to default" );
        values.resize( nb_properties, default_value );
    }

    void complete_prop_header( geode::internal::PropHeaderData& header )
    {
        const auto nb_properties = header.names.size();
        for( const auto property : geode::Indices{ header.names } )
        {
            auto& name = header.names[property];
            if( name.empty() )
            {
                name = absl::StrCat( "property_", property );
                geode::Logger::warning(
                    "[GocadInput] Empty property name renamed to ", name );
            }
        }
        complete_property_values( header.prop_legal_ranges, nb_properties,
            "PROP_LEGAL_RANGES",
            std::pair< std::string, std::string >{ "**none**", "**none**" } );
        complete_property_values(
            header.no_data_values, nb_properties, "NO_DATA_VALUES", -99999. );
        if( header.property_classes.size() < nb_properties )
        {
            geode::Logger::warning( "[GocadInput] PROPERTY_CLASSES has ",
                header.property_classes.size(), " values for ", nb_properties,
                " properties, missing values are set to property names" );
            for( const auto property :
                geode::Range{ header.property_classes.size(), nb_properties } )
            {
                header.property_classes.push_back( header.names[property] );
            }
        }
        complete_property_values( header.property_classes, nb_properties,
            "PROPERTY_CLASSES", std::string{} );
        complete_property_values( header.kinds, nb_properties, "PROPERTY_KINDS",
            std::string{ "Real Number" } );
        complete_property_values( header.property_subclass, nb_properties,
            "PROPERTY_SUBCLASSES",
            std::pair< std::string, std::string >{ "QUANTITY", "Float" } );
        complete_property_values(
            header.esizes, nb_properties, "ESIZES", geode::index_t{ 1 } );
        for( auto& esize : header.esizes )
        {
            if( esize == 0 )
            {
                geode::Logger::warning(
                    "[GocadInput] ESIZES of 0 replaced by 1" );
                esize = 1;
            }
        }
        complete_property_values(
            header.units, nb_properties, "UNITS", std::string{ "unitless" } );
    }

    bool read_prop_header_line( geode::internal::PropHeaderData& header,
        std::string_view line,
        std::string_view prefix )
    {
        const auto tokens = geode::string_split( line );
        if( tokens.empty() )
        {
            return false;
        }
        auto keyword = tokens.front();
        if( !absl::ConsumePrefix( &keyword, prefix ) )
        {
            return false;
        }
        if( keyword == "PROPERTIES" )
        {
            header.names = read_property_strings( line );
        }
        else if( keyword == "PROP_LEGAL_RANGES" )
        {
            header.prop_legal_ranges = read_property_string_pairs( line );
        }
        else if( keyword == "NO_DATA_VALUES" )
        {
            header.no_data_values = read_property_numbers< double >(
                line, geode::string_to_double );
        }
        else if( keyword == "PROPERTY_CLASSES" )
        {
            header.property_classes = read_property_strings( line );
        }
        else if( keyword == "PROPERTY_KINDS" )
        {
            header.kinds = read_property_strings( line );
        }
        else if( keyword == "PROPERTY_SUBCLASSES" )
        {
            header.property_subclass = read_property_string_pairs( line );
        }
        else if( keyword == "ESIZES" )
        {
            header.esizes = read_property_numbers< geode::index_t >(
                line, geode::string_to_index );
        }
        else if( keyword == "UNITS" )
        {
            header.units = read_property_strings( line );
        }
        else
        {
            return false;
        }
        return true;
    }

    struct ObjectHeaderData
    {
        geode::internal::HeaderData header;
        geode::internal::CRSData crs;
        geode::internal::PropHeaderData vertices_properties_header;
    };

    bool is_object_body_keyword( std::string_view keyword )
    {
        static constexpr std::array< std::string_view, 12 > BODY_KEYWORDS{
            "GOCAD", "TFACE", "ILINE", "VRTX", "PVRTX", "ATOM", "PATOM", "TRGL",
            "SEG", "BSTONE", "BORDER", "END"
        };
        return absl::c_linear_search( BODY_KEYWORDS, keyword );
    }

    void skip_block( std::ifstream& file )
    {
        std::string line;
        while( std::getline( file, line ) )
        {
            if( absl::StartsWith(
                    absl::StripLeadingAsciiWhitespace( line ), "}" ) )
            {
                return;
            }
        }
    }

    ObjectHeaderData read_object_header( std::ifstream& file )
    {
        ObjectHeaderData data;
        data.header = geode::internal::read_header( file );
        std::string line;
        auto position = file.tellg();
        while( std::getline( file, line ) )
        {
            const auto tokens = geode::string_split( line );
            if( !is_empty_or_comment( tokens ) )
            {
                if( is_object_body_keyword( tokens.front() ) )
                {
                    file.seekg( position );
                    break;
                }
                if( tokens.front() == "GOCAD_ORIGINAL_COORDINATE_SYSTEM" )
                {
                    file.seekg( position );
                    data.crs = geode::internal::read_CRS( file );
                }
                else if( !read_prop_header_line(
                             data.vertices_properties_header, line, "" )
                         && tokens.back().back() == '{' )
                {
                    skip_block( file );
                }
            }
            position = file.tellg();
        }
        file.clear();
        complete_prop_header( data.vertices_properties_header );
        return data;
    }

    template < typename Container >
    void add_vertices_container_attribute( std::string_view attribute_name,
        absl::Span< const double > attribute_values,
        geode::AttributeManager& attribute_manager,
        geode::index_t nb_vertices,
        absl::Span< const geode::index_t > inverse_mapping,
        Container value_array )
    {
        geode::AttributeValues< Container > default_attribute_values;
        default_attribute_values.default_value = value_array;
        default_attribute_values.no_value = value_array;
        geode::AttributeProperties attribute_properties;
        attribute_properties.assignable = false;
        attribute_properties.interpolable = false;
        attribute_properties.transferable = true;
        const auto attribute_id = attribute_manager.template create_attribute<
            geode::VariableAttribute, Container >(
            attribute_name, default_attribute_values, attribute_properties );
        auto attribute =
            attribute_manager
                .find_attribute< geode::VariableAttribute, Container >(
                    attribute_id );
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

namespace geode::internal
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
        throw geode::OpenGeodeGeosciencesIOMeshException{ nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[read_header] Cannot find the end of \"HEADER\" section" };
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
            if( string_starts_with( line, "END_ORIGINAL_COORDINATE_SYSTEM" ) )
            {
                return crs;
            }
            const auto tokens = split_string_considering_quotes( line );
            if( tokens.size() < 2 )
            {
                continue;
            }
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
        throw geode::OpenGeodeGeosciencesIOMeshException{ nullptr,
            geode::OpenGeodeException::TYPE::data,
            "Cannot find the end of CRS section" };
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
        file << "ZPOSITIVE " << ( data.z_sign_positive ? "Elevation" : "Depth" )
             << EOL;
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
        read_prop_header_line( header, opt_line.value(), prefix );
        std::string line;
        auto position = file.tellg();
        while( std::getline( file, line ) )
        {
            if( !read_prop_header_line( header, line, prefix )
                && !is_empty_or_comment( geode::string_split( line ) ) )
            {
                file.seekg( position );
                break;
            }
            position = file.tellg();
        }
        file.clear();
        complete_prop_header( header );
        return header;
    }

    void read_properties( const PropHeaderData& properties_header,
        std::vector< std::vector< double > >& attribute_values,
        absl::Span< const std::string_view > tokens,
        geode::index_t line_properties_position )
    {
        for( const auto attr_id : geode::Indices{ properties_header.names } )
        {
            for( const auto item :
                geode::LRange{ properties_header.esizes[attr_id] } )
            {
                geode_unused( item );
                OpenGeodeGeosciencesIOMeshException::check_exception(
                    line_properties_position < tokens.size(), nullptr,
                    geode::OpenGeodeException::TYPE::data,
                    "[GocadInput::read_point_properties] Cannot read "
                    "properties: number of property items is higher than "
                    "number of tokens." );
                attribute_values[attr_id].push_back( geode::string_to_double(
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
        for( const auto attr_id : geode::Indices{ attributes_header.names } )
        {
            const auto nb_attribute_items = attributes_header.esizes[attr_id];
            if( nb_attribute_items == 1 )
            {
                geode::AttributeValues< double > attribute_values;
                attribute_values.default_value =
                    attributes_header.no_data_values[attr_id];
                attribute_values.no_value =
                    attributes_header.no_data_values[attr_id];
                geode::AttributeProperties attribute_properties;
                attribute_properties.assignable = false;
                attribute_properties.interpolable = false;
                attribute_properties.transferable = true;
                const auto attribute_id =
                    attribute_manager
                        .create_attribute< geode::VariableAttribute, double >(
                            attributes_header.names[attr_id], attribute_values,
                            attribute_properties );
                auto attribute =
                    attribute_manager
                        .find_attribute< geode::VariableAttribute, double >(
                            attribute_id );
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
                    attributes_values[attr_id], attribute_manager, nb_vertices,
                    inverse_vertex_mapping, container );
            }
            else if( nb_attribute_items == 3 )
            {
                std::array< double, 3 > container;
                container.fill( attributes_header.no_data_values[attr_id] );
                add_vertices_container_attribute(
                    attributes_header.names[attr_id],
                    attributes_values[attr_id], attribute_manager, nb_vertices,
                    inverse_vertex_mapping, container );
            }
            else
            {
                std::vector< double > container( nb_attribute_items,
                    attributes_header.no_data_values[attr_id] );
                add_vertices_container_attribute<>(
                    attributes_header.names[attr_id],
                    attributes_values[attr_id], attribute_manager, nb_vertices,
                    inverse_vertex_mapping, container );
            }
        }
    }

    void write_prop_header( std::ofstream& file, const PropHeaderData& data )
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
        file << "PROPERTY_CLASS_HEADER" << SPACE << data.name << SPACE << "{"
             << EOL;
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
        auto object_header = read_object_header( file );
        tsurf.header = std::move( object_header.header );
        tsurf.crs = std::move( object_header.crs );
        tsurf.vertices_properties_header =
            std::move( object_header.vertices_properties_header );
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
        auto object_header = read_object_header( file );
        ecurve.header = std::move( object_header.header );
        ecurve.crs = std::move( object_header.crs );
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
        auto object_header = read_object_header( file );
        vertex_set.header = std::move( object_header.header );
        vertex_set.crs = std::move( object_header.crs );
        vertex_set.vertices_properties_header =
            std::move( object_header.vertices_properties_header );
        vertex_set.vertices_attribute_values.resize(
            vertex_set.vertices_properties_header.names.size() );
        read_VSet_vertices( file, vertex_set );
        return vertex_set;
    }
} // namespace geode::internal
