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

#include <geode/geosciences_io/model/internal/geos_xml.hpp>

#include <algorithm>
#include <cctype>

#include <pugixml.hpp>

#include <absl/strings/str_cat.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_join.h>

namespace geode::internal
{
    std::string geos_real( double value )
    {
        return absl::StrFormat( "%.15g", value );
    }

    std::string geos_array( absl::Span< const double > values )
    {
        return absl::StrCat( "{ ",
            absl::StrJoin( values, ", ",
                []( std::string* out, double value ) {
                    absl::StrAppend( out, geos_real( value ) );
                } ),
            " }" );
    }

    std::string geos_names( absl::Span< const std::string > names )
    {
        return absl::StrCat( "{ ", absl::StrJoin( names, ", " ), " }" );
    }

    void set_attribute( pugi::xml_node& node,
        std::string_view attribute,
        std::string_view value )
    {
        const std::string attribute_name{ attribute };
        auto xml_attribute = node.attribute( attribute_name.c_str() );
        if( !xml_attribute )
        {
            xml_attribute = node.append_attribute( attribute_name.c_str() );
        }
        xml_attribute.set_value( std::string{ value }.c_str() );
    }

    void set_real_attribute(
        pugi::xml_node& node, std::string_view attribute, double value )
    {
        set_attribute( node, attribute, geos_real( value ) );
    }

    void add_comment( pugi::xml_node& node, std::string_view comment )
    {
        node.append_child( pugi::node_comment )
            .set_value( absl::StrCat( " ", comment, " " ).c_str() );
    }

    bool is_geos_name( std::string_view name )
    {
        return !name.empty()
               && std::all_of( name.begin(), name.end(), []( char character ) {
                      return std::isalnum(
                                 static_cast< unsigned char >( character ) )
                                 != 0
                             || character == '_' || character == '-';
                  } );
    }

    bool is_xml_attribute_name( std::string_view name )
    {
        const auto is_letter = []( char character ) {
            return std::isalpha( static_cast< unsigned char >( character ) )
                       != 0
                   || character == '_';
        };
        return !name.empty() && is_letter( name.front() )
               && std::all_of(
                   name.begin(), name.end(), [&is_letter]( char character ) {
                       return is_letter( character )
                              || std::isdigit(
                                     static_cast< unsigned char >( character ) )
                                     != 0;
                   } );
    }
} // namespace geode::internal
