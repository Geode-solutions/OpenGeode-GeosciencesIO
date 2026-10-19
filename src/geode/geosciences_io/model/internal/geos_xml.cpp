#include <geode/geosciences_io/model/internal/geos_xml.hpp>

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
} // namespace geode::internal
