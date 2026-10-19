#pragma once

#include <string>
#include <string_view>

#include <absl/types/span.h>

#include <geode/geosciences_io/model/common.hpp>

namespace pugi
{
    class xml_node;
} // namespace pugi

namespace geode::internal
{
    /// Real value written with 15 significant digits
    [[nodiscard]] std::string geos_real( double value );

    /// GEOS array of reals: "{ a, b, c }"
    [[nodiscard]] std::string geos_array( absl::Span< const double > values );

    /// GEOS array of names: "{ a, b, c }"
    [[nodiscard]] std::string geos_names(
        absl::Span< const std::string > names );

    void set_attribute( pugi::xml_node& node,
        std::string_view attribute,
        std::string_view value );

    void set_real_attribute(
        pugi::xml_node& node, std::string_view attribute, double value );

    void add_comment( pugi::xml_node& node, std::string_view comment );
} // namespace geode::internal
