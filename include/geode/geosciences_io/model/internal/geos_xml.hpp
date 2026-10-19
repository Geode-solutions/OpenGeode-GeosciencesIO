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
