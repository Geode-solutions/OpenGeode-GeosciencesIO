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

#pragma once

#include <string>
#include <vector>

#include <geode/mesh/io/point_set_output.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSet );
    ALIAS_3D( PointSet );
} // namespace geode

namespace geode
{
    namespace internal
    {
        class VSOutput final : public PointSetOutput< 3 >
        {
        public:
            explicit VSOutput( std::string_view filename )
                : PointSetOutput< 3 >( filename )
            {
            }

            static std::string_view extension()
            {
                static constexpr auto EXT = "vs";
                return EXT;
            }

            std::vector< std::string > write(
                const PointSet3D& surface ) const final;
        };
    } // namespace internal
} // namespace geode