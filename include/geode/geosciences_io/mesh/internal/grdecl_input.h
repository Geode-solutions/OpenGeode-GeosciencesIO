/*
 * Copyright (c) 2019 - 2024 Geode-solutions
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
#include <geode/mesh/io/hybrid_solid_input.h>

#include <geode/geosciences_io/mesh/common.h>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( HybridSolid );
    ALIAS_3D( HybridSolid );
} // namespace geode

namespace geode
{
    namespace internal
    {
        class GRDECLInput : public HybridSolidInput< 3 >
        {
        public:
            explicit GRDECLInput( std::string_view filename )
                : HybridSolidInput< 3 >( filename )
            {
            }

            static std::string_view extension()
            {
                static constexpr auto EXT = "grdecl";
                return EXT;
            }

            std::unique_ptr< HybridSolid3D > read( const MeshImpl& impl ) final;
        };
    } // namespace internal
} // namespace geode