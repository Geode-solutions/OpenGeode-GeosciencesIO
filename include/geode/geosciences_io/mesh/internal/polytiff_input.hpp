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

#include <geode/geosciences_io/mesh/common.hpp>

#include <geode/mesh/io/polygonal_surface_input.hpp>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PolygonalSurface );
    ALIAS_3D( PolygonalSurface );
} // namespace geode

namespace geode
{
    namespace internal
    {
        class PolyTIFFInput final : public PolygonalSurfaceInput3D
        {
        public:
            explicit PolyTIFFInput( std::string_view filename )
                : PolygonalSurfaceInput3D( filename )
            {
            }

            static std::vector< std::string > extensions()
            {
                static const std::vector< std::string > extensions{ "tiff",
                    "tif" };
                return extensions;
            }

            std::unique_ptr< PolygonalSurface3D > read(
                const MeshImpl& impl ) final;
        };
    } // namespace internal
} // namespace geode