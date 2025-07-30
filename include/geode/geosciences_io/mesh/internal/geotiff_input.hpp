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

#include <geode/mesh/io/light_regular_grid_input.hpp>

namespace geode
{
    namespace internal
    {
        class GEOTIFFInput final : public LightRegularGridInput2D
        {
        public:
            explicit GEOTIFFInput( std::string_view filename )
                : LightRegularGridInput2D( filename )
            {
            }

            static std::vector< std::string > extensions()
            {
                static const std::vector< std::string > extensions{ "tiff",
                    "tif" };
                return extensions;
            }

            LightRegularGrid< 2 > read() final;

            // TODO: implement missing files
            // A tiff can use completementary file to open... but actually I
            // donnot use it:
            //* Sidecar Metadata Files: TIFF files may be accompanied by XMP
            //(.xmp) or other sidecar files containing metadata, especially in
            // workflows that use Adobe software or other image management
            // tools.
            //* Associated Color Profiles: TIFF images may rely on ICC color
            // profiles (.icc or .icm) to ensure color consistency across
            // different devices.
            //* Layered or Multi-Page TIFFs: Some TIFF files store multiple
            // images (e.g., multi-page scans). In such cases, associated files
            // might be used to store extracted or processed versions of
            // individual pages.
            //* Auxiliary Data Files: In geospatial imaging (GeoTIFF),
            // additional files such as .tfw (world file) or .prj may provide
            // spatial reference data.
            //* Thumbnails or Proxy Files: Some workflows generate
            // low-resolution preview images (.jpg or .png) for quick browsing
            // without opening large TIFF files.
            //
            // LightRegularGridInput2D::AdditionalFiles additional_files()
            // const final;

            bool is_loadable() const final;

            AdditionalFiles additional_files() const final;

            index_t object_priority() const final
            {
                return 1;
            }
        };
    } // namespace internal
} // namespace geode