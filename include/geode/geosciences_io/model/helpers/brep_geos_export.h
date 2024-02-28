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

#include <geode/geosciences_io/model/common.h>

#include <absl/strings/string_view.h>

#include <geode/basic/pimpl.h>

namespace geode
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSet );
    ALIAS_3D( PointSet );
    class BRep;
} // namespace geode

namespace geode
{
    class opengeode_geosciencesio_model_api BRepGeosExporter
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( BRepGeosExporter );

    public:
        BRepGeosExporter() = delete;
        BRepGeosExporter( const BRep& brep, absl::string_view files_directory );
        ~BRepGeosExporter();

        void add_well_perforations( const PointSet3D& well_perforations );
        void add_cell_property_1d( absl::string_view name );
        void add_cell_property_2d( absl::string_view name );
        void add_cell_property_3d( absl::string_view name );

        void run();

    private:
        IMPLEMENTATION_MEMBER( impl_ );
    };
} // namespace geode