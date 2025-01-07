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

#include <geode/basic/logger.hpp>
#include <geode/tests_config.hpp>

#include <absl/strings/str_cat.h>

#include <geode/geosciences_io/model/helpers/brep_geos_export.hpp>

#include <geode/geometry/point.hpp>
#include <geode/io/mesh/common.hpp>
#include <geode/io/model/common.hpp>
#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/core/point_set.hpp>

#include <geode/mesh/core/geode/geode_point_set.hpp>

#include <geode/model/representation/core/brep.hpp>
#include <geode/model/representation/io/brep_input.hpp>

void test_picasso()
{
    // Load structural model
    auto model =
        geode::load_brep( absl::StrCat( geode::DATA_PATH, "picasso.og_brep" ) );
    geode::BRepGeosExporter exporter( model, "picasso" );
    exporter.run();
}
void toy_model()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "adaptive_brep_perm_and_poro.og_brep" ) );
    geode::BRepGeosExporter exporter( model, "toy_model" );
    exporter.add_cell_property_1d( "permeability" );
    exporter.add_cell_property_1d( "porosity" );
    auto point_set = geode::PointSet3D::create(
        geode::OpenGeodePointSet3D::impl_name_static() );
    auto builder = geode::PointSetBuilder3D::create( *point_set );
    builder->create_point( geode::Point3D{ { 20., 20., 10. } } );
    exporter.add_well_perforations( *point_set );
    exporter.run();
}
int main()
{
    try
    {
        geode::GeosciencesIOModelLibrary::initialize();
        geode::IOMeshLibrary::initialize();
        geode::IOModelLibrary::initialize();

        test_picasso();
        toy_model();
        geode::Logger::info( "TEST SUCCESS" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
