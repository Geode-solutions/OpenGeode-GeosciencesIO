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

#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/variable_attribute.hpp>
#include <geode/tests_config.hpp>

#include <absl/strings/str_cat.h>
#include <absl/types/span.h>

#include <geode/geosciences_io/model/helpers/brep_geos_export.hpp>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>
#include <geode/io/mesh/common.hpp>
#include <geode/io/model/common.hpp>
#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/core/point_set.hpp>

#include <geode/mesh/core/geode/geode_point_set.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/builder/brep_builder.hpp>
#include <geode/model/representation/core/brep.hpp>
#include <geode/model/representation/io/brep_input.hpp>

std::vector< geode::uuid > east_surfaces( const geode::BRep& model )
{
    static constexpr double TOLERANCE{ 1e-3 };
    const auto east_x = model.bounding_box().max().value( 0 );
    std::vector< geode::uuid > surfaces;
    for( const auto& surface : model.surfaces() )
    {
        const auto box = surface.mesh().bounding_box();
        if( std::fabs( box.min().value( 0 ) - east_x ) < TOLERANCE
            && std::fabs( box.max().value( 0 ) - east_x ) < TOLERANCE )
        {
            surfaces.push_back( surface.id() );
        }
    }
    return surfaces;
}

template < template < typename > class Attribute, typename T >
geode::uuid add_boundary_condition( geode::BRep& model,
    absl::Span< const geode::uuid > surfaces,
    geode::PHYSICAL_PROPERTY_NAME property,
    std::string_view name,
    T value )
{
    const geode::uuid attribute_id;
    geode::AttributeValues< T > values;
    values.default_value = value;
    for( const auto& surface_id : surfaces )
    {
        model.surface( surface_id )
            .mesh()
            .polygon_attribute_manager()
            .create_attribute< Attribute, T >( name, attribute_id, values, {} );
    }
    geode::BRepBuilder{ model }.set_physical_property(
        property, geode::Surface3D::component_type_static(), attribute_id );
    return attribute_id;
}

void test_picasso()
{
    // Load structural model
    auto model =
        geode::load_brep( absl::StrCat( geode::DATA_PATH, "picasso.og_brep" ) );
    const auto surfaces = east_surfaces( model );
    geode::OpenGeodeGeosciencesIOModelException::test(
        !surfaces.empty(), "[Test] No east Surface found" );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "pressure", 1e7 );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_temperature, "temperature",
        350. );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_oil_fraction, "oil_fraction",
        0.9995 );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction,
        "water_fraction", 0.0005 );
    geode::BRepGeosExporter exporter( model, "picasso" );
    exporter.run();
}

void test_variable_boundary_condition()
{
    auto model =
        geode::load_brep( absl::StrCat( geode::DATA_PATH, "picasso.og_brep" ) );
    const auto surfaces = east_surfaces( model );
    add_boundary_condition< geode::VariableAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "pressure", 1e7 );
    geode::BRepGeosExporter exporter( model, "picasso_variable" );
    try
    {
        exporter.run();
    }
    catch( const geode::OpenGeodeException& exception )
    {
        geode::Logger::info( "[Test] Expected error: ", exception.what() );
        return;
    }
    throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
        geode::OpenGeodeException::TYPE::internal,
        "[Test] Variable boundary condition should not be exported" };
}
void toy_model()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "adaptive_brep_perm_and_poro.og_brep" ) );
    geode::BRepGeosExporter exporter( model, "toy_model" );
    auto point_set = geode::PointSet3D::create(
        geode::OpenGeodePointSet3D::impl_name_static() );
    auto builder = geode::PointSetBuilder3D::create( *point_set );
    builder->create_point( geode::Point3D{ { 20., 20., 10. } } );
    exporter.add_well_perforations( *point_set, "well" );
    exporter.run();
}

void add_vertical_well( geode::BRepGeosExporter& exporter,
    std::string_view name,
    double x,
    double y )
{
    static constexpr std::array< double, 2 > LAYER_CENTERS_Z{ 0.305, 0.915 };
    auto point_set = geode::PointSet3D::create(
        geode::OpenGeodePointSet3D::impl_name_static() );
    auto builder = geode::PointSetBuilder3D::create( *point_set );
    for( const auto z : LAYER_CENTERS_Z )
    {
        builder->create_point( geode::Point3D{ { x, y, z } } );
    }
    exporter.add_well_perforations( *point_set, name );
}

void add_spe10_wells( geode::BRepGeosExporter& exporter )
{
    add_vertical_well( exporter, "source", 185.93, 336.8 );
    add_vertical_well( exporter, "sink1", 3.048, 1.524 );
    add_vertical_well( exporter, "sink2", 3.048, 669.036 );
    add_vertical_well( exporter, "sink3", 362.712, 1.524 );
    add_vertical_well( exporter, "sink4", 362.712, 669.036 );
}

void test_grid_geos()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) );
    geode::BRepGeosExporter exporter( model, "grid_geos" );
    add_spe10_wells( exporter );
    exporter.run();
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        geode::OpenGeodeIOMeshLibrary::initialize();
        geode::OpenGeodeIOModelLibrary::initialize();
        test_picasso();
        test_variable_boundary_condition();
        toy_model();
        test_grid_geos();
        geode::Logger::info( "TEST SUCCESS" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
