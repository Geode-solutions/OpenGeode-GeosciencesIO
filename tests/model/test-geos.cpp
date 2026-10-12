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

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>
#include <geode/tests_config.hpp>

#include <absl/algorithm/container.h>
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

std::string read_file( const std::string& filename )
{
    std::ifstream file{ filename };
    geode::OpenGeodeGeosciencesIOModelException::test(
        file.good(), "[Test] Cannot open file ", filename );
    std::stringstream content;
    content << file.rdbuf();
    return content.str();
}

std::vector< double > read_table_values( const std::string& filename )
{
    std::istringstream content{ read_file( filename ) };
    std::vector< double > values;
    double value{ 0. };
    while( content >> value )
    {
        values.push_back( value );
    }
    return values;
}

geode::index_t count_occurrences(
    std::string_view content, std::string_view pattern )
{
    geode::index_t nb_occurrences{ 0 };
    for( auto position = content.find( pattern );
        position != std::string_view::npos;
        position = content.find( pattern, position + pattern.size() ) )
    {
        nb_occurrences++;
    }
    return nb_occurrences;
}

constexpr std::array< std::string_view, 3 > AXIS_NAMES{ "x", "y", "z" };
constexpr double PRESSURE_TOLERANCE{ 1e-3 };

// Mimics GEOS TableFunction nearest lookup along one axis, clamped outside
// the table bounds, ties going to the lower node. See
// https://github.com/GEOS-DEV/GEOS/blob/develop/src/coreComponents/functions/TableFunction.hpp
std::size_t nearest_index(
    absl::Span< const double > coordinates, double input )
{
    if( input <= coordinates.front() )
    {
        return 0;
    }
    if( input >= coordinates.back() )
    {
        return coordinates.size() - 1;
    }
    const auto upper = static_cast< std::size_t >(
        std::lower_bound( coordinates.begin(), coordinates.end(), input )
        - coordinates.begin() );
    if( input - coordinates[upper - 1] <= coordinates[upper] - input )
    {
        return upper - 1;
    }
    return upper;
}

void test_nearest_index()
{
    const std::array< double, 3 > coordinates{ 0., 1., 3. };
    const std::array< std::pair< double, std::size_t >, 9 > expected_indices{
        { { -1., 0 }, { 0., 0 }, { 0.5, 0 }, { 0.6, 1 }, { 1., 1 }, { 2., 1 },
            { 2.5, 2 }, { 3., 2 }, { 4., 2 } }
    };
    for( const auto& [input, expected_index] : expected_indices )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            nearest_index( coordinates, input ) == expected_index,
            "[Test] Wrong nearest table index for input ", input );
    }
}

// Mimics GEOS TableFunction nearest interpolation, values in Fortran order.
// The picasso Surfaces are triangulated, so GEOS face centers are the polygon
// barycenters.
double table_value( const std::array< std::vector< double >, 3 >& coordinates,
    absl::Span< const double > values,
    const geode::Point3D& point )
{
    const auto nb_x = coordinates[0].size();
    const auto nb_xy = nb_x * coordinates[1].size();
    return values[nearest_index( coordinates[0], point.value( 0 ) )
                  + nearest_index( coordinates[1], point.value( 1 ) ) * nb_x
                  + nearest_index( coordinates[2], point.value( 2 ) ) * nb_xy];
}

void test_uniform_variable_boundary_condition()
{
    auto model =
        geode::load_brep( absl::StrCat( geode::DATA_PATH, "picasso.og_brep" ) );
    const auto surfaces = east_surfaces( model );
    add_boundary_condition< geode::VariableAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "pressure", 1e7 );
    geode::BRepGeosExporter exporter( model, "picasso_uniform" );
    exporter.run();
    const auto xml =
        read_file( "picasso_uniform/picasso_uniform_simulation.xml" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_occurrences( xml, "<TableFunction" ) == 0
            && count_occurrences( xml, "objectPath=\"faceManager\"" ) == 1,
        "[Test] Uniform VariableAttribute should be exported as a single "
        "scaled boundary condition" );
}

std::pair< double, double > set_depth_pressures( geode::BRep& model,
    absl::Span< const geode::uuid > surfaces,
    const geode::uuid& attribute_id )
{
    static constexpr double TOP_PRESSURE{ 1e7 };
    static constexpr double PRESSURE_GRADIENT{ 1e4 };
    const auto top_z = model.bounding_box().max().value( 2 );
    double min_pressure{ std::numeric_limits< double >::max() };
    double max_pressure{ std::numeric_limits< double >::lowest() };
    for( const auto& surface_id : surfaces )
    {
        const auto& mesh = model.surface( surface_id ).mesh();
        auto attribute =
            mesh.polygon_attribute_manager()
                .find_attribute< geode::VariableAttribute, double >(
                    attribute_id );
        for( const auto polygon_id : geode::Range{ mesh.nb_polygons() } )
        {
            const auto depth =
                top_z - mesh.polygon_barycenter( polygon_id ).value( 2 );
            const auto pressure = TOP_PRESSURE + PRESSURE_GRADIENT * depth;
            attribute->set_value( polygon_id, pressure );
            min_pressure = std::min( min_pressure, pressure );
            max_pressure = std::max( max_pressure, pressure );
        }
    }
    return { min_pressure, max_pressure };
}

void check_surface_pressure_table( const geode::Surface3D& surface,
    geode::index_t region_id,
    const geode::uuid& attribute_id,
    std::pair< double, double > pressure_range )
{
    const auto table_prefix =
        absl::StrCat( "picasso_variable/boundary_tables/surface", region_id );
    std::array< std::vector< double >, 3 > coordinates;
    std::size_t nb_grid_nodes{ 1 };
    for( const auto axis : geode::LRange{ 3 } )
    {
        coordinates[axis] = read_table_values(
            absl::StrCat( table_prefix, "_", AXIS_NAMES[axis], "lin.geos" ) );
        nb_grid_nodes *= coordinates[axis].size();
    }
    const auto values = read_table_values(
        absl::StrCat( "picasso_variable/boundary_tables/boundaryPressure_"
                      "surface",
            region_id, ".geos" ) );
    geode::OpenGeodeGeosciencesIOModelException::test(
        values.size() == nb_grid_nodes,
        "[Test] Wrong number of table values for Surface region ", region_id );
    const auto [min_pressure, max_pressure] = pressure_range;
    for( const auto value : values )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            value >= min_pressure - PRESSURE_TOLERANCE
                && value <= max_pressure + PRESSURE_TOLERANCE,
            "[Test] Table value ", value, " is out of the Surface range" );
    }
    const auto& mesh = surface.mesh();
    const auto attribute =
        mesh.polygon_attribute_manager().find_read_only_attribute< double >(
            attribute_id );
    for( const auto polygon_id : geode::Range{ mesh.nb_polygons() } )
    {
        const auto table_pressure = table_value(
            coordinates, values, mesh.polygon_barycenter( polygon_id ) );
        geode::OpenGeodeGeosciencesIOModelException::test(
            table_pressure == attribute->value( polygon_id ),
            "[Test] Boundary pressure table of Surface region ", region_id,
            " does not give the exact value of polygon ", polygon_id );
    }
}

void test_variable_boundary_condition()
{
    auto model =
        geode::load_brep( absl::StrCat( geode::DATA_PATH, "picasso.og_brep" ) );
    const auto surfaces = east_surfaces( model );
    const auto attribute_id =
        add_boundary_condition< geode::VariableAttribute >( model, surfaces,
            geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "pressure", 0. );
    const auto pressure_range =
        set_depth_pressures( model, surfaces, attribute_id );
    geode::BRepGeosExporter exporter( model, "picasso_variable" );
    exporter.run();

    const auto xml =
        read_file( "picasso_variable/picasso_variable_simulation.xml" );
    const auto nb_tables = count_occurrences( xml, "<TableFunction" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        nb_tables == surfaces.size()
            && count_occurrences( xml, "faceCenter" ) == nb_tables
            && count_occurrences( xml, "functionName=" ) == nb_tables,
        "[Test] Wrong number of boundary condition tables: ", nb_tables );
    // Surface regions are numbered after the Block regions, following the
    // Surface order of the model
    auto region_id = model.nb_blocks();
    for( const auto& surface : model.surfaces() )
    {
        const auto current_region = region_id++;
        if( absl::c_find( surfaces, surface.id() ) != surfaces.end() )
        {
            check_surface_pressure_table(
                surface, current_region, attribute_id, pressure_range );
        }
    }
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
        test_nearest_index();
        test_uniform_variable_boundary_condition();
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
