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
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/logger_client.hpp>
#include <geode/basic/logger_manager.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>
#include <geode/tests_config.hpp>

#include <absl/algorithm/container.h>
#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>
#include <absl/types/span.h>

#include <geode/geosciences/simulation/builder/physical_problem_builder.hpp>
#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/schedule.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>
#include <geode/geosciences/simulation/helpers/physical_problem_simulation.hpp>

#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.hpp>
#include <geode/geosciences/explicit/representation/builder/structural_model_builder.hpp>
#include <geode/geosciences/explicit/representation/core/structural_model.hpp>

#include <geode/geosciences_io/model/helpers/brep_geos_export.hpp>
#include <geode/geosciences_io/model/helpers/structural_model_geos_export.hpp>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>
#include <geode/io/mesh/common.hpp>
#include <geode/io/model/common.hpp>
#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/builder/regular_grid_solid_builder.hpp>
#include <geode/mesh/builder/surface_mesh_builder.hpp>
#include <geode/mesh/core/point_set.hpp>

#include <geode/mesh/core/geode/geode_point_set.hpp>
#include <geode/mesh/core/geode/geode_regular_grid_solid.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/builder/brep_builder.hpp>
#include <geode/model/representation/core/brep.hpp>
#include <geode/model/representation/io/brep_input.hpp>

std::string read_file( const std::string& filename );
geode::index_t count_occurrences(
    std::string_view content, std::string_view pattern );
geode::PhysicalProblem spe10_reference_problem( const geode::BRep& model );
geode::PhysicalProblem toy_problem( bool thermal );

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
    // Oil and water fluid, without the source term well
    auto problem = spe10_reference_problem( model );
    geode::PhysicalProblemBuilder{ problem }.remove_source_term( 0 );
    geode::BRepGeosExporter exporter( model, problem, "picasso" );
    exporter.run();
    const auto deck = read_file( "picasso/picasso_simulation.xml" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_occurrences( deck, "boundaryComposition_oil" ) == 1
            && count_occurrences( deck, "boundaryComposition_water" ) == 1
            && !absl::StrContains( deck, "component=\"2\"" ),
        "[Test] Boundary fractions must use the components of the problem "
        "fluid" );
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
    const auto problem = toy_problem( false );
    geode::BRepGeosExporter exporter( model, problem, "picasso_uniform" );
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
    const auto problem = toy_problem( false );
    geode::BRepGeosExporter exporter( model, problem, "picasso_variable" );
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

std::shared_ptr< geode::DeadOilFluid > three_phase_dead_oil()
{
    geode::HydrocarbonPvtTable oil_table;
    oil_table.phase = geode::PHASE::oil;
    oil_table.pressures = { 1e5, 2e7, 6e7 };
    oil_table.formation_volume_factors = { 1.02, 1.01, 0.99 };
    oil_table.viscosities = { 2.9e-3, 3e-3, 3.1e-3 };
    geode::HydrocarbonPvtTable gas_table;
    gas_table.phase = geode::PHASE::gas;
    gas_table.pressures = { 1e5, 2e7, 6e7 };
    gas_table.formation_volume_factors = { 0.9, 0.006, 0.0025 };
    gas_table.viscosities = { 1e-5, 1.5e-5, 2.5e-5 };
    return std::make_shared< geode::DeadOilFluid >(
        std::vector< geode::PHASE >{
            geode::PHASE::oil, geode::PHASE::gas, geode::PHASE::water },
        std::vector< geode::FluidComponent >{
            { "oil", 114e-3 }, { "gas", 16e-3 }, { "water", 18e-3 } },
        std::vector< double >{ 800.0, 0.9907, 1022.0 },
        std::vector< geode::HydrocarbonPvtTable >{ oil_table, gas_table },
        geode::WaterPvt{ 30600000.1, 1.03, 4.1e-10, 3e-4 } );
}

geode::Material spe10_material( std::string name, double compressibility )
{
    auto relperm = std::make_shared< geode::BrooksCoreyRelativePermeability >(
        std::vector< geode::PHASE >{
            geode::PHASE::oil, geode::PHASE::gas, geode::PHASE::water },
        std::vector< double >{ 0., 0., 0. },
        std::vector< double >{ 2., 2., 2. },
        std::vector< double >{ 1., 1., 1. } );
    geode::Material material;
    material.name = std::move( name );
    material.pore_compressibility = compressibility;
    material.reference_pressure = 1e7;
    material.relative_permeability = relperm;
    return material;
}

geode::PhysicalProblem spe10_problem( const geode::BRep& model )
{
    geode::PhysicalProblem problem;
    geode::PhysicalProblemBuilder builder{ problem };
    builder.set_name( "spe10" );
    builder.set_model_id( model.id() );
    auto flow = std::make_shared< geode::MultiphaseFlow >(
        geode::COMPOSITION_BASIS::mass );
    flow->set_fluid( three_phase_dead_oil() );
    builder.set_flow( flow );
    builder.add_material( spe10_material( "reservoir", 1e-10 ) );
    builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{}, 4.1369e7, 300.,
            std::vector< double >{ 0.9995, 0., 0.0005 } ) );
    builder.set_schedule( { 0., 2e6, 1e3, 5e5 } );
    geode::OutputSettings output;
    output.frequency = 5e5;
    builder.set_output( output );
    builder.set_solver_option( "geos", "nonlinear/newtonTol", "1e-4" );
    return problem;
}

std::vector< geode::uuid > all_surfaces( const geode::BRep& model )
{
    std::vector< geode::uuid > surfaces;
    for( const auto& surface : model.surfaces() )
    {
        surfaces.push_back( surface.id() );
    }
    return surfaces;
}

// The only Surface of the SPE10 grid is its west face
void add_spe10_boundary_conditions( geode::BRep& model )
{
    const auto surfaces = all_surfaces( model );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure",
        4.5e7 );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_oil_fraction, "bc_oil", 0.1 );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction, "bc_water",
        0.9 );
}

void add_spe10_gas_fraction( geode::BRep& model )
{
    add_boundary_condition< geode::ConstantAttribute >( model,
        all_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_gas_fraction, "bc_gas", 0. );
}

// Values of the GEOS reference deck
// inputFiles/compositionalMultiphaseFlow/benchmarks/SPE10/
// deadOilSpe10Layers84_85_base_iterative.xml, without its wells
std::shared_ptr< geode::DeadOilFluid > spe10_reference_fluid()
{
    // pres_pvdo.txt, B_o_pvdo.txt and visc_pvdo.txt
    geode::HydrocarbonPvtTable oil_table;
    oil_table.phase = geode::PHASE::oil;
    oil_table.pressures = { 2.068e6, 5.516e6, 5.516e7 };
    oil_table.formation_volume_factors = { 1.05, 1.02, 1.01 };
    oil_table.viscosities = { 0.00285, 0.00299, 0.003 };
    return std::make_shared< geode::DeadOilFluid >(
        std::vector< geode::PHASE >{ geode::PHASE::oil, geode::PHASE::water },
        std::vector< geode::FluidComponent >{
            { "oil", 114e-3 }, { "water", 18e-3 } },
        std::vector< double >{ 800.0, 1022.0 },
        std::vector< geode::HydrocarbonPvtTable >{ oil_table },
        geode::WaterPvt{ 30600000.1, 1.03, 4.1e-10, 3e-4 } );
}

geode::Material spe10_reference_material()
{
    auto relperm = std::make_shared< geode::BrooksCoreyRelativePermeability >(
        std::vector< geode::PHASE >{ geode::PHASE::oil, geode::PHASE::water },
        std::vector< double >{ 0., 0. }, std::vector< double >{ 2., 2. },
        std::vector< double >{ 1., 1. } );
    geode::Material material;
    material.name = "rock";
    material.pore_compressibility = 1e-10;
    material.reference_pressure = 1e7;
    material.porosity = 0.1;
    material.permeability = std::array< double, 3 >{ 1e-14, 1e-14, 1e-18 };
    material.relative_permeability = relperm;
    return material;
}

geode::PhysicalProblem spe10_reference_problem( const geode::BRep& model )
{
    geode::PhysicalProblem problem;
    geode::PhysicalProblemBuilder builder{ problem };
    builder.set_name( "spe10_reference" );
    builder.set_model_id( model.id() );
    auto flow = std::make_shared< geode::MultiphaseFlow >(
        geode::COMPOSITION_BASIS::mass );
    flow->set_fluid( spe10_reference_fluid() );
    builder.set_flow( flow );
    builder.add_material( spe10_reference_material() );
    builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{}, 4.1369e7, 300.,
            std::vector< double >{ 0.9995, 0.0005 } ) );
    builder.set_schedule( { 0., 2e6, 1e3, 5e5 } );
    geode::OutputSettings output;
    output.frequency = 5e5;
    output.restart_frequency = 1e6;
    builder.set_output( output );
    builder.set_solver_option( "geos", "flow/maxCompFractionChange", "0.1" );
    builder.set_solver_option( "geos", "nonlinear/newtonTol", "1e-4" );
    builder.set_solver_option( "geos", "nonlinear/newtonMaxIter", "40" );
    builder.set_solver_option( "geos", "nonlinear/lineSearchAction", "None" );
    // Water injected in the central well, 10 kg/s
    builder.add_source_term( { "source", "water", 10. } );
    return problem;
}

// Boundary conditions of the reference grid_geos_simulation.xml
void add_spe10_reference_boundary_conditions( geode::BRep& model )
{
    const auto surfaces = all_surfaces( model );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure",
        2.7579e7 );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_temperature, "bc_temperature",
        300. );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_oil_fraction, "bc_oil",
        0.9995 );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction, "bc_water",
        0.0005 );
}

void test_spe10_reference()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) );
    add_spe10_reference_boundary_conditions( model );
    const auto problem = spe10_reference_problem( model );
    geode::BRepGeosExporter exporter( model, problem, "spe10_reference" );
    add_spe10_wells( exporter );
    exporter.run();
    const auto deck =
        read_file( "spe10_reference/spe10_reference_simulation.xml" );
    for( const auto* pattern :
        { "phaseNames=\"{ oil, water }\"", "\"pvt_oil_fvf\"",
            "maxCompFractionChange=\"0.1\"", "newtonMaxIter=\"40\"",
            "lineSearchAction=\"None\"", "newtonTol=\"1e-4\"",
            "defaultReferencePorosity=\"0.1\"",
            "name=\"boundaryTemperature_0\"",
            "permeabilityComponents=\"{ 1e-14, 1e-14, 1e-18 }\"",
            "name=\"restarts\"", "timeFrequency=\"1000000\"",
            "targetExactTimestep=\"0\"", "<Restart", "<SourceFlux",
            "name=\"sourceFlux_source_water\"", "objectPath=\"ElementRegions\"",
            "setNames=\"{ source }\"", "scale=\"-10\"" } )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::StrContains( deck, pattern ), "[Test] ", pattern,
            " missing in the SPE10 reference deck" );
    }
    // Water is the second component, in the initial and boundary conditions
    // and in the source term
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_occurrences( deck, "boundaryComposition_" ) == 2
            && count_occurrences( deck, "component=\"1\"" ) == 3
            && !absl::StrContains( deck, "component=\"2\"" )
            && !absl::StrContains( deck, "pvt_gas" )
            && !absl::StrContains( deck, "isothermal" ),
        "[Test] Wrong components in the SPE10 reference deck" );
    // Perforation boxes of the reference grid_geos_simulation.xml
    for( const auto* pattern : { "name=\"source\"", "name=\"sink1\"",
             "name=\"sink2\"", "name=\"sink3\"", "name=\"sink4\"",
             "xMin=\"{182.88, 335.28, -0.0001}\"",
             "xMax=\"{188.976, 338.328, 1.2201}\"" } )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::StrContains( deck, pattern ), "[Test] ", pattern,
            " missing in the SPE10 reference deck" );
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_occurrences( deck, "<Box" ) == 5,
        "[Test] Wrong number of wells in the SPE10 reference deck" );
}

void test_toy_model_thermal()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "adaptive_brep_perm_and_poro.og_brep" ) );
    const auto surfaces = east_surfaces( model );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure", 2e7 );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_temperature, "bc_temperature",
        330. );

    geode::PhysicalProblem problem;
    geode::PhysicalProblemBuilder builder{ problem };
    auto fluid = std::make_shared< geode::SinglePhaseFluid >(
        std::vector< geode::PHASE >{ geode::PHASE::water },
        std::vector< geode::FluidComponent >{ { "water", 18e-3 } }, 1000., 1e-3,
        5e-10, 0., 1e7, geode::FluidThermalProperties{ 4180., 2e-4, 0. } );
    auto flow = std::make_shared< geode::SinglePhaseFlow >();
    flow->set_fluid( fluid );
    builder.set_flow( flow );
    builder.set_thermal( geode::ThermalParameters{ 300. } );
    geode::Material rock;
    rock.name = "rock";
    rock.pore_compressibility = 1e-9;
    rock.reference_pressure = 1e7;
    rock.porosity = 0.2;
    rock.permeability = std::array< double, 3 >{ 1e-13, 1e-13, 1e-14 };
    rock.thermal = geode::RockThermalProperties{ 2e6, { 2.5, 2.5, 2.5 } };
    builder.add_material( rock );
    builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{}, 1e7, 300. ) );
    // Small first step: the boundary temperature jump is stiff on small cells
    builder.set_schedule( { 0., 1e5, 10., 1e4 } );
    geode::OutputSettings output;
    output.frequency = 5e4;
    builder.set_output( output );

    geode::BRepGeosExporter exporter( model, problem, "toy_model_thermal" );
    exporter.run();
    const auto deck =
        read_file( "toy_model_thermal/toy_model_thermal_simulation.xml" );
    for( const auto* pattern : { "<SinglePhaseFVM",
             "<ThermalCompressibleSinglePhaseFluid", "isThermal=\"1\"",
             "<SolidInternalEnergy", "<SinglePhaseThermalConductivity",
             "thermalCond_rock", "boundaryTemperature" } )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::StrContains( deck, pattern ), "[Test] Thermal deck misses ",
            pattern );
    }
}

void check_export_failure( const geode::BRep& model,
    const geode::PhysicalProblem& problem,
    std::string_view keyword )
{
    std::filesystem::remove_all( "grid_geos_failure" );
    geode::BRepGeosExporter exporter( model, problem, "grid_geos_failure" );
    try
    {
        exporter.run();
    }
    catch( const geode::OpenGeodeException& exception )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::StrContains( exception.what(), keyword ),
            "[Test] Export error does not mention \"", keyword,
            "\": ", exception.what() );
        geode::OpenGeodeGeosciencesIOModelException::test(
            !std::filesystem::exists( "grid_geos_failure" ),
            "[Test] A failed export should not create its directory" );
        return;
    }
    throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
        geode::OpenGeodeException::TYPE::internal, "[Test] Export with ",
        keyword, " should fail" };
}

class TestFluid final : public geode::FluidModel
{
public:
    TestFluid( std::vector< geode::PHASE > phases,
        std::vector< geode::FluidComponent > components )
        : geode::FluidModel( std::move( phases ), std::move( components ) )
    {
    }

    [[nodiscard]] std::string_view type_name() const final
    {
        return "TestFluid";
    }
};

bool has_issue(
    absl::Span< const std::string > issues, std::string_view keyword )
{
    return absl::c_any_of( issues, [keyword]( const std::string& issue ) {
        return absl::StrContains( issue, keyword );
    } );
}

void check_geos_error( const geode::PhysicalProblem& problem,
    const geode::BRep& model,
    std::string_view keyword )
{
    const auto report =
        geode::check_physical_problem_simulation( problem, model, "geos" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        has_issue( report.errors, keyword ), "[Test] No GEOS error with \"",
        keyword, "\": ", report.string() );
}

void check_no_geos_error( const geode::PhysicalProblem& problem,
    const geode::BRep& model,
    std::string_view keyword )
{
    const auto report =
        geode::check_physical_problem_simulation( problem, model, "geos" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        !has_issue( report.errors, keyword ),
        "[Test] Unexpected GEOS error "
        "with \"",
        keyword, "\": ", report.string() );
}

void check_geos_warning( const geode::PhysicalProblem& problem,
    const geode::BRep& model,
    std::string_view keyword )
{
    const auto report =
        geode::check_physical_problem_simulation( problem, model, "geos" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        has_issue( report.warnings, keyword ), "[Test] No GEOS warning with \"",
        keyword, "\": ", report.string() );
}

void check_no_geos_warning( const geode::PhysicalProblem& problem,
    const geode::BRep& model,
    std::string_view keyword )
{
    const auto report =
        geode::check_physical_problem_simulation( problem, model, "geos" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        !has_issue( report.warnings, keyword ),
        "[Test] Unexpected GEOS warning with \"", keyword,
        "\": ", report.string() );
}

void test_geos_rules()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) );
    add_spe10_boundary_conditions( model );
    add_spe10_gas_fraction( model );
    {
        auto problem = spe10_problem( model );
        auto flow = std::make_shared< geode::MultiphaseFlow >();
        auto fluid = std::make_shared< TestFluid >(
            std::vector< geode::PHASE >{
                geode::PHASE::oil, geode::PHASE::gas, geode::PHASE::water },
            std::vector< geode::FluidComponent >{
                { "oil", 0.1 }, { "gas", 0.01 }, { "water", 0.018 } } );
        flow->set_fluid( fluid );
        geode::PhysicalProblemBuilder{ problem }.set_flow( flow );
        check_export_failure( model, problem, "TestFluid is not supported" );
    }
    {
        auto problem = spe10_problem( model );
        auto material = problem.materials()[0];
        material.regions = { geode::uuid{} };
        geode::PhysicalProblemBuilder{ problem }.set_material( 0, material );
        check_export_failure( model, problem, "is not a Block" );
    }
    {
        auto problem = spe10_problem( model );
        geode::PhysicalProblemBuilder{ problem }.add_source_term(
            { "missing", "water", 1. } );
        check_export_failure(
            model, problem, "has no perforations given to the exporter" );
    }
    {
        auto problem = spe10_problem( model );
        geode::PhysicalProblemBuilder builder{ problem };
        builder.set_thermal( geode::ThermalParameters{ 300. } );
        builder.add_source_term( { "missing", "water", 1. } );
        check_export_failure( model, problem, "carry no energy" );
    }
    {
        // Without its source term, which needs wells given to the exporter.
        // The gas fraction of the model is ignored by an oil-water fluid.
        auto problem = spe10_reference_problem( model );
        geode::PhysicalProblemBuilder{ problem }.remove_source_term( 0 );
        check_geos_warning(
            problem, model, "boundary_gas_fraction is ignored by GEOS" );
        check_no_geos_error( problem, model, "has no gas phase" );
    }
}

const geode::Block3D& first_block( const geode::BRep& model )
{
    for( const auto& block : model.blocks() )
    {
        return block;
    }
    throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
        geode::OpenGeodeException::TYPE::internal, "[Test] No Block" };
}

void test_structural_model_regions()
{
    geode::StructuralModel model{ geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "adaptive_brep_perm_and_poro.og_brep" ) ) };
    add_spe10_boundary_conditions( model );
    add_spe10_gas_fraction( model );
    geode::StructuralModelBuilder model_builder{ model };
    const auto& unit_id = model_builder.add_stratigraphic_unit();
    const auto& unit = model.stratigraphic_unit( unit_id );
    for( const auto& block : model.blocks() )
    {
        model_builder.add_block_in_stratigraphic_unit( block, unit );
    }
    auto problem = spe10_problem( model );
    // The model has no porosity and permeability properties
    auto reservoir = problem.materials()[0];
    reservoir.porosity = 0.2;
    reservoir.permeability = std::array< double, 3 >{ 1e-13, 1e-13, 1e-14 };
    geode::PhysicalProblemBuilder{ problem }.set_material( 0, reservoir );
    auto unit_material = spe10_material( "unit", 1e-10 );
    unit_material.porosity = reservoir.porosity;
    unit_material.permeability = reservoir.permeability;
    unit_material.regions = { first_block( model ).id() };
    const auto unit_material_id =
        geode::PhysicalProblemBuilder{ problem }.add_material( unit_material );
    const auto report =
        geode::check_physical_problem_simulation( problem, model, "geos" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        report.is_simulable()
            && absl::c_any_of( report.warnings,
                []( const std::string& warning ) {
                    return absl::StrContains(
                        warning, "covers only part of the GEOS region" );
                } ),
        "[Test] The public check should only warn about a Block merged in a "
        "StratigraphicUnit region: ",
        report.string() );
    bool has_thrown{ false };
    try
    {
        geode::StructuralModelGeosExporter exporter(
            model, problem, "structural_model_failure" );
        exporter.run();
    }
    catch( const geode::OpenGeodeException& exception )
    {
        has_thrown = true;
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::StrContains(
                exception.what(), "covers only part of the GEOS region" ),
            "[Test] A Block merged in a StratigraphicUnit region should be "
            "reported: ",
            exception.what() );
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        has_thrown, "[Test] Export with a merged Block region should fail" );

    unit_material.regions = { unit_id };
    geode::PhysicalProblemBuilder{ problem }.set_material(
        unit_material_id, unit_material );
    geode::StructuralModelGeosExporter exporter(
        model, problem, "structural_model_regions" );
    exporter.run();
    const auto deck = read_file(
        "structural_model_regions/structural_model_regions_simulation.xml" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        absl::StrContains( deck, "rock_unit" )
            && count_occurrences( deck, "<CellElementRegion" ) == 1,
        "[Test] StructuralModel deck should merge the unit Blocks" );
}

geode::BRep spe10_reference_model()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) );
    add_spe10_reference_boundary_conditions( model );
    return model;
}

template < typename T >
geode::uuid add_block_property( geode::BRep& model,
    geode::PHYSICAL_PROPERTY_NAME property,
    std::string_view name,
    T value,
    geode::AttributeProperties properties = {} )
{
    geode::AttributeValues< T > values;
    values.default_value = value;
    const geode::uuid attribute_id;
    for( const auto& block : model.blocks() )
    {
        block.mesh()
            .polyhedron_attribute_manager()
            .create_attribute< geode::VariableAttribute, T >(
                name, attribute_id, values, properties );
    }
    geode::BRepBuilder{ model }.set_physical_property(
        property, geode::Block3D::component_type_static(), attribute_id );
    return attribute_id;
}

void test_geos_checker_problem()
{
    const auto model = spe10_reference_model();
    auto problem = spe10_reference_problem( model );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::SimulationCheckerFactory::has_creator( "geos" ),
        "[Test] GEOS checker should be registered" );
    const auto report =
        geode::check_physical_problem_simulation( problem, model, "GEOS" );
    geode::OpenGeodeGeosciencesIOModelException::test( report.is_simulable(),
        "[Test] SPE10 reference should be simulable: ", report.string() );
    for( const auto* keyword : { "Source term well source is not checked",
             "porosity is only a GEOS default" } )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            has_issue( report.warnings, keyword ), "[Test] Missing warning ",
            keyword, ": ", report.string() );
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        !geode::is_physical_problem_simulable( problem, model, "unknown" ),
        "[Test] Unknown simulator should not be simulable" );

    geode::PhysicalProblemBuilder builder{ problem };
    geode::OutputSettings restart_only;
    restart_only.restart = true;
    builder.set_output( restart_only );
    builder.set_solver_option( "geos", "vtk/plotLevel", "2" );
    builder.set_solver_option( "geos", "nonlinear/bad attribute", "1" );
    check_geos_warning( problem, model, "the problem has no VTK output" );
    check_geos_error( problem, model, "is not an XML attribute name" );
    builder.set_output( geode::OutputSettings{} );
    check_geos_error( problem, model, "The problem has no output" );

    auto named = spe10_reference_problem( model );
    auto fluid = spe10_reference_fluid();
    auto renamed = std::make_shared< geode::DeadOilFluid >(
        std::vector< geode::PHASE >{ geode::PHASE::oil, geode::PHASE::water },
        std::vector< geode::FluidComponent >{
            { "dead oil", 114e-3 }, { "water", 18e-3 } },
        std::vector< double >{ 800.0, 1022.0 },
        std::vector< geode::HydrocarbonPvtTable >{
            fluid->hydrocarbon_pvt().begin(), fluid->hydrocarbon_pvt().end() },
        fluid->water_pvt() );
    auto flow = std::make_shared< geode::MultiphaseFlow >();
    flow->set_fluid( renamed );
    geode::PhysicalProblemBuilder{ named }.set_flow( flow );
    check_geos_error( named, model, "Fluid component name \"dead oil\"" );
}

void test_geos_checker_model()
{
    const auto problem = spe10_reference_problem( spe10_reference_model() );
    {
        auto model = spe10_reference_model();
        add_block_property( model, geode::PHYSICAL_PROPERTY_NAME::permeability,
            "scalar_permeability", 1e-14 );
        check_geos_error(
            problem, model, "must store std::array< double, 3 > values" );
    }
    {
        auto model = spe10_reference_model();
        add_block_property(
            model, geode::PHYSICAL_PROPERTY_NAME::porosity, "attribute", 0.2 );
        check_geos_error( problem, model, "is already used in the GEOS mesh" );
    }
    {
        auto model = spe10_reference_model();
        add_block_property( model, geode::PHYSICAL_PROPERTY_NAME::porosity,
            "bad_porosity", 2. );
        check_geos_error( problem, model, "values outside (0, 1]" );
    }
    {
        auto model = spe10_reference_model();
        static_cast< void >( geode::BRepBuilder{ model }.add_block() );
        check_geos_error( problem, model, "has no polyhedron" );
    }
    {
        auto model = spe10_reference_model();
        geode::BRepBuilder{ model }.set_physical_property(
            geode::PHYSICAL_PROPERTY_NAME::boundary_gas_fraction,
            geode::Surface3D::component_type_static(), geode::uuid{} );
        check_geos_warning(
            problem, model, "No Surface carries a boundary_gas_fraction" );
    }
    {
        auto model = spe10_reference_model();
        auto two_models = spe10_reference_problem( model );
        auto other = spe10_reference_material();
        other.name = "other";
        other.pore_compressibility = 2e-10;
        other.regions = { first_block( model ).id() };
        static_cast< void >(
            geode::PhysicalProblemBuilder{ two_models }.add_material( other ) );
        const auto& porosity_id =
            model
                .physical_property_info(
                    geode::PHYSICAL_PROPERTY_NAME::porosity )
                .attribute_id;
        const auto& block = first_block( model );
        const auto porosity_name = block.mesh()
                                       .polyhedron_attribute_manager()
                                       .find_generic_attribute( porosity_id )
                                       ->name()
                                       .value();
        geode::AttributeValues< double > values;
        static_cast< void >(
            block.mesh()
                .polyhedron_attribute_manager()
                .create_attribute< geode::VariableAttribute, double >(
                    absl::StrCat( porosity_name, "_0" ), values, {} ) );
        // Only the imported arrays are written in the GEOS mesh
        check_no_geos_error(
            two_models, model, "is needed for a copy of the porosity" );
        add_block_property( model, geode::PHYSICAL_PROPERTY_NAME::permeability,
            absl::StrCat( porosity_name, "_1" ),
            std::array< double, 3 >{ 1e-13, 1e-13, 1e-13 } );
        check_geos_error(
            two_models, model, "is needed for a copy of the porosity" );
    }
    {
        auto model = spe10_reference_model();
        geode::AttributeProperties time_step;
        time_step.time = 1.;
        add_block_property( model, geode::PHYSICAL_PROPERTY_NAME::porosity,
            "porosity_step", 0.2, time_step );
        check_geos_error(
            problem, model, "porosity attribute is a time series step" );
    }
    {
        auto model = spe10_reference_model();
        add_block_property(
            model, geode::PHYSICAL_PROPERTY_NAME::porosity, "bad name!", 0.2 );
        check_geos_error( problem, model,
            "Block attribute name \"bad name!\" must only contain" );
    }
    {
        auto model = geode::load_brep( absl::StrCat(
            geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) );
        const auto surfaces = all_surfaces( model );
        add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
            geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure",
            2.7579e7 );
        add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
            geode::PHYSICAL_PROPERTY_NAME::boundary_oil_fraction, "bc_oil",
            0.5 );
        add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
            geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction, "bc_water",
            0.2 );
        check_geos_error(
            problem, model, "boundary fluid fractions do not sum to 1" );
    }
}

geode::PhysicalProblem toy_problem( bool thermal )
{
    geode::PhysicalProblem problem;
    geode::PhysicalProblemBuilder builder{ problem };
    auto fluid = std::make_shared< geode::SinglePhaseFluid >(
        std::vector< geode::PHASE >{ geode::PHASE::water },
        std::vector< geode::FluidComponent >{ { "water", 18e-3 } }, 1000., 1e-3,
        5e-10, 0., 1e7, geode::FluidThermalProperties{ 4180., 2e-4, 0. } );
    auto flow = std::make_shared< geode::SinglePhaseFlow >();
    flow->set_fluid( fluid );
    builder.set_flow( flow );
    geode::Material rock;
    rock.name = "rock";
    rock.porosity = 0.2;
    rock.permeability = std::array< double, 3 >{ 1e-13, 1e-13, 1e-14 };
    if( thermal )
    {
        builder.set_thermal( geode::ThermalParameters{ 300. } );
        rock.thermal = geode::RockThermalProperties{ 2e6, { 2.5, 2.5, 2.5 } };
    }
    builder.add_material( rock );
    builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{}, 1e7, 300. ) );
    builder.set_schedule( { 0., 1e5, 10., 1e4 } );
    geode::OutputSettings output;
    output.frequency = 5e4;
    builder.set_output( output );
    return problem;
}

geode::BRep load_toy_model()
{
    return geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "adaptive_brep_perm_and_poro.og_brep" ) );
}

// Surface outside the Blocks, its vertices have no unique vertex
geode::uuid add_floating_surface( geode::BRep& model,
    absl::Span< const geode::Point3D > points,
    absl::Span< const std::array< geode::index_t, 3 > > triangles )
{
    geode::BRepBuilder builder{ model };
    const auto surface_id = builder.add_surface();
    auto mesh_builder =
        builder.surface_mesh_builder( model.surface( surface_id ) );
    for( const auto& point : points )
    {
        mesh_builder->create_point( point );
    }
    for( const auto& triangle : triangles )
    {
        mesh_builder->create_polygon( triangle );
    }
    return surface_id;
}

void test_geos_checker_mesh()
{
    {
        geode::BRep model;
        geode::BRepBuilder builder{ model };
        const auto& block_id = builder.add_block(
            geode::OpenGeodeRegularGrid3D::impl_name_static() );
        builder
            .block_mesh_builder< geode::RegularGrid3D >(
                model.block( block_id ) )
            ->initialize_grid(
                geode::Point3D{ { 0., 0., 0. } }, { 2, 2, 2 }, 1. );
        check_geos_error( toy_problem( false ), model,
            "The Block meshes merge into a RegularGrid3D" );
    }
    {
        auto model = load_toy_model();
        const auto problem = toy_problem( false );
        check_no_geos_error( problem, model, "are not Block" );
        const std::array< geode::Point3D, 3 > points{ geode::Point3D{
                                                          { 1e4, 1e4, 1e4 } },
            geode::Point3D{ { 1e4 + 1, 1e4, 1e4 } },
            geode::Point3D{ { 1e4, 1e4 + 1, 1e4 } } };
        const std::array< std::array< geode::index_t, 3 >, 1 > triangles{
            { { 0, 1, 2 } }
        };
        const auto surface_id =
            add_floating_surface( model, points, triangles );
        check_geos_error( problem, model,
            "polygon vertices which are not Block polyhedron vertices" );
        check_no_geos_error(
            problem, model, "unique vertices are not Block vertices" );
        geode::BRepBuilder builder{ model };
        const auto first_unique_vertex = builder.create_unique_vertices( 3 );
        const auto& surface = model.surface( surface_id );
        for( const auto vertex : geode::Range{ 3 } )
        {
            builder.set_unique_vertex( { surface.component_id(), vertex },
                first_unique_vertex + vertex );
        }
        check_geos_error(
            problem, model, "3 unique vertices are not Block vertices" );
        check_geos_error( problem, model,
            "polygon vertices which are not Block polyhedron vertices" );
    }
}

// Pressure and oil fraction on the given Surfaces, water fraction on the
// water Surfaces
void add_multiphase_boundary_conditions( geode::BRep& model,
    absl::Span< const geode::uuid > surfaces,
    absl::Span< const geode::uuid > water_surfaces )
{
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure",
        2.7579e7 );
    add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_oil_fraction, "bc_oil",
        0.9995 );
    add_boundary_condition< geode::ConstantAttribute >( model, water_surfaces,
        geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction, "bc_water",
        0.0005 );
}

std::shared_ptr< geode::Co2BrineFluid > co2_brine_fluid(
    geode::CO2_BRINE_CORRELATION correlation,
    std::vector< double > density_coefficients,
    std::vector< double > viscosity_coefficients,
    geode::ValidityRange pressures = { 1e6, 5e7, 1e6 },
    geode::ValidityRange temperatures = { 290., 400., 10. } )
{
    return std::make_shared< geode::Co2BrineFluid >(
        std::vector< geode::PHASE >{ geode::PHASE::gas, geode::PHASE::water },
        std::vector< geode::FluidComponent >{
            { "co2", 44e-3 }, { "water", 18e-3 } },
        correlation, 1., pressures, temperatures,
        std::move( density_coefficients ),
        std::move( viscosity_coefficients ) );
}

geode::PhysicalProblem co2_brine_problem(
    const geode::BRep& model, std::shared_ptr< geode::Co2BrineFluid > fluid )
{
    auto problem = spe10_reference_problem( model );
    auto flow = std::make_shared< geode::MultiphaseFlow >(
        geode::COMPOSITION_BASIS::mass );
    flow->set_fluid( std::move( fluid ) );
    geode::PhysicalProblemBuilder{ problem }.set_flow( flow );
    return problem;
}

// GEOS writes the materials covering no Block too
void test_geos_checker_unused_material()
{
    // No rock property in the model: materials give every value
    const auto model =
        geode::load_brep( absl::StrCat( geode::DATA_PATH, "picasso.og_brep" ) );
    auto problem = toy_problem( false );
    auto rock = problem.materials()[0];
    auto unused = rock;
    for( const auto& block : model.blocks() )
    {
        rock.regions.push_back( block.id() );
    }
    unused.name = "unused";
    unused.permeability.reset();
    geode::PhysicalProblemBuilder builder{ problem };
    builder.set_material( 0, rock );
    static_cast< void >( builder.add_material( unused ) );
    check_geos_warning( problem, model, "Material unused covers no Block" );
    const auto report =
        geode::check_physical_problem_simulation( problem, model, "geos" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        report.errors.size() == 1
            && has_issue(
                report.errors, "Material unused needs a permeability" ),
        "[Test] The unused material permeability should be the only error: ",
        report.string() );
}

// Values of properties GEOS does not write are not checked
void test_geos_checker_ignored_boundary_values()
{
    auto model = load_toy_model();
    add_boundary_condition< geode::ConstantAttribute >( model,
        all_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_temperature, "celsius", -5. );
    const auto isothermal = toy_problem( false );
    check_geos_warning( isothermal, model,
        "boundary_temperature is ignored by GEOS single phase" );
    check_no_geos_error( isothermal, model, "values outside" );
    check_geos_error( toy_problem( true ), model, "values outside ]0, +inf[" );
}

void test_geos_checker_problem_warnings()
{
    const auto model = spe10_reference_model();
    {
        // The SPE10 values are inside the dead oil PVT pressures
        auto problem = spe10_reference_problem( model );
        check_no_geos_warning(
            problem, model, "outside the fluid validity range" );
        geode::PhysicalProblemBuilder builder{ problem };
        builder.remove_initial_condition( 0 );
        builder.add_initial_condition(
            std::make_shared< geode::UniformInitialCondition >(
                std::vector< geode::uuid >{}, 1e8, 300.,
                std::vector< double >{ 0.9995, 0.0005 } ) );
        check_geos_warning( problem, model,
            "Initial condition #0 pressure is outside the fluid validity "
            "range" );
        auto high_pressure_model = geode::load_brep( absl::StrCat(
            geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) );
        const auto surfaces = all_surfaces( high_pressure_model );
        add_multiphase_boundary_conditions(
            high_pressure_model, surfaces, surfaces );
        check_no_geos_warning( spe10_reference_problem( high_pressure_model ),
            high_pressure_model, "Boundary pressures are outside" );
        add_boundary_condition< geode::ConstantAttribute >( high_pressure_model,
            surfaces, geode::PHYSICAL_PROPERTY_NAME::boundary_pressure,
            "bc_pressure_2", 1e8 );
        check_geos_warning( spe10_reference_problem( high_pressure_model ),
            high_pressure_model,
            "Boundary pressures are outside the fluid validity range" );
    }
    {
        const auto wide = co2_brine_problem( model,
            co2_brine_fluid( geode::CO2_BRINE_CORRELATION::phillips, {}, {} ) );
        check_no_geos_warning(
            wide, model, "outside the fluid validity range" );
        // The SPE10 pressures and temperatures are outside these ranges
        const auto narrow = co2_brine_problem(
            model, co2_brine_fluid( geode::CO2_BRINE_CORRELATION::phillips, {},
                       {}, { 1e6, 2e7, 1e6 }, { 310., 400., 10. } ) );
        for( const auto* keyword :
            { "Boundary pressures are outside the fluid validity range",
                "Boundary temperatures are outside the fluid validity range",
                "Initial condition #0 pressure is outside the fluid validity "
                "range",
                "Initial condition #0 temperature is outside the fluid "
                "validity range" } )
        {
            check_geos_warning( narrow, model, keyword );
        }
    }
    {
        check_geos_warning(
            co2_brine_problem(
                model, co2_brine_fluid( geode::CO2_BRINE_CORRELATION::phillips,
                           { 1., 2. }, {} ) ),
            model,
            "Co2BrineFluid Ezrokhi coefficients are ignored by the Phillips "
            "correlation" );
        check_geos_warning(
            co2_brine_problem(
                model, co2_brine_fluid( geode::CO2_BRINE_CORRELATION::ezrokhi,
                           { 1., 2. }, {} ) ),
            model,
            "Co2BrineFluid Ezrokhi correlation without density or viscosity "
            "coefficients" );
        check_no_geos_warning(
            co2_brine_problem(
                model, co2_brine_fluid( geode::CO2_BRINE_CORRELATION::ezrokhi,
                           { 1., 2. }, { 3., 4. } ) ),
            model, "Ezrokhi" );
    }
    {
        auto toy_model = load_toy_model();
        add_boundary_condition< geode::ConstantAttribute >( toy_model,
            east_surfaces( toy_model ),
            geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure",
            2e7 );
        auto problem = toy_problem( false );
        check_no_geos_warning( problem, toy_model, "Material rock" );
        auto rock = problem.materials()[0];
        rock.relative_permeability =
            spe10_reference_material().relative_permeability;
        rock.thermal = geode::RockThermalProperties{ 2e6, { 2.5, 2.5, 2.5 } };
        geode::PhysicalProblemBuilder{ problem }.set_material( 0, rock );
        check_geos_warning( problem, toy_model,
            "Material rock saturation functions are ignored by single phase "
            "GEOS simulations" );
        check_geos_warning( problem, toy_model,
            "Material rock thermal properties are ignored by isothermal GEOS "
            "simulations" );
        check_no_geos_warning(
            toy_problem( true ), toy_model, "thermal properties are ignored" );
        check_no_geos_warning(
            problem, toy_model, "nothing drives the flow in the GEOS deck" );
        check_geos_warning( problem, load_toy_model(),
            "The model has no boundary_pressure PhysicalProperty: nothing "
            "drives the flow in the GEOS deck" );
    }
}

void test_geos_checker_boundaries()
{
    {
        auto model = load_toy_model();
        const auto surfaces = east_surfaces( model );
        add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
            geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure",
            2e7 );
        add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
            geode::PHYSICAL_PROPERTY_NAME::boundary_oil_fraction, "bc_oil",
            1. );
        check_geos_warning( toy_problem( false ), model,
            "boundary_oil_fraction is ignored by GEOS single phase" );
        check_geos_error( toy_problem( true ), model,
            "need both a boundary pressure and a boundary temperature" );
        add_boundary_condition< geode::ConstantAttribute >( model,
            all_surfaces( model ),
            geode::PHYSICAL_PROPERTY_NAME::boundary_temperature,
            "bc_temperature", 330. );
        check_geos_error(
            toy_problem( true ), model, "but not boundary_pressure" );
        // Single phase isothermal decks ignore the boundary temperatures
        check_no_geos_error(
            toy_problem( false ), model, "these properties must be given" );
        check_geos_warning( toy_problem( false ), model,
            "boundary_temperature is ignored by GEOS single phase" );
    }
    {
        auto model = load_toy_model();
        const auto surfaces = all_surfaces( model );
        const std::vector< geode::uuid > other_surfaces{ surfaces.begin() + 1,
            surfaces.end() };
        add_multiphase_boundary_conditions( model, surfaces, other_surfaces );
        check_geos_error( spe10_reference_problem( model ), model,
            "carries boundary_pressure, boundary_oil_fraction but not "
            "boundary_water_fraction" );
    }
    {
        // Isothermal multiphase decks take the boundary temperatures given
        // on the model, so they must be on every pressure Surface
        auto model = load_toy_model();
        const auto surfaces = all_surfaces( model );
        add_multiphase_boundary_conditions( model, surfaces, surfaces );
        const auto problem = spe10_reference_problem( model );
        check_no_geos_error( problem, model, "these properties must be given" );
        add_boundary_condition< geode::ConstantAttribute >( model,
            std::vector< geode::uuid >{ surfaces.front() },
            geode::PHYSICAL_PROPERTY_NAME::boundary_temperature,
            "bc_temperature", 300. );
        check_geos_error( problem, model,
            "carries boundary_pressure, boundary_oil_fraction, "
            "boundary_water_fraction but not boundary_temperature" );
    }
    {
        auto model = load_toy_model();
        const auto problem = toy_problem( false );
        // Two triangles sharing the (1, 1, 0) center
        const std::array< geode::Point3D, 6 > points{
            geode::Point3D{ { 0., 0., 0. } }, geode::Point3D{ { 3., 0., 0. } },
            geode::Point3D{ { 0., 3., 0. } }, geode::Point3D{ { 2., 2., 0. } },
            geode::Point3D{ { 1., -1., 0. } }, geode::Point3D{ { 0., 2., 0. } }
        };
        const std::array< std::array< geode::index_t, 3 >, 2 > triangles{
            { { 0, 1, 2 }, { 3, 4, 5 } }
        };
        const auto surface_id =
            add_floating_surface( model, points, triangles );
        const auto attribute_id =
            add_boundary_condition< geode::VariableAttribute >( model,
                std::vector< geode::uuid >{ surface_id },
                geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure",
                2e7 );
        check_no_geos_error(
            problem, model, "GEOS tables cannot tell them apart" );
        model.surface( surface_id )
            .mesh()
            .polygon_attribute_manager()
            .find_attribute< geode::VariableAttribute, double >( attribute_id )
            ->set_value( 1, 3e7 );
        check_geos_error(
            problem, model, "GEOS tables cannot tell them apart" );
    }
}

void check_well_failure( std::string_view name,
    const geode::Point3D& point,
    std::string_view keyword )
{
    const auto model = spe10_reference_model();
    auto problem = spe10_reference_problem( model );
    geode::PhysicalProblemBuilder{ problem }.remove_source_term( 0 );
    geode::BRepGeosExporter exporter(
        model, problem, "grid_geos_well_failure" );
    auto point_set = geode::PointSet3D::create(
        geode::OpenGeodePointSet3D::impl_name_static() );
    geode::PointSetBuilder3D::create( *point_set )->create_point( point );
    exporter.add_well_perforations( *point_set, name );
    try
    {
        exporter.run();
    }
    catch( const geode::OpenGeodeException& exception )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::StrContains( exception.what(), keyword ),
            "[Test] Well error does not mention \"", keyword,
            "\": ", exception.what() );
        return;
    }
    throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
        geode::OpenGeodeException::TYPE::internal, "[Test] Export with ",
        keyword, " should fail" };
}

class WarningRecorder final : public geode::LoggerClient
{
public:
    explicit WarningRecorder( std::vector< std::string >& warnings )
        : warnings_( warnings )
    {
    }

    void trace( const std::string& /*unused*/ ) override {}

    void debug( const std::string& /*unused*/ ) override {}

    void info( const std::string& /*unused*/ ) override {}

    void warning( const std::string& message ) override
    {
        warnings_.push_back( message );
    }

    void error( const std::string& /*unused*/ ) override {}

    void critical( const std::string& /*unused*/ ) override {}

private:
    std::vector< std::string >& warnings_;
};

std::vector< std::string >& logged_warnings()
{
    static std::vector< std::string > warnings;
    static const auto& client = geode::LoggerManager::register_client(
        std::make_unique< WarningRecorder >( warnings ) );
    geode_unused( client );
    return warnings;
}

void test_geos_checker_wells()
{
    {
        const auto model = spe10_reference_model();
        auto problem = spe10_reference_problem( model );
        geode::PhysicalProblemBuilder{ problem }.remove_source_term( 0 );
        geode::BRepGeosExporter exporter(
            model, problem, "grid_geos_partial_well" );
        auto point_set = geode::PointSet3D::create(
            geode::OpenGeodePointSet3D::impl_name_static() );
        auto builder = geode::PointSetBuilder3D::create( *point_set );
        builder->create_point( geode::Point3D{ { 185.93, 336.8, 0.305 } } );
        builder->create_point( geode::Point3D{ { 1e5, 1e5, 1e5 } } );
        exporter.add_well_perforations( *point_set, "partial" );
        auto& warnings = logged_warnings();
        warnings.clear();
        exporter.run();
        geode::OpenGeodeGeosciencesIOModelException::test(
            has_issue( warnings, "Well partial: 1 of 2 perforations are "
                                 "outside the Blocks and ignored" ),
            "[Test] Missing warning on the perforation outside the Blocks" );
    }
    check_well_failure( "bad name", geode::Point3D{ { 185.93, 336.8, 0.305 } },
        "Well name \"bad name\"" );
    check_well_failure( "outside", geode::Point3D{ { 1e5, 1e5, 1e5 } },
        "has no perforation in a Block" );
}

void test_geos_deck_boundary_temperature()
{
    {
        auto model = geode::load_brep( absl::StrCat(
            geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) );
        const auto surfaces = all_surfaces( model );
        add_multiphase_boundary_conditions( model, surfaces, surfaces );
        auto problem = spe10_reference_problem( model );
        geode::PhysicalProblemBuilder{ problem }.remove_source_term( 0 );
        geode::BRepGeosExporter exporter( model, problem, "spe10_isothermal" );
        exporter.run();
        const auto deck =
            read_file( "spe10_isothermal/spe10_isothermal_simulation.xml" );
        const auto node_start =
            deck.find( R"(name="boundaryTemperature_isothermal")" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            node_start != std::string::npos,
            "[Test] Isothermal multiphase deck misses the boundary "
            "temperature of its initial conditions" );
        // Initial temperature on the faces of the boundary pressure
        const auto node = deck.substr(
            node_start, deck.find( "/>", node_start ) - node_start );
        for( const auto* pattern : { R"(fieldName="temperature")",
                 R"(scale="300")", R"(setNames="{ 1 }")" } )
        {
            geode::OpenGeodeGeosciencesIOModelException::test(
                absl::StrContains( node, pattern ),
                "[Test] Isothermal boundary temperature misses ", pattern );
        }
    }
    {
        auto model = load_toy_model();
        const auto surfaces = east_surfaces( model );
        add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
            geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, "bc_pressure",
            2e7 );
        add_boundary_condition< geode::ConstantAttribute >( model, surfaces,
            geode::PHYSICAL_PROPERTY_NAME::boundary_temperature,
            "bc_temperature", 330. );
        const auto problem = toy_problem( false );
        geode::BRepGeosExporter exporter(
            model, problem, "toy_model_isothermal" );
        exporter.run();
        const auto deck = read_file(
            "toy_model_isothermal/toy_model_isothermal_simulation.xml" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::StrContains( deck, "name=\"boundaryPressure_0\"" )
                && !absl::StrContains( deck, "boundaryTemperature" )
                && count_occurrences( deck, "objectPath=\"faceManager\"" ) == 1,
            "[Test] Single phase isothermal deck should not write the "
            "boundary temperature" );
    }
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        geode::OpenGeodeIOMeshLibrary::initialize();
        geode::OpenGeodeIOModelLibrary::initialize();
        test_picasso();
        test_uniform_variable_boundary_condition();
        test_variable_boundary_condition();
        test_spe10_reference();
        test_toy_model_thermal();
        test_geos_rules();
        test_structural_model_regions();
        test_geos_checker_problem();
        test_geos_checker_model();
        test_geos_checker_mesh();
        test_geos_checker_unused_material();
        test_geos_checker_ignored_boundary_values();
        test_geos_checker_problem_warnings();
        test_geos_checker_boundaries();
        test_geos_checker_wells();
        test_geos_deck_boundary_temperature();
        geode::Logger::info( "TEST SUCCESS" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
