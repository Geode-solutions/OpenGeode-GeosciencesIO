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
#include <geode/mesh/builder/tetrahedral_solid_builder.hpp>
#include <geode/mesh/core/geode/geode_tetrahedral_solid.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/io/tetrahedral_solid_input.hpp>

#include <geode/mesh/core/geode/geode_point_set.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/builder/brep_builder.hpp>
#include <geode/model/representation/core/brep.hpp>
#include <geode/model/representation/io/brep_input.hpp>

std::string read_file( const std::string& filename );

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
    const auto deck = read_file( "picasso/picasso_simulation.xml" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        absl::StrContains( deck, "component=\"0\"" )
            && absl::StrContains( deck, "component=\"2\"" )
            && !absl::StrContains( deck, "component=\"1\"" ),
        "[Test] Boundary fractions must use the fixed oil=0, water=2 "
        "component indices" );
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
    geode::BRepGeosExporter exporter( model, "spe10_reference" );
    exporter.set_physical_problem( problem );
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

    geode::BRepGeosExporter exporter( model, "toy_model_thermal" );
    exporter.set_physical_problem( problem );
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
    geode::BRepGeosExporter exporter( model, "grid_geos_failure" );
    exporter.set_physical_problem( problem );
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
        // Without its source term, which needs wells given to the exporter
        auto problem = spe10_reference_problem( model );
        geode::PhysicalProblemBuilder{ problem }.remove_source_term( 0 );
        check_export_failure( model, problem, "the fluid has no gas phase" );
    }
}

void test_structural_model_regions()
{
    geode::StructuralModel model{ geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) ) };
    add_spe10_boundary_conditions( model );
    add_spe10_gas_fraction( model );
    geode::StructuralModelBuilder model_builder{ model };
    const auto& unit_id = model_builder.add_stratigraphic_unit();
    const auto& unit = model.stratigraphic_unit( unit_id );
    std::optional< geode::uuid > unit_block;
    for( const auto& block : model.blocks() )
    {
        model_builder.add_block_in_stratigraphic_unit( block, unit );
        unit_block = block.id();
        break;
    }
    auto problem = spe10_problem( model );
    auto unit_material = spe10_material( "unit", 1e-10 );
    unit_material.regions = { unit_block.value() };
    const auto unit_material_id =
        geode::PhysicalProblemBuilder{ problem }.add_material( unit_material );
    bool has_thrown{ false };
    try
    {
        geode::StructuralModelGeosExporter exporter(
            model, "structural_model_failure" );
        exporter.set_physical_problem( problem );
        exporter.run();
    }
    catch( const geode::OpenGeodeException& exception )
    {
        has_thrown = true;
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::StrContains( exception.what(), "is part of the GEOS region" ),
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
        model, "structural_model_regions" );
    exporter.set_physical_problem( problem );
    exporter.run();
    const auto deck = read_file(
        "structural_model_regions/structural_model_regions_simulation.xml" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        absl::StrContains( deck, "rock_unit" )
            && absl::StrContains( deck, "rock_reservoir" ),
        "[Test] StructuralModel deck should use both materials" );
}

geode::BRep tetrahedra_brep(
    geode::index_t nb_blocks, const geode::Point3D& origin )
{
    geode::BRep model;
    geode::BRepBuilder builder{ model };
    builder.create_unique_vertices( 4 * nb_blocks );
    for( const auto block_id : geode::Range{ nb_blocks } )
    {
        const auto& block = model.block( builder.add_block(
            geode::OpenGeodeTetrahedralSolid3D::impl_name_static() ) );
        auto mesh_builder =
            builder.block_mesh_builder< geode::TetrahedralSolid3D >( block );
        auto x = origin.value( 0 ) + 2. * block_id;
        mesh_builder->create_point(
            geode::Point3D{ { x, origin.value( 1 ), origin.value( 2 ) } } );
        mesh_builder->create_point( geode::Point3D{
            { x + 1., origin.value( 1 ), origin.value( 2 ) } } );
        mesh_builder->create_point( geode::Point3D{
            { x, origin.value( 1 ) + 1., origin.value( 2 ) } } );
        mesh_builder->create_point( geode::Point3D{
            { x, origin.value( 1 ), origin.value( 2 ) + 1. } } );
        mesh_builder->create_tetrahedron( { 0, 1, 2, 3 } );
        for( const auto vertex : geode::LRange{ 4 } )
        {
            builder.set_unique_vertex(
                { block.component_id(), vertex }, 4 * block_id + vertex );
        }
    }
    return model;
}

template < template < typename > class Attribute, typename T >
geode::uuid add_block_property( geode::BRep& model,
    geode::PHYSICAL_PROPERTY_NAME property,
    std::string_view name,
    absl::Span< const T > block_values )
{
    const geode::uuid attribute_id;
    geode::index_t block_id{ 0 };
    for( const auto& block : model.blocks() )
    {
        geode::AttributeValues< T > values;
        values.default_value = block_values[block_id++];
        block.mesh()
            .polyhedron_attribute_manager()
            .create_attribute< Attribute, T >( name, attribute_id, values, {} );
    }
    geode::BRepBuilder{ model }.set_physical_property(
        property, geode::Block3D::component_type_static(), attribute_id );
    return attribute_id;
}

template < typename T >
std::vector< T > vtu_cell_values(
    std::string_view filename, std::string_view attribute_name )
{
    const auto solid = geode::load_tetrahedral_solid< 3 >( filename );
    const auto& manager = solid->polyhedron_attribute_manager();
    const auto ids = manager.attribute_ids_matching_name( attribute_name );
    geode::OpenGeodeGeosciencesIOModelException::test(
        ids.has_value() && ids->size() == 1, "[Test] No attribute ",
        attribute_name, " in ", filename );
    const auto attribute =
        manager.find_read_only_attribute< T >( ids->front() );
    std::vector< T > values;
    for( const auto polyhedron : geode::Range{ solid->nb_polyhedra() } )
    {
        values.push_back( attribute->value( polyhedron ) );
    }
    return values;
}

void test_constant_block_properties()
{
    auto model = tetrahedra_brep( 3, geode::Point3D{ { 0., 0., 0. } } );
    const std::array< double, 3 > porosities{ 0.1, 0.2, 0.3 };
    add_block_property< geode::ConstantAttribute, double >(
        model, geode::PHYSICAL_PROPERTY_NAME::porosity, "poro", porosities );
    geode::BRepGeosExporter exporter( model, "constant_blocks" );
    exporter.run();
    const auto values = vtu_cell_values< double >(
        "constant_blocks/constant_blocks.vtu", "poro" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        values.size() == porosities.size(),
        "[Test] Wrong number of exported cells" );
    for( const auto block : geode::Indices{ porosities } )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            values[block] == porosities[block],
            "[Test] Wrong porosity exported for Block ", block, ": ",
            values[block] );
    }
}

void test_structural_model_block_without_unit()
{
    geode::StructuralModel model{ tetrahedra_brep(
        3, geode::Point3D{ { 0., 0., 0. } } ) };
    geode::StructuralModelBuilder builder{ model };
    const auto& unit =
        model.stratigraphic_unit( builder.add_stratigraphic_unit() );
    std::vector< geode::uuid > block_ids;
    for( const auto& block : model.blocks() )
    {
        block_ids.push_back( block.id() );
    }
    builder.add_block_in_stratigraphic_unit(
        model.block( block_ids[0] ), unit );
    builder.add_block_in_stratigraphic_unit(
        model.block( block_ids[2] ), unit );
    geode::StructuralModelGeosExporter exporter( model, "unit_regions" );
    exporter.run();
    const auto regions = vtu_cell_values< geode::index_t >(
        "unit_regions/unit_regions.vtu", "attribute" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        regions == std::vector< geode::index_t >{ 0, 1, 0 },
        "[Test] Wrong regions for a StructuralModel Block in no "
        "StratigraphicUnit" );
}

void test_structural_model_block_in_two_units()
{
    geode::StructuralModel model{ tetrahedra_brep(
        2, geode::Point3D{ { 0., 0., 0. } } ) };
    geode::StructuralModelBuilder builder{ model };
    const auto& first_unit =
        model.stratigraphic_unit( builder.add_stratigraphic_unit() );
    const auto& second_unit =
        model.stratigraphic_unit( builder.add_stratigraphic_unit() );
    std::vector< geode::uuid > block_ids;
    for( const auto& block : model.blocks() )
    {
        block_ids.push_back( block.id() );
    }
    builder.add_block_in_stratigraphic_unit(
        model.block( block_ids[0] ), first_unit );
    builder.add_block_in_stratigraphic_unit(
        model.block( block_ids[0] ), second_unit );
    builder.add_block_in_stratigraphic_unit(
        model.block( block_ids[1] ), second_unit );
    std::vector< geode::uuid > unit_ids;
    for( const auto& unit : model.stratigraphic_units() )
    {
        unit_ids.push_back( unit.id() );
    }
    geode::StructuralModelGeosExporter exporter( model, "two_units" );
    exporter.run();
    const auto regions = vtu_cell_values< geode::index_t >(
        "two_units/two_units.vtu", "attribute" );
    const auto first_region = unit_ids.front() == first_unit.id() ? 0u : 1u;
    geode::OpenGeodeGeosciencesIOModelException::test(
        regions.size() == 2 && regions[1] == ( first_region == 0 ? 1u : 0u )
            && ( regions[0] == 0 || regions[0] == 1 ),
        "[Test] Wrong regions for a Block in two StratigraphicUnits" );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        geode::OpenGeodeIOMeshLibrary::initialize();
        geode::OpenGeodeIOModelLibrary::initialize();
        test_constant_block_properties();
        test_structural_model_block_without_unit();
        test_structural_model_block_in_two_units();
        test_picasso();
        test_uniform_variable_boundary_condition();
        test_variable_boundary_condition();
        test_spe10_reference();
        test_toy_model_thermal();
        test_geos_rules();
        test_structural_model_regions();
        geode::Logger::info( "TEST SUCCESS" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
