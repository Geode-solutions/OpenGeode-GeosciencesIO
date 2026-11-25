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
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <pugixml.hpp>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_set.h>
#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/types.hpp>
#include <geode/basic/variable_attribute.hpp>
#include <geode/tests_config.hpp>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/builder/hybrid_solid_builder.hpp>
#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/builder/polyhedral_solid_builder.hpp>
#include <geode/mesh/builder/regular_grid_solid_builder.hpp>
#include <geode/mesh/core/geode/geode_hybrid_solid.hpp>
#include <geode/mesh/core/geode/geode_point_set.hpp>
#include <geode/mesh/core/geode/geode_polyhedral_solid.hpp>
#include <geode/mesh/core/geode/geode_regular_grid_solid.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/helpers/component_mesh_polygons.hpp>
#include <geode/model/helpers/convert_to_mesh.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/builder/brep_builder.hpp>
#include <geode/model/representation/io/brep_input.hpp>

#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.hpp>
#include <geode/geosciences/explicit/representation/builder/structural_model_builder.hpp>
#include <geode/geosciences/explicit/representation/core/structural_model.hpp>

#include <geode/geosciences/simulation/builder/physical_problem_builder.hpp>
#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/schedule.hpp>
#include <geode/geosciences/simulation/core/solid_mechanics.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>
#include <geode/geosciences/simulation/core/thermal_parameters.hpp>
#include <geode/geosciences/simulation/helpers/simulation_checker.hpp>

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/helpers/brep_ogs_export.hpp>
#include <geode/geosciences_io/model/helpers/structural_model_ogs_export.hpp>
#include <geode/geosciences_io/model/internal/ogs_simulation_rules.hpp>
#include <geode/geosciences_io/model/internal/ogs_vtu.hpp>
#include <geode/geosciences_io/model/internal/simulation_regions.hpp>

namespace
{
    pugi::xml_document load_xml( const std::string& filename )
    {
        pugi::xml_document document;
        const auto result = document.load_file( filename.c_str() );
        geode::OpenGeodeGeosciencesIOModelException::test(
            static_cast< bool >( result ), "[Test] Cannot parse ", filename );
        return document;
    }

    pugi::xml_node data_array( const pugi::xml_document& document,
        std::string_view section,
        std::string_view name )
    {
        const auto xpath =
            absl::StrCat( "//", section, "/DataArray[@Name='", name, "']" );
        return document.select_node( xpath.c_str() ).node();
    }

    std::string text( const pugi::xml_node& node )
    {
        return node.text().as_string();
    }

    void check( bool condition, std::string_view message )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            condition, "[Test] ", message );
    }
} // namespace

void test_ogs_mesh_file()
{
    geode::internal::OgsMeshFile mesh;
    mesh.points = { geode::Point3D{ { 0, 0, 0 } },
        geode::Point3D{ { 1, 0, 0 } }, geode::Point3D{ { 1, 1, 0 } },
        geode::Point3D{ { 0, 1, 0 } } };
    mesh.cells = { { 0, 1, 2 }, { 0, 2, 3 } };
    mesh.bulk_node_ids = { 7, 3, 9, 4 };
    mesh.point_arrays.emplace_back(
        "values", std::vector< double >{ 1.5, 2.5, 3.5, 4.5 } );
    geode::internal::save_ogs_mesh_file( "ogs_mesh_file.vtu", mesh );
    const auto document = load_xml( "ogs_mesh_file.vtu" );
    const auto piece = document.select_node( "//Piece" ).node();
    check( piece.attribute( "NumberOfPoints" ).as_uint() == 4,
        "Wrong number of points" );
    check( piece.attribute( "NumberOfCells" ).as_uint() == 2,
        "Wrong number of cells" );
    const auto ids = data_array( document, "PointData", "bulk_node_ids" );
    check( std::string{ ids.attribute( "type" ).value() } == "UInt64",
        "bulk_node_ids must be UInt64" );
    check( text( ids ) == "7 3 9 4", "Wrong bulk_node_ids" );
    const auto values = data_array( document, "PointData", "values" );
    check( std::string{ values.attribute( "type" ).value() } == "Float64",
        "Point arrays must be Float64" );
    check( text( data_array( document, "Cells", "types" ) ) == "5 5",
        "Triangles must be VTK type 5" );
    check( !data_array( document, "CellData", "bulk_element_ids" ),
        "Empty cell id arrays must not be written" );

    geode::internal::OgsMeshFile vertices;
    vertices.points = { geode::Point3D{ { 0.5, 0.5, 0.5 } } };
    vertices.cells = { { 0 } };
    vertices.bulk_node_ids = { 12 };
    geode::internal::save_ogs_mesh_file( "ogs_vertex_file.vtu", vertices );
    const auto vertex_document = load_xml( "ogs_vertex_file.vtu" );
    check( text( data_array( vertex_document, "Cells", "types" ) ) == "1",
        "Single-vertex cells must be VTK type 1" );
}

void test_ogs_structural_model_regions()
{
    const auto model_filename =
        absl::StrCat( geode::DATA_PATH, "picasso.og_brep" );
    geode::StructuralModel model{ geode::load_one(
        geode::load_brep( model_filename ), model_filename ) };
    check( model.nb_blocks() >= 3, "picasso should have 3 Blocks or more" );
    geode::StructuralModelBuilder builder{ model };
    const auto& first_id = builder.add_stratigraphic_unit();
    const auto& second_id = builder.add_stratigraphic_unit();
    std::vector< geode::uuid > first_blocks;
    std::vector< geode::uuid > second_blocks;
    for( const auto& block : model.blocks() )
    {
        if( first_blocks.size() < 2 )
        {
            builder.add_block_in_stratigraphic_unit(
                block, model.stratigraphic_unit( first_id ) );
            first_blocks.push_back( block.id() );
        }
        else if( second_blocks.empty() )
        {
            builder.add_block_in_stratigraphic_unit(
                block, model.stratigraphic_unit( second_id ) );
            second_blocks.push_back( block.id() );
        }
    }
    const auto regions =
        geode::internal::structural_model_block_regions( model );
    check( regions.solid_regions.size() == model.nb_blocks(),
        "One region per Block expected" );
    check( regions.merged_blocks.empty(), "No Block should be merged" );
    const auto check_unit = [&regions]( const geode::uuid& unit,
                                absl::Span< const geode::uuid > blocks ) {
        std::vector< geode::index_t > expected;
        for( const auto& block : blocks )
        {
            expected.push_back( regions.block_regions.at( block ) );
        }
        check( regions.regions_by_uuid.at( unit ) == expected,
            "A unit should map to the regions of its Blocks" );
    };
    check_unit( first_id, first_blocks );
    check_unit( second_id, second_blocks );
}

std::shared_ptr< geode::SinglePhaseFluid > water_fluid( double compressibility )
{
    return std::make_shared< geode::SinglePhaseFluid >(
        std::vector< geode::PHASE >{ geode::PHASE::water },
        std::vector< geode::FluidComponent >{ { "water", 0.018 } }, 1000., 1e-3,
        compressibility, 0., 1e5 );
}

std::shared_ptr< geode::DeadOilFluid > water_oil_fluid(
    std::vector< geode::PHASE > phases )
{
    std::vector< geode::FluidComponent > components;
    std::vector< double > densities;
    std::vector< geode::HydrocarbonPvtTable > tables;
    std::optional< geode::WaterPvt > water;
    for( const auto phase : phases )
    {
        components.push_back(
            { std::string{ geode::phase_name( phase ) }, 0.1 } );
        densities.push_back( 1000. );
        if( phase == geode::PHASE::water )
        {
            water = geode::WaterPvt{ 1e5, 1., 0., 1e-3 };
            continue;
        }
        geode::HydrocarbonPvtTable table;
        table.phase = phase;
        table.pressures = { 1e5, 1e7 };
        table.formation_volume_factors = { 1., 1. };
        table.viscosities = { 1e-3, 1e-3 };
        tables.push_back( std::move( table ) );
    }
    return std::make_shared< geode::DeadOilFluid >( std::move( phases ),
        std::move( components ), std::move( densities ), std::move( tables ),
        water );
}

geode::PhysicalProblem flow_problem( std::shared_ptr< geode::FlowModel > flow,
    std::shared_ptr< geode::FluidModel > fluid )
{
    geode::PhysicalProblem problem;
    geode::PhysicalProblemBuilder builder{ problem };
    flow->set_fluid( std::move( fluid ) );
    builder.set_flow( std::move( flow ) );
    return problem;
}

void test_ogs_rules_process()
{
    using geode::internal::OGS_PROCESS;
    auto single = flow_problem(
        std::make_shared< geode::SinglePhaseFlow >(), water_fluid( 0. ) );
    check( geode::internal::ogs_process( single ) == OGS_PROCESS::liquid_flow,
        "Single-phase flow is LIQUID_FLOW" );
    geode::PhysicalProblemBuilder{ single }.set_thermal(
        geode::ThermalParameters{ 300. } );
    check( !geode::internal::ogs_process( single ),
        "Thermal single-phase flow is not supported yet" );
    const auto two_phase =
        flow_problem( std::make_shared< geode::MultiphaseFlow >(),
            water_oil_fluid( { geode::PHASE::water, geode::PHASE::oil } ) );
    check( geode::internal::ogs_process( two_phase )
               == OGS_PROCESS::two_phase_flow_pp,
        "Water-oil dead oil is TWOPHASE_FLOW_PP" );
    const auto three_phase = flow_problem(
        std::make_shared< geode::MultiphaseFlow >(),
        water_oil_fluid(
            { geode::PHASE::water, geode::PHASE::oil, geode::PHASE::gas } ) );
    check( !geode::internal::ogs_process( three_phase ),
        "Three-phase dead oil is not supported" );
}

void test_ogs_rules_storage_and_temperature()
{
    check(
        std::fabs( geode::internal::ogs_storage( 0.2, 1e-9 ) - 2e-10 ) < 1e-25,
        "Storage is porosity times pore compressibility" );
    auto problem = flow_problem(
        std::make_shared< geode::SinglePhaseFlow >(), water_fluid( 0. ) );
    geode::PhysicalProblemBuilder builder{ problem };
    builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{}, 1e5 ) );
    check( geode::internal::ogs_reference_temperature( problem ) == 293.15,
        "Default reference temperature is 293.15 K" );
    auto with_temperature = flow_problem(
        std::make_shared< geode::SinglePhaseFlow >(), water_fluid( 0. ) );
    geode::PhysicalProblemBuilder temperature_builder{ with_temperature };
    temperature_builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{ geode::uuid{} }, 1e5, 350. ) );
    temperature_builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{ geode::uuid{} }, 2e5, 350. ) );
    check(
        geode::internal::ogs_reference_temperature( with_temperature ) == 350.,
        "Shared initial temperature is the reference temperature" );
}

void test_ogs_rules_output_times()
{
    const geode::Schedule schedule{ 0., 100., 1., std::nullopt };
    geode::OutputSettings frequency;
    frequency.frequency = 25.;
    check( geode::internal::ogs_output_times( schedule, frequency )
               == std::vector< double >{ 25., 50., 75., 100. },
        "Frequency output times" );
    geode::OutputSettings times;
    times.times = { 10., 30. };
    check( geode::internal::ogs_output_times( schedule, times )
               == std::vector< double >{ 10., 30. },
        "Explicit output times" );
}

double interpolate( const geode::internal::OgsCurve& curve, double coord )
{
    for( const auto index : geode::Range{ 1, curve.coords.size() } )
    {
        if( coord <= curve.coords[index] )
        {
            const auto ratio =
                ( coord - curve.coords[index - 1] )
                / ( curve.coords[index] - curve.coords[index - 1] );
            return curve.values[index - 1]
                   + ratio * ( curve.values[index] - curve.values[index - 1] );
        }
    }
    return curve.values.back();
}

void test_ogs_rules_saturation()
{
    const auto fluid =
        water_oil_fluid( { geode::PHASE::water, geode::PHASE::oil } );
    const std::array< double, 2 > composition{ 0.3, 0.7 };
    check( std::fabs( geode::internal::water_saturation( *fluid,
                          geode::COMPOSITION_BASIS::mass, composition, 1e6 )
                      - 0.3 )
               < 1e-12,
        "Equal densities: Sw is the water mass fraction" );

    const geode::BrooksCoreyCapillaryPressure brooks_corey{
        { geode::PHASE::water, geode::PHASE::oil }, { 0.1, 0.2 }, { 5000., 0. },
        { 2., 2. }
    };
    const auto curve = geode::internal::saturation_curve( brooks_corey );
    check(
        curve.coords.size() == curve.values.size(), "Saturation curve sizes" );
    check( curve.coords.front() == 0. && curve.values.front() == 0.8,
        "Saturation curve starts at (0, 1 - Snr)" );
    for( const auto index : geode::Range{ 1, curve.coords.size() } )
    {
        check( curve.coords[index] > curve.coords[index - 1],
            "Capillary pressures must strictly increase" );
        check( curve.values[index] <= curve.values[index - 1],
            "Saturations must not increase" );
    }
    const auto expected = 0.1 + ( 1. - 0.1 - 0.2 ) * std::pow( 4., -2. );
    check( std::fabs( interpolate( curve, 4. * 5000. ) - expected ) < 1e-3,
        "Brooks-Corey saturation at 4 entry pressures" );

    geode::TwoPhaseCapillaryPressureTable table;
    table.pressure.saturations = { 0.2, 0.5, 0.9 };
    table.pressure.values = { 1e4, 2e4, 0. };
    const geode::TableCapillaryPressure non_monotone{
        { geode::PHASE::water, geode::PHASE::oil }, { table }
    };
    bool has_thrown{ false };
    try
    {
        geode_unused( geode::internal::saturation_curve( non_monotone ) );
    }
    catch( const geode::OpenGeodeException& )
    {
        has_thrown = true;
    }
    check( has_thrown, "A non-monotone capillary table must throw" );
}

void test_ogs_rules_relperm()
{
    const geode::BrooksCoreyRelativePermeability brooks_corey{
        { geode::PHASE::water, geode::PHASE::oil }, { 0.1, 0.2 }, { 2., 3. },
        { 0.8, 0.9 }
    };
    const auto curves =
        geode::internal::relative_permeability_curves( brooks_corey );
    const auto normalized_water = ( 0.5 - 0.1 ) / 0.7;
    const auto normalized_oil = ( 0.5 - 0.2 ) / 0.7;
    check( std::fabs( interpolate( curves[0], 0.5 )
                      - 0.8 * std::pow( normalized_water, 2. ) )
               < 1e-6,
        "Brooks-Corey wetting relative permeability" );
    check( std::fabs( interpolate( curves[1], 0.5 )
                      - 0.9 * std::pow( normalized_oil, 3. ) )
               < 1e-6,
        "Brooks-Corey non-wetting relative permeability" );

    geode::TwoPhaseRelativePermeabilityTable table;
    table.wetting_kr.saturations = { 0.1, 0.5, 1. };
    table.wetting_kr.values = { 0., 0.2, 1. };
    table.non_wetting_kr.saturations = { 0., 0.4, 0.9 };
    table.non_wetting_kr.values = { 0., 0.3, 1. };
    const geode::TableRelativePermeability tabulated{
        { geode::PHASE::water, geode::PHASE::oil }, { table }
    };
    const auto table_curves =
        geode::internal::relative_permeability_curves( tabulated );
    check( table_curves[0].coords == table.wetting_kr.saturations
               && table_curves[0].values == table.wetting_kr.values,
        "Wetting table written as given" );
    check( table_curves[1].values == std::vector< double >{ 1., 0.3, 0. },
        "Non-wetting table values versus the water saturation" );
    check( std::fabs( table_curves[1].coords[0] - 0.1 ) < 1e-12
               && std::fabs( table_curves[1].coords[1] - 0.6 ) < 1e-12
               && table_curves[1].coords[2] == 1.,
        "Non-wetting table coordinates are 1 - So" );
}

geode::BRep load_picasso()
{
    const auto model_filename =
        absl::StrCat( geode::DATA_PATH, "picasso.og_brep" );
    return geode::load_one(
        geode::load_brep( model_filename ), model_filename );
}

// Surfaces lying on the maximum x face of the model box
std::vector< geode::uuid > east_surfaces( const geode::BRep& model )
{
    static constexpr double TOLERANCE{ 1e-3 };
    const auto side = model.bounding_box().max().value( 0 );
    std::vector< geode::uuid > surfaces;
    for( const auto& surface : model.surfaces() )
    {
        const auto box = surface.mesh().bounding_box();
        if( std::fabs( box.min().value( 0 ) - side ) < TOLERANCE
            && std::fabs( box.max().value( 0 ) - side ) < TOLERANCE )
        {
            surfaces.push_back( surface.id() );
        }
    }
    return surfaces;
}

void add_constant_boundary( geode::BRep& model,
    absl::Span< const geode::uuid > surfaces,
    geode::PHYSICAL_PROPERTY_NAME property,
    double value )
{
    const geode::uuid attribute_id;
    geode::AttributeValues< double > values;
    values.default_value = value;
    const auto name = geode::physical_property_string( property );
    for( const auto& surface_id : surfaces )
    {
        model.surface( surface_id )
            .mesh()
            .polygon_attribute_manager()
            .template create_attribute< geode::ConstantAttribute, double >(
                name, attribute_id, values, {} );
    }
    geode::BRepBuilder{ model }.set_physical_property(
        property, geode::Surface3D::component_type_static(), attribute_id );
}

geode::Material rock_material( std::string name )
{
    geode::Material rock;
    rock.name = std::move( name );
    rock.pore_compressibility = 1e-9;
    rock.reference_pressure = 1e7;
    rock.porosity = 0.2;
    rock.permeability = std::array< double, 3 >{ 1e-13, 1e-13, 1e-14 };
    return rock;
}

geode::PhysicalProblem single_phase_problem( double compressibility )
{
    auto problem = flow_problem( std::make_shared< geode::SinglePhaseFlow >(),
        water_fluid( compressibility ) );
    geode::PhysicalProblemBuilder builder{ problem };
    builder.add_material( rock_material( "rock" ) );
    builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{}, 1e7 ) );
    builder.set_schedule( { 0., 1e5, 1e3, 1e4 } );
    geode::OutputSettings output;
    output.frequency = 5e4;
    builder.set_output( output );
    return problem;
}

std::shared_ptr< geode::BrooksCoreyCapillaryPressure > brooks_corey_pc()
{
    return std::make_shared< geode::BrooksCoreyCapillaryPressure >(
        std::vector< geode::PHASE >{ geode::PHASE::water, geode::PHASE::oil },
        std::vector< double >{ 0.1, 0.2 }, std::vector< double >{ 5000., 0. },
        std::vector< double >{ 2., 2. } );
}

geode::PhysicalProblem two_phase_problem()
{
    auto problem = flow_problem( std::make_shared< geode::MultiphaseFlow >(
                                     geode::COMPOSITION_BASIS::mass ),
        water_oil_fluid( { geode::PHASE::water, geode::PHASE::oil } ) );
    geode::PhysicalProblemBuilder builder{ problem };
    auto rock = rock_material( "rock" );
    rock.permeability = std::array< double, 3 >{ 1e-13, 1e-13, 1e-13 };
    rock.relative_permeability =
        std::make_shared< geode::BrooksCoreyRelativePermeability >(
            std::vector< geode::PHASE >{
                geode::PHASE::water, geode::PHASE::oil },
            std::vector< double >{ 0.1, 0.2 }, std::vector< double >{ 2., 2. },
            std::vector< double >{ 1., 1. } );
    rock.capillary_pressure = brooks_corey_pc();
    builder.add_material( rock );
    builder.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{}, 1e7, std::nullopt,
            std::vector< double >{ 0.3, 0.7 } ) );
    builder.set_schedule( { 0., 1e3, 10., 100. } );
    return problem;
}

geode::SimulationReport ogs_report(
    const geode::PhysicalProblem& problem, const geode::BRep& model )
{
    return geode::SimulationCheckerFactory::create( "ogs", problem, model )
        ->check();
}

void check_refused( const geode::PhysicalProblem& problem,
    const geode::BRep& model,
    std::string_view expected )
{
    const auto report = ogs_report( problem, model );
    check( !report.is_simulable(),
        absl::StrCat( "The problem should be refused: ", expected ) );
    check( absl::StrContains( report.string(), expected ),
        absl::StrCat(
            "Missing issue \"", expected, "\" in\n", report.string() ) );
}

void test_ogs_registered_checker()
{
    check( geode::SimulationCheckerFactory::has_creator( "ogs" ),
        "The ogs checker should be registered" );
    auto model = load_picasso();
    add_constant_boundary( model, east_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, 1e7 );
    const auto report = ogs_report( single_phase_problem( 4.5e-10 ), model );
    check( report.is_simulable(),
        absl::StrCat(
            "Picasso single-phase should be simulable:\n", report.string() ) );
    check( report.simulator() == "OGS", "Report simulator should be OGS" );
}

void test_ogs_checker_refusals()
{
    auto model = load_picasso();
    add_constant_boundary( model, east_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, 1e7 );
    {
        auto problem = single_phase_problem( 0. );
        geode::PhysicalProblemBuilder{ problem }.set_thermal(
            geode::ThermalParameters{ 300. } );
        check_refused(
            problem, model, "not supported by the OGS exporter yet" );
    }
    {
        auto problem = single_phase_problem( 0. );
        geode::PhysicalProblemBuilder{ problem }.set_solid_mechanics(
            geode::SolidMechanics{} );
        check_refused(
            problem, model, "not supported by the OGS exporter yet" );
    }
    {
        auto problem = two_phase_problem();
        auto& fluid = const_cast< geode::FlowModel& >( problem.flow() );
        fluid.set_fluid( water_oil_fluid(
            { geode::PHASE::water, geode::PHASE::oil, geode::PHASE::gas } ) );
        check_refused( problem, model, "gas phase" );
    }
    {
        auto problem = two_phase_problem();
        auto material = problem.materials()[0];
        material.capillary_pressure.reset();
        geode::PhysicalProblemBuilder{ problem }.set_material( 0, material );
        check_refused( problem, model, "needs a capillary pressure" );
    }
    {
        auto problem = two_phase_problem();
        geode::PhysicalProblemBuilder{ problem }.add_source_term(
            { "well", "water", 1. } );
        check_refused( problem, model, "source terms" );
    }
    {
        auto problem = two_phase_problem();
        auto material = problem.materials()[0];
        geode::TwoPhaseCapillaryPressureTable table;
        table.pressure.saturations = { 0.2, 0.5, 0.9 };
        table.pressure.values = { 1e4, 2e4, 0. };
        material.capillary_pressure =
            std::make_shared< geode::TableCapillaryPressure >(
                std::vector< geode::PHASE >{
                    geode::PHASE::water, geode::PHASE::oil },
                std::vector< geode::TwoPhaseCapillaryPressureTable >{ table } );
        geode::PhysicalProblemBuilder{ problem }.set_material( 0, material );
        check_refused( problem, model, "must strictly decrease" );
    }
    {
        auto problem = single_phase_problem( 0. );
        geode::PhysicalProblemBuilder{ problem }.set_solver_option(
            "ogs", "unknown/solver_type", "BiCGSTAB" );
        check_refused( problem, model, "unknown" );
    }
    auto gas_model = load_picasso();
    add_constant_boundary( gas_model, east_surfaces( gas_model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, 1e7 );
    add_constant_boundary( gas_model, east_surfaces( gas_model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_gas_fraction, 0.1 );
    check_refused( two_phase_problem(), gas_model, "boundary_gas_fraction" );
}

void test_ogs_unsupported_cells()
{
    geode::BRep model;
    geode::BRepBuilder builder{ model };
    const auto& block = model.block( builder.add_block(
        geode::OpenGeodePolyhedralSolid3D::impl_name_static() ) );
    builder.set_block_name( block, "pentagonal_prism" );
    auto mesh_builder =
        builder.block_mesh_builder< geode::PolyhedralSolid3D >( block );
    constexpr geode::index_t NB_SIDES{ 5 };
    for( const auto level : geode::Range{ 2 } )
    {
        for( const auto side : geode::Range{ NB_SIDES } )
        {
            const auto angle = 2. * M_PI * side / NB_SIDES;
            mesh_builder->create_point( geode::Point3D{
                { std::cos( angle ), std::sin( angle ), 1. * level } } );
        }
    }
    std::vector< std::vector< geode::local_index_t > > facets{
        { 4, 3, 2, 1, 0 }, { 5, 6, 7, 8, 9 }
    };
    for( const auto side : geode::LRange{ NB_SIDES } )
    {
        const auto next =
            static_cast< geode::local_index_t >( ( side + 1 ) % NB_SIDES );
        facets.push_back( { side, next,
            static_cast< geode::local_index_t >( next + NB_SIDES ),
            static_cast< geode::local_index_t >( side + NB_SIDES ) } );
    }
    mesh_builder->create_polyhedron( { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }, facets );
    builder.create_unique_vertices( 2 * NB_SIDES );
    for( const auto vertex : geode::Range{ 2 * NB_SIDES } )
    {
        builder.set_unique_vertex( { block.component_id(), vertex }, vertex );
    }
    auto problem = single_phase_problem( 0. );
    check_refused( problem, model, "pentagonal_prism" );
    check_refused(
        problem, model, "tetrahedra, hexahedra, prisms or pyramids" );
}

void test_ogs_invalid_names()
{
    auto model = load_picasso();
    add_constant_boundary( model, east_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, 1e7 );
    auto problem = single_phase_problem( 0. );
    geode::PhysicalProblemBuilder builder{ problem };
    auto material = problem.materials()[0];
    material.name = "a b";
    builder.set_material( 0, material );
    check_refused( problem, model, "a b" );
    material.name = "rock";
    builder.set_material( 0, material );
    builder.add_source_term( { "w/1", "water", 1. } );
    check_refused( problem, model, "w/1" );
}

std::string read_file( const std::string& filename )
{
    std::ifstream file{ filename, std::ios::binary };
    check( file.good(), absl::StrCat( "Cannot open file ", filename ) );
    return { std::istreambuf_iterator< char >{ file },
        std::istreambuf_iterator< char >{} };
}

// DataArray tag of the named array in a VTU file
std::string vtu_array_tag( const std::string& filename, std::string_view name )
{
    const auto content = read_file( filename );
    const auto name_position =
        content.find( absl::StrCat( "Name=\"", name, "\"" ) );
    check( name_position != std::string::npos,
        absl::StrCat( filename, " misses ", name ) );
    const auto start = content.rfind( '<', name_position );
    return content.substr( start, content.find( '>', name_position ) - start );
}

// Values of a named array of a VTU file in appended raw format
template < typename T >
std::vector< T > read_vtu_values(
    const std::string& filename, std::string_view name )
{
    const auto content = read_file( filename );
    const auto name_position =
        content.find( absl::StrCat( "Name=\"", name, "\"" ) );
    check( name_position != std::string::npos,
        absl::StrCat( filename, " misses ", name ) );
    constexpr std::string_view OFFSET{ "offset=\"" };
    const auto offset_position =
        content.find( OFFSET, name_position ) + OFFSET.size();
    const auto offset = std::stoull( content.substr( offset_position,
        content.find( '"', offset_position ) - offset_position ) );
    const auto data_position =
        content.find( '_', content.find( "<AppendedData" ) ) + 1 + offset;
    std::uint64_t nb_bytes{ 0 };
    std::memcpy( &nb_bytes, content.data() + data_position, sizeof nb_bytes );
    std::vector< T > values( nb_bytes / sizeof( T ) );
    std::memcpy( values.data(),
        content.data() + data_position + sizeof nb_bytes, nb_bytes );
    return values;
}

template < typename T >
std::vector< T > ascii_values( const pugi::xml_node& array )
{
    std::istringstream stream{ text( array ) };
    std::vector< T > values;
    T value{};
    while( stream >> value )
    {
        values.push_back( value );
    }
    return values;
}

std::vector< geode::Point3D > ascii_points( const pugi::xml_document& mesh )
{
    const auto coordinates = ascii_values< double >(
        mesh.select_node( "//Points/DataArray" ).node() );
    std::vector< geode::Point3D > points;
    for( geode::index_t index = 0; index + 2 < coordinates.size(); index += 3 )
    {
        points.push_back( geode::Point3D{ { coordinates[index],
            coordinates[index + 1], coordinates[index + 2] } } );
    }
    return points;
}

std::vector< std::string > boundary_files(
    std::string_view directory, std::string_view property )
{
    std::vector< std::string > files;
    for( const auto& entry :
        std::filesystem::directory_iterator{ std::string{ directory } } )
    {
        const auto name = entry.path().filename().string();
        if( absl::StrContains( name, property ) )
        {
            files.push_back( entry.path().string() );
        }
    }
    std::sort( files.begin(), files.end() );
    return files;
}

void test_ogs_mesh_only()
{
    auto model = load_picasso();
    add_constant_boundary( model, east_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, 1e7 );
    std::filesystem::remove_all( "ogs_mesh_only" );
    geode::BRepOgsExporter exporter{ model, "ogs_mesh_only" };
    exporter.run();
    const std::string bulk{ "ogs_mesh_only/ogs_mesh_only.vtu" };
    check( absl::StrContains(
               vtu_array_tag( bulk, "MaterialIDs" ), "type=\"Int32\"" ),
        "MaterialIDs must be Int32" );
    const auto material_ids =
        read_vtu_values< std::int32_t >( bulk, "MaterialIDs" );
    std::vector< std::int32_t > expected;
    std::int32_t region{ 0 };
    for( const auto& block : model.blocks() )
    {
        expected.insert(
            expected.end(), block.mesh().nb_polyhedra(), region++ );
    }
    check( material_ids == expected,
        "Without problem, MaterialIDs are the Block regions" );
    check( boundary_files( "ogs_mesh_only", "boundary_pressure" ).size() == 1,
        "One boundary mesh per uniform boundary_pressure group" );
    check( !std::filesystem::exists( "ogs_mesh_only/ogs_mesh_only.prj" ),
        "No project file without problem" );
}

void test_ogs_boundary_meshes()
{
    auto model = load_picasso();
    add_constant_boundary( model, east_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, 1e7 );
    const auto& property_info = model.physical_property_info(
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure );
    // A west Surface with one value per polygon
    std::optional< geode::uuid > west;
    const auto side = model.bounding_box().min().value( 0 );
    for( const auto& surface : model.surfaces() )
    {
        const auto box = surface.mesh().bounding_box();
        if( std::fabs( box.max().value( 0 ) - side ) < 1e-3 )
        {
            west = surface.id();
            break;
        }
    }
    check( west.has_value(), "picasso should have a west Surface" );
    const auto& west_surface = model.surface( west.value() );
    {
        geode::AttributeValues< double > values;
        values.default_value = 0.;
        auto& manager = west_surface.mesh().polygon_attribute_manager();
        manager.create_attribute< geode::VariableAttribute, double >(
            geode::physical_property_string(
                geode::PHYSICAL_PROPERTY_NAME::boundary_pressure ),
            property_info.attribute_id, values, {} );
        auto attribute =
            manager.find_attribute< geode::VariableAttribute, double >(
                property_info.attribute_id );
        for( const auto polygon :
            geode::Range{ west_surface.mesh().nb_polygons() } )
        {
            attribute->set_value( polygon, 1e6 + 10. * polygon );
        }
    }
    std::filesystem::remove_all( "ogs_boundaries" );
    geode::BRepOgsExporter exporter{ model, "ogs_boundaries" };
    exporter.run();
    auto [solid, mappings] = geode::convert_brep_into_solid( model );
    const auto files = boundary_files( "ogs_boundaries", "boundary_pressure" );
    check( files.size() == 2, "One uniform and one variable boundary mesh" );
    bool has_variable{ false };
    for( const auto& file : files )
    {
        const auto mesh = load_xml( file );
        const auto ids = ascii_values< std::uint64_t >(
            data_array( mesh, "PointData", "bulk_node_ids" ) );
        const auto points = ascii_points( mesh );
        check( ids.size() == points.size(), "One bulk node id per point" );
        for( const auto index : geode::Indices{ ids } )
        {
            check( solid->point( ids[index] ) == points[index],
                "bulk_node_ids must point to the same bulk vertex" );
        }
        const auto values_node = data_array( mesh, "PointData", "values" );
        if( !values_node )
        {
            for( const auto& point : points )
            {
                check( std::fabs( point.value( 0 )
                                  - model.bounding_box().max().value( 0 ) )
                           < 1e-3,
                    "Uniform group points lie on the east Surfaces" );
            }
            continue;
        }
        has_variable = true;
        check(
            std::string{ values_node.attribute( "type" ).value() } == "Float64",
            "Variable values must be Float64" );
        const auto values = ascii_values< double >( values_node );
        const auto& surface_mesh = west_surface.mesh();
        for( const auto vertex : geode::Range{ surface_mesh.nb_vertices() } )
        {
            const auto& polygons =
                surface_mesh.polygons_around_vertex( vertex );
            double sum{ 0 };
            for( const auto& polygon : polygons )
            {
                sum += 1e6 + 10. * polygon.polygon_id;
            }
            const auto expected = sum / polygons.size();
            const auto position = std::find(
                points.begin(), points.end(), surface_mesh.point( vertex ) );
            check( position != points.end(),
                "Every west Surface vertex is in its boundary mesh" );
            check( std::fabs( values[position - points.begin()] - expected )
                       < 1e-6,
                "A vertex value is the mean of its polygon values" );
        }
    }
    check( has_variable, "The west Surface should have a values array" );
}

void test_ogs_well_mesh()
{
    auto model = load_picasso();
    auto [solid, mappings] = geode::convert_brep_into_solid( model );
    geode::OpenGeodePointSet3D well;
    auto builder = geode::PointSetBuilder3D::create( well );
    std::vector< std::uint64_t > expected;
    for( const geode::index_t polyhedron : { 0, 100, 200 } )
    {
        const auto vertex = solid->polyhedron_vertex( { polyhedron, 0 } );
        const auto& point = solid->point( vertex );
        const auto center = solid->polyhedron_barycenter( polyhedron );
        builder->create_point( point + ( center - point ) * 0.01 );
        expected.push_back( vertex );
    }
    std::filesystem::remove_all( "ogs_well" );
    geode::BRepOgsExporter exporter{ model, "ogs_well" };
    exporter.add_well_perforations( well, "producer" );
    exporter.run();
    const auto mesh = load_xml( "ogs_well/ogs_well_well_producer.vtu" );
    check( ascii_values< std::uint64_t >(
               data_array( mesh, "PointData", "bulk_node_ids" ) )
               == expected,
        "Well vertices are the nearest vertices of the containing cells" );
    check( text( data_array( mesh, "Cells", "types" ) ) == "1 1 1",
        "Well cells are VTK vertices" );
}

void test_ogs_regular_grid()
{
    geode::BRep model;
    geode::BRepBuilder builder{ model };
    const auto& block = model.block( builder.add_block(
        geode::OpenGeodeRegularGrid3D::impl_name_static() ) );
    builder.block_mesh_builder< geode::RegularGrid3D >( block )
        ->initialize_grid( geode::Point3D{ { 0., 0., 0. } }, { 2, 2, 2 }, 1. );
    const auto nb_vertices = block.mesh().nb_vertices();
    builder.create_unique_vertices( nb_vertices );
    for( const auto vertex : geode::Range{ nb_vertices } )
    {
        builder.set_unique_vertex( { block.component_id(), vertex }, vertex );
    }
    std::filesystem::remove_all( "ogs_grid" );
    geode::BRepOgsExporter exporter{ model, "ogs_grid" };
    exporter.run();
    const auto types =
        read_vtu_values< std::uint8_t >( "ogs_grid/ogs_grid.vtu", "types" );
    check( types == std::vector< std::uint8_t >( 8, 12 ),
        "RegularGrid Blocks are exported as hexahedra" );
}

void test_ogs_output_directories()
{
    const auto model = load_picasso();
    bool has_thrown{ false };
    try
    {
        geode::BRepOgsExporter exporter{ model, "" };
    }
    catch( const geode::OpenGeodeException& )
    {
        has_thrown = true;
    }
    check( has_thrown, "An empty directory must throw" );
    std::filesystem::remove_all( "ogs_nested" );
    geode::BRepOgsExporter nested{ model, "ogs_nested/./directory/" };
    nested.run();
    check( std::filesystem::exists( "ogs_nested/directory/directory.vtu" ),
        "Nested directories must be created" );
}

pugi::xml_node select(
    const pugi::xml_document& document, std::string_view xpath )
{
    return document.select_node( std::string{ xpath }.c_str() ).node();
}

std::string select_text(
    const pugi::xml_document& document, std::string_view xpath )
{
    return text( select( document, xpath ) );
}

// Parameter node of the named parameter
pugi::xml_node parameter(
    const pugi::xml_document& project, std::string_view name )
{
    return select(
        project, absl::StrCat( "/OpenGeoSysProject/parameters/parameter[name='",
                     name, "']" ) );
}

geode::BRep picasso_with_boundary()
{
    auto model = load_picasso();
    add_constant_boundary( model, east_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure, 1e7 );
    return model;
}

pugi::xml_document export_project( const geode::BRep& model,
    const geode::PhysicalProblem& problem,
    std::string_view name )
{
    std::filesystem::remove_all( std::string{ name } );
    geode::BRepOgsExporter exporter{ model, name };
    exporter.set_physical_problem( problem );
    exporter.run();
    return load_xml( absl::StrCat( name, "/", name, ".prj" ) );
}

void test_ogs_liquid_flow_picasso()
{
    const auto model = picasso_with_boundary();
    auto problem = single_phase_problem( 4.5e-10 );
    geode::PhysicalProblemBuilder builder{ problem };
    auto other = rock_material( "other" );
    other.regions = { ( *model.blocks().begin() ).id() };
    builder.add_material( other );
    const auto project = export_project( model, problem, "ogs_picasso" );
    check( select_text( project, "/OpenGeoSysProject/meshes/mesh[1]" )
               == "ogs_picasso.vtu",
        "The bulk mesh comes first" );
    check( select_text( project, "//processes/process/type" ) == "LIQUID_FLOW",
        "LIQUID_FLOW process" );
    check( select_text( project, "//processes/process/equation_balance_type" )
               == "mass",
        "Mass balance" );
    check( project.select_nodes( "//media/medium" ).size() == 2,
        "One medium per material" );
    const auto permeability = select( project,
        "//media/medium[@id='0']/properties/property[name='permeability']" );
    check( text( permeability.child( "type" ) ) == "Constant"
               || text( permeability.child( "type" ) ) == "Parameter",
        "Permeability property" );
    for( const auto& mesh_node :
        project.select_nodes( "//boundary_condition/mesh" ) )
    {
        const auto stem = text( mesh_node.node() );
        check( std::filesystem::exists(
                   absl::StrCat( "ogs_picasso/", stem, ".vtu" ) ),
            absl::StrCat( "Boundary mesh ", stem, " must exist" ) );
        check( static_cast< bool >( select( project,
                   absl::StrCat( "/OpenGeoSysProject/meshes/mesh[text()='",
                       stem, ".vtu']" ) ) ),
            "Boundary meshes are listed in <meshes>" );
    }
    check( project.select_nodes( "//boundary_condition" ).size() == 1,
        "One Dirichlet condition for the east Surfaces" );
}

void test_ogs_permeability_field()
{
    auto model = picasso_with_boundary();
    const geode::uuid attribute_id;
    geode::AttributeValues< std::array< double, 3 > > values;
    values.default_value = { 1e-12, 2e-12, 3e-12 };
    for( const auto& block : model.blocks() )
    {
        block.mesh()
            .polyhedron_attribute_manager()
            .create_attribute< geode::VariableAttribute,
                std::array< double, 3 > >( "perm", attribute_id, values, {} );
    }
    geode::BRepBuilder{ model }.set_physical_property(
        geode::PHYSICAL_PROPERTY_NAME::permeability,
        geode::Block3D::component_type_static(), attribute_id );
    const auto project = export_project(
        model, single_phase_problem( 0. ), "ogs_permeability_field" );
    const auto property = select(
        project, "//media/medium/properties/property[name='permeability']" );
    check( text( property.child( "type" ) ) == "Parameter",
        "A permeability field is a Parameter property" );
    const auto field =
        parameter( project, text( property.child( "parameter_name" ) ) );
    check( text( field.child( "type" ) ) == "MeshElement"
               && text( field.child( "field_name" ) ) == "permeability",
        "The permeability parameter is MeshElement on the permeability "
        "array" );
    const auto cells = read_vtu_values< double >(
        "ogs_permeability_field/ogs_permeability_field.vtu", "permeability" );
    check( cells.size() >= 3 && cells[0] == 1e-12 && cells[1] == 2e-12
               && cells[2] == 3e-12,
        "The permeability array carries the Block values" );
}

void test_ogs_compressible_fluid()
{
    const auto model = picasso_with_boundary();
    const auto project = export_project(
        model, single_phase_problem( 4.5e-10 ), "ogs_compressible" );
    const auto density =
        select( project, "//phase[type='AqueousLiquid']/properties/"
                         "property[name='density']" );
    check( text( density.child( "type" ) ) == "Exponential",
        "A compressible fluid has an Exponential density" );
    check( std::stod( text( density.child( "exponent" ).child( "factor" ) ) )
               == 4.5e-10,
        "The density factor is the fluid compressibility" );
    const auto storage =
        select( project, "//medium/properties/property[name='storage']" );
    check( text( storage.child( "type" ) ) == "Constant"
               && std::fabs(
                      std::stod( text( storage.child( "value" ) ) ) - 2e-10 )
                      < 1e-24,
        "Storage is porosity times pore compressibility" );
    const auto incompressible = export_project(
        model, single_phase_problem( 0. ), "ogs_incompressible" );
    check( select_text( incompressible, "//phase[type='AqueousLiquid']/"
                                        "properties/property[name='density']/"
                                        "type" )
               == "Constant",
        "An incompressible fluid has a Constant density" );
}

void test_ogs_wells_rates()
{
    const auto model = picasso_with_boundary();
    auto [solid, mappings] = geode::convert_brep_into_solid( model );
    geode::OpenGeodePointSet3D well;
    auto builder = geode::PointSetBuilder3D::create( well );
    for( const geode::index_t polyhedron : { 0, 100, 200 } )
    {
        const auto vertex = solid->polyhedron_vertex( { polyhedron, 0 } );
        const auto& point = solid->point( vertex );
        builder->create_point(
            point
            + ( solid->polyhedron_barycenter( polyhedron ) - point ) * 0.01 );
    }
    auto problem = single_phase_problem( 0. );
    geode::PhysicalProblemBuilder{ problem }.add_source_term(
        { "injector", "water", 3. } );
    std::filesystem::remove_all( "ogs_wells" );
    geode::BRepOgsExporter exporter{ model, "ogs_wells" };
    exporter.add_well_perforations( well, "injector" );
    exporter.set_physical_problem( problem );
    exporter.run();
    const auto project = load_xml( "ogs_wells/ogs_wells.prj" );
    const auto source = select( project, "//source_terms/source_term" );
    check( text( source.child( "type" ) ) == "Nodal"
               && text( source.child( "mesh" ) ) == "ogs_wells_well_injector",
        "A Nodal source term on the well mesh" );
    const auto rate = parameter( project, text( source.child( "parameter" ) ) );
    check( text( rate.child( "type" ) ) == "Constant"
               && std::stod( text( rate.child( "value" ) ) ) == 1.,
        "The rate is split evenly over the well vertices" );
}

// Two hexahedra Blocks sharing the face x = 1, top face first if mirrored
geode::BRep two_hexahedra_brep( bool mirrored = false )
{
    geode::BRep model;
    geode::BRepBuilder builder{ model };
    builder.create_unique_vertices( 12 );
    const auto unique_vertex = []( geode::index_t x, geode::index_t y,
                                   geode::index_t z ) {
        return x + 3 * ( y + 2 * z );
    };
    for( const geode::index_t block_x : { 0, 1 } )
    {
        const auto& block = model.block( builder.add_block(
            geode::OpenGeodeHybridSolid3D::impl_name_static() ) );
        auto mesh_builder =
            builder.block_mesh_builder< geode::HybridSolid3D >( block );
        geode::index_t vertex{ 0 };
        for( const geode::index_t z : { 0, 1 } )
        {
            for( const auto& [dx, y] :
                std::array< std::pair< geode::index_t, geode::index_t >, 4 >{
                    { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } } } )
            {
                mesh_builder->create_point( geode::Point3D{
                    { 1. * ( block_x + dx ), 1. * y, 1. * z } } );
                builder.set_unique_vertex( { block.component_id(), vertex++ },
                    unique_vertex( block_x + dx, y, z ) );
            }
        }
        if( mirrored )
        {
            mesh_builder->create_hexahedron( { 4, 5, 6, 7, 0, 1, 2, 3 } );
        }
        else
        {
            mesh_builder->create_hexahedron( { 0, 1, 2, 3, 4, 5, 6, 7 } );
        }
    }
    return model;
}

void test_ogs_initial_pressure()
{
    const auto model = two_hexahedra_brep();
    auto problem = single_phase_problem( 0. );
    const auto uniform = export_project( model, problem, "ogs_uniform_p0" );
    const auto initial = select_text(
        uniform, "//process_variables/process_variable/initial_condition" );
    check( text( parameter( uniform, initial ).child( "type" ) ) == "Constant",
        "One initial pressure gives a Constant parameter" );
    geode::PhysicalProblemBuilder{ problem }.add_initial_condition(
        std::make_shared< geode::UniformInitialCondition >(
            std::vector< geode::uuid >{ ( *model.blocks().begin() ).id() },
            2e7 ) );
    const auto varying = export_project( model, problem, "ogs_varying_p0" );
    const auto field = parameter( varying,
        select_text( varying,
            "//process_variables/process_variable/initial_condition" ) );
    check( text( field.child( "type" ) ) == "MeshNode"
               && text( field.child( "field_name" ) ) == "initial_pressure",
        "Several initial pressures give a MeshNode parameter" );
    const auto values = read_vtu_values< double >(
        "ogs_varying_p0/ogs_varying_p0.vtu", "initial_pressure" );
    check( values.size() == 12, "One initial pressure per vertex" );
    check( values[0] == 2e7 && values[1] == 1.5e7 && values[2] == 1e7,
        "Shared vertices get the mean of their cells" );
}

void test_ogs_solver_options()
{
    const auto model = picasso_with_boundary();
    auto problem = single_phase_problem( 0. );
    geode::PhysicalProblemBuilder{ problem }.set_solver_option(
        "ogs", "linear/solver_type", "BiCGSTAB" );
    geode::PhysicalProblemBuilder{ problem }.set_solver_option(
        "ogs", "nonlinear/max_iter", "7" );
    const auto project = export_project( model, problem, "ogs_options" );
    check( select_text( project, "//linear_solvers/linear_solver/eigen/"
                                 "solver_type" )
               == "BiCGSTAB",
        "A linear option sets the element text" );
    check(
        select_text( project, "//nonlinear_solvers/nonlinear_solver/max_iter" )
            == "7",
        "A nonlinear option sets the element text" );
}

void test_ogs_conflicting_dirichlet()
{
    auto model = load_picasso();
    const auto first = ( *model.surfaces().begin() ).id();
    std::optional< geode::uuid > neighbour;
    absl::flat_hash_set< geode::index_t > first_vertices;
    for( const auto vertex :
        geode::Range{ model.surface( first ).mesh().nb_vertices() } )
    {
        first_vertices.insert( model.unique_vertex(
            { model.surface( first ).component_id(), vertex } ) );
    }
    for( const auto& surface : model.surfaces() )
    {
        if( surface.id() == first )
        {
            continue;
        }
        for( const auto vertex : geode::Range{ surface.mesh().nb_vertices() } )
        {
            if( first_vertices.contains( model.unique_vertex(
                    { surface.component_id(), vertex } ) ) )
            {
                neighbour = surface.id();
                break;
            }
        }
        if( neighbour )
        {
            break;
        }
    }
    check( neighbour.has_value(), "A Surface should touch the first one" );
    const geode::uuid attribute_id;
    for( const auto& [surface_id, value] :
        std::array< std::pair< geode::uuid, double >, 2 >{
            { { first, 1e5 }, { neighbour.value(), 2e5 } } } )
    {
        geode::AttributeValues< double > values;
        values.default_value = value;
        model.surface( surface_id )
            .mesh()
            .polygon_attribute_manager()
            .create_attribute< geode::ConstantAttribute, double >(
                "boundary_pressure", attribute_id, values, {} );
    }
    geode::BRepBuilder{ model }.set_physical_property(
        geode::PHYSICAL_PROPERTY_NAME::boundary_pressure,
        geode::Surface3D::component_type_static(), attribute_id );
    const auto report = ogs_report( single_phase_problem( 0. ), model );
    geode::index_t nb_conflicts{ 0 };
    for( const auto& issue : report.issues() )
    {
        if( !issue.is_error()
            && absl::StrContains( issue.message, "share vertices" ) )
        {
            nb_conflicts++;
        }
    }
    check( report.is_simulable() && nb_conflicts == 1,
        absl::StrCat( "One warning per pair of conflicting Surfaces:\n",
            report.string() ) );
}

void test_ogs_run_twice()
{
    const auto model = picasso_with_boundary();
    const auto problem = single_phase_problem( 4.5e-10 );
    std::filesystem::remove_all( "ogs_twice" );
    geode::BRepOgsExporter exporter{ model, "ogs_twice" };
    exporter.set_physical_problem( problem );
    exporter.run();
    std::vector< std::string > first;
    std::vector< std::string > names;
    for( const auto& entry :
        std::filesystem::directory_iterator{ "ogs_twice" } )
    {
        names.push_back( entry.path().string() );
    }
    std::sort( names.begin(), names.end() );
    for( const auto& name : names )
    {
        first.push_back( read_file( name ) );
    }
    exporter.run();
    for( const auto index : geode::Indices{ names } )
    {
        check( read_file( names[index] ) == first[index],
            absl::StrCat( names[index], " changed on the second run" ) );
    }
}

void test_ogs_structural_model_materials()
{
    geode::StructuralModel model{ picasso_with_boundary() };
    geode::StructuralModelBuilder builder{ model };
    const auto& unit_id = builder.add_stratigraphic_unit();
    std::vector< geode::uuid > unit_blocks;
    for( const auto& block : model.blocks() )
    {
        builder.add_block_in_stratigraphic_unit(
            block, model.stratigraphic_unit( unit_id ) );
        unit_blocks.push_back( block.id() );
        if( unit_blocks.size() == 2 )
        {
            break;
        }
    }
    auto problem = single_phase_problem( 0. );
    auto unit_material = rock_material( "unit" );
    unit_material.regions = { unit_id };
    const auto unit_index =
        geode::PhysicalProblemBuilder{ problem }.add_material( unit_material );
    std::filesystem::remove_all( "ogs_units" );
    geode::StructuralModelOgsExporter exporter{ model, "ogs_units" };
    exporter.set_physical_problem( problem );
    exporter.run();
    const auto material_ids = read_vtu_values< std::int32_t >(
        "ogs_units/ogs_units.vtu", "MaterialIDs" );
    geode::index_t cell{ 0 };
    for( const auto& block : model.blocks() )
    {
        const auto in_unit = absl::c_linear_search( unit_blocks, block.id() );
        for( const auto polyhedron :
            geode::Range{ block.mesh().nb_polyhedra() } )
        {
            geode_unused( polyhedron );
            check( ( material_ids[cell]
                       == static_cast< std::int32_t >( unit_index ) )
                       == in_unit,
                "Cells of the unit Blocks get the unit material" );
            cell++;
        }
    }
}

std::shared_ptr< geode::DeadOilFluid > varying_oil_fluid()
{
    geode::HydrocarbonPvtTable oil;
    oil.phase = geode::PHASE::oil;
    oil.pressures = { 1e5, 1e7, 2e7 };
    oil.formation_volume_factors = { 1.2, 1.1, 1.05 };
    oil.viscosities = { 2e-3, 1.5e-3, 1.4e-3 };
    return std::make_shared< geode::DeadOilFluid >(
        std::vector< geode::PHASE >{ geode::PHASE::water, geode::PHASE::oil },
        std::vector< geode::FluidComponent >{
            { "water", 0.018 }, { "oil", 0.2 } },
        std::vector< double >{ 1000., 800. },
        std::vector< geode::HydrocarbonPvtTable >{ oil },
        geode::WaterPvt{ 1e5, 1.01, 0., 1e-3 } );
}

geode::PhysicalProblem varying_oil_problem()
{
    auto problem = two_phase_problem();
    const_cast< geode::FlowModel& >( problem.flow() )
        .set_fluid( varying_oil_fluid() );
    return problem;
}

double expected_capillary_pressure(
    const geode::PhysicalProblem& problem, double water_fraction )
{
    const auto& fluid =
        dynamic_cast< const geode::DeadOilFluid& >( problem.fluid() );
    const std::array< double, 2 > composition{ water_fraction,
        1. - water_fraction };
    const auto saturation = geode::internal::water_saturation(
        fluid, geode::COMPOSITION_BASIS::mass, composition, 1e7 );
    return geode::internal::capillary_pressure(
        *problem.materials()[0].capillary_pressure, saturation );
}

void test_ogs_two_phase_project()
{
    const auto model = two_hexahedra_brep();
    const auto problem = varying_oil_problem();
    const auto project = export_project( model, problem, "ogs_two_phase" );
    check( select_text( project, "//processes/process/type" )
               == "TWOPHASE_FLOW_PP",
        "TWOPHASE_FLOW_PP process" );
    for( const auto& property_name : { "saturation", "relative_permeability",
             "relative_permeability_nonwetting_phase" } )
    {
        const auto property = select(
            project, absl::StrCat( "//medium/properties/property[name='",
                         property_name, "']" ) );
        check( text( property.child( "type" ) ) == "Curve",
            absl::StrCat( property_name, " is a Curve" ) );
        const auto curve = select(
            project, absl::StrCat( "/OpenGeoSysProject/curves/curve[name='",
                         text( property.child( "curve" ) ), "']" ) );
        check( static_cast< bool >( curve ),
            absl::StrCat( "Curve of ", property_name, " must exist" ) );
        const auto coords = ascii_values< double >( curve.child( "coords" ) );
        const auto values = ascii_values< double >( curve.child( "values" ) );
        check( coords.size() == values.size() && coords.size() > 1,
            "As many coords as values" );
        for( const auto index : geode::Range{ 1, coords.size() } )
        {
            check( coords[index] > coords[index - 1],
                "Curve coords must strictly increase" );
        }
    }
    const auto density_curve = select(
        project, absl::StrCat( "/OpenGeoSysProject/curves/curve[name='",
                     select_text( project, "//phase[type='Gas']/properties/"
                                           "property[name='density']/curve" ),
                     "']" ) );
    const auto densities =
        ascii_values< double >( density_curve.child( "values" ) );
    check( ascii_values< double >( density_curve.child( "coords" ) )
                   == std::vector< double >{ 1e5, 1e7, 2e7 }
               && std::fabs( densities[0] - 800. / 1.2 ) < 1e-9
               && std::fabs( densities[2] - 800. / 1.05 ) < 1e-9,
        "Oil density is surface density / Bo at the PVT pressures" );
    const auto initial = parameter( project,
        select_text( project, "//process_variable[name='capillary_pressure']/"
                              "initial_condition" ) );
    check( text( initial.child( "type" ) ) == "Constant",
        "Uniform initial capillary pressure" );
    check( std::fabs( std::stod( text( initial.child( "value" ) ) )
                      - expected_capillary_pressure( problem, 0.3 ) )
               < 1e-6,
        "Initial capillary pressure is pc(Sw(composition))" );
}

void test_ogs_two_phase_boundary()
{
    auto model = picasso_with_boundary();
    add_constant_boundary( model, east_surfaces( model ),
        geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction, 0.8 );
    const auto problem = varying_oil_problem();
    const auto project =
        export_project( model, problem, "ogs_two_phase_boundary" );
    const auto capillary = select( project,
        "//process_variable[name='capillary_pressure']/boundary_conditions/"
        "boundary_condition" );
    check( text( capillary.child( "type" ) ) == "Dirichlet",
        "Water fraction gives a Dirichlet on capillary_pressure" );
    const auto value =
        parameter( project, text( capillary.child( "parameter" ) ) );
    check( text( value.child( "type" ) ) == "Constant"
               && std::fabs( std::stod( text( value.child( "value" ) ) )
                             - expected_capillary_pressure( problem, 0.8 ) )
                      < 1e-6,
        "The capillary pressure comes from the water fraction" );
    const auto gas = select( project,
        "//process_variable[name='gas_pressure']/boundary_conditions/"
        "boundary_condition" );
    check( text( gas.child( "type" ) ) == "Dirichlet"
               && std::stod( text(
                      parameter( project, text( gas.child( "parameter" ) ) )
                          .child( "value" ) ) )
                      == 1e7,
        "Boundary pressure gives a Dirichlet on gas_pressure" );
}

void test_ogs_two_phase_permeability()
{
    const auto model = two_hexahedra_brep();
    auto problem = varying_oil_problem();
    auto material = problem.materials()[0];
    material.permeability = std::array< double, 3 >{ 1e-13, 1e-13, 1e-13 };
    geode::PhysicalProblemBuilder{ problem }.set_material( 0, material );
    const auto project =
        export_project( model, problem, "ogs_two_phase_permeability" );
    const auto property =
        select( project, "//medium/properties/property[name='permeability']" );
    const auto values =
        parameter( project, text( property.child( "parameter_name" ) ) );
    check( ascii_values< double >( values.child( "values" ) )
               == std::vector< double >{ 1e-13 },
        "Two-phase permeability is a scalar" );
    material.permeability = std::array< double, 3 >{ 1e-13, 1e-13, 1e-14 };
    geode::PhysicalProblemBuilder{ problem }.set_material( 0, material );
    check_refused( problem, model, "isotropic permeability" );

    auto field_model = two_hexahedra_brep();
    const geode::uuid attribute_id;
    geode::AttributeValues< std::array< double, 3 > > field_values;
    field_values.default_value = { 1e-12, 1e-12, 1e-13 };
    for( const auto& block : field_model.blocks() )
    {
        block.mesh()
            .polyhedron_attribute_manager()
            .create_attribute< geode::VariableAttribute,
                std::array< double, 3 > >(
                "perm", attribute_id, field_values, {} );
    }
    geode::BRepBuilder{ field_model }.set_physical_property(
        geode::PHYSICAL_PROPERTY_NAME::permeability,
        geode::Block3D::component_type_static(), attribute_id );
    bool has_thrown{ false };
    try
    {
        geode_unused( export_project(
            field_model, varying_oil_problem(), "ogs_anisotropic_field" ) );
    }
    catch( const geode::OpenGeodeException& exception )
    {
        has_thrown = absl::StrContains( exception.what(), "isotropic" );
    }
    check( has_thrown,
        "An anisotropic permeability field must be refused in two-phase" );
}

// VTK cells must be positively oriented for OGS: det(v1 - v0, v3 - v0,
// v4 - v0) > 0 for hexahedra
void test_ogs_cell_orientation()
{
    const auto model = two_hexahedra_brep( true );
    std::filesystem::remove_all( "ogs_orientation" );
    geode::BRepOgsExporter exporter{ model, "ogs_orientation" };
    exporter.run();
    auto [solid, mappings] = geode::convert_brep_into_solid( model );
    const auto connectivity = read_vtu_values< std::int64_t >(
        "ogs_orientation/ogs_orientation.vtu", "connectivity" );
    check( connectivity.size() == 16, "Two hexahedra expected" );
    for( const geode::index_t cell : { 0, 1 } )
    {
        const auto point = [&]( geode::index_t local ) {
            return solid->point(
                static_cast< geode::index_t >( connectivity[8 * cell + local] ) );
        };
        const geode::Vector3D first{ point( 0 ), point( 1 ) };
        const geode::Vector3D second{ point( 0 ), point( 3 ) };
        const geode::Vector3D third{ point( 0 ), point( 4 ) };
        check( first.cross( second ).dot( third ) > 0,
            "Hexahedra must be positively oriented" );
    }
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        test_ogs_mesh_file();
        test_ogs_structural_model_regions();
        test_ogs_rules_process();
        test_ogs_rules_storage_and_temperature();
        test_ogs_rules_output_times();
        test_ogs_rules_saturation();
        test_ogs_rules_relperm();
        test_ogs_registered_checker();
        test_ogs_checker_refusals();
        test_ogs_unsupported_cells();
        test_ogs_invalid_names();
        test_ogs_mesh_only();
        test_ogs_boundary_meshes();
        test_ogs_well_mesh();
        test_ogs_regular_grid();
        test_ogs_output_directories();
        test_ogs_liquid_flow_picasso();
        test_ogs_permeability_field();
        test_ogs_compressible_fluid();
        test_ogs_wells_rates();
        test_ogs_initial_pressure();
        test_ogs_solver_options();
        test_ogs_conflicting_dirichlet();
        test_ogs_run_twice();
        test_ogs_structural_model_materials();
        test_ogs_two_phase_project();
        test_ogs_two_phase_boundary();
        test_ogs_two_phase_permeability();
        test_ogs_cell_orientation();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
