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
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <pugixml.hpp>

#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/types.hpp>
#include <geode/tests_config.hpp>

#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/polyhedral_solid_builder.hpp>
#include <geode/mesh/core/geode/geode_polyhedral_solid.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

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
#include <geode/geosciences/simulation/helpers/simulation_checker.hpp>
#include <geode/geosciences/simulation/core/thermal_parameters.hpp>

#include <geode/geosciences_io/model/common.hpp>
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
        const auto xpath = absl::StrCat( "//", section, "/DataArray[@Name='",
            name, "']" );
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
    const auto regions = geode::internal::structural_model_block_regions(
        model );
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
        std::vector< geode::FluidComponent >{ { "water", 0.018 } }, 1000.,
        1e-3, compressibility, 0., 1e5 );
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
    const auto three_phase =
        flow_problem( std::make_shared< geode::MultiphaseFlow >(),
            water_oil_fluid( { geode::PHASE::water, geode::PHASE::oil,
                geode::PHASE::gas } ) );
    check( !geode::internal::ogs_process( three_phase ),
        "Three-phase dead oil is not supported" );
}

void test_ogs_rules_storage_and_temperature()
{
    check( std::fabs( geode::internal::ogs_storage( 0.2, 1e-9 ) - 2e-10 )
               < 1e-25,
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
            const auto ratio = ( coord - curve.coords[index - 1] )
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
    check( curve.coords.size() == curve.values.size(),
        "Saturation curve sizes" );
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
        geode_unused(
            geode::internal::saturation_curve( non_monotone ) );
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
    return geode::load_one( geode::load_brep( model_filename ), model_filename );
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
    auto problem = flow_problem(
        std::make_shared< geode::MultiphaseFlow >(
            geode::COMPOSITION_BASIS::mass ),
        water_oil_fluid( { geode::PHASE::water, geode::PHASE::oil } ) );
    geode::PhysicalProblemBuilder builder{ problem };
    auto rock = rock_material( "rock" );
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
        absl::StrCat( "Missing issue \"", expected, "\" in\n",
            report.string() ) );
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
        absl::StrCat( "Picasso single-phase should be simulable:\n",
            report.string() ) );
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
        fluid.set_fluid( water_oil_fluid( { geode::PHASE::water,
            geode::PHASE::oil, geode::PHASE::gas } ) );
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
                std::vector< geode::TwoPhaseCapillaryPressureTable >{
                    table } );
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
        const auto next = static_cast< geode::local_index_t >(
            ( side + 1 ) % NB_SIDES );
        facets.push_back( { side, next,
            static_cast< geode::local_index_t >( next + NB_SIDES ),
            static_cast< geode::local_index_t >( side + NB_SIDES ) } );
    }
    mesh_builder->create_polyhedron(
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }, facets );
    builder.create_unique_vertices( 2 * NB_SIDES );
    for( const auto vertex : geode::Range{ 2 * NB_SIDES } )
    {
        builder.set_unique_vertex( { block.component_id(), vertex }, vertex );
    }
    auto problem = single_phase_problem( 0. );
    check_refused( problem, model, "pentagonal_prism" );
    check_refused( problem, model, "tetrahedra, hexahedra, prisms or pyramids" );
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
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
