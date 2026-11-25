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

#include <absl/strings/str_cat.h>

#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/types.hpp>
#include <geode/tests_config.hpp>

#include <geode/geometry/point.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/io/brep_input.hpp>

#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.hpp>
#include <geode/geosciences/explicit/representation/builder/structural_model_builder.hpp>
#include <geode/geosciences/explicit/representation/core/structural_model.hpp>

#include <geode/geosciences/simulation/builder/physical_problem_builder.hpp>
#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/schedule.hpp>
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
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
