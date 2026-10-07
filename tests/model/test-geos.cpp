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
#include <optional>
#include <string>
#include <string_view>

#include <pugixml.hpp>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/variable_attribute.hpp>
#include <geode/tests_config.hpp>

#include <absl/strings/str_cat.h>

#include <geode/geosciences_io/model/helpers/brep_geos_export.hpp>

#include <geode/geometry/point.hpp>
#include <geode/io/mesh/common.hpp>
#include <geode/io/model/common.hpp>
#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/core/point_set.hpp>

#include <geode/mesh/core/geode/geode_point_set.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>
#include <geode/model/mixin/core/surface.hpp>
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

// GEOS imports each surface as a node set named after its region attribute:
// the exporter numbers the Blocks first, then the Surfaces
std::string geos_surface_set_name(
    const geode::BRep& model, const geode::Surface3D& surface )
{
    auto region_id = model.nb_blocks();
    for( const auto& model_surface : model.surfaces() )
    {
        if( model_surface.id() == surface.id() )
        {
            break;
        }
        region_id++;
    }
    return absl::StrCat( region_id );
}

struct FaceBoundaryCondition
{
    std::string_view name;
    std::string_view field_name;
    std::optional< int > component;
    std::string_view scale;
};

// CompositionalMultiphaseFVM requires pressure, temperature and the whole
// composition for face Dirichlet boundary conditions
static constexpr std::array< FaceBoundaryCondition, 4 > SINK_FACE_CONDITIONS{ {
    { "boundaryPressure", "pressure", std::nullopt, "2.7579e+7" },
    { "boundaryTemperature", "temperature", std::nullopt, "300" },
    { "boundaryComposition_oil", "globalCompFraction", 0, "0.9995" },
    { "boundaryComposition_water", "globalCompFraction", 1, "0.0005" },
} };

void add_surface_boundary_condition( std::string_view xml_file,
    const geode::BRep& model,
    const geode::Surface3D& surface )
{
    pugi::xml_document document;
    geode::OpenGeodeGeosciencesIOModelException::test(
        static_cast< bool >(
            document.load_file( geode::to_string( xml_file ).c_str() ) ),
        "[Test] Cannot load ", xml_file );
    auto field_specifications =
        document.child( "Problem" ).append_child( "FieldSpecifications" );
    const auto set_name = geos_surface_set_name( model, surface );
    const auto set_names = absl::StrCat( "{ ", set_name, " }" );
    for( const auto& condition : SINK_FACE_CONDITIONS )
    {
        auto field = field_specifications.append_child( "FieldSpecification" );
        field.append_attribute( "name" ).set_value(
            absl::StrCat( condition.name, "_", set_name ).c_str() );
        field.append_attribute( "setNames" ).set_value( set_names.c_str() );
        field.append_attribute( "objectPath" ).set_value( "faceManager" );
        field.append_attribute( "fieldName" )
            .set_value( geode::to_string( condition.field_name ).c_str() );
        if( condition.component )
        {
            field.append_attribute( "component" )
                .set_value( *condition.component );
        }
        field.append_attribute( "scale" ).set_value(
            geode::to_string( condition.scale ).c_str() );
    }
    document.save_file( geode::to_string( xml_file ).c_str(),
        PUGIXML_TEXT( "    " ), pugi::format_indent_attributes );
}

void test_grid_geos()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "grid_geos_with_physical_properties.og_brep" ) );
    // The model only has one Surface: its border at minimum X
    const auto& border_surface = *model.surfaces().begin();
    geode::BRepGeosExporter exporter( model, "grid_geos" );
    add_spe10_wells( exporter );
    exporter.run();
    add_surface_boundary_condition(
        "grid_geos/grid_geos_simulation.xml", model, border_surface );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        geode::OpenGeodeIOMeshLibrary::initialize();
        geode::OpenGeodeIOModelLibrary::initialize();
        test_picasso();
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
