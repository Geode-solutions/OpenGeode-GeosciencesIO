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

#include <geode/tests_config.hpp>

#include <fstream>
#include <string>

#include <absl/container/flat_hash_set.h>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/string.hpp>

#include <geode/geometry/bounding_box.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/line.hpp>
#include <geode/model/mixin/core/model_boundary.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/mixin/core/surface_collection.hpp>

#include <geode/geosciences/explicit/mixin/core/fault_block.hpp>
#include <geode/geosciences/explicit/mixin/core/horizon.hpp>
#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.hpp>
#include <geode/geosciences/explicit/representation/builder/structural_model_builder.hpp>
#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences/explicit/representation/io/structural_model_output.hpp>
#include <geode/geosciences_io/model/internal/ml_input.hpp>
#include <geode/geosciences_io/model/internal/ml_output_structural_model.hpp>

void check_model( const geode::StructuralModel& model,
    geode::index_t nb_corners,
    geode::index_t nb_lines,
    geode::index_t nb_surfaces,
    geode::index_t nb_blocks,
    geode::index_t nb_faults,
    geode::index_t nb_horizons,
    geode::index_t nb_model_boundaries )
{
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_corners() == nb_corners, "Number of Corners in the loaded "
                                          "StructuralModel is not correct" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_lines() == nb_lines, "Number of Lines in the loaded "
                                      "StructuralModel is not correct" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_surfaces() == nb_surfaces,
        "Number of Surfaces in the loaded StructuralModel is not "
        "correct" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_blocks() == nb_blocks, "Number of Blocks in the loaded "
                                        "StructuralModel is not correct" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_faults() == nb_faults, "Number of Faults in the loaded "
                                        "StructuralModel is not correct" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_horizons() == nb_horizons, "Number of Horizons in the loaded "
                                            "StructuralModel is not correct" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_model_boundaries() == nb_model_boundaries,
        "Number of ModelBoundary in the loaded StructuralModel is "
        "not correct" );
}

void test_modelA4()
{
    // Load structural model
    auto model = geode::load_structural_model( absl::StrCat( geode::DATA_PATH,
        "/modelA4.", geode::internal::MLInput::extension() ) );
    check_model( model, 52, 98, 55, 8, 2, 3, 6 );

    geode::index_t nb_block_internals{ 0 };
    for( const auto& block : model.blocks() )
    {
        const auto nb_internals = model.nb_internals( block.id() );
        if( nb_internals )
        {
            auto token = block.name()->substr( block.name()->size() - 3 );
            geode::OpenGeodeGeosciencesIOModelException::test(
                token == "b_2", "Block name should end by b_2" );
        }
        nb_block_internals += nb_internals;
    }
    geode::OpenGeodeGeosciencesIOModelException::test( nb_block_internals == 4,
        "Number of Block internals in the "
        "loaded StructuralModel is not correct" );

    geode::index_t nb_surface_internals{ 0 };
    for( const auto& surface : model.surfaces() )
    {
        const auto nb_internals = model.nb_internals( surface.id() );
        if( nb_internals )
        {
            for( const auto& collection : model.collections( surface.id() ) )
            {
                const auto& name = model.model_boundary( collection.id ).name();
                geode::OpenGeodeGeosciencesIOModelException::test(
                    name == "voi_top_boundary" || name == "voi_bottom_boundary",
                    "ModelBoundary name is not correct" );
            }
        }
        nb_surface_internals += model.nb_internals( surface.id() );
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        nb_surface_internals == 2,
        "Number of Surface internals in the loaded StructuralModel "
        "is not correct" );

    // Save structural model
    geode::save_structural_model(
        model, absl::StrCat( "modelA4.", model.native_extension() ) );

    geode::save_structural_model( model, "modelA4_saved.ml" );
    auto reload = geode::load_structural_model( "modelA4_saved.ml" );
    check_model( reload, 52, 98, 55, 8, 2, 3, 6 );

    /// Test save from a brep
    geode::save_brep( model, "modelA4_saved.og_brep" );
    auto reload_brep = geode::load_brep( "modelA4_saved.og_brep" );
    geode::save_brep( reload_brep, "modelA4_saved_from_brep.ml" );
    // Unclassified Surfaces share their names: TSURF names made unique
    auto reload_from_brep =
        geode::load_structural_model( "modelA4_saved_from_brep.ml" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        reload_from_brep.nb_surfaces() == 55,
        "Number of Surfaces in the StructuralModel reloaded from BRep "
        "export is not correct" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        reload_from_brep.nb_model_boundaries() == 6,
        "Number of ModelBoundaries in the StructuralModel reloaded from BRep "
        "export is not correct" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        reload_from_brep.nb_blocks() == 8,
        "Number of Blocks in the StructuralModel reloaded from BRep "
        "export is not correct" );
}

void test_shared_surfaces_and_contact_types()
{
    auto model = geode::load_structural_model( absl::StrCat( geode::DATA_PATH,
        "/modelA4.", geode::internal::MLInput::extension() ) );
    geode::StructuralModelBuilder builder{ model };

    // Erosion contact type has no dedicated .ml type
    const auto& eroded = *model.horizons().begin();
    const auto eroded_name = eroded.name().value();
    builder.set_horizon_contact_type(
        eroded, geode::Horizon3D::CONTACT_TYPE::erosion );

    // Topographic Horizon made only of ModelBoundary Surfaces
    const auto& top_boundary = *model.model_boundaries().begin();
    const auto& topography_id =
        builder.add_horizon( geode::Horizon3D::CONTACT_TYPE::topography );
    const auto& topography = model.horizon( topography_id );
    builder.set_horizon_name( topography, "topography" );
    for( const auto& surface : model.model_boundary_items( top_boundary ) )
    {
        builder.add_surface_in_horizon( surface, topography );
    }

    geode::save_structural_model( model, "modelA4_shared_surfaces.ml" );
    auto reload = geode::load_structural_model( "modelA4_shared_surfaces.ml" );
    // The topographic Horizon has no Surface of its own: not exported
    check_model( reload, 52, 98, 55, 8, 2, 3, 6 );
    bool eroded_found{ false };
    for( const auto& horizon : reload.horizons() )
    {
        if( horizon.name() == eroded_name )
        {
            eroded_found = true;
            geode::OpenGeodeGeosciencesIOModelException::test(
                horizon.contact_type()
                    == geode::Horizon3D::CONTACT_TYPE::discontinuity,
                "Erosion Horizon should be reloaded as a discontinuity" );
        }
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        eroded_found, "Erosion Horizon should be found after reload" );
}

void test_shared_model_boundary_surfaces()
{
    auto model = geode::load_structural_model( absl::StrCat( geode::DATA_PATH,
        "/modelA4.", geode::internal::MLInput::extension() ) );
    geode::StructuralModelBuilder builder{ model };

    // One Surface of the second ModelBoundary also put in the first one
    auto boundary_it = model.model_boundaries().begin();
    const auto& first_boundary = *boundary_it;
    ++boundary_it;
    const auto& second_boundary = *boundary_it;
    const auto& shared_surface =
        *model.model_boundary_items( second_boundary ).begin();
    builder.add_surface_in_model_boundary( shared_surface, first_boundary );

    // ModelBoundary made only of Surfaces of the first one
    const auto& duplicate_id = builder.add_model_boundary();
    const auto& duplicate = model.model_boundary( duplicate_id );
    builder.set_model_boundary_name( duplicate, "duplicate_boundary" );
    for( const auto& surface : model.model_boundary_items( first_boundary ) )
    {
        builder.add_surface_in_model_boundary( surface, duplicate );
    }

    geode::save_structural_model(
        model, "modelA4_shared_boundary_surfaces.ml" );
    auto reload =
        geode::load_structural_model( "modelA4_shared_boundary_surfaces.ml" );
    // The duplicate ModelBoundary has no Surface of its own: not exported
    check_model( reload, 52, 98, 55, 8, 2, 3, 6 );

    // Same ModelBoundaries exported through the BRep writer
    geode::save_brep( model, "modelA4_shared_boundary_surfaces_brep.ml" );
}

void test_duplicate_tsurf_names()
{
    auto model = geode::load_structural_model( absl::StrCat( geode::DATA_PATH,
        "/modelA4.", geode::internal::MLInput::extension() ) );
    geode::StructuralModelBuilder builder{ model };

    // Same name, and a name only differing by a space
    auto horizon_it = model.horizons().begin();
    builder.set_horizon_name( *horizon_it, "same_name" );
    ++horizon_it;
    builder.set_horizon_name( *horizon_it, "same_name" );
    ++horizon_it;
    builder.set_horizon_name( *horizon_it, "same name" );

    geode::save_structural_model( model, "modelA4_duplicate_names.ml" );
    auto reload = geode::load_structural_model( "modelA4_duplicate_names.ml" );
    check_model( reload, 52, 98, 55, 8, 2, 3, 6 );
    absl::flat_hash_set< std::string > names;
    for( const auto& horizon : reload.horizons() )
    {
        names.insert( horizon.name().value() );
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        names.size() == 3, "Reloaded Horizon names should be unique" );
}

void write_file( std::string_view filename, std::string_view content )
{
    std::ofstream file{ geode::to_string( filename ) };
    file << content;
}

// Tetrahedron: 4 Corners, 6 Lines, 4 Surfaces, 1 Block, 1 ModelBoundary
std::string tetra_ml()
{
    return "GOCAD Model3d 1\n"
           "HEADER {\n"
           "name: tetra\n"
           "}\n"
           "TSURF bnd\n"
           "TFACE 1 boundary bnd\n"
           "  0 0 0\n"
           "  1 0 0\n"
           "  0 1 0\n"
           "TFACE 2 boundary bnd\n"
           "  0 0 0\n"
           "  1 0 0\n"
           "  0 0 1\n"
           "TFACE 3 boundary bnd\n"
           "  0 0 0\n"
           "  0 1 0\n"
           "  0 0 1\n"
           "TFACE 4 boundary bnd\n"
           "  1 0 0\n"
           "  0 1 0\n"
           "  0 0 1\n"
           "REGION 5  Universe \n"
           "  +1  +2  +3  +4  0\n"
           "REGION 6 block\n"
           "  -1  -2  -3  -4  0\n"
           "END\n"
           "GOCAD TSurf 1\n"
           "HEADER {\n"
           "name: bnd\n"
           "}\n"
           "TFACE\n"
           "VRTX 1 0 0 0\n"
           "VRTX 2 1 0 0\n"
           "VRTX 3 0 1 0\n"
           "TRGL 1 2 3\n"
           "TFACE\n"
           "VRTX 4 0 0 0\n"
           "VRTX 5 1 0 0\n"
           "VRTX 6 0 0 1\n"
           "TRGL 4 5 6\n"
           "TFACE\n"
           "VRTX 7 0 0 0\n"
           "VRTX 8 0 1 0\n"
           "VRTX 9 0 0 1\n"
           "TRGL 7 8 9\n"
           "TFACE\n"
           "VRTX 10 1 0 0\n"
           "VRTX 11 0 1 0\n"
           "VRTX 12 0 0 1\n"
           "TRGL 10 11 12\n"
           "BSTONE 1\nBSTONE 2\nBSTONE 3\nBSTONE 4\nBSTONE 5\nBSTONE 6\n"
           "BSTONE 7\nBSTONE 8\nBSTONE 9\nBSTONE 10\nBSTONE 11\nBSTONE 12\n"
           "BORDER 13 1 2\nBORDER 14 2 3\nBORDER 15 3 1\n"
           "BORDER 16 4 5\nBORDER 17 5 6\nBORDER 18 6 4\n"
           "BORDER 19 7 8\nBORDER 20 8 9\nBORDER 21 9 7\n"
           "BORDER 22 10 11\nBORDER 23 11 12\nBORDER 24 12 10\n"
           "END\n";
}

std::string replace( std::string_view content,
    std::string_view old_text,
    std::string_view new_text )
{
    const auto position = content.find( old_text );
    geode::OpenGeodeGeosciencesIOModelException::test(
        position != std::string_view::npos,
        "Text to replace not found: ", old_text );
    return absl::StrCat( content.substr( 0, position ), new_text,
        content.substr( position + old_text.size() ) );
}

geode::StructuralModel load_ml(
    std::string_view filename, std::string_view content )
{
    write_file( filename, content );
    return geode::load_structural_model( filename );
}

void check_load_throws( std::string_view filename, std::string_view content )
{
    write_file( filename, content );
    try
    {
        geode::load_structural_model( filename );
    }
    catch( const geode::OpenGeodeException& exception )
    {
        geode::Logger::info( "Expected exception: ", exception.what() );
        return;
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        false, "Loading ", filename, " should throw an OpenGeodeException" );
}

void test_tetra()
{
    const auto model = load_ml( "tetra.ml", tetra_ml() );
    check_model( model, 4, 6, 4, 1, 0, 0, 1 );
}

void test_truncated_regions()
{
    const auto tetra = tetra_ml();
    const auto end_of_model3d = tetra.find( "END\nGOCAD TSurf" );
    // REGION boundary lists without their 0 terminator, up to end of file
    check_load_throws( "truncated_region.ml",
        replace( tetra.substr( 0, end_of_model3d ), "-4  0", "-4" ) );
    check_load_throws( "truncated_universe.ml",
        replace( tetra.substr( 0, tetra.find( "REGION 6" ) ), "+4  0", "+4" ) );
    check_load_throws( "truncated_layer.ml",
        absl::StrCat( tetra.substr( 0, end_of_model3d ), "LAYER layer\n  6" ) );
    check_load_throws( "truncated_fault_block.ml",
        absl::StrCat(
            tetra.substr( 0, end_of_model3d ), "FAULT_BLOCK fb\n  6" ) );
}

void test_empty_lines()
{
    auto content = replace( tetra_ml(), "TSURF bnd\n", "TSURF bnd\n\n" );
    content = replace( content, "REGION 6", "   \nREGION 6" );
    const auto model = load_ml( "empty_lines.ml", content );
    check_model( model, 4, 6, 4, 1, 0, 0, 1 );
}

void test_end_prefixed_keyword()
{
    // Not read by the Model3d section reader, which must not stop there
    const auto content = replace( tetra_ml(), "TSURF bnd\n",
        "TSURF bnd\n"
        "END_UNKNOWN_KEYWORD\n" );
    const auto model = load_ml( "end_prefixed_keyword.ml", content );
    check_model( model, 4, 6, 4, 1, 0, 0, 1 );
}

void test_invalid_tfaces()
{
    check_load_throws( "short_tface.ml",
        replace( tetra_ml(), "TFACE 2 boundary bnd", "TFACE 2 boundary" ) );
    check_load_throws( "undeclared_tface_tsurf.ml",
        replace(
            tetra_ml(), "TFACE 2 boundary bnd", "TFACE 2 boundary other" ) );
}

void test_region_ids()
{
    auto content = tetra_ml();
    for( const auto tface : { "1", "2", "3", "4" } )
    {
        content = replace( content, absl::StrCat( "TFACE ", tface, " " ),
            absl::StrCat( "TFACE 1", tface, " " ) );
    }
    // Universe not first, unknown TFACE and REGION ids
    content = replace( content,
        "REGION 5  Universe \n"
        "  +1  +2  +3  +4  0\n"
        "REGION 6 block\n"
        "  -1  -2  -3  -4  0\n",
        "REGION 22 block\n"
        "  -11  -12  -13  -14  -9  0\n"
        "REGION 21  Universe \n"
        "  +11  +12  +13  +14  +7  0\n"
        "LAYER layer\n"
        "  22  21  0\n"
        "FAULT_BLOCK fault_block\n"
        "  22  30  0\n" );
    const auto model = load_ml( "region_ids.ml", content );
    check_model( model, 4, 6, 4, 1, 0, 0, 1 );
    const auto& block = *model.blocks().begin();
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_boundaries( block.id() ) == 4,
        "Block should have 4 boundary Surfaces" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_stratigraphic_units() == 1
            && model.nb_items( ( *model.stratigraphic_units().begin() ).id() )
                   == 1,
        "StratigraphicUnit should contain the Block" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_fault_blocks() == 1
            && model.nb_items( ( *model.fault_blocks().begin() ).id() ) == 1,
        "FaultBlock should contain the Block" );

    check_load_throws(
        "duplicate_tface_ids.ml", replace( tetra_ml(), "TFACE 2", "TFACE 1" ) );
}

void test_small_regions()
{
    const auto one_surface = load_ml( "one_surface_region.ml",
        replace( tetra_ml(), "  -1  -2  -3  -4  0", "  -1  0" ) );
    geode::OpenGeodeGeosciencesIOModelException::test(
        one_surface.nb_boundaries( ( *one_surface.blocks().begin() ).id() )
            == 1,
        "Block should have 1 boundary Surface" );
    const auto no_surface = load_ml( "no_surface_region.ml",
        replace( tetra_ml(), "  -1  -2  -3  -4  0", "  0" ) );
    geode::OpenGeodeGeosciencesIOModelException::test(
        no_surface.nb_boundaries( ( *no_surface.blocks().begin() ).id() ) == 0,
        "Block should have no boundary Surface" );
    const auto internal = load_ml(
        "internal_region.ml", replace( tetra_ml(), "  -1  -2  -3  -4  0",
                                  "  -1  -2  -3  -4  +4  0" ) );
    const auto& block = *internal.blocks().begin();
    geode::OpenGeodeGeosciencesIOModelException::test(
        internal.nb_boundaries( block.id() ) == 3
            && internal.nb_internals( block.id() ) == 1,
        "Block should have 3 boundary Surfaces and 1 internal Surface" );
}

// Square bordered by a single closed Line
std::string square_ml()
{
    return "GOCAD Model3d 1\n"
           "HEADER {\n"
           "name: square\n"
           "}\n"
           "TSURF square\n"
           "TFACE 1 boundary square\n"
           "  0 0 0\n"
           "  1 0 0\n"
           "  1 1 0\n"
           "REGION 2  Universe \n"
           "  +1  0\n"
           "END\n"
           "GOCAD TSurf 1\n"
           "HEADER {\n"
           "name: square\n"
           "}\n"
           "TFACE\n"
           "VRTX 1 0 0 0\n"
           "VRTX 2 1 0 0\n"
           "VRTX 3 1 1 0\n"
           "VRTX 4 0 1 0\n"
           "TRGL 1 2 3\n"
           "TRGL 1 3 4\n"
           "BSTONE 1\n"
           "BORDER 5 1 2\n"
           "END\n";
}

void test_single_line_surface()
{
    const auto model = load_ml( "single_line_surface.ml", square_ml() );
    check_model( model, 1, 1, 1, 0, 0, 0, 1 );
    const auto& surface = *model.surfaces().begin();
    const auto& line = *model.lines().begin();
    geode::OpenGeodeGeosciencesIOModelException::test(
        line.mesh().nb_vertices() == 5, "Closed Line should have 5 vertices" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.is_boundary( line, surface )
            && !model.is_internal( line, surface ),
        "Single Line should be a boundary of the Surface" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        surface.mesh().nb_vertices() == 4,
        "Surface should not have been split along its Line" );
}

void test_missing_tsurf_section()
{
    const auto content =
        replace( tetra_ml(), "TSURF bnd\n", "TSURF bnd\nTSURF other\n" );
    check_load_throws( "missing_tsurf_section_load.ml",
        replace( content, "REGION 5", "TFACE 7 top other\nREGION 5" ) );
}

using TetraPoints = std::array< std::array< double, 3 >, 4 >;

std::string tsurf_section( std::string_view name, const TetraPoints& points )
{
    std::string section =
        absl::StrCat( "GOCAD TSurf 1\nHEADER {\nname: ", name, "\n}\n" );
    constexpr std::array< std::array< geode::index_t, 3 >, 4 > faces{ {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 },
    } };
    geode::index_t vertex{ 1 };
    for( const auto& face : faces )
    {
        absl::StrAppend( &section, "TFACE\n" );
        for( const auto v : face )
        {
            absl::StrAppend( &section, "VRTX ", vertex++, " ", points[v][0],
                " ", points[v][1], " ", points[v][2], "\n" );
        }
        absl::StrAppend( &section, "TRGL ", vertex - 3, " ", vertex - 2, " ",
            vertex - 1, "\n" );
    }
    for( const auto v : geode::Range{ 1, 13 } )
    {
        absl::StrAppend( &section, "BSTONE ", v, "\n" );
    }
    for( const auto f : geode::Range{ 4 } )
    {
        const auto first = 3 * f + 1;
        absl::StrAppend( &section, "BORDER ", 13 + 3 * f, " ", first, " ",
            first + 1, "\n", "BORDER ", 14 + 3 * f, " ", first + 1, " ",
            first + 2, "\n", "BORDER ", 15 + 3 * f, " ", first + 2, " ", first,
            "\n" );
    }
    absl::StrAppend( &section, "END\n" );
    return section;
}

// Small tetrahedron Block (inner) inside a tetrahedron Block (outer)
std::string nested_tetra_ml( bool inner_section_first )
{
    std::string content = "GOCAD Model3d 1\n"
                          "HEADER {\n"
                          "name: nested\n"
                          "}\n"
                          "TSURF outer\n"
                          "TSURF inner\n";
    for( const auto f : geode::Range{ 1, 5 } )
    {
        absl::StrAppend( &content, "TFACE ", f,
            " boundary outer\n  0 0 0\n  0 0 0\n  0 0 0\n" );
    }
    for( const auto f : geode::Range{ 5, 9 } )
    {
        absl::StrAppend(
            &content, "TFACE ", f, " top inner\n  0 0 0\n  0 0 0\n  0 0 0\n" );
    }
    absl::StrAppend( &content, "REGION 9  Universe \n  +1  +2  +3  +4  0\n",
        "REGION 10 outer_block\n  -1  -2  -3  -4  +5\n  +6  +7  +8  0\n",
        "REGION 11 inner_block\n  -5  -6  -7  -8  0\n", "END\n" );
    const auto outer = tsurf_section(
        "outer", { { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } } );
    const auto inner =
        tsurf_section( "inner", { { { 0.2, 0.2, 0.2 }, { 0.3, 0.2, 0.2 },
                                    { 0.2, 0.3, 0.2 }, { 0.2, 0.2, 0.3 } } } );
    if( inner_section_first )
    {
        absl::StrAppend( &content, inner, outer );
    }
    else
    {
        absl::StrAppend( &content, outer, inner );
    }
    return content;
}

void test_tsurf_section_names()
{
    const auto model =
        load_ml( "tsurf_section_names.ml", nested_tetra_ml( true ) );
    check_model( model, 8, 12, 8, 2, 0, 1, 1 );
    for( const auto& boundary : model.model_boundaries() )
    {
        for( const auto& surface : model.model_boundary_items( boundary ) )
        {
            const auto max = surface.mesh().bounding_box().max();
            geode::OpenGeodeGeosciencesIOModelException::test(
                std::max( { max.value( 0 ), max.value( 1 ), max.value( 2 ) } )
                    > 0.5,
                "ModelBoundary Surfaces should be the outer tetrahedron" );
        }
    }
}

// REGION name -> TFACE id -> signs (several for internal Surfaces)
absl::flat_hash_map< std::string,
    absl::flat_hash_map< geode::index_t, std::vector< bool > > >
    read_region_signs( std::string_view filename )
{
    absl::flat_hash_map< std::string,
        absl::flat_hash_map< geode::index_t, std::vector< bool > > >
        regions;
    std::ifstream file{ geode::to_string( filename ) };
    std::string line;
    while( std::getline( file, line ) )
    {
        const auto tokens = geode::string_split( line );
        if( tokens.empty() || tokens.front() != "REGION" )
        {
            continue;
        }
        auto& region = regions[tokens.back()];
        bool end{ false };
        while( !end && std::getline( file, line ) )
        {
            for( const auto token : geode::string_split( line ) )
            {
                const auto item = geode::string_to_int( token );
                if( item == 0 )
                {
                    end = true;
                    break;
                }
                region[static_cast< geode::index_t >( std::abs( item ) )]
                    .push_back( item > 0 );
            }
        }
    }
    return regions;
}

void test_nested_blocks_output()
{
    const auto model = load_ml( "nested_blocks.ml", nested_tetra_ml( false ) );
    geode::save_structural_model( model, "nested_blocks_saved.ml" );
    const auto reload =
        geode::load_structural_model( "nested_blocks_saved.ml" );
    check_model( reload, 8, 12, 8, 2, 0, 1, 1 );
    const auto regions = read_region_signs( "nested_blocks_saved.ml" );
    const auto& universe = regions.at( "Universe" );
    const auto& outer = regions.at( "outer_block" );
    const auto& inner = regions.at( "inner_block" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        universe.size() == 4 && outer.size() == 8 && inner.size() == 4,
        "Wrong number of Surfaces in exported REGIONs" );
    // Each Surface is seen with opposite signs from its two sides
    for( const auto& [tface, signs] : universe )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            outer.at( tface ).front() != signs.front(),
            "Universe and outer Block should see the outer Surfaces from "
            "opposite sides" );
    }
    for( const auto& [tface, signs] : inner )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            outer.at( tface ).front() != signs.front(),
            "Outer and inner Blocks should see the inclusion Surfaces from "
            "opposite sides" );
    }
}

void test_block_without_boundary_output()
{
    auto model = load_ml( "empty_block.ml", tetra_ml() );
    geode::StructuralModelBuilder builder{ model };
    builder.set_block_name( model.block( builder.add_block() ), "empty_block" );
    geode::save_structural_model( model, "empty_block_saved.ml" );
    const auto reload = geode::load_structural_model( "empty_block_saved.ml" );
    check_model( reload, 4, 6, 4, 2, 0, 0, 1 );
}

void test_empty_line_output()
{
    auto model = load_ml( "empty_line.ml", tetra_ml() );
    geode::StructuralModelBuilder builder{ model };
    builder.add_line();
    geode::save_structural_model( model, "empty_line_saved.ml" );
    const auto reload = geode::load_structural_model( "empty_line_saved.ml" );
    check_model( reload, 4, 6, 4, 1, 0, 0, 1 );
}

void test_surface_collection_output()
{
    auto model = load_ml( "surface_collection.ml", nested_tetra_ml( false ) );
    geode::StructuralModelBuilder builder{ model };
    builder.remove_horizon( *model.horizons().begin() );
    const auto& collection_id = builder.add_surface_collection();
    const auto& collection = model.surface_collection( collection_id );
    for( const auto& surface : model.surfaces() )
    {
        if( model.nb_collections( surface.id() ) == 0 )
        {
            builder.add_surface_in_surface_collection( surface, collection );
        }
    }
    geode::save_structural_model( model, "surface_collection_saved.ml" );
    const auto reload =
        geode::load_structural_model( "surface_collection_saved.ml" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        reload.nb_surfaces() == 8 && reload.nb_blocks() == 2,
        "Surfaces only in a SurfaceCollection should be exported" );
}

// Square with a square hole: 2 closed border Lines
std::string annulus_ml()
{
    return "GOCAD Model3d 1\n"
           "HEADER {\n"
           "name: annulus\n"
           "}\n"
           "TSURF annulus\n"
           "TFACE 1 boundary annulus\n"
           "  0 0 0\n"
           "  3 0 0\n"
           "  2 1 1\n"
           "REGION 2  Universe \n"
           "  +1  0\n"
           "END\n"
           "GOCAD TSurf 1\n"
           "HEADER {\n"
           "name: annulus\n"
           "}\n"
           "TFACE\n"
           "VRTX 1 0 0 0\n"
           "VRTX 2 3 0 0\n"
           "VRTX 3 3 3 0\n"
           "VRTX 4 0 3 0\n"
           "VRTX 5 1 1 1\n"
           "VRTX 6 2 1 1\n"
           "VRTX 7 2 2 1\n"
           "VRTX 8 1 2 1\n"
           "TRGL 1 2 6\nTRGL 1 6 5\nTRGL 2 3 7\nTRGL 2 7 6\n"
           "TRGL 3 4 8\nTRGL 3 8 7\nTRGL 4 1 5\nTRGL 4 5 8\n"
           "BSTONE 1\n"
           "BSTONE 5\n"
           "BORDER 9 1 2\n"
           "BORDER 10 5 6\n"
           "END\n";
}

void check_annulus( const geode::StructuralModel& model )
{
    check_model( model, 2, 2, 1, 0, 0, 0, 1 );
    const auto& surface = *model.surfaces().begin();
    for( const auto& line : model.lines() )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            model.is_boundary( line, surface ),
            "Both Lines should be boundaries of the Surface" );
    }
}

void test_several_border_loops_output()
{
    const auto model = load_ml( "annulus.ml", annulus_ml() );
    check_annulus( model );
    geode::save_structural_model( model, "annulus_saved.ml" );
    check_annulus( geode::load_structural_model( "annulus_saved.ml" ) );
}

void test_empty_surface_output()
{
    auto model = load_ml( "empty_surface.ml", tetra_ml() );
    geode::StructuralModelBuilder builder{ model };
    const auto& surface = model.surface( builder.add_surface() );
    builder.add_surface_in_model_boundary(
        surface, *model.model_boundaries().begin() );
    builder.add_surface_block_boundary_relationship(
        surface, *model.blocks().begin() );
    // Unclassified empty Surface
    builder.add_surface();
    geode::save_structural_model( model, "empty_surface_saved.ml" );
    const auto reload =
        geode::load_structural_model( "empty_surface_saved.ml" );
    check_model( reload, 4, 6, 4, 1, 0, 0, 1 );
}

void test_coordinates_precision()
{
    constexpr auto COORDINATE = "0.12345678901234567";
    auto content = tetra_ml();
    for( const auto vertex : { "VRTX 2 1 0 ", "VRTX 5 1 0 ", "VRTX 10 1 0 " } )
    {
        content = replace( content, absl::StrCat( vertex, "0\n" ),
            absl::StrCat( vertex, COORDINATE, "\n" ) );
    }
    const auto model = load_ml( "precision.ml", content );
    geode::save_structural_model( model, "precision_saved.ml" );
    const auto reload = geode::load_structural_model( "precision_saved.ml" );
    const auto expected = geode::string_to_double( COORDINATE );
    geode::index_t nb_found{ 0 };
    for( const auto& surface : reload.surfaces() )
    {
        const auto& mesh = surface.mesh();
        for( const auto v : geode::Range{ mesh.nb_vertices() } )
        {
            if( mesh.point( v ).value( 0 ) == 1 )
            {
                nb_found++;
                geode::OpenGeodeGeosciencesIOModelException::test(
                    mesh.point( v ).value( 2 ) == expected,
                    "Coordinates should be exported without precision loss" );
            }
        }
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        nb_found == 3, "Modified vertices not found" );
}

void test_border_without_bstone()
{
    // BORDER loop without any BSTONE
    const auto no_bstone = load_ml(
        "border_without_bstone.ml", replace( square_ml(), "BSTONE 1\n", "" ) );
    check_model( no_bstone, 0, 0, 1, 0, 0, 0, 1 );
}

void test_bstone_without_line()
{
    // BORDER ignored: the BSTONE ends no Line
    const auto model = load_ml( "bstone_without_line.ml",
        replace( square_ml(), "BORDER 5 1 2", "BORDER 5 2 3" ) );
    check_model( model, 1, 0, 1, 0, 0, 0, 1 );
}

void test_is_loadable()
{
    const auto model_a4 = absl::StrCat(
        geode::DATA_PATH, "/modelA4.", geode::internal::MLInput::extension() );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_structural_model_loadable( model_a4 ).value() == 1,
        "modelA4 should be fully loadable" );

    write_file( "undeclared_tsurf.ml", "GOCAD Model3d 1\n"
                                       "TSURF A\n"
                                       "TFACE 1 top B\n"
                                       "END\n"
                                       "GOCAD TSurf 1\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_structural_model_loadable( "undeclared_tsurf.ml" ).value()
            == 0,
        "A TFACE referring to an undeclared TSURF should not be loadable" );

    write_file( "unknown_type.ml", "GOCAD Model3d 1\n"
                                   "TSURF A\n"
                                   "TFACE 1 erosion A\n"
                                   "END\n"
                                   "GOCAD TSurf 1\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_structural_model_loadable( "unknown_type.ml" ).value() == 0,
        "A TFACE with an unknown type should not be loadable" );

    write_file( "crs_unknown_type.ml", "GOCAD Model3d 1\n"
                                       "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
                                       "NAME Default\n"
                                       "END_ORIGINAL_COORDINATE_SYSTEM\n"
                                       "TSURF A\n"
                                       "TFACE 1 erosion A\n"
                                       "END\n"
                                       "GOCAD TSurf 1\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_structural_model_loadable( "crs_unknown_type.ml" ).value()
            == 0,
        "The Model3d section should be checked after a coordinate system "
        "block" );

    write_file( "missing_tsurf_section.ml", "GOCAD Model3d 1\n"
                                            "TSURF A\n"
                                            "TFACE 1 top A\n"
                                            "END\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_structural_model_loadable( "missing_tsurf_section.ml" )
                .value()
            == 0,
        "A Model3d without TSurf section should not be loadable" );

    write_file( "empty_tsurf.ml", "GOCAD Model3d 1\n"
                                  "TSURF A\n"
                                  "TSURF B\n"
                                  "TFACE 1 top A\n"
                                  "END\n"
                                  "GOCAD TSurf 1\n"
                                  "GOCAD TSurf 1\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_structural_model_loadable( "empty_tsurf.ml" ).value() == 1,
        "A TSURF without TFACE is ignored by the reader: the file should be "
        "loadable" );

    write_file( "duplicate_tsurf.ml", "GOCAD Model3d 1\n"
                                      "TSURF A\n"
                                      "TSURF A\n"
                                      "TFACE 1 top A\n"
                                      "END\n"
                                      "GOCAD TSurf 1\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_structural_model_loadable( "duplicate_tsurf.ml" ).value()
            == 0,
        "Duplicate TSURF names should not be loadable" );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        test_modelA4();
        test_shared_surfaces_and_contact_types();
        test_shared_model_boundary_surfaces();
        test_duplicate_tsurf_names();
        test_is_loadable();
        test_tetra();
        test_truncated_regions();
        test_empty_lines();
        test_end_prefixed_keyword();
        test_invalid_tfaces();
        test_region_ids();
        test_small_regions();
        test_single_line_surface();
        test_missing_tsurf_section();
        test_tsurf_section_names();
        test_nested_blocks_output();
        test_block_without_boundary_output();
        test_empty_line_output();
        test_surface_collection_output();
        test_several_border_loops_output();
        test_empty_surface_output();
        test_coordinates_precision();
        test_border_without_bstone();
        test_bstone_without_line();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
