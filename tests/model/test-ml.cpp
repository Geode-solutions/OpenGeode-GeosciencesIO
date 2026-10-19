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

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/model_boundary.hpp>
#include <geode/model/mixin/core/surface.hpp>

#include <geode/geosciences/explicit/mixin/core/horizon.hpp>
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

void test_is_loadable()
{
    const auto model_a4 = absl::StrCat(
        geode::DATA_PATH, "/modelA4.", geode::internal::MLInput::extension() );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::internal::MLInput{ model_a4 }.is_loadable().value() == 1,
        "modelA4 should be fully loadable" );

    write_file( "undeclared_tsurf.ml", "GOCAD Model3d 1\n"
                                       "TSURF A\n"
                                       "TFACE 1 top B\n"
                                       "END\n"
                                       "GOCAD TSurf 1\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::internal::MLInput{ "undeclared_tsurf.ml" }.is_loadable().value()
            == 0,
        "A TFACE referring to an undeclared TSURF should not be loadable" );

    write_file( "unknown_type.ml", "GOCAD Model3d 1\n"
                                   "TSURF A\n"
                                   "TFACE 1 erosion A\n"
                                   "END\n"
                                   "GOCAD TSurf 1\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::internal::MLInput{ "unknown_type.ml" }.is_loadable().value()
            == 0,
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
        geode::internal::MLInput{ "crs_unknown_type.ml" }.is_loadable().value()
            == 0,
        "The Model3d section should be checked after a coordinate system "
        "block" );

    write_file( "missing_tsurf_section.ml", "GOCAD Model3d 1\n"
                                            "TSURF A\n"
                                            "TFACE 1 top A\n"
                                            "END\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::internal::MLInput{ "missing_tsurf_section.ml" }
                .is_loadable()
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
        geode::internal::MLInput{ "empty_tsurf.ml" }.is_loadable().value()
            == 0.5,
        "A TSURF without TFACE should not be counted as loadable" );

    write_file( "duplicate_tsurf.ml", "GOCAD Model3d 1\n"
                                      "TSURF A\n"
                                      "TSURF A\n"
                                      "TFACE 1 top A\n"
                                      "END\n"
                                      "GOCAD TSurf 1\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::internal::MLInput{ "duplicate_tsurf.ml" }.is_loadable().value()
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

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
