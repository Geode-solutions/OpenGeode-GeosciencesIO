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

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/model_boundary.hpp>
#include <geode/model/mixin/core/surface.hpp>

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
        geode::internal::MLInput{ "crs_unknown_type.ml" }
                .is_loadable()
                .value()
            == 0,
        "The Model3d section should be checked after a coordinate system "
        "block" );

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
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        test_modelA4();
        test_is_loadable();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
