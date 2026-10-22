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

#include <filesystem>
#include <fstream>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/geosciences_io/mesh/internal/grdecl_input.hpp>
#include <geode/mesh/core/geode/geode_hybrid_solid.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/io/hybrid_solid_input.hpp>
#include <geode/mesh/io/hybrid_solid_output.hpp>

void check_solid( const geode::HybridSolid3D& solid,
    geode::index_t nb_polyhedra,
    geode::index_t nb_vertices )
{
    geode::OpenGeodeGeosciencesIOMeshException::test(
        solid.nb_polyhedra() == nb_polyhedra,
        "Number of polyhedra in the GrdeclHybridSolid is not correct" );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        solid.nb_vertices() == nb_vertices,
        "Number of vertices in the GrdeclHybridSolid is not correct " );
    for( const auto polyhedra : geode::Range{ solid.nb_polyhedra() } )
    {
        geode::OpenGeodeGeosciencesIOMeshException::test(
            solid.polyhedron_volume( polyhedra ) > 0.,
            "Found negative volume of polyhedron is not correct" );
    }
}

void check_file( std::string_view filename,
    geode::index_t nb_polyhedra,
    geode::index_t nb_vertices )
{
    const auto solid = geode::load_hybrid_solid< 3 >( filename );
    check_solid( *solid, nb_polyhedra, nb_vertices );
}

std::string write_grdecl( std::string_view name, std::string_view content )
{
    auto filename = absl::StrCat( name, ".grdecl" );
    std::ofstream file{ filename };
    file << content;
    return filename;
}

void check_load_throws( std::string_view filename )
{
    bool thrown{ false };
    try
    {
        const auto solid = geode::load_hybrid_solid< 3 >( filename );
    }
    catch( const geode::OpenGeodeException& )
    {
        thrown = true;
    }
    geode::OpenGeodeGeosciencesIOMeshException::test(
        thrown, "Loading ", filename, " should have thrown" );
}

constexpr std::string_view SPECGRID_2X1X1{ "SPECGRID\n2 1 1 1 F\n/\n" };
constexpr std::string_view COORD_2X1X1{
    "COORD\n"
    "0 0 1000 0 0 1100\n10 0 1000 10 0 1100\n20 0 1000 20 0 1100\n"
    "0 10 1000 0 10 1100\n10 10 1000 10 10 1100\n20 10 1000 20 10 1100\n"
    "/\n"
};
constexpr std::string_view ZCORN_2X1X1{
    "ZCORN\n"
    "1000 1000 1000 1000 1000 1000 1000 1000\n"
    "1010 1010 1010 1010 1010 1010 1010 1010\n"
    "/\n"
};

void check_grid_2x1x1( std::string_view name, std::string_view content )
{
    check_file( write_grdecl( name, content ), 2, 12 );
}

void test_specgrid_variants()
{
    check_grid_2x1x1( "grdecl_specgrid_no_flag",
        absl::StrCat(
            "SPECGRID\n-- From Petrel\n2 1 1 /\n", COORD_2X1X1, ZCORN_2X1X1 ) );
    check_grid_2x1x1( "grdecl_specgrid_true_flag",
        absl::StrCat( "SPECGRID\n\n2 1 1 1 T /\n", COORD_2X1X1, ZCORN_2X1X1 ) );
    check_load_throws( write_grdecl( "grdecl_specgrid_short",
        absl::StrCat( "SPECGRID\n2 1 /\n", COORD_2X1X1, ZCORN_2X1X1 ) ) );
    check_load_throws( write_grdecl( "grdecl_specgrid_eof", "SPECGRID\n" ) );
}

void test_truncated_data()
{
    check_load_throws( write_grdecl( "grdecl_truncated_coord",
        absl::StrCat( SPECGRID_2X1X1,
            "COORD\n0 0 1000 0 0 1100\n10 0 1000 10 0 1100\n" ) ) );
    check_load_throws( write_grdecl(
        "grdecl_truncated_zcorn", absl::StrCat( SPECGRID_2X1X1, COORD_2X1X1,
                                      "ZCORN\n1000 1000 1000 1000\n" ) ) );
    check_grid_2x1x1(
        "grdecl_unterminated_zcorn", absl::StrCat( SPECGRID_2X1X1, COORD_2X1X1,
                                         "ZCORN\n8*1000\n8*1010\n" ) );
}

void test_wrong_number_of_values()
{
    check_grid_2x1x1( "grdecl_extra_zcorn",
        absl::StrCat( SPECGRID_2X1X1, COORD_2X1X1,
            "ZCORN\n8*1000 8*1010 1020 1020 1020\n/\n" ) );
    check_grid_2x1x1( "grdecl_extra_coord",
        absl::StrCat( SPECGRID_2X1X1,
            "COORD\n"
            "0 0 1000 0 0 1100\n10 0 1000 10 0 1100\n20 0 1000 20 0 1100\n"
            "0 10 1000 0 10 1100\n10 10 1000 10 10 1100\n"
            "20 10 1000 20 10 1100 30 10 1000 30 10 1100\n/\n",
            ZCORN_2X1X1 ) );
    check_load_throws( write_grdecl(
        "grdecl_missing_zcorn", absl::StrCat( SPECGRID_2X1X1, COORD_2X1X1,
                                    "ZCORN\n8*1000 1010\n/\n" ) ) );
    check_load_throws( write_grdecl( "grdecl_huge_specgrid",
        absl::StrCat(
            "SPECGRID\n1024 1024 1024 1 F /\n", COORD_2X1X1, ZCORN_2X1X1 ) ) );
    check_load_throws( write_grdecl(
        "grdecl_huge_repeat", absl::StrCat( SPECGRID_2X1X1, COORD_2X1X1,
                                  "ZCORN\n99999999999999999999*1000\n/\n" ) ) );
}

void test_eclipse_syntax()
{
    check_grid_2x1x1( "grdecl_eclipse_syntax",
        absl::StrCat( "-- Header comment\n", SPECGRID_2X1X1,
            "COORD 0 0 1000 0 0 1100 10 0 -- trailing comment\n"
            "  1000 10 0 1100 20 0 1000 20 0 1100 0 10 1000\n"
            "0 10 1100 10 10\n1000 10 10 1100\n20 10 1000\n20 10 1100\n"
            " /\n"
            "ZCORN\n"
            "-- top\n"
            "4*1000 2*1000 1000 1000\n"
            "8*1010/\n" ) );
}

void test_keyword_order()
{
    check_grid_2x1x1( "grdecl_keyword_order",
        absl::StrCat( "COORDSYS\n1 1 /\n", SPECGRID_2X1X1, ZCORN_2X1X1,
            "  MAPAXES\n0 1 0 0 1 0 /\n", COORD_2X1X1 ) );
    check_load_throws( write_grdecl(
        "grdecl_no_zcorn", absl::StrCat( SPECGRID_2X1X1, COORD_2X1X1 ) ) );
}

void test_include()
{
    {
        std::ofstream grid{ "grdecl include grid.inc" };
        grid << "-- Grid geometry\n" << SPECGRID_2X1X1 << COORD_2X1X1;
        std::ofstream zcorn{ "grdecl_include_zcorn.inc" };
        zcorn << ZCORN_2X1X1;
    }
    const auto zcorn_path =
        std::filesystem::absolute( "grdecl_include_zcorn.inc" ).string();
    check_grid_2x1x1( "grdecl_include",
        absl::StrCat( "INCLUDE\n'grdecl include grid.inc' /\n", "INCLUDE\n  \"",
            zcorn_path, "\"\n/\n" ) );
    check_load_throws( write_grdecl( "grdecl_include_missing",
        absl::StrCat( SPECGRID_2X1X1, COORD_2X1X1,
            "INCLUDE\n'grdecl_does_not_exist.inc' /\n" ) ) );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::trace );
        check_file( absl::StrCat( geode::DATA_PATH, "Simple20x20x5_Fault.",
                        geode::internal::GRDECLInput::extension() ),

            20 * 20 * 5, 21 * 6 * ( 21 + 1 ) );
        check_file( absl::StrCat( geode::DATA_PATH, "EclipseGridTest.",
                        geode::internal::GRDECLInput::extension() ),

            24, 60 );
        test_specgrid_variants();
        test_truncated_data();
        test_wrong_number_of_values();
        test_eclipse_syntax();
        test_keyword_order();
        test_include();
        geode::Logger::info( "[TEST SUCCESS]" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
