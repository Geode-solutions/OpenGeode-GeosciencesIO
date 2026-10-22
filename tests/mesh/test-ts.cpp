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

#include <geode/mesh/core/geode/geode_triangulated_surface.hpp>
#include <geode/mesh/core/triangulated_surface.hpp>
#include <geode/mesh/io/triangulated_surface_input.hpp>
#include <geode/mesh/io/triangulated_surface_output.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>
#include <geode/geosciences_io/mesh/internal/ts_input.hpp>

namespace
{
    std::unique_ptr< geode::TriangulatedSurface3D > load_test_file(
        std::string_view filename, std::string_view content )
    {
        {
            std::ofstream file{ geode::to_string( filename ),
                std::ios::binary };
            file << content;
        }
        return geode::load_triangulated_surface< 3 >( filename );
    }

    void test_blank_lines()
    {
        const auto surface = load_test_file( "test_ts_blank_lines.ts",
            "GOCAD TSurf 1\nHEADER {\nname: blank\n}\n"
            "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n\nNAME Default\nZPOSITIVE\n"
            "END_ORIGINAL_COORDINATE_SYSTEM\nTFACE\n\n# comment\n"
            "VRTX 1 0 0 0\nVRTX 2 1 0 0\n   \nVRTX 3 0 1 0\nTRGL 1 2\n"
            "TRGL 1 2 3\nBSTONE\nBORDER 1 2\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 3,
            "Wrong number of vertices with blank lines" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_polygons() == 1,
            "Wrong number of triangles with blank lines" );
    }

    void test_invalid_indices()
    {
        constexpr auto content =
            "GOCAD TSurf 1\nHEADER {\nname: invalid\n}\nTFACE\n"
            "VRTX 1 0 0 0\nVRTX 2 1 0 0\nVRTX 3 0 1 0\nTRGL 1 2 3\n"
            "TRGL 1 2 9\nTRGL 0 1 2\nTFACE\nVRTX 4 1 1 0\nVRTX 5 2 1 0\n"
            "VRTX 6 1 2 0\nTRGL 4 5 6\nTRGL 3 5 6\nBSTONE 1\nBSTONE 0\n"
            "BSTONE 42\nBORDER 7 1 2\nBORDER 8 3 4\nBORDER 9 1 10\n"
            "BORDER 10 4 5\nEND\n";
        const auto surface =
            load_test_file( "test_ts_invalid_indices.ts", content );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 6,
            "Wrong number of vertices with invalid indices" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_polygons() == 2,
            "Wrong number of triangles with invalid indices" );

        std::ifstream file{ "test_ts_invalid_indices.ts", std::ios::binary };
        const auto tsurf = geode::internal::read_tsurf( file ).value();
        geode::OpenGeodeGeosciencesIOMeshException::test(
            tsurf.tface_triangles_offset.size() == 3
                && tsurf.tface_triangles_offset[1] == 1
                && tsurf.tface_triangles_offset[2] == 2,
            "Wrong TFACE triangle offsets with invalid indices" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            tsurf.triangles[1] == std::array< geode::index_t, 3 >{ 3, 4, 5 },
            "Wrong triangle with invalid indices" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            tsurf.bstones.size() == 1 && tsurf.bstones[0] == 0,
            "Wrong BSTONE with invalid indices" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            tsurf.borders.size() == 2 && tsurf.borders[0].corner_id == 0
                && tsurf.borders[0].next_id == 1
                && tsurf.borders[1].corner_id == 3
                && tsurf.borders[1].next_id == 4,
            "Wrong BORDER with invalid indices" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            tsurf.tface_id( 2 ) == 0 && tsurf.tface_id( 3 ) == 1
                && tsurf.tface_id( 5 ) == 1,
            "Wrong TFACE id" );
    }

    double vertex_value( const geode::SurfaceMesh3D& surface,
        std::string_view name,
        geode::index_t vertex )
    {
        const auto& manager = surface.vertex_attribute_manager();
        const auto ids = manager.attribute_ids_matching_name( name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            ids.has_value() && ids->size() == 1,
            "Attribute should exist once: ", name );
        return manager.find_read_only_attribute< double >( ids->front() )
            ->value( vertex );
    }

    void test_partial_property_header()
    {
        const auto surface = load_test_file( "test_ts_partial_properties.ts",
            "GOCAD TSurf 1\nHEADER {\nname: s1\n}\nPROPERTIES p\n"
            "PROP_LEGAL_RANGES **none** **none**\nNO_DATA_VALUES -99999\n"
            "PROPERTY_CLASSES p\nPROPERTY_KINDS \"Real Number\"\nESIZES 1\n"
            "UNITS none\nTFACE\nPVRTX 1 0 0 1 11\nPVRTX 2 1 0 1 12\n"
            "PVRTX 3 0 1 1 13\nTRGL 1 2 3\nEND\n"
            "GOCAD TSurf 1\nHEADER {\nname: s2\n}\nPROPERTIES a b\n"
            "UNITS m\nNO_DATA_VALUES -1 -2 -3\nESIZES 1\n\nTFACE\n"
            "PVRTX 1 0 0 2 21 31\nPVRTX 2 1 0 2 22 32\nPVRTX 3 0 1 2 23 33\n"
            "TRGL 1 2 3\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 6,
            "Wrong number of vertices with partial property header" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_polygons() == 2,
            "Wrong number of triangles with partial property header" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_value( *surface, "p", 2 ) == 13,
            "Wrong property value with partial property header" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_value( *surface, "b", 5 ) == 33
                && vertex_value( *surface, "b", 0 ) == -2
                && vertex_value( *surface, "p", 4 ) == -99999,
            "Wrong property value in second TSurf" );
    }

    void test_multiple_tsurf_properties()
    {
        const auto surface = load_test_file( "test_ts_multi_properties.ts",
            "GOCAD TSurf 1\nHEADER {\nname: s1\n}\nPROPERTIES p\nESIZES 1\n"
            "TFACE\nPVRTX 1 0 0 1 11\nPVRTX 2 1 0 1 12\nPVRTX 3 0 1 1 13\n"
            "TRGL 1 2 3\nEND\nGOCAD TSurf 1\nHEADER {\nname: s2\n}\n"
            "PROPERTIES p\nESIZES 1\nTFACE\nPVRTX 1 0 0 2 21\n"
            "PVRTX 2 1 0 2 22\nPVRTX 3 0 1 2 23\nTRGL 1 2 3\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_value( *surface, "p", 1 ) == 12
                && vertex_value( *surface, "p", 4 ) == 22,
            "Wrong property values with several TSurf" );
    }

    void test_vertex_lines_with_properties()
    {
        const auto surface = load_test_file( "test_ts_vrtx_atom.ts",
            "GOCAD TSurf 1\nHEADER {\nname: s1\n}\nPROPERTIES p\n"
            "NO_DATA_VALUES -5\nTFACE\nPVRTX 1 0 0 0 11\nVRTX 2 1 0 0\n"
            "PVRTX 3 0 1 0 13\nTRGL 1 2 3\nTFACE\nATOM 4 1\nATOM 5 3\n"
            "PVRTX 6 1 1 0 16\nTRGL 4 5 6\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 6 && surface->nb_polygons() == 2,
            "Wrong surface with VRTX and ATOM" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_value( *surface, "p", 1 ) == -5
                && vertex_value( *surface, "p", 3 ) == 11
                && vertex_value( *surface, "p", 4 ) == 13
                && vertex_value( *surface, "p", 5 ) == 16,
            "Wrong property values with VRTX and ATOM" );
        auto has_thrown = false;
        try
        {
            load_test_file( "test_ts_bad_atom.ts",
                "GOCAD TSurf 1\nHEADER {\nname: s1\n}\nTFACE\n"
                "VRTX 1 0 0 0\nATOM 2 7\nEND\n" );
        }
        catch( const geode::OpenGeodeException& exception )
        {
            geode::Logger::info( "Expected exception: ", exception.what() );
            has_thrown = true;
        }
        geode::OpenGeodeGeosciencesIOMeshException::test(
            has_thrown, "ATOM referencing an unknown vertex should throw" );
    }

    void test_object_header_order()
    {
        const auto surface = load_test_file( "test_ts_header_order.ts",
            "GOCAD TSurf 1\nHEADER {\nname: s1\n}\nGEOLOGICAL_FEATURE f\n"
            "PROPERTY_CLASS_HEADER Z {\nkind: Z\nis_z: on\n}\n"
            "GOCAD_ORIGINAL_COORDINATE_SYSTEM\nNAME Default\nZPOSITIVE Depth\n"
            "END_ORIGINAL_COORDINATE_SYSTEM\nGEOLOGICAL_TYPE top\n"
            "PROPERTIES p\nESIZES 1\nPVRTX 1 0 0 1 11\nPVRTX 2 1 0 1 12\n"
            "PVRTX 3 0 1 1 13\nTRGL 1 2 3\n"
            "GOCAD TSurf 1\nHEADER {\nname: s2\n}\nTFACE\nVRTX 1 0 0 2\n"
            "VRTX 2 1 0 2\nVRTX 3 0 1 2\nTRGL 1 2 3\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->nb_vertices() == 6 && surface->nb_polygons() == 2,
            "Wrong surface with unordered object header" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            surface->point( 0 ).value( 2 ) == -1
                && surface->point( 3 ).value( 2 ) == 2,
            "ZPOSITIVE not applied with unordered object header" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_value( *surface, "p", 1 ) == 12,
            "Wrong property with unordered object header" );
    }
} // namespace

void check_surface( const geode::SurfaceMesh3D& surface,
    geode::index_t nb_vertices,
    geode::index_t nb_polygons,
    std::string name )
{
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface.nb_vertices() == nb_vertices,
        "Number of vertices in the TSurf 3D is not correct" );
    geode::OpenGeodeGeosciencesIOMeshException::test(
        surface.nb_polygons() == nb_polygons,
        "Number of polygons in the TSurf 3D is not correct" );
    geode::OpenGeodeGeosciencesIOMeshException::test( surface.name() == name,
        "Wrong name: ", surface.name().value(), " should be ", name );
}

void check_file( std::string file,
    geode::index_t nb_vertices,
    geode::index_t nb_polygons,
    std::string name )
{
    // Load file
    auto surface = geode::load_triangulated_surface< 3 >( file );
    check_surface( *surface, nb_vertices, nb_polygons, name );

    // Save triangulated tsurf
    const auto output_file_native =
        absl::StrCat( "test_output.", surface->native_extension() );
    geode::save_triangulated_surface( *surface, output_file_native );
    const auto output_file_ts = "test_output.ts";
    geode::save_triangulated_surface( *surface, output_file_ts );

    // Load file
    auto reloaded_surface =
        geode::load_triangulated_surface< 3 >( output_file_native );
    check_surface( *reloaded_surface, nb_vertices, nb_polygons, name );
    auto reloaded_surface_ts =
        geode::load_triangulated_surface< 3 >( output_file_ts );
    check_surface( *reloaded_surface_ts, nb_vertices, nb_polygons, name );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        geode::Logger::set_level( geode::Logger::LEVEL::trace );
        check_file( absl::StrCat( geode::DATA_PATH, "/surf2d_multi.",
                        geode::internal::TSInput::extension() ),
            92, 92, "section1" );
        check_file( absl::StrCat( geode::DATA_PATH, "/surf2d.",
                        geode::internal::TSInput::extension() ),
            46, 46, "section1" );
        check_file( absl::StrCat( geode::DATA_PATH, "/2triangles.",
                        geode::internal::TSInput::extension() ),
            4, 2, "2triangles.ts" );
        check_file( absl::StrCat( geode::DATA_PATH, "/sgrid_tsurf.",
                        geode::internal::TSInput::extension() ),
            4, 2, "FractureSet_N20_1fract_visualisation" );
        check_file( absl::StrCat( geode::DATA_PATH, "Fault_without_crs.",
                        geode::internal::TSInput::extension() ),
            189, 324, "Fault" );
        check_file( absl::StrCat( geode::DATA_PATH, "ts-2props.",
                        geode::internal::TSInput::extension() ),
            4, 2, "test" );
        test_blank_lines();
        test_invalid_indices();
        test_partial_property_header();
        test_object_header_order();
        test_multiple_tsurf_properties();
        test_vertex_lines_with_properties();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
