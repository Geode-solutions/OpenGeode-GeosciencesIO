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

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_map.h>
#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_replace.h>
#include <absl/strings/str_split.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <geode/basic/assert.hpp>
#include <geode/basic/attribute.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/string.hpp>
#include <geode/basic/types.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/mesh_factory.hpp>
#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/surface_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/core/triangulated_surface.hpp>
#include <geode/mesh/io/tetrahedral_solid_input.hpp>
#include <geode/mesh/io/tetrahedral_solid_output.hpp>

#include <geode/mesh/builder/point_set_builder.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/corner.hpp>
#include <geode/model/mixin/core/model_boundary.hpp>
#include <geode/model/mixin/core/surface.hpp>

#include <geode/geosciences/explicit/representation/builder/structural_model_builder.hpp>
#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences/explicit/representation/io/structural_model_output.hpp>
#include <geode/geosciences_io/model/internal/lso_input.hpp>

constexpr auto nb_mandatory_attributes = 2;

void check_model( const geode::StructuralModel& model,
    geode::index_t nb_corners,
    geode::index_t nb_lines,
    geode::index_t nb_surfaces,
    geode::index_t nb_blocks,
    geode::index_t nb_horizons,
    geode::index_t nb_block_internals,
    geode::index_t nb_vertices_attributes )
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
        model.nb_horizons() == nb_horizons, "Number of Horizons in the loaded "
                                            "StructuralModel is not correct" );

    geode::index_t count_block_internals{ 0 };
    for( const auto& block : model.blocks() )
    {
        count_block_internals += model.nb_internals( block.id() );
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_block_internals == nb_block_internals,
        "Number of Block internals in the "
        "loaded StructuralModel is not correct" );
    for( const auto& block : model.blocks() )
    {
        const auto block_attribute_ids =
            block.mesh().vertex_attribute_manager().attribute_ids();
        geode::OpenGeodeGeosciencesIOModelException::test(
            block_attribute_ids.size()
                == nb_vertices_attributes + nb_mandatory_attributes,
            "Number of Block attributes in the loaded "
            "StructuralModel is not correct: expected ",
            nb_vertices_attributes + nb_mandatory_attributes,
            " attributes, got ", block_attribute_ids.size(), " attributes." );
    }
}

void test_file( std::string file,
    geode::index_t nb_corners,
    geode::index_t nb_lines,
    geode::index_t nb_surfaces,
    geode::index_t nb_blocks,
    geode::index_t nb_horizons,
    geode::index_t nb_block_internals,
    geode::index_t nb_vertices_attributes )
{
    const auto model = geode::load_structural_model( file );
    check_model( model, nb_corners, nb_lines, nb_surfaces, nb_blocks,
        nb_horizons, nb_block_internals, nb_vertices_attributes );
    geode::index_t counter{ 0 };
    for( const auto& block : model.blocks() )
    {
        const auto& block_mesh = block.mesh< geode::TetrahedralSolid3D >();
        const auto block_mesh_filename =
            absl::StrCat( geode::filename_without_extension( file ).string(),
                "_b", counter, ".og_tso3d" );
        geode::save_tetrahedral_solid( block_mesh, block_mesh_filename );
        const auto reload_block_mesh =
            geode::load_tetrahedral_solid< 3 >( block_mesh_filename );

        for( const auto ids :
            block_mesh.vertex_attribute_manager().attribute_ids() )
        {
            geode::OpenGeodeGeosciencesIOModelException::test(
                reload_block_mesh->vertex_attribute_manager().attribute_exists(
                    ids ) );
            geode::OpenGeodeGeosciencesIOModelException::test(
                reload_block_mesh->vertex_attribute_manager().attribute_type(
                    ids )
                    == block_mesh.vertex_attribute_manager().attribute_type(
                        ids ),
                "Wrong attribute type for reloaded mesh" );
        }

        counter++;
    }

    geode::save_structural_model( model, "test.lso" );
    const auto reload_model = geode::load_structural_model( "test.lso" );
    check_model( reload_model, nb_corners, nb_lines, nb_surfaces, nb_blocks,
        nb_horizons, nb_block_internals, 0 );
}

std::vector< std::string > read_lines( std::string_view file )
{
    std::ifstream input{ geode::to_string( file ) };
    std::vector< std::string > lines;
    std::string line;
    while( std::getline( input, line ) )
    {
        if( !line.empty() && line.back() == '\r' )
        {
            line.pop_back();
        }
        lines.push_back( line );
    }
    return lines;
}

std::string write_modified_vri( std::string_view name,
    const std::function< void( std::vector< std::string >& ) >& modifier )
{
    auto lines = read_lines( absl::StrCat( geode::DATA_PATH, "vri.", "lso" ) );
    modifier( lines );
    const auto filename = absl::StrCat( name, ".lso" );
    std::ofstream output{ filename };
    for( const auto& line : lines )
    {
        output << line << "\n";
    }
    return filename;
}

geode::index_t find_line( const std::vector< std::string >& lines,
    std::string_view start,
    geode::index_t from = 0 )
{
    for( const auto i : geode::Range{ from, lines.size() } )
    {
        if( absl::StartsWith( lines[i], start ) )
        {
            return i;
        }
    }
    throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
        geode::OpenGeodeException::TYPE::data, "Line not found: ", start };
}

geode::index_t count_tetrahedra( const geode::StructuralModel& model )
{
    geode::index_t count{ 0 };
    for( const auto& block : model.blocks() )
    {
        count += block.mesh().nb_polyhedra();
    }
    return count;
}

void check_load_throws( std::string_view file )
{
    try
    {
        const auto model = geode::load_structural_model( file );
    }
    catch( const geode::OpenGeodeException& exception )
    {
        geode::Logger::info( "Expected exception: ", exception.what() );
        return;
    }
    throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
        geode::OpenGeodeException::TYPE::data,
        "Loading should have thrown an OpenGeodeException: ", file };
}

void test_malformed_tokens()
{
    check_load_throws( write_modified_vri(
        "lso_truncated_after_tetra", []( std::vector< std::string >& lines ) {
            lines.resize( find_line( lines, "TETRA " ) + 1 );
        } ) );
    check_load_throws( write_modified_vri(
        "lso_short_tetra", []( std::vector< std::string >& lines ) {
            lines[find_line( lines, "TETRA " )] = "TETRA 1 2 3";
        } ) );
    check_load_throws( write_modified_vri(
        "lso_short_vertex", []( std::vector< std::string >& lines ) {
            lines[find_line( lines, "PVRTX 5 " )] = "PVRTX 5 0.5";
        } ) );
    check_load_throws( write_modified_vri(
        "lso_short_shared", []( std::vector< std::string >& lines ) {
            lines[find_line( lines, "SHAREDPVRTX" )] = "SHAREDPVRTX 104";
        } ) );
    check_load_throws( write_modified_vri(
        "lso_short_triangle", []( std::vector< std::string >& lines ) {
            lines[find_line( lines, "TRGL" )] = "TRGL 1 2";
        } ) );
    check_load_throws( write_modified_vri(
        "lso_short_region", []( std::vector< std::string >& lines ) {
            lines[find_line( lines, "MODEL_REGION" )] = "MODEL_REGION";
        } ) );
    const auto model = geode::load_structural_model( write_modified_vri(
        "lso_missing_ctetra", []( std::vector< std::string >& lines ) {
            lines.erase( lines.begin() + find_line( lines, "# CTETRA" ) );
        } ) );
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_tetrahedra( model ) == 291,
        "Wrong number of tetrahedra with missing CTETRA line" );
}

std::vector< std::string > split( std::string_view line )
{
    return absl::StrSplit( line, ' ', absl::SkipEmpty() );
}

void renumber_vertices( std::vector< std::string >& lines,
    const std::function< geode::index_t( geode::index_t ) >& renumber )
{
    for( auto& line : lines )
    {
        auto tokens = split( line );
        if( tokens.empty() )
        {
            continue;
        }
        geode::index_t nb_ids{ 0 };
        if( tokens[0] == "PVRTX" )
        {
            nb_ids = 1;
        }
        else if( tokens[0] == "SHAREDPVRTX" )
        {
            nb_ids = 2;
        }
        else if( tokens[0] == "TETRA" )
        {
            nb_ids = 4;
        }
        else if( tokens[0] == "TRGL" || tokens[0] == "KEYVERTICES" )
        {
            nb_ids = 3;
        }
        for( const auto i : geode::Range{ 1, nb_ids + 1 } )
        {
            tokens[i] =
                absl::StrCat( renumber( geode::string_to_index( tokens[i] ) ) );
        }
        line = absl::StrJoin( tokens, " " );
    }
}

void check_vri_model( const geode::StructuralModel& model,
    geode::index_t nb_vertices_attributes = 9 )
{
    check_model( model, 12, 20, 11, 2, 7, 0, nb_vertices_attributes );
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_tetrahedra( model ) == 292, "Wrong number of tetrahedra" );
}

void add_shared_vertex( std::vector< std::string >& lines,
    std::string_view id,
    std::string_view reference )
{
    const auto shared = split( lines[find_line( lines, "SHAREDPVRTX 104 " )] );
    std::vector< std::string > new_line{ "SHAREDPVRTX", geode::to_string( id ),
        geode::to_string( reference ) };
    new_line.insert( new_line.end(), shared.begin() + 3, shared.end() );
    lines.insert(
        lines.begin() + find_line( lines, "BEGIN_VERTEX_REGION_INDICATORS" ),
        absl::StrJoin( new_line, " " ) );
}

void test_shared_vertices()
{
    check_load_throws( write_modified_vri(
        "lso_shared_out_of_range", []( std::vector< std::string >& lines ) {
            add_shared_vertex( lines, "124", "999" );
        } ) );
    check_load_throws( write_modified_vri(
        "lso_shared_zero", []( std::vector< std::string >& lines ) {
            add_shared_vertex( lines, "124", "0" );
        } ) );
    check_load_throws( write_modified_vri(
        "lso_tetra_out_of_range", []( std::vector< std::string >& lines ) {
            lines[find_line( lines, "TETRA " )] = "TETRA 1 2 3 999";
        } ) );
    check_load_throws( write_modified_vri(
        "lso_triangle_zero", []( std::vector< std::string >& lines ) {
            lines[find_line( lines, "TRGL" )] = "TRGL 0 1 2";
        } ) );
    // A shared vertex referring to another shared vertex
    check_vri_model( geode::load_structural_model( write_modified_vri(
        "lso_shared_chain", []( std::vector< std::string >& lines ) {
            add_shared_vertex( lines, "124", "104" );
        } ) ) );
    // Shared vertex 104 (copy of 13) moved right after vertex 13
    check_vri_model( geode::load_structural_model( write_modified_vri(
        "lso_shared_interleaved", []( std::vector< std::string >& lines ) {
            renumber_vertices( lines, []( geode::index_t id ) {
                if( id == 104 )
                {
                    return geode::index_t{ 14 };
                }
                if( id >= 14 && id <= 103 )
                {
                    return id + 1;
                }
                return id;
            } );
            const auto shared = find_line( lines, "SHAREDPVRTX 14 " );
            auto line = lines[shared];
            lines.erase( lines.begin() + shared );
            lines.insert(
                lines.begin() + find_line( lines, "PVRTX 15 " ), line );
        } ) ) );
    // Triangles referring to shared vertices instead of original ones
    check_vri_model( geode::load_structural_model( write_modified_vri(
        "lso_triangle_on_shared", []( std::vector< std::string >& lines ) {
            absl::flat_hash_map< std::string, std::string > copies;
            for( const auto& line : lines )
            {
                const auto tokens = split( line );
                if( !tokens.empty() && tokens[0] == "SHAREDPVRTX" )
                {
                    copies.emplace( tokens[2], tokens[1] );
                }
            }
            for( auto& line : lines )
            {
                auto tokens = split( line );
                if( tokens.empty() || tokens[0] != "TRGL" )
                {
                    continue;
                }
                for( const auto i : geode::LRange{ 1, 4 } )
                {
                    const auto it = copies.find( tokens[i] );
                    if( it != copies.end() )
                    {
                        tokens[i] = it->second;
                    }
                }
                line = absl::StrJoin( tokens, " " );
            }
        } ) ) );
}

void test_vertex_ids()
{
    check_vri_model( geode::load_structural_model( write_modified_vri(
        "lso_vertex_ids_with_gaps", []( std::vector< std::string >& lines ) {
            renumber_vertices( lines, []( geode::index_t id ) {
                return 2 * id + 10;
            } );
        } ) ) );
    check_load_throws( write_modified_vri(
        "lso_duplicated_vertex_id", []( std::vector< std::string >& lines ) {
            auto tokens = split( lines[find_line( lines, "PVRTX 5 " )] );
            tokens[1] = "3";
            lines[find_line( lines, "PVRTX 5 " )] =
                absl::StrJoin( tokens, " " );
        } ) );
}

geode::index_t count_triangles( const geode::StructuralModel& model )
{
    geode::index_t count{ 0 };
    for( const auto& surface : model.surfaces() )
    {
        count += surface.mesh().nb_polygons();
    }
    return count;
}

void test_model_section()
{
    const auto nb_triangles = count_triangles( geode::load_structural_model(
        absl::StrCat( geode::DATA_PATH, "vri.lso" ) ) );
    const auto without_keyvertices =
        geode::load_structural_model( write_modified_vri(
            "lso_without_keyvertices", []( std::vector< std::string >& lines ) {
                lines.erase( std::remove_if( lines.begin(), lines.end(),
                                 []( const std::string& line ) {
                                     return absl::StartsWith(
                                         line, "KEYVERTICES" );
                                 } ),
                    lines.end() );
            } ) );
    check_vri_model( without_keyvertices );
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_triangles( without_keyvertices ) == nb_triangles,
        "Wrong number of triangles without KEYVERTICES" );
    const auto with_unknown_lines =
        geode::load_structural_model( write_modified_vri(
            "lso_with_unknown_lines", []( std::vector< std::string >& lines ) {
                lines.insert( lines.begin() + find_line( lines, "MODEL" ) + 1,
                    "UNKNOWN 1" );
                lines.insert( lines.begin() + find_line( lines, "SURFACE" ) + 1,
                    "UNKNOWN 2" );
                lines.insert( lines.begin()
                                  + find_line( lines, "SURFACE",
                                      find_line( lines, "SURFACE" ) + 1 ),
                    "UNKNOWN 3" );
                lines.insert(
                    lines.begin() + find_line( lines, "MODEL_REGION" ) + 1,
                    "UNKNOWN 4" );
            } ) );
    check_vri_model( with_unknown_lines );
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_triangles( with_unknown_lines ) == nb_triangles,
        "Wrong number of triangles with unknown lines" );
    const auto unnamed_surface =
        geode::load_structural_model( write_modified_vri(
            "lso_unnamed_surface", []( std::vector< std::string >& lines ) {
                lines[find_line( lines, "SURFACE" )] = "SURFACE";
            } ) );
    check_vri_model( unnamed_surface );
}

void test_region_names()
{
    const auto rename = []( std::vector< std::string >& lines ) {
        for( auto& line : lines )
        {
            if( absl::StartsWith( line, "# CTETRA Region_2 " )
                || absl::StartsWith( line, "MODEL_REGION Region_2 " ) )
            {
                line = absl::StrReplaceAll(
                    line, { { "Region_2", "\"Region 2\"" } } );
            }
        }
    };
    const auto check_names = []( const geode::StructuralModel& model,
                                 geode::index_t nb_vertices_attributes ) {
        check_vri_model( model, nb_vertices_attributes );
        for( const auto& block : model.blocks() )
        {
            geode::OpenGeodeGeosciencesIOModelException::test(
                block.name() == "Region 2" || block.name() == "Region_1",
                "Wrong Block name: ", block.name().value() );
        }
    };
    const auto model = geode::load_structural_model(
        write_modified_vri( "lso_region_names_with_spaces", rename ) );
    check_names( model, 9 );
    geode::save_structural_model( model, "lso_region_names_saved.lso" );
    const auto reloaded =
        geode::load_structural_model( "lso_region_names_saved.lso" );
    check_names( reloaded, 0 );

    const auto unknown_region =
        geode::load_structural_model( write_modified_vri(
            "lso_unknown_region", []( std::vector< std::string >& lines ) {
                lines[find_line( lines, "# CTETRA" )] =
                    "# CTETRA Unknown none none none none";
            } ) );
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_tetrahedra( unknown_region ) == 291,
        "Wrong number of tetrahedra with unknown region" );
}

void check_save_throws(
    const geode::StructuralModel& model, std::string_view file )
{
    try
    {
        geode::save_structural_model( model, file );
    }
    catch( const geode::OpenGeodeException& exception )
    {
        geode::Logger::info( "Expected exception: ", exception.what() );
        return;
    }
    throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
        geode::OpenGeodeException::TYPE::data,
        "Saving should have thrown an OpenGeodeException: ", file };
}

void test_output_vertices_outside_blocks()
{
    auto model = geode::load_structural_model(
        absl::StrCat( geode::DATA_PATH, "vri.lso" ) );
    geode::StructuralModelBuilder builder{ model };
    const auto& corner = model.corner( builder.add_corner() );
    builder.corner_mesh_builder( corner )->create_point(
        geode::Point3D{ { 10, 10, 10 } } );
    builder.set_unique_vertex(
        { corner.component_id(), 0 }, builder.create_unique_vertex() );
    geode::save_structural_model( model, "lso_free_corner.lso" );
    check_vri_model( geode::load_structural_model( "lso_free_corner.lso" ), 0 );

    builder.create_unique_vertex();
    check_save_throws( model, "lso_orphan_unique_vertex.lso" );
}

void test_output_surfaces_outside_collections()
{
    auto model = geode::load_structural_model(
        absl::StrCat( geode::DATA_PATH, "vri.lso" ) );
    geode::StructuralModelBuilder builder{ model };
    std::vector< geode::uuid > horizons;
    for( const auto& horizon : model.horizons() )
    {
        horizons.push_back( horizon.id() );
    }
    for( const auto& horizon : horizons )
    {
        builder.remove_horizon( model.horizon( horizon ) );
    }
    geode_unused( builder.add_surface( geode::MeshFactory::default_impl(
        geode::TriangulatedSurface3D::type_name_static() ) ) );
    geode::save_structural_model( model, "lso_surfaces_without_horizon.lso" );
    const auto reloaded =
        geode::load_structural_model( "lso_surfaces_without_horizon.lso" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        reloaded.nb_surfaces() == 11 && reloaded.nb_blocks() == 2
            && count_tetrahedra( reloaded ) == 292,
        "Wrong reloaded model with Surfaces outside any collection" );
}

void test_border_without_corner()
{
    auto lines = read_lines( absl::StrCat( geode::DATA_PATH, "test.lso" ) );
    const auto first_triangle =
        find_line( lines, "TRGL", find_line( lines, "TFACE 4" ) );
    lines.insert( lines.begin() + first_triangle + 1, lines[first_triangle] );
    std::ofstream output{ "lso_duplicated_triangle.lso" };
    for( const auto& line : lines )
    {
        output << line << "\n";
    }
    output.close();
    check_load_throws( "lso_duplicated_triangle.lso" );
}

geode::index_t count_block_boundaries( const geode::StructuralModel& model )
{
    geode::index_t count{ 0 };
    for( const auto& block : model.blocks() )
    {
        count += model.nb_boundaries( block.id() );
    }
    return count;
}

void test_partial_conformity()
{
    const auto nb_boundaries =
        count_block_boundaries( geode::load_structural_model(
            absl::StrCat( geode::DATA_PATH, "vri.lso" ) ) );
    // Remove a tetrahedron having facets on two surfaces
    const auto model = geode::load_structural_model( write_modified_vri(
        "lso_partial_conformity", []( std::vector< std::string >& lines ) {
            const auto tetra = find_line( lines, "TETRA 112 59 110 111" );
            lines.erase( lines.begin() + tetra, lines.begin() + tetra + 2 );
        } ) );
    geode::OpenGeodeGeosciencesIOModelException::test(
        count_block_boundaries( model ) == nb_boundaries,
        "Wrong number of Block boundaries with partial conformity" );
}

void test_border_loop_without_corner()
{
    const auto filename = "lso_border_loop_without_corner.lso";
    {
        std::ofstream output{ filename };
        output << "GOCAD LightTSolid 1\n"
                  "HEADER {\nname: two_tetrahedra\n}\n"
                  "VRTX 1 0 0 0\nVRTX 2 1 0 0\nVRTX 3 0 1 0\n"
                  "VRTX 4 0 0 1\nVRTX 5 0 0 -1\n"
                  "TETRA 1 2 3 4\n# CTETRA A none none none none\n"
                  "TETRA 1 2 3 5\n# CTETRA B none none none none\n"
                  "MODEL\n"
                  "SURFACE top\nTFACE 1\nKEYVERTICES 1 2 4\n"
                  "TRGL 1 2 4\nTRGL 2 3 4\nTRGL 3 1 4\n"
                  "SURFACE middle\nTFACE 2\nKEYVERTICES 1 2 3\n"
                  "TRGL 1 2 3\n"
                  "SURFACE bottom\nTFACE 3\nKEYVERTICES 1 2 5\n"
                  "TRGL 1 2 5\nTRGL 2 3 5\nTRGL 3 1 5\n"
                  "MODEL_REGION A +1\nMODEL_REGION B -3\nEND\n";
    }
    const auto model = geode::load_structural_model( filename );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_blocks() == 2 && model.nb_surfaces() == 3,
        "Wrong number of Blocks or Surfaces" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        model.nb_corners() == 1 && model.nb_lines() == 1,
        "A border loop without Corner should get one Corner and one Line" );
    for( const auto& surface : model.surfaces() )
    {
        geode::OpenGeodeGeosciencesIOModelException::test(
            model.nb_boundaries( surface.id() ) == 1,
            "Surface should be bounded by the Line" );
    }
}

std::vector< std::array< double, 3 > > sorted_block_points(
    const geode::StructuralModel& model )
{
    std::vector< std::array< double, 3 > > points;
    for( const auto& block : model.blocks() )
    {
        const auto& mesh = block.mesh();
        for( const auto v : geode::Range{ mesh.nb_vertices() } )
        {
            const auto& point = mesh.point( v );
            points.push_back(
                { point.value( 0 ), point.value( 1 ), point.value( 2 ) } );
        }
    }
    absl::c_sort( points );
    return points;
}

void test_output_precision()
{
    const auto model = geode::load_structural_model(
        absl::StrCat( geode::DATA_PATH, "vri.lso" ) );
    geode::save_structural_model( model, "lso_precision.lso" );
    const auto reloaded = geode::load_structural_model( "lso_precision.lso" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        sorted_block_points( model ) == sorted_block_points( reloaded ),
        "Coordinates are not exactly preserved by LSO output" );
}

void test_tetra_properties()
{
    const auto model = geode::load_structural_model(
        absl::StrCat( geode::DATA_PATH, "lso_tetra_properties.lso" ) );
    geode::index_t nb_tetrahedra{ 0 };
    for( const auto& block : model.blocks() )
    {
        const auto& mesh = block.mesh();
        const auto& manager = mesh.polyhedron_attribute_manager();
        const auto x_ids = manager.attribute_ids_matching_name( "centroid_x" );
        const auto y_ids = manager.attribute_ids_matching_name( "centroid_y" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            x_ids && x_ids->size() == 1 && y_ids && y_ids->size() == 1,
            "Tetrahedra properties not found" );
        const auto centroid_x =
            manager.find_read_only_attribute< double >( x_ids->front() );
        const auto centroid_y =
            manager.find_read_only_attribute< std::array< double, 2 > >(
                y_ids->front() );
        for( const auto t : geode::Range{ mesh.nb_polyhedra() } )
        {
            const auto barycenter = mesh.polyhedron_barycenter( t );
            geode::OpenGeodeGeosciencesIOModelException::test(
                std::fabs( centroid_x->value( t ) - barycenter.value( 0 ) )
                    < geode::GLOBAL_EPSILON,
                "Wrong tetrahedron property value on tetrahedron ", t,
                " of Block ", block.name().value() );
            for( const auto value : centroid_y->value( t ) )
            {
                geode::OpenGeodeGeosciencesIOModelException::test(
                    std::fabs( value - barycenter.value( 1 ) )
                        < geode::GLOBAL_EPSILON,
                    "Wrong tetrahedron array property value on "
                    "tetrahedron ",
                    t, " of Block ", block.name().value() );
            }
            nb_tetrahedra++;
        }
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        nb_tetrahedra == 292, "Wrong number of tetrahedra" );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();

        geode::Logger::info( "Reading the test.lso file" );
        test_file( absl::StrCat( geode::DATA_PATH, "test.",
                       geode::internal::LSOInput::extension() ),
            22, 39, 23, 4, 4, 2, 1 );
        geode::Logger::info( "Reading the vri.lso file" );
        test_file( absl::StrCat( geode::DATA_PATH, "vri.",
                       geode::internal::LSOInput::extension() ),
            12, 20, 11, 2, 7, 0, 9 );
        test_tetra_properties();
        test_malformed_tokens();
        test_shared_vertices();
        test_vertex_ids();
        test_model_section();
        test_region_names();
        test_output_vertices_outside_blocks();
        test_output_surfaces_outside_collections();
        test_border_without_corner();
        test_partial_conformity();
        test_border_loop_without_corner();
        test_output_precision();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
