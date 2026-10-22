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

#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>

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

#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/io/tetrahedral_solid_input.hpp>
#include <geode/mesh/io/tetrahedral_solid_output.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/model_boundary.hpp>
#include <geode/model/mixin/core/surface.hpp>

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

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
