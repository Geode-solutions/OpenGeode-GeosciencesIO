/*
 * Copyright (c) 2019 - 2022 Geode-solutions
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

#include <geode/tests_config.h>

#include <geode/basic/assert.h>
#include <geode/basic/attribute_manager.h>
#include <geode/basic/filename.h>
#include <geode/basic/logger.h>

#include <geode/mesh/core/solid_mesh.h>
#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/mesh/io/tetrahedral_solid_input.h>
#include <geode/mesh/io/tetrahedral_solid_output.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>

#include <geode/geosciences/private/lso_input.h>
#include <geode/geosciences/representation/core/structural_model.h>
#include <geode/geosciences/representation/io/structural_model_output.h>

void check_model( const geode::StructuralModel& model,
    geode::index_t nb_corners,
    geode::index_t nb_lines,
    geode::index_t nb_surfaces,
    geode::index_t nb_blocks,
    geode::index_t nb_horizons,
    geode::index_t nb_block_internals,
    geode::index_t nb_vertices_attributes )
{
    OPENGEODE_EXCEPTION( model.nb_corners() == nb_corners,
        "[Test] Number of Corners in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_lines() == nb_lines,
        "[Test] Number of Lines in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_surfaces() == nb_surfaces,
        "[Test] Number of Surfaces in the loaded StructuralModel is not "
        "correct" );
    OPENGEODE_EXCEPTION( model.nb_blocks() == nb_blocks,
        "[Test] Number of Blocks in the loaded "
        "StructuralModel is not correct" );
    OPENGEODE_EXCEPTION( model.nb_horizons() == nb_horizons,
        "[Test] Number of Horizons in the loaded "
        "StructuralModel is not correct" );

    geode::index_t count_block_internals{ 0 };
    for( const auto& block : model.blocks() )
    {
        count_block_internals += model.nb_internals( block.id() );
    }
    OPENGEODE_EXCEPTION( count_block_internals == nb_block_internals,
        "[Test] Number of Block internals in the "
        "loaded StructuralModel is not correct" );
    for( const auto& block : model.blocks() )
    {
        const auto block_attribute_names =
            block.mesh().vertex_attribute_manager().attribute_names();
        OPENGEODE_EXCEPTION(
            block_attribute_names.size() == nb_vertices_attributes + 3,
            "[Test] Number of Block attributes in the loaded "
            "StructuralModel "
            "is not correct: expected ",
            nb_vertices_attributes + 3, " attributes, got ",
            block_attribute_names.size(), " attributes." );
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
            absl::StrCat( geode::filename_without_extension( file ), "_b",
                counter, ".og_tso3d" );
        geode::save_tetrahedral_solid( block_mesh, block_mesh_filename );
        const auto reload_block_mesh =
            geode::load_tetrahedral_solid< 3 >( block_mesh_filename );

        for( const auto name :
            block_mesh.vertex_attribute_manager().attribute_names() )
        {
            OPENGEODE_EXCEPTION(
                reload_block_mesh->vertex_attribute_manager().attribute_exists(
                    name ) );
            OPENGEODE_EXCEPTION(
                reload_block_mesh->vertex_attribute_manager().attribute_type(
                    name )
                    == block_mesh.vertex_attribute_manager().attribute_type(
                        name ),
                "[Test] Wrong attribute type for reloaded mesh" );
        }

        counter++;
    }

    geode::save_structural_model( model, "test.lso" );
    const auto reload_model = geode::load_structural_model( "test.lso" );
    check_model( reload_model, nb_corners, nb_lines, nb_surfaces, nb_blocks,
        nb_horizons, nb_block_internals, 0 );
}

int main()
{
    try
    {
        geode::detail::initialize_geosciences_io();

        geode::Logger::info( "Reading the test.lso file" );
        test_file( absl::StrCat( geode::data_path, "test.",
                       geode::detail::LSOInput::extension() ),
            22, 39, 23, 4, 4, 2, 1 );
        geode::Logger::info( "Reading the vri.lso file" );
        test_file( absl::StrCat( geode::data_path, "vri.",
                       geode::detail::LSOInput::extension() ),
            12, 20, 11, 2, 7, 0, 9 );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
