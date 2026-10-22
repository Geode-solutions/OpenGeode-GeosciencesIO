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

#include <geode/basic/assert.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/hybrid_solid_builder.hpp>
#include <geode/mesh/core/geode/geode_hybrid_solid.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>

#include <geode/geosciences_io/model/internal/brep_fem_output.hpp>

#include <geode/mesh/io/tetrahedral_solid_input.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/builder/brep_builder.hpp>
#include <geode/model/representation/io/brep_input.hpp>
#include <geode/model/representation/io/brep_output.hpp>

namespace
{

    void test_brep_fem_output()
    {
        auto brep = geode::load_brep(
            absl::StrCat( geode::DATA_PATH, "remeshed_yom-faults.og_brep" ) );
        geode::save_brep( brep, "test_brep.fem" );
    }

    bool write_throws( const geode::BRep& brep, std::string_view filename )
    {
        try
        {
            geode::internal::BRepFemOutput{ filename }.write( brep );
        }
        catch( const geode::OpenGeodeException& )
        {
            return true;
        }
        return false;
    }

    void test_non_tetrahedral_brep()
    {
        geode::BRep empty_brep;
        geode::internal::BRepFemOutput output{ "test_empty_brep.fem" };
        geode::OpenGeodeGeosciencesIOModelException::test(
            !output.is_saveable( empty_brep ),
            "[Test] BRep without Block should not be saveable" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            write_throws( empty_brep, "test_empty_brep.fem" ),
            "[Test] Saving a BRep without Block should throw" );

        geode::BRep brep;
        geode::BRepBuilder builder{ brep };
        const auto& block_id = builder.add_block(
            geode::OpenGeodeHybridSolid3D::impl_name_static() );
        const auto& block = brep.block( block_id );
        auto mesh_builder =
            builder.block_mesh_builder< geode::HybridSolid3D >( block );
        mesh_builder->create_point( geode::Point3D{ { 0, 0, 0 } } );
        mesh_builder->create_point( geode::Point3D{ { 1, 0, 0 } } );
        mesh_builder->create_point( geode::Point3D{ { 0, 1, 0 } } );
        mesh_builder->create_point( geode::Point3D{ { 0, 0, 1 } } );
        mesh_builder->create_tetrahedron( { 0, 1, 2, 3 } );
        builder.create_unique_vertices( 4 );
        for( const auto v : geode::Range{ 4 } )
        {
            builder.set_unique_vertex( { block.component_id(), v }, v );
        }
        geode::OpenGeodeGeosciencesIOModelException::test(
            !output.is_saveable( brep ),
            "[Test] BRep with non tetrahedral Block should not be saveable" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            write_throws( brep, "test_hybrid_brep.fem" ),
            "[Test] Saving a BRep with non tetrahedral Block should throw" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        test_brep_fem_output();
        test_non_tetrahedral_brep();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}