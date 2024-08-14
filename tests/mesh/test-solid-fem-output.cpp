/*
 * Copyright (c) 2019 - 2024 Geode-solutions
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

#include <geode/mesh/core/tetrahedral_solid.hpp>

#include <geode/geosciences_io/mesh/internal/fem_output.hpp>

#include <geode/mesh/io/tetrahedral_solid_input.hpp>

namespace
{
    void test_solid_fem_output()
    {
        // auto tet_solid = geode::TetrahedralSolid3D::create();
        auto tet_solid = geode::load_tetrahedral_solid< 3 >(
            absl::StrCat( geode::DATA_PATH, "bmsh_342.og_tso3d" ) );
        tet_solid->polyhedron_attribute_manager()
            .find_or_create_attribute< geode::VariableAttribute,
                geode::index_t >( "geode_aspect_ratio", 10 );
        geode::internal::SolidFemOutput fem_output( "test.fem" );
        fem_output.write( *tet_solid );
    }
} // namespace

int main()
{
    try
    {
        geode::GeosciencesIOMeshLibrary::initialize();
        test_solid_fem_output();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}