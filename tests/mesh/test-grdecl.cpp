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
#include <geode/basic/logger.hpp>

#include <geode/geosciences_io/mesh/internal/grdecl_input.hpp>
#include <geode/mesh/core/geode/geode_hybrid_solid.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/io/hybrid_solid_input.hpp>
#include <geode/mesh/io/hybrid_solid_output.hpp>

void check_solid( const geode::HybridSolid3D& solid,
    geode::index_t nb_polyhedra,
    geode::index_t nb_vertices )
{
    OPENGEODE_EXCEPTION( solid.nb_polyhedra() == nb_polyhedra,
        "Number of polyhedra in the GrdeclHybridSolid is not correct" );
    OPENGEODE_EXCEPTION( solid.nb_vertices() == nb_vertices,
        "Number of vertices in the GrdeclHybridSolid is not correct " );
    for( const auto polyhedra : geode::Range{ solid.nb_polyhedra() } )
    {
        OPENGEODE_EXCEPTION( solid.polyhedron_volume( polyhedra ) > 0.,
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

int main()
{
    try
    {
        geode::GeosciencesIOMeshLibrary::initialize();

        check_file( absl::StrCat( geode::DATA_PATH, "Simple20x20x5_Fault.",
                        geode::internal::GRDECLInput::extension() ),

            20 * 20 * 5, 21 * 6 * ( 21 + 1 ) );
        geode::Logger::info( "[TEST SUCCESS]" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
