/*
 * Copyright (c) 2019 - 2023 Geode-solutions
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
#include <geode/basic/logger.h>

#include <geode/geosciences_io/mesh/private/grdecl_input.h>
#include <geode/mesh/core/geode/geode_hybrid_solid.h>
#include <geode/mesh/core/hybrid_solid.h>
#include <geode/mesh/io/hybrid_solid_input.h>
#include <geode/mesh/io/hybrid_solid_output.h>

#include <geode/io/mesh/common.h>

void check_solid( geode::HybridSolid3D& solid, geode::index_t nb_polyhedra )
{
    OPENGEODE_EXCEPTION( solid.nb_polyhedra() == nb_polyhedra,
        "Number of polyhedra in the GrdeclHybridSolid is not correct" );
}

void check_file( std::string filename, geode::index_t nb_polyhedra )
{
    // Load File
    auto solid = geode::load_hybrid_solid< 3 >( filename );
    check_solid( *solid, nb_polyhedra );
    geode::save_hybrid_solid( *solid, "toto.vtu" );
}

int main()
{
    try
    {
        geode::IOMeshLibrary::initialize();
        geode::GeosciencesIOMeshLibrary::initialize();
        check_file( absl::StrCat( geode::data_path, "Johansen.",
                        geode::detail::GRDECLInput::extension() ),
            100 * 100 * 11 );
        geode::Logger::info( "[TEST SUCCESS]" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
