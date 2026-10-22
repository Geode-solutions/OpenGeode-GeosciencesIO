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

#include <absl/container/flat_hash_map.h>
#include <absl/strings/numbers.h>
#include <absl/strings/str_split.h>

#include <geode/basic/assert.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/tetrahedral_solid_builder.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>

#include <geode/geosciences_io/mesh/internal/fem_output.hpp>

#include <geode/mesh/io/tetrahedral_solid_input.hpp>

namespace
{
    std::vector< std::string > read_lines( std::string_view filename )
    {
        std::ifstream file{ geode::to_string( filename ) };
        geode::OpenGeodeGeosciencesIOMeshException::test(
            file.good(), "[Test] Cannot open file ", filename );
        std::vector< std::string > lines;
        std::string line;
        while( std::getline( file, line ) )
        {
            lines.push_back( line );
        }
        return lines;
    }

    std::optional< geode::index_t > find_line(
        absl::Span< const std::string > lines,
        std::string_view content,
        geode::index_t start = 0 )
    {
        for( const auto l : geode::Range{ start, lines.size() } )
        {
            if( lines[l] == content )
            {
                return l;
            }
        }
        return std::nullopt;
    }

    // Reads the "value id id-id ..." lines following the line named header
    // inside the given section
    absl::flat_hash_map< geode::index_t, double > read_distribution(
        absl::Span< const std::string > lines,
        std::string_view section,
        std::string_view header )
    {
        const auto section_line = find_line( lines, section );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            section_line.has_value(), "[Test] Missing section ", section );
        const auto header_line =
            find_line( lines, header, section_line.value() );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            header_line.has_value(), "[Test] Missing header ", header );
        absl::flat_hash_map< geode::index_t, double > distribution;
        for( const auto l :
            geode::Range{ header_line.value() + 1, lines.size() } )
        {
            const std::vector< std::string_view > tokens =
                absl::StrSplit( lines[l], ' ', absl::SkipEmpty() );
            double value;
            if( tokens.empty() || !absl::SimpleAtod( tokens[0], &value ) )
            {
                break;
            }
            for( const auto t : geode::Range{ 1, tokens.size() } )
            {
                const std::vector< std::string_view > range =
                    absl::StrSplit( tokens[t], '-' );
                geode::index_t first;
                geode::index_t last;
                geode::OpenGeodeGeosciencesIOMeshException::test(
                    absl::SimpleAtoi( range.front(), &first )
                        && absl::SimpleAtoi( range.back(), &last ),
                    "[Test] Wrong id token ", tokens[t] );
                for( const auto id : geode::Range{ first, last + 1 } )
                {
                    geode::OpenGeodeGeosciencesIOMeshException::test(
                        distribution.emplace( id, value ).second, "[Test] Id ",
                        id, " listed twice in ", header );
                }
            }
        }
        return distribution;
    }

    template < typename T >
    std::shared_ptr< geode::VariableAttribute< T > > create_attribute(
        geode::AttributeManager& manager,
        std::string_view name,
        T default_value )
    {
        geode::AttributeValues< T > values;
        values.default_value = default_value;
        const auto id = manager.create_attribute< geode::VariableAttribute, T >(
            name, values, geode::AttributeProperties{} );
        return manager.find_attribute< geode::VariableAttribute, T >( id );
    }

    // Values of the lines following the line named header, in file order
    std::vector< double > read_distribution_values(
        absl::Span< const std::string > lines,
        std::string_view section,
        std::string_view header )
    {
        const auto section_line = find_line( lines, section );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            section_line.has_value(), "[Test] Missing section ", section );
        const auto header_line =
            find_line( lines, header, section_line.value() );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            header_line.has_value(), "[Test] Missing header ", header );
        std::vector< double > values;
        for( const auto l :
            geode::Range{ header_line.value() + 1, lines.size() } )
        {
            const std::vector< std::string_view > tokens =
                absl::StrSplit( lines[l], ' ', absl::SkipEmpty() );
            double value;
            if( tokens.empty() || !absl::SimpleAtod( tokens[0], &value ) )
            {
                break;
            }
            values.push_back( value );
        }
        return values;
    }

    std::unique_ptr< geode::TetrahedralSolid3D > create_solid()
    {
        auto solid = geode::TetrahedralSolid3D::create();
        auto builder = geode::TetrahedralSolidBuilder3D::create( *solid );
        builder->create_point( geode::Point3D{ { 0, 0, 0 } } );
        builder->create_point( geode::Point3D{ { 1, 0, 0 } } );
        builder->create_point( geode::Point3D{ { 0, 1, 0 } } );
        builder->create_point( geode::Point3D{ { 0, 0, 1 } } );
        builder->create_point( geode::Point3D{ { 0, 0, -1 } } );
        builder->create_tetrahedron( { 0, 1, 2, 3 } );
        builder->create_tetrahedron( { 0, 2, 1, 4 } );
        builder->compute_polyhedron_adjacencies();
        return solid;
    }

    void test_vertex_attribute_distribution()
    {
        auto solid = create_solid();
        auto attribute = create_attribute< double >(
            solid->vertex_attribute_manager(), "vertex_prop", 0 );
        for( const auto v : geode::Range{ solid->nb_vertices() } )
        {
            attribute->set_value( v, 10. * v );
        }
        geode::save_tetrahedral_solid( *solid, "test_vertex_attribute.fem" );
        const auto lines = read_lines( "test_vertex_attribute.fem" );
        const auto distribution =
            read_distribution( lines, "REF_DIS_I", "vertex_prop" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            distribution.size() == solid->nb_vertices(),
            "[Test] Wrong number of nodes in vertex_prop distribution" );
        for( const auto v : geode::Range{ solid->nb_vertices() } )
        {
            geode::OpenGeodeGeosciencesIOMeshException::test(
                distribution.at( v + 1 ) == 10. * v,
                "[Test] Wrong vertex_prop value for node ", v + 1 );
        }
    }

    void test_nan_attribute_values()
    {
        auto solid = create_solid();
        const auto nan = std::numeric_limits< double >::quiet_NaN();
        auto vertex_attribute = create_attribute< double >(
            solid->vertex_attribute_manager(), "vertex_nan", nan );
        vertex_attribute->set_value( 1, 2. );
        vertex_attribute->set_value( 3, 1. );
        auto polyhedron_attribute = create_attribute< double >(
            solid->polyhedron_attribute_manager(), "polyhedron_nan", nan );
        polyhedron_attribute->set_value( 1, 3. );
        geode::save_tetrahedral_solid( *solid, "test_nan_attribute.fem" );
        const auto lines = read_lines( "test_nan_attribute.fem" );
        const auto vertex_distribution =
            read_distribution( lines, "REF_DIS_I", "vertex_nan" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            vertex_distribution.size() == 2 && vertex_distribution.at( 2 ) == 2.
                && vertex_distribution.at( 4 ) == 1.,
            "[Test] Wrong vertex_nan distribution" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            read_distribution_values( lines, "REF_DIS_I", "vertex_nan" )
                == std::vector< double >{ 1., 2. },
            "[Test] vertex_nan values should be sorted" );
        const auto polyhedron_distribution =
            read_distribution( lines, "REF_DISE_I", "polyhedron_nan" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            polyhedron_distribution.size() == 1
                && polyhedron_distribution.at( 2 ) == 3.,
            "[Test] Wrong polyhedron_nan distribution" );
    }

    void test_solid_fem_output()
    {
        auto tet_solid = geode::load_tetrahedral_solid< 3 >(
            absl::StrCat( geode::DATA_PATH, "bmsh_342.og_tso3d" ) );
        geode::AttributeValues< geode::index_t > attribute_values;
        attribute_values.default_value = 10;
        attribute_values.no_value = geode::NO_ID;
        geode::AttributeProperties attribute_properties;
        attribute_properties.assignable = false;
        attribute_properties.interpolable = false;
        attribute_properties.transferable = true;
        const auto attribute_id =
            tet_solid->polyhedron_attribute_manager()
                .create_attribute< geode::VariableAttribute, geode::index_t >(
                    "geode_aspect_ratio", attribute_values,
                    attribute_properties );
        auto attribute =
            tet_solid->polyhedron_attribute_manager()
                .find_attribute< geode::VariableAttribute, geode::index_t >(
                    attribute_id );
        geode_unused( attribute );
        geode::save_tetrahedral_solid( *tet_solid, "test.fem" );
    }
} // namespace

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        test_solid_fem_output();
        test_vertex_attribute_distribution();
        test_nan_attribute_values();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}