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
#include <geode/mesh/core/solid_edges.hpp>
#include <geode/mesh/core/solid_facets.hpp>
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

    // Reads the "value id id-id ..." lines of the given DFE material
    absl::flat_hash_map< geode::index_t, double > read_material(
        absl::Span< const std::string > lines, std::string_view material )
    {
        const auto material_line = find_line( lines,
            absl::StrCat( "          <material id=\"", material, "\">" ) );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            material_line.has_value(), "[Test] Missing material ", material );
        absl::flat_hash_map< geode::index_t, double > values;
        for( const auto l :
            geode::Range{ material_line.value() + 2, lines.size() } )
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
                if( tokens[t] == "]]>" )
                {
                    return values;
                }
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
                        values.emplace( id, value ).second, "[Test] Feature ",
                        id, " listed twice in ", material );
                }
            }
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

    void test_facet_conductivity()
    {
        auto solid = create_solid();
        solid->enable_facets();
        auto& facet_manager = solid->facets().facet_attribute_manager();
        auto aperture = create_attribute< double >(
            facet_manager, "diagres_discontinuity_aperture", -1 );
        auto conductivity = create_attribute< double >(
            facet_manager, "diagres_conductivity", -1 );
        const auto facet0 =
            solid->facets().facet_from_vertices( { 0, 1, 2 } ).value();
        const auto facet1 =
            solid->facets().facet_from_vertices( { 0, 1, 3 } ).value();
        aperture->set_value( facet0, 0.5 );
        aperture->set_value( facet1, 0.25 );
        conductivity->set_value( facet0, 5 );
        conductivity->set_value( facet1, 7 );
        const auto nb_attributes = facet_manager.attribute_ids().size();
        geode::save_tetrahedral_solid( *solid, "test_facet_conductivity.fem" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            facet_manager.attribute_ids().size() == nb_attributes,
            "[Test] Saving should not create attributes on the solid" );
        const auto lines = read_lines( "test_facet_conductivity.fem" );
        const geode::index_t feature0 = facet0 < facet1 ? 1 : 2;
        const geode::index_t feature1 = facet0 < facet1 ? 2 : 1;
        const auto area = read_material( lines, "AREA" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            area.size() == 2 && area.at( feature0 ) == 0.5
                && area.at( feature1 ) == 0.25,
            "[Test] Wrong AREA material" );
        const auto cond = read_material( lines, "COND" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            cond.size() == 2 && cond.at( feature0 ) == 5
                && cond.at( feature1 ) == 7,
            "[Test] Wrong COND material" );
    }

    void test_wrong_feature_attribute_type()
    {
        auto solid = create_solid();
        solid->enable_edges();
        auto conduit_area =
            create_attribute< int >( solid->edges().edge_attribute_manager(),
                "diagres_conduit_area", 1 );
        auto vertex_block = create_attribute< int >(
            solid->vertex_attribute_manager(), "Block_ID_vertex", 1 );
        auto polyhedron_block = create_attribute< double >(
            solid->polyhedron_attribute_manager(), "Block_ID_polyhedron", 1 );
        geode::save_tetrahedral_solid( *solid, "test_wrong_type.fem" );
        const auto lines = read_lines( "test_wrong_type.fem" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !find_line( lines, "DFE" ),
            "[Test] Non double conduit area should be ignored" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !find_line( lines, "NODALSETS" ),
            "[Test] Non string vector Block_ID_vertex should be ignored" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !find_line( lines, "ELEMENTALSETS" ),
            "[Test] Non string Block_ID_polyhedron should be ignored" );
    }

    void test_double_precision()
    {
        auto solid = create_solid();
        auto builder = geode::TetrahedralSolidBuilder3D::create( *solid );
        const geode::Point3D point{ { 512345.678, 4567890.123, -1e-9 } };
        builder->set_point( 4, point );
        auto vertex_attribute = create_attribute< double >(
            solid->vertex_attribute_manager(), "vertex_precision", 1e-9 );
        vertex_attribute->set_value( 0, 1.0000001 );
        vertex_attribute->set_value( 1, 1.0000002 );
        auto polyhedron_attribute =
            create_attribute< double >( solid->polyhedron_attribute_manager(),
                "polyhedron_precision", 0.1 );
        polyhedron_attribute->set_value( 1, 123456789.123 );
        solid->enable_facets();
        auto aperture = create_attribute< double >(
            solid->facets().facet_attribute_manager(),
            "diagres_discontinuity_aperture", -1 );
        const auto facet =
            solid->facets().facet_from_vertices( { 0, 1, 2 } ).value();
        aperture->set_value( facet, 1.23456789e-5 );
        geode::save_tetrahedral_solid( *solid, "test_precision.fem" );
        const auto lines = read_lines( "test_precision.fem" );

        const auto coordinates = find_line( lines, "XYZCOOR" ).value();
        const std::vector< std::string_view > tokens =
            absl::StrSplit( lines[coordinates + 5], absl::ByAnyChar( ", " ),
                absl::SkipEmpty() );
        for( const auto d : geode::LRange{ 3 } )
        {
            double value;
            geode::OpenGeodeGeosciencesIOMeshException::test(
                tokens.size() == 3 && absl::SimpleAtod( tokens[d], &value )
                    && value == point.value( d ),
                "[Test] Wrong coordinate ", d, " of node 5" );
        }
        geode::OpenGeodeGeosciencesIOMeshException::test(
            read_distribution( lines, "REF_DIS_I", "X" ).at( 5 )
                    == point.value( 0 )
                && read_distribution( lines, "REF_DIS_I", "Y" ).at( 5 )
                       == point.value( 1 )
                && read_distribution( lines, "REF_DIS_I", "Z" ).at( 5 )
                       == point.value( 2 ),
            "[Test] Wrong coordinate distributions of node 5" );
        const auto vertex_distribution =
            read_distribution( lines, "REF_DIS_I", "vertex_precision" );
        for( const auto v : geode::Range{ solid->nb_vertices() } )
        {
            geode::OpenGeodeGeosciencesIOMeshException::test(
                vertex_distribution.at( v + 1 ) == vertex_attribute->value( v ),
                "[Test] Wrong vertex_precision value for node ", v + 1 );
        }
        const auto polyhedron_distribution =
            read_distribution( lines, "REF_DISE_I", "polyhedron_precision" );
        for( const auto p : geode::Range{ solid->nb_polyhedra() } )
        {
            geode::OpenGeodeGeosciencesIOMeshException::test(
                polyhedron_distribution.at( p + 1 )
                    == polyhedron_attribute->value( p ),
                "[Test] Wrong polyhedron_precision value for element ", p + 1 );
        }
        geode::OpenGeodeGeosciencesIOMeshException::test(
            read_material( lines, "AREA" ).at( 1 ) == 1.23456789e-5,
            "[Test] Wrong aperture value" );
    }

    void test_distribution_counts()
    {
        auto solid = create_solid();
        geode::save_tetrahedral_solid( *solid, "test_counts_empty.fem" );
        auto lines = read_lines( "test_counts_empty.fem" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            lines.at( find_line( lines, "REF_DIS_I" ).value() + 1 )
                == " 3, 5, 0",
            "[Test] Wrong REF_DIS_I count without attribute" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            lines.at( find_line( lines, "REF_DISE_I" ).value() + 1 )
                == " 0, 2, 0",
            "[Test] Wrong REF_DISE_I count without attribute" );

        create_attribute< double >(
            solid->vertex_attribute_manager(), "vertex_count", 1 );
        create_attribute< double >(
            solid->polyhedron_attribute_manager(), "polyhedron_count", 1 );
        create_attribute< std::string >(
            solid->polyhedron_attribute_manager(), "polyhedron_string", "a" );
        geode::save_tetrahedral_solid( *solid, "test_counts.fem" );
        lines = read_lines( "test_counts.fem" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            lines.at( find_line( lines, "REF_DIS_I" ).value() + 1 )
                == " 4, 5, 0",
            "[Test] Wrong REF_DIS_I count" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            lines.at( find_line( lines, "REF_DISE_I" ).value() + 1 )
                == " 1, 2, 0",
            "[Test] Wrong REF_DISE_I count" );
    }

    void test_multi_chunk_material()
    {
        auto solid = geode::TetrahedralSolid3D::create();
        auto builder = geode::TetrahedralSolidBuilder3D::create( *solid );
        for( const auto t : geode::Range{ 1000 } )
        {
            const auto x = static_cast< double >( t );
            const auto v0 =
                builder->create_point( geode::Point3D{ { x, 0, 0 } } );
            builder->create_point( geode::Point3D{ { x + 1, 0, 0 } } );
            builder->create_point( geode::Point3D{ { x, 1, 0 } } );
            builder->create_point( geode::Point3D{ { x, 0, 1 } } );
            builder->create_tetrahedron( { v0, v0 + 1, v0 + 2, v0 + 3 } );
        }
        solid->enable_facets();
        auto aperture = create_attribute< double >(
            solid->facets().facet_attribute_manager(),
            "diagres_discontinuity_aperture", 1 );
        for( const auto facet : geode::Range{ solid->facets().nb_facets() } )
        {
            aperture->set_value( facet, facet % 2 == 0 ? 0.5 : 0.25 );
        }
        geode::save_tetrahedral_solid( *solid, "test_multi_chunk.fem" );
        const auto lines = read_lines( "test_multi_chunk.fem" );
        const auto area = read_material( lines, "AREA" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            area.size() == solid->facets().nb_facets(),
            "[Test] Wrong number of features in multi chunk AREA" );
        for( const auto facet : geode::Range{ solid->facets().nb_facets() } )
        {
            geode::OpenGeodeGeosciencesIOMeshException::test(
                area.at( facet + 1 ) == aperture->value( facet ),
                "[Test] Wrong multi chunk AREA value for feature ", facet + 1 );
        }
    }

    void test_facet_and_edge_features()
    {
        auto solid = create_solid();
        solid->enable_facets();
        solid->enable_edges();
        auto& facet_manager = solid->facets().facet_attribute_manager();
        auto& edge_manager = solid->edges().edge_attribute_manager();
        auto aperture = create_attribute< double >(
            facet_manager, "diagres_discontinuity_aperture", -1 );
        auto facet_conductivity = create_attribute< double >(
            facet_manager, "diagres_conductivity", -1 );
        auto conduit_area = create_attribute< double >(
            edge_manager, "diagres_conduit_area", -1 );
        auto edge_conductivity = create_attribute< double >(
            edge_manager, "diagres_conductivity", -1 );
        const auto facet =
            solid->facets().facet_from_vertices( { 0, 1, 2 } ).value();
        aperture->set_value( facet, 0.5 );
        facet_conductivity->set_value( facet, 5 );
        const auto edge = solid->edges().edge_from_vertices( { 0, 3 } ).value();
        conduit_area->set_value( edge, 0.75 );
        edge_conductivity->set_value( edge, 9 );
        geode::save_tetrahedral_solid( *solid, "test_facet_edge_features.fem" );
        const auto lines = read_lines( "test_facet_edge_features.fem" );
        const auto area = read_material( lines, "AREA" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            area.size() == 2 && area.at( 1 ) == 0.5 && area.at( 2 ) == 0.75,
            "[Test] Wrong AREA material with facet and edge features" );
        const auto cond = read_material( lines, "COND" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            cond.size() == 2 && cond.at( 1 ) == 5 && cond.at( 2 ) == 9,
            "[Test] Wrong COND material with facet and edge features" );
    }

    void test_features_without_valid_values()
    {
        auto solid = create_solid();
        solid->enable_facets();
        solid->enable_edges();
        auto aperture = create_attribute< double >(
            solid->facets().facet_attribute_manager(),
            "diagres_discontinuity_aperture", -1 );
        create_attribute< double >( solid->edges().edge_attribute_manager(),
            "diagres_conduit_area", -1 );
        geode::save_tetrahedral_solid( *solid, "test_no_feature.fem" );
        auto lines = read_lines( "test_no_feature.fem" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !find_line( lines, "DFE" ),
            "[Test] No discrete feature should be written" );

        aperture->set_value(
            solid->facets().facet_from_vertices( { 0, 1, 2 } ).value(), 0.5 );
        aperture->set_value(
            solid->facets().facet_from_vertices( { 0, 1, 3 } ).value(),
            std::numeric_limits< double >::quiet_NaN() );
        geode::save_tetrahedral_solid( *solid, "test_nan_feature.fem" );
        lines = read_lines( "test_nan_feature.fem" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            find_line( lines, "    <fep count=\"1\">" ).has_value()
                && find_line( lines, "    <groups count=\"1\">" ).has_value(),
            "[Test] Only one facet feature should be written" );
    }

    // First tokens of the lines following the given line
    std::vector< std::string > first_tokens(
        absl::Span< const std::string > lines,
        std::string_view start,
        geode::index_t nb_lines )
    {
        const auto start_line = find_line( lines, start ).value();
        std::vector< std::string > tokens;
        for( const auto l : geode::Range{ nb_lines } )
        {
            const std::vector< std::string > line_tokens = absl::StrSplit(
                lines.at( start_line + 1 + l ), ' ', absl::SkipEmpty() );
            tokens.push_back( line_tokens.front() );
        }
        return tokens;
    }

    void test_element_set_names()
    {
        auto solid = create_solid();
        auto polyhedron_sets = create_attribute< std::string >(
            solid->polyhedron_attribute_manager(), "Block_ID_polyhedron", "" );
        polyhedron_sets->set_value( 0, "first_block_with_a_long_name" );
        polyhedron_sets->set_value( 1, "second_block_with_a_long_name" );
        geode::save_tetrahedral_solid( *solid, "test_element_sets.fem" );
        const auto lines = read_lines( "test_element_sets.fem" );
        auto element_sets = first_tokens( lines, "ELEMENTALSETS", 2 );
        absl::c_sort( element_sets );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            element_sets
                == std::vector< std::string >{ "first_block_with_a_long_name",
                    "second_block_with_a_long_name" },
            "[Test] Wrong elemental set names" );
    }

    void test_deterministic_order()
    {
        auto solid = create_solid();
        const std::vector< std::string > names{ "zeta", "alpha", "mu", "beta",
            "omega", "gamma", "delta", "kappa" };
        auto vertex_sets = create_attribute< std::vector< std::string > >(
            solid->vertex_attribute_manager(), "Block_ID_vertex", {} );
        for( const auto v : geode::Range{ solid->nb_vertices() } )
        {
            vertex_sets->set_value( v, { names[v], names[v + 3] } );
        }
        auto polyhedron_sets = create_attribute< std::string >(
            solid->polyhedron_attribute_manager(), "Block_ID_polyhedron", "" );
        polyhedron_sets->set_value( 0, "zeta" );
        polyhedron_sets->set_value( 1, "alpha" );
        solid->enable_facets();
        auto aperture = create_attribute< double >(
            solid->facets().facet_attribute_manager(),
            "diagres_discontinuity_aperture", -1 );
        for( const auto facet : geode::Range{ solid->facets().nb_facets() } )
        {
            aperture->set_value( facet, 10. - facet );
        }
        geode::save_tetrahedral_solid( *solid, "test_order.fem" );
        const auto lines = read_lines( "test_order.fem" );
        const auto nodal_sets = first_tokens( lines, "NODALSETS", 8 );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            absl::c_is_sorted( nodal_sets ), "[Test] Unsorted nodal sets" );
        const auto element_sets = first_tokens( lines, "ELEMENTALSETS", 2 );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            element_sets == std::vector< std::string >{ "alpha", "zeta" },
            "[Test] Unsorted elemental sets" );
        const auto nb_facets = solid->facets().nb_facets();
        const auto material_values =
            first_tokens( lines, "            <![CDATA[", nb_facets );
        std::vector< double > values;
        for( const auto& token : material_values )
        {
            double value;
            geode::OpenGeodeGeosciencesIOMeshException::test(
                absl::SimpleAtod( token, &value ), "[Test] Wrong value" );
            values.push_back( value );
        }
        geode::OpenGeodeGeosciencesIOMeshException::test(
            values.size() == nb_facets && absl::c_is_sorted( values ),
            "[Test] Unsorted material values" );
    }

    void test_vector_attributes()
    {
        auto solid = create_solid();
        create_attribute< std::array< double, 3 > >(
            solid->vertex_attribute_manager(), "vertex_vector",
            { 1., 2., 3. } );
        create_attribute< geode::Point3D >(
            solid->polyhedron_attribute_manager(), "polyhedron_point",
            geode::Point3D{ { 1., 2., 3. } } );
        geode::save_tetrahedral_solid( *solid, "test_vector_attributes.fem" );
        const auto lines = read_lines( "test_vector_attributes.fem" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            !find_line( lines, "vertex_vector" )
                && !find_line( lines, "polyhedron_point" ),
            "[Test] Vector attributes should not be exported" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            lines.at( find_line( lines, "REF_DIS_I" ).value() + 1 )
                    == " 3, 5, 0"
                && lines.at( find_line( lines, "REF_DISE_I" ).value() + 1 )
                       == " 0, 2, 0",
            "[Test] Wrong distribution counts with vector attributes" );
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
        test_facet_conductivity();
        test_wrong_feature_attribute_type();
        test_double_precision();
        test_distribution_counts();
        test_multi_chunk_material();
        test_facet_and_edge_features();
        test_features_without_valid_values();
        test_element_set_names();
        test_deterministic_order();
        test_vector_attributes();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}