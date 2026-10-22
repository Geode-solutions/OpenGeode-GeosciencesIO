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

#include <geode/mesh/builder/edged_curve_builder.hpp>
#include <geode/mesh/builder/hybrid_solid_builder.hpp>
#include <geode/mesh/builder/surface_mesh_builder.hpp>
#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/geode/geode_hybrid_solid.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/core/surface_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>

#include <geode/geosciences_io/model/internal/brep_fem_output.hpp>

#include <geode/mesh/io/tetrahedral_solid_input.hpp>

#include <geode/model/helpers/convert_to_mesh.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/line.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/builder/brep_builder.hpp>
#include <geode/model/representation/io/brep_input.hpp>
#include <geode/model/representation/io/brep_output.hpp>

namespace
{
    std::vector< std::string > read_lines( std::string_view filename )
    {
        std::ifstream file{ geode::to_string( filename ) };
        geode::OpenGeodeGeosciencesIOModelException::test(
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
        geode::OpenGeodeGeosciencesIOModelException::test(
            section_line.has_value(), "[Test] Missing section ", section );
        const auto header_line =
            find_line( lines, header, section_line.value() );
        geode::OpenGeodeGeosciencesIOModelException::test(
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
                geode::OpenGeodeGeosciencesIOModelException::test(
                    absl::SimpleAtoi( range.front(), &first )
                        && absl::SimpleAtoi( range.back(), &last ),
                    "[Test] Wrong id token ", tokens[t] );
                for( const auto id : geode::Range{ first, last + 1 } )
                {
                    geode::OpenGeodeGeosciencesIOModelException::test(
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

    void add_constant_attribute(
        geode::AttributeManager& manager, std::string_view name, double value )
    {
        geode::AttributeValues< double > values;
        values.default_value = value;
        const auto id =
            manager.create_attribute< geode::VariableAttribute, double >(
                name, values, geode::AttributeProperties{} );
        geode_unused( id );
    }

    void test_unmapped_surface_and_line_elements()
    {
        auto brep = geode::load_brep(
            absl::StrCat( geode::DATA_PATH, "remeshed_yom-faults.og_brep" ) );
        geode::BRepBuilder builder{ brep };
        const auto& block = *brep.blocks().begin();
        const auto nb_block_vertices = block.mesh().nb_vertices();
        const std::array< geode::index_t, 3 > block_vertices{ 0,
            nb_block_vertices / 2, nb_block_vertices - 1 };
        const auto& surface = brep.surface( builder.add_surface() );
        auto surface_builder = builder.surface_mesh_builder( surface );
        for( const auto v : geode::LRange{ 3 } )
        {
            surface_builder->create_point(
                block.mesh().point( block_vertices[v] ) );
            builder.set_unique_vertex( { surface.component_id(), v },
                brep.unique_vertex(
                    { block.component_id(), block_vertices[v] } ) );
        }
        surface_builder->create_polygon( { 0, 1, 2 } );
        add_constant_attribute( surface.mesh().polygon_attribute_manager(),
            "diagres_conductivity", 2.5 );
        const auto& line = brep.line( builder.add_line() );
        auto line_builder = builder.line_mesh_builder( line );
        for( const auto v : geode::LRange{ 2 } )
        {
            line_builder->create_point(
                block.mesh().point( block_vertices[v] ) );
            builder.set_unique_vertex( { line.component_id(), v },
                brep.unique_vertex(
                    { block.component_id(), block_vertices[v] } ) );
        }
        line_builder->create_edge( 0, 1 );
        add_constant_attribute(
            line.mesh().edge_attribute_manager(), "diagres_conductivity", 2.5 );
        std::optional< geode::ModelToMeshMappings > converted_mappings;
        try
        {
            converted_mappings.emplace(
                std::get< 1 >( geode::convert_brep_into_solid( brep ) ) );
        }
        catch( const std::exception& )
        {
            // Older OpenGeode versions cannot convert BReps with Surface
            // polygons or Line edges that are not in the converted Solid
            geode::Logger::warning( "[Test] Skipping unmapped elements test: "
                                    "BRep conversion not supported" );
            return;
        }
        const auto& mappings = converted_mappings.value();
        geode::OpenGeodeGeosciencesIOModelException::test(
            !mappings.surface_polygons_mapping.has_mapping_input(
                { surface.id(), 0 } ),
            "[Test] New polygon should not be a solid facet" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            !mappings.line_edges_mapping.has_mapping_input( { line.id(), 0 } ),
            "[Test] New edge should not be a solid edge" );
        geode::save_brep( brep, "test_brep_unmapped.fem" );
    }

    void test_block_attributes()
    {
        auto brep = geode::load_brep(
            absl::StrCat( geode::DATA_PATH, "remeshed_yom-faults.og_brep" ) );
        absl::flat_hash_map< geode::uuid, double > block_values;
        double value{ 1.5 };
        for( const auto& block : brep.blocks() )
        {
            if( block_values.size() == brep.nb_blocks() - 1 )
            {
                break;
            }
            geode::AttributeValues< double > values;
            values.default_value = value;
            const auto id =
                block.mesh()
                    .polyhedron_attribute_manager()
                    .create_attribute< geode::VariableAttribute, double >(
                        "permeability", values, geode::AttributeProperties{} );
            geode_unused( id );
            block_values.emplace( block.id(), value );
            value += 1.5;
        }
        geode::save_brep( brep, "test_brep_block_attributes.fem" );
        const auto lines = read_lines( "test_brep_block_attributes.fem" );
        geode::OpenGeodeGeosciencesIOModelException::test(
            absl::c_count( lines, "permeability" ) == 1,
            "[Test] permeability should be written once" );
        const auto distribution =
            read_distribution( lines, "REF_DISE_I", "permeability" );
        const auto [solid, mappings] = geode::convert_brep_into_solid( brep );
        geode::index_t nb_values{ 0 };
        for( const auto& block : brep.blocks() )
        {
            const auto block_value = block_values.find( block.id() );
            for( const auto polyhedron :
                geode::Range{ block.mesh().nb_polyhedra() } )
            {
                const auto element = mappings.solid_polyhedra_mapping.in2out(
                                         { block.id(), polyhedron } )
                                     + 1;
                if( block_value == block_values.end() )
                {
                    geode::OpenGeodeGeosciencesIOModelException::test(
                        !distribution.contains( element ), "[Test] Element ",
                        element, " should not have a permeability" );
                    continue;
                }
                geode::OpenGeodeGeosciencesIOModelException::test(
                    distribution.at( element ) == block_value->second,
                    "[Test] Wrong permeability for element ", element );
                nb_values++;
            }
        }
        geode::OpenGeodeGeosciencesIOModelException::test(
            distribution.size() == nb_values,
            "[Test] Wrong number of permeability values" );
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
        test_unmapped_surface_and_line_elements();
        test_block_attributes();
        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}