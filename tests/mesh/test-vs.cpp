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

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/io/point_set_input.hpp>
#include <geode/mesh/io/point_set_output.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>
#include <geode/geosciences_io/mesh/internal/vs_input.hpp>

namespace
{
    std::string write_test_file(
        std::string_view filename, std::string_view content )
    {
        std::ofstream file{ geode::to_string( filename ), std::ios::binary };
        file << content;
        return geode::to_string( filename );
    }

    std::unique_ptr< geode::PointSet3D > load_test_file(
        std::string_view filename, std::string_view content )
    {
        return geode::load_point_set< 3 >(
            write_test_file( filename, content ) );
    }

    void test_points_normal()
    {
        const auto pointset = geode::load_point_set< 3 >(
            absl::StrCat( geode::DATA_PATH, "points.vs" ) );
        const auto& manager = pointset->vertex_attribute_manager();
        const auto ids = manager.attribute_ids_matching_name( "normal" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            ids.has_value() && ids->size() == 1,
            "points.vs normal property is missing" );
        const auto normal =
            manager.find_read_only_attribute< std::array< double, 3 > >(
                ids->front() );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            normal->value( 5 ) == std::array< double, 3 >{ 0, 0, 0 },
            "Wrong points.vs normal value" );
    }

    void test_multiple_vset_properties()
    {
        const auto pointset = load_test_file( "test_vs_multi_properties.vs",
            "GOCAD VSet 1\nHEADER {\nname: a\n}\nPROPERTIES p\n"
            "PVRTX 1 0 0 0 1\nPVRTX 2 1 0 0 2\nEND\nGOCAD VSet 1\n"
            "HEADER {\nname: b\n}\nPROPERTIES p\nPVRTX 1 0 1 0 3\nEND\n" );
        const auto& manager = pointset->vertex_attribute_manager();
        const auto ids = manager.attribute_ids_matching_name( "p" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            pointset->nb_vertices() == 3 && ids.has_value() && ids->size() == 1,
            "Wrong point set with several VSet" );
        const auto attribute =
            manager.find_read_only_attribute< double >( ids->front() );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            attribute->value( 0 ) == 1 && attribute->value( 2 ) == 3,
            "Wrong property values with several VSet" );
    }

    template < typename T >
    std::shared_ptr< geode::VariableAttribute< T > > create_attribute(
        geode::PointSet3D& pointset, std::string_view name, T default_value )
    {
        auto& manager = pointset.vertex_attribute_manager();
        geode::AttributeValues< T > values;
        values.default_value = default_value;
        values.no_value = default_value;
        geode::AttributeProperties properties;
        properties.transferable = true;
        return manager.find_attribute< geode::VariableAttribute, T >(
            manager.create_attribute< geode::VariableAttribute, T >(
                name, values, properties ) );
    }

    template < typename T >
    T read_attribute_value( const geode::PointSet3D& pointset,
        std::string_view name,
        geode::index_t vertex )
    {
        const auto& manager = pointset.vertex_attribute_manager();
        const auto ids = manager.attribute_ids_matching_name( name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            ids.has_value() && ids->size() == 1,
            "Attribute should exist once: ", name );
        return manager.find_read_only_attribute< T >( ids->front() )
            ->value( vertex );
    }

    void test_output_properties()
    {
        auto pointset = geode::PointSet3D::create();
        auto builder = geode::PointSetBuilder3D::create( *pointset );
        builder->create_point( geode::Point3D{ { 1, 2, 3 } } );
        builder->create_point( geode::Point3D{ { 4, 5, 6 } } );
        const std::array< double, 3 > normal{ 0.1, 1e-9, 123456.789012345 };
        create_attribute< std::array< double, 3 > >(
            *pointset, "normal", { 0, 0, 0 } )
            ->set_value( 1, normal );
        create_attribute< double >( *pointset, "third", 0 )
            ->set_value( 0, 1. / 3. );
        geode::save_point_set( *pointset, "test_vs_properties.vs" );
        const auto reloaded =
            geode::load_point_set< 3 >( "test_vs_properties.vs" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            read_attribute_value< std::array< double, 3 > >(
                *reloaded, "normal", 1 )
                == normal,
            "Wrong saved vector property" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            read_attribute_value< double >( *reloaded, "third", 0 ) == 1. / 3.,
            "Wrong saved double property" );
    }

    void test_blank_lines()
    {
        const auto pointset = load_test_file( "test_vs_blank_lines.vs",
            "GOCAD VSet 1\nHEADER {\nname: a\n}\nVRTX 1 0 0 0\n\n"
            "# comment\nVRTX 2 1 1 1\n\t\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            pointset->nb_vertices() == 2,
            "Wrong number of vertices with blank lines" );
    }

    void test_quotes()
    {
        const auto filename = write_test_file( "test_vs_quotes.txt",
            "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
            "NAME \"SKUA Local\"\n"
            "DATUM \" a  b \"\n"
            "PROJECTION \"\n"
            "END_ORIGINAL_COORDINATE_SYSTEM\n" );
        std::ifstream file{ filename, std::ios::binary };
        const auto crs = geode::internal::read_CRS( file );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            crs.name == "SKUA Local", "Wrong CRS name: ", crs.name );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            crs.datum == "a b", "Wrong CRS datum: ", crs.datum );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            crs.projection.empty(), "Wrong CRS projection: ", crs.projection );

        const auto pointset = load_test_file( "test_vs_unclosed_quote.vs",
            "GOCAD VSet 1\nHEADER {\nname: a\n}\nPROPERTIES \"\n"
            "PVRTX 1 0 0 0 4\nEND\n" );
        geode::OpenGeodeGeosciencesIOMeshException::test(
            pointset->nb_vertices() == 1
                && pointset->vertex_attribute_manager()
                    .attribute_ids_matching_name( "property_0" )
                    .has_value(),
            "Wrong point set with unclosed quote" );
    }
} // namespace

void check_pointset(
    const geode::PointSet3D& pointset, geode::index_t nb_vertices )
{
    geode::OpenGeodeGeosciencesIOMeshException::test(
        pointset.nb_vertices() == nb_vertices,
        "Number of vertices in the VSet 3D is not correct" );
}

void check_file( std::string file, geode::index_t nb_vertices )
{
    // Load file
    auto pointset = geode::load_point_set< 3 >( file );
    check_pointset( *pointset, nb_vertices );

    // Save point set
    const auto output_file_native =
        absl::StrCat( "test_output.", pointset->native_extension() );
    geode::save_point_set( *pointset, output_file_native );
    const auto output_file_vs = "test_output.vs";
    geode::save_point_set( *pointset, output_file_vs );

    // Load file
    auto reloaded_pointset = geode::load_point_set< 3 >( output_file_native );
    check_pointset( *reloaded_pointset, nb_vertices );
    auto reloaded_pointset_vs = geode::load_point_set< 3 >( output_file_vs );
    check_pointset( *reloaded_pointset_vs, nb_vertices );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOMeshLibrary::initialize();
        check_file( absl::StrCat( geode::DATA_PATH, "points.",
                        geode::internal::VSInput::extension() ),
            6 );
        test_quotes();
        test_blank_lines();
        test_points_normal();
        test_multiple_vset_properties();
        test_output_properties();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
