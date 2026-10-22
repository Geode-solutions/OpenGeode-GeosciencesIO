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
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geosciences/explicit/mixin/core/horizon.hpp>
#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.hpp>
#include <geode/geosciences/implicit/representation/core/horizons_stack.hpp>
#include <geode/geosciences/implicit/representation/io/horizons_stack_input.hpp>
#include <geode/geosciences/implicit/representation/io/horizons_stack_output.hpp>
#include <geode/geosciences_io/model/internal/horizons_stack_skua_input.hpp>

void test_file()
{
    const auto horizons_stack = geode::load_horizons_stack< 3 >(
        absl::StrCat( geode::DATA_PATH, "test_skua_horizons_stack.xml" ) );

    geode::OpenGeodeGeosciencesIOModelException::test(
        horizons_stack.name() == "skua_model",
        "[TEST] HorizonsStack should be named 'skua_model'" );

    geode::OpenGeodeGeosciencesIOModelException::test(
        horizons_stack.nb_horizons() == 4,
        "[TEST] Wrong number of horizons in the "
        "loaded HorizonsStack." );
    geode::OpenGeodeGeosciencesIOModelException::test(
        horizons_stack.nb_stratigraphic_units() == 5,
        "[TEST] Wrong number of units in the loaded HorizonsStack." );
    geode::uuid erosion_horizon_uuid;
    geode::uuid h1_horizon_uuid;
    for( const auto& horizon : horizons_stack.horizons() )
    {
        if( horizon.name() == "model_erosion" )
        {
            erosion_horizon_uuid = horizon.id();
        }
        if( horizon.name() == "model_horizon_h1" )
        {
            h1_horizon_uuid = horizon.id();
        }
    }
    geode::uuid eroded_unit_uuid;
    for( const auto& unit : horizons_stack.stratigraphic_units() )
    {
        if( unit.name() == "eroded_unit" )
        {
            eroded_unit_uuid = unit.id();
        }
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        horizons_stack.is_above( erosion_horizon_uuid, eroded_unit_uuid ),
        "[TEST] Horizon 'model_erosion' should be above unit 'eroded_unit'" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        horizons_stack.is_above( eroded_unit_uuid, h1_horizon_uuid ),
        "[TEST] Unit 'eroded_unit' should be above horizon "
        "'model_horizon_h1'" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        horizons_stack.horizon( erosion_horizon_uuid ).contact_type()
            == geode::Horizon3D::CONTACT_TYPE::erosion,
        "[TEST] Horizon 'model_erosion' should erode unit 'eroded_unit'" );

    geode::save_horizons_stack(
        horizons_stack, "test_horizons_stack_import.og_hst3d" );
}

void write_file( std::string_view filename, std::string_view content )
{
    std::ofstream file{ geode::to_string( filename ) };
    file << content;
}

void test_is_loadable()
{
    write_file( "skua_no_indent.xml", "<UserObjects>\n"
                                      "<LocalStratigraphicColumn>\n"
                                      "</LocalStratigraphicColumn>\n"
                                      "</UserObjects>\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_horizons_stack_loadable< 3 >( "skua_no_indent.xml" ).value()
            == 1,
        "[TEST] Column without indentation should be loadable" );
    write_file( "skua_one_line.xml",
        "<UserObjects><LocalStratigraphicColumn version=\"2.0\">"
        "</LocalStratigraphicColumn></UserObjects>\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_horizons_stack_loadable< 3 >( "skua_one_line.xml" ).value()
            == 1,
        "[TEST] Column on a single line should be loadable" );
    write_file( "skua_no_column.xml", "<UserObjects>\n</UserObjects>\n" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        geode::is_horizons_stack_loadable< 3 >( "skua_no_column.xml" ).value()
            == 0,
        "[TEST] File without column should not be loadable" );
}

void test_optional_attributes()
{
    write_file( "skua_no_base_names.xml", "<UserObjects>\n"
                                          " <LocalStratigraphicColumn>\n"
                                          "  <name>column</name>\n"
                                          "  <units>\n"
                                          "   <unit>\n"
                                          "    <name>unit</name>\n"
                                          "    <top><name>top</name></top>\n"
                                          "    <base><name>base</name></base>\n"
                                          "   </unit>\n"
                                          "  </units>\n"
                                          " </LocalStratigraphicColumn>\n"
                                          "</UserObjects>\n" );
    const auto horizons_stack =
        geode::load_horizons_stack< 3 >( "skua_no_base_names.xml" );
    geode::OpenGeodeGeosciencesIOModelException::test(
        horizons_stack.nb_horizons() == 2
            && horizons_stack.nb_stratigraphic_units() >= 1,
        "[TEST] Column without use_base_names should be loaded" );
}

void test_missing_column()
{
    write_file( "skua_missing_column.xml",
        "<UserObjects>\n <Other/>\n</UserObjects>\n" );
    try
    {
        geode::load_horizons_stack< 3 >( "skua_missing_column.xml" );
    }
    catch( const geode::OpenGeodeException& )
    {
        return;
    }
    geode::OpenGeodeGeosciencesIOModelException::test(
        false, "[TEST] File without column should throw" );
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();

        geode::Logger::info( "Reading stratigraphic column file." );
        test_file();
        test_is_loadable();
        test_optional_attributes();
        test_missing_column();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
