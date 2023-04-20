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

#include <geode/geosciences_io/model/private/su_stack_skua_input.h>

#include <fstream>

#include <pugixml.hpp>

#include <absl/strings/str_split.h>

#include <geode/basic/attribute_manager.h>

#include <geode/geosciences/explicit/mixin/core/horizon.h>
#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.h>
#include <geode/geosciences/explicit/representation/io/structural_model_input.h>

#include <geode/geosciences/implicit/representation/builder/stratigraphic_units_stack_builder.h>
#include <geode/geosciences/implicit/representation/core/stratigraphic_units_stack.h>

namespace
{
    template < geode::index_t dimension >
    class SUStackSKUAInputImpl
    {
    public:
        SUStackSKUAInputImpl( absl::string_view filename )
            : filename_( filename )
        {
        }

        geode::StratigraphicUnitsStack< dimension > read_file()
        {
            load_file();
            geode::StratigraphicUnitsStack< dimension > su_stack;
            geode::StratigraphicUnitsStackBuilder< dimension > builder{
                su_stack
            };
            builder.set_name( absl::StripAsciiWhitespace(
                root_.child( "name" ).child_value() ) );
            const auto& units = root_.child( "units" );
            bool use_base_names;
            const auto ok = absl::SimpleAtob(
                units.attribute( "use_base_names" ).value(), &use_base_names );
            OPENGEODE_EXCEPTION( ok, "[SUStackSKUAInput::read_units] Failed to "
                                     "read use_base_names attribute." );
            absl::flat_hash_map< std::string, geode::uuid > name_map;
            for( const auto& unit : units.children( "unit" ) )
            {
                const auto unit_name = absl::StripAsciiWhitespace(
                    unit.child( "name" ).child_value() );
                const auto& unit_uuid = builder.add_stratigraphic_unit();
                const auto& strati_unit =
                    su_stack.stratigraphic_unit( unit_uuid );
                builder.set_stratigraphic_unit_name( unit_uuid, unit_name );
                for( const auto& unit_top : unit.children( "top" ) )
                {
                    const auto top_horizon_name = absl::StripAsciiWhitespace(
                        unit_top.child( "name" ).child_value() );
                    const auto top_relation = absl::StripAsciiWhitespace(
                        unit_top.child( "relation" ).child_value() );
                    if( !name_map.contains( top_horizon_name ) )
                    {
                        const auto& top_horizon_uuid = builder.add_horizon();
                        builder.set_horizon_name(
                            top_horizon_uuid, top_horizon_name );
                        name_map[top_horizon_name] = top_horizon_uuid;
                    }
                    const auto& top_horizon =
                        su_stack.horizon( name_map.at( top_horizon_name ) );
                    builder.add_horizon_above( top_horizon, strati_unit );
                    if( top_relation == "structural" )
                    {
                        builder.add_erosion_above( top_horizon, strati_unit );
                    }
                }
                for( const auto& unit_base : unit.children( "base" ) )
                {
                    const auto base_horizon_name = absl::StripAsciiWhitespace(
                        unit_base.child( "name" ).child_value() );
                    const auto base_relation = absl::StripAsciiWhitespace(
                        unit_base.child( "relation" ).child_value() );
                    if( !name_map.contains( base_horizon_name ) )
                    {
                        const auto& base_horizon_uuid = builder.add_horizon();
                        builder.set_horizon_name(
                            base_horizon_uuid, base_horizon_name );
                        name_map[base_horizon_name] = base_horizon_uuid;
                    }
                    const auto& base_horizon =
                        su_stack.horizon( name_map.at( base_horizon_name ) );
                    builder.add_horizon_under( base_horizon, strati_unit );
                    if( base_relation == "structural" )
                    {
                        builder.add_baselap_under( base_horizon, strati_unit );
                    }
                }
            }
            return su_stack;
        }

    private:
        void load_file()
        {
            std::ifstream file{ geode::to_string( filename_ ) };
            OPENGEODE_EXCEPTION( file.good(),
                "[SUStackSKUAInput] Error while opening file: ", filename_ );
            const auto status =
                document_.load_file( geode::to_string( filename_ ).c_str() );
            OPENGEODE_EXCEPTION( status, "[SUStackSKUAInput] Error ",
                status.description(), " while parsing file: ", filename_ );
            root_ = document_.child( "UserObjects" )
                        .child( "LocalStratigraphicColumn" );
        }

    private:
        std::string filename_;
        pugi::xml_document document_;
        pugi::xml_node root_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        template < index_t dimension >
        StratigraphicUnitsStack< dimension >
            SUStackSKUAInput< dimension >::read()
        {
            SUStackSKUAInputImpl< dimension > impl{ this->filename() };
            return impl.read_file();
        }

        template class SUStackSKUAInput< 2 >;
        template class SUStackSKUAInput< 3 >;
    } // namespace detail
} // namespace geode
