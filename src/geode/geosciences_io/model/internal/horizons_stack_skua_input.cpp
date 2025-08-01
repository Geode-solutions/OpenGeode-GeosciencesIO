/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

#include <geode/geosciences_io/model/internal/horizons_stack_skua_input.hpp>

#include <fstream>

#include <pugixml.hpp>

#include <absl/strings/str_split.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/file.hpp>

#include <geode/geosciences/explicit/mixin/core/horizon.hpp>
#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.hpp>
#include <geode/geosciences/explicit/representation/io/structural_model_input.hpp>

#include <geode/geosciences/implicit/representation/builder/horizons_stack_builder.hpp>
#include <geode/geosciences/implicit/representation/core/detail/helpers.hpp>
#include <geode/geosciences/implicit/representation/core/horizons_stack.hpp>

namespace
{
    template < geode::index_t dimension >
    class HorizonStackSKUAInputImpl
    {
        using CONTACT_TYPE = typename geode::Horizon< dimension >::CONTACT_TYPE;

    public:
        HorizonStackSKUAInputImpl( std::string_view filename )
            : filename_( filename )
        {
        }

        geode::HorizonsStack< dimension > read_file()
        {
            load_file();
            geode::HorizonsStack< dimension > horizons_stack;
            geode::HorizonsStackBuilder< dimension > builder{ horizons_stack };
            builder.set_name( absl::StripAsciiWhitespace(
                root_.child( "name" ).child_value() ) );
            const auto& units = root_.child( "units" );
            bool use_base_names;
            const auto ok = absl::SimpleAtob(
                units.attribute( "use_base_names" ).value(), &use_base_names );
            OPENGEODE_EXCEPTION( ok,
                "[HorizonStackSKUAInput::read_units] Failed to "
                "read use_base_names attribute." );
            absl::flat_hash_map< std::string, geode::uuid > name_map;
            for( const auto& unit : units.children( "unit" ) )
            {
                const auto unit_name = absl::StripAsciiWhitespace(
                    unit.child( "name" ).child_value() );
                const auto& unit_uuid = builder.add_stratigraphic_unit();
                const auto& strati_unit =
                    horizons_stack.stratigraphic_unit( unit_uuid );
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
                    const auto& top_horizon = horizons_stack.horizon(
                        name_map.at( top_horizon_name ) );
                    builder.set_horizon_above( top_horizon, strati_unit );
                    if( top_relation != "structural" )
                    {
                        continue;
                    }
                    if( top_horizon.contact_type() == CONTACT_TYPE::conformal )
                    {
                        builder.set_horizon_contact_type(
                            top_horizon.id(), CONTACT_TYPE::erosion );
                    }
                    else
                    {
                        builder.set_horizon_contact_type(
                            top_horizon.id(), CONTACT_TYPE::discontinuity );
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
                    const auto& base_horizon = horizons_stack.horizon(
                        name_map.at( base_horizon_name ) );
                    builder.set_horizon_under( base_horizon, strati_unit );
                    if( base_relation != "structural" )
                    {
                        continue;
                    }
                    if( base_horizon.contact_type() == CONTACT_TYPE::conformal )
                    {
                        builder.set_horizon_contact_type(
                            base_horizon.id(), CONTACT_TYPE::baselap );
                    }
                    else
                    {
                        builder.set_horizon_contact_type(
                            base_horizon.id(), CONTACT_TYPE::discontinuity );
                    }
                }
            }
            geode::detail::repair_horizon_stack_if_possible(
                horizons_stack, builder );
            return horizons_stack;
        }

    private:
        void load_file()
        {
            std::ifstream file{ geode::to_string( filename_ ) };
            OPENGEODE_EXCEPTION( file.good(),
                "[HorizonStackSKUAInput] Error while opening file: ",
                filename_ );
            const auto status =
                document_.load_file( geode::to_string( filename_ ).c_str() );
            OPENGEODE_EXCEPTION( status, "[HorizonStackSKUAInput] Error ",
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
    namespace internal
    {
        template < index_t dimension >
        HorizonsStack< dimension > HorizonStackSKUAInput< dimension >::read()
        {
            HorizonStackSKUAInputImpl< dimension > impl{ this->filename() };
            return impl.read_file();
        }

        template < index_t dimension >
        Percentage HorizonStackSKUAInput< dimension >::is_loadable() const
        {
            std::ifstream file{ geode::to_string( this->filename() ) };
            if( goto_keyword_if_it_exists(
                    file, " <LocalStratigraphicColumn" ) )
            {
                return Percentage{ 1 };
            }
            return Percentage{ 0 };
        }

        template class HorizonStackSKUAInput< 2 >;
        template class HorizonStackSKUAInput< 3 >;
    } // namespace internal
} // namespace geode
