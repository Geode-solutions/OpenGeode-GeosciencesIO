/*
 * Copyright (c) 2019 - 2020 Geode-solutions
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

#include <geode/geosciences/private/pl_output.h>

#include <fstream>
#include <memory>
#include <vector>

#include <geode/basic/attribute_manager.h>

#include <geode/basic/filename.h>
#include <geode/mesh/core/edged_curve.h>

#include <geode/geosciences/private/gocad_common.h>

namespace
{
    class PLOutputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };
        PLOutputImpl(
            absl::string_view filename, const geode::EdgedCurve3D& edged_curve )
            : file_( filename.data() ), edged_curve_( edged_curve )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void write_prop_header()
        {
            const auto names =
                edged_curve_.vertex_attribute_manager().attribute_names();
            geode::detail::PropHeaderData prop_header;
            std::vector< geode::detail::PropClassHeaderData > prop_class_header;

            for( const auto& name : names )
            {
                const auto attribute = edged_curve_.vertex_attribute_manager()
                                           .find_generic_attribute( name );
                if( !attribute || !attribute->is_genericable() )
                {
                    continue;
                }
                generic_att_.push_back( attribute );
                prop_header.names.push_back( name.data() );
                prop_header.prop_legal_ranges.emplace_back(
                     "**none**", "**none**" );
                prop_header.no_data_values.push_back( -99999. );
                prop_header.property_classes.push_back( name.data() );
                prop_header.kinds.push_back( " Real Number" );
                prop_header.property_subclass.emplace_back(
                    "QUANTITY", "Float" );
                prop_header.esizes.push_back( 1 );
                prop_header.units.push_back( "unitless" );

                geode::detail::PropClassHeaderData tmp_prop_class_header;
                tmp_prop_class_header.name = name.data();
                tmp_prop_class_header.kind = " Real Number";
                tmp_prop_class_header.unit = "unitless";
                tmp_prop_class_header.is_z = false;

                prop_class_header.push_back( tmp_prop_class_header );
            }
            if( !prop_header.empty() )
            {
                geode::detail::write_prop_header( file_, prop_header );
            }
            write_XYZ_prop_class_header();
            for( const auto& ch : prop_class_header )
            {
                write_property_class_header( file_, ch );
            }
        }

        void write_XYZ_prop_class_header()
        {
            geode::detail::PropClassHeaderData x_prop_header;
            x_prop_header.name = "X";
            x_prop_header.kind = "X";
            x_prop_header.unit = "m";
            x_prop_header.is_z = false;
            write_property_class_header( file_, x_prop_header );

            geode::detail::PropClassHeaderData y_prop_header;
            y_prop_header.name = "Y";
            y_prop_header.kind = "Y";
            y_prop_header.unit = "m";
            y_prop_header.is_z = false;
            write_property_class_header( file_, y_prop_header );

            geode::detail::PropClassHeaderData z_prop_header;
            z_prop_header.name = "Z";
            z_prop_header.kind = "Z";
            z_prop_header.unit = "m";
            z_prop_header.is_z = true;
            write_property_class_header( file_, z_prop_header );
        }

        void write_iline()
        {
            geode::index_t current_offset = OFFSET_START;
            file_ << "ILINE" << EOL;
            for( const auto v : geode::Range{ edged_curve_.nb_vertices() } )
            {
                file_ << "PVRTX" << SPACE << current_offset + v << SPACE
                      << edged_curve_.point( v ).value( 0 ) << SPACE
                      << edged_curve_.point( v ).value( 1 ) << SPACE
                      << edged_curve_.point( v ).value( 2 );
                for( const auto att : generic_att_ )
                {
                    file_ << SPACE << att->generic_value( v );
                }
                file_ << EOL;
            }

            for( const auto s : geode::Range{ edged_curve_.nb_edges() } )
            {
                file_ << "SEG" << SPACE
                      << OFFSET_START + edged_curve_.edge_vertex( { s, 0 } )
                      << SPACE
                      << OFFSET_START + edged_curve_.edge_vertex( { s, 1 } )
                      << EOL;
            }
        }

        void write_file()
        {
            geode::Logger::info( "[PLOutput::write] Do what have to be done." );

            file_ << "GOCAD PLine 1" << EOL;
            geode::detail::HeaderData header;
            header.name = "edged_curve_name";
            geode::detail::write_header( file_, header );
            geode::detail::write_CRS( file_, {} );
            write_prop_header();
            write_iline();
            file_ << "END" << EOL;
            return;
        }

    private:
        std::ofstream file_;
        const geode::EdgedCurve3D& edged_curve_;
        std::vector< std::shared_ptr< geode::AttributeBase > > generic_att_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        void PLOutput::write() const
        {
            PLOutputImpl impl{ filename(), edged_curve() };
            impl.write_file();
        }
    } // namespace detail
} // namespace geode
