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

#include <geode/geosciences_io/mesh/private/well_dev_input.h>

#include <fstream>

#include <geode/basic/attribute_manager.h>
#include <geode/basic/string.h>

#include <geode/geometry/point.h>

#include <geode/mesh/builder/edged_curve_builder.h>
#include <geode/mesh/core/edged_curve.h>

#include <geode/geosciences_io/mesh/private/utils.h>

namespace
{
    struct HeaderData
    {
        absl::string_view name;
        std::vector< absl::string_view > attribute_names;
        std::array< geode::index_t, 3 > xyz_attributes_position;
    };

    class WellDevInputImpl
    {
    public:
        WellDevInputImpl(
            absl::string_view filename, geode::EdgedCurve3D& curve )
            : file_{ geode::to_string( filename ) },
              curve_( curve ),
              builder_( geode::EdgedCurveBuilder3D::create( curve ) )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void read_file()
        {
            read_header();
            builder_->set_name( header_.name );
            create_attributes();
            std::string line;
            while( std::getline( file_, line ) )
            {
                const auto split_line = geode::string_split( line );
                OPENGEODE_ASSERT(
                    split_line.size() == header_.attribute_names.size() + 3,
                    "[WellDevInput::read_coord_and_attributes] Wrong number of "
                    "split_line: ",
                    split_line.size(), " - ", header_.attribute_names.size() );
                const auto point_id = create_point( split_line );
                assign_point_attributes( split_line, point_id );
            }
            for( const auto pt_id : geode::Range{ curve_.nb_vertices() - 1 } )
            {
                builder_->create_edge( pt_id, pt_id + 1 );
            }
        }

    private:
        void read_header()
        {
            geode::detail::check_keyword( file_, "# WELL TRACE FROM PETREL" );
            std::string line;
            auto first_pass = true;
            while( std::getline( file_, line ) )
            {
                if( geode::detail::string_starts_with( line, "#===" ) )
                {
                    if( first_pass )
                    {
                        first_pass = false;
                        continue;
                    }
                    else
                    {
                        return;
                    }
                }
                if( geode::detail::string_starts_with( line, "# WELL NAME:" ) )
                {
                    const auto split_line = geode::string_split( line );
                    header_.name = split_line[3];
                }
                else if( !first_pass )
                {
                    const auto split_line = geode::string_split( line );
                    OPENGEODE_EXCEPTION( split_line.size() >= 3,
                        "[WellDevInut::read_header] There are less than 3 "
                        "attributes given for the well" );
                    header_.attribute_names.reserve( split_line.size() - 3 );
                    for( const auto i : geode::Indices{ split_line } )
                    {
                        if( split_line[i] == "X" )
                        {
                            header_.xyz_attributes_position[0] = i;
                            continue;
                        }
                        else if( split_line[i] == "Y" )
                        {
                            header_.xyz_attributes_position[1] = i;
                            continue;
                        }
                        else if( split_line[i] == "Z" )
                        {
                            header_.xyz_attributes_position[2] = i;
                            continue;
                        }
                        header_.attribute_names.push_back( split_line[i] );
                    }
                    OPENGEODE_EXCEPTION(
                        header_.xyz_attributes_position[0]
                                != header_.xyz_attributes_position[1]
                            && header_.xyz_attributes_position[0]
                                   != header_.xyz_attributes_position[2]
                            && header_.xyz_attributes_position[1]
                                   != header_.xyz_attributes_position[2],
                        "[WellDevInput::header] Cannot find the X, Y and Z "
                        "point "
                        "position attributes in the header." );
                }
            }
            throw geode::OpenGeodeException{
                "[read_header] Cannot find the end of \"HEADER\" section"
            };
        }

        void create_attributes()
        {
            attributes_.reserve( header_.attribute_names.size() );
            for( const auto attr_name : header_.attribute_names )
            {
                attributes_.push_back(
                    curve_.vertex_attribute_manager()
                        .template find_or_create_attribute<
                            geode::VariableAttribute, double >(
                            attr_name, 0 ) );
            }
        }

        geode::index_t create_point(
            absl::Span< const absl::string_view > split_line )
        {
            return builder_->create_point(
                { { geode::string_to_double(
                        split_line[header_.xyz_attributes_position[0]] ),
                    geode::string_to_double(
                        split_line[header_.xyz_attributes_position[1]] ),
                    geode::string_to_double(
                        split_line[header_.xyz_attributes_position[2]] ) } } );
        }

        void assign_point_attributes(
            absl::Span< const absl::string_view > split_line,
            geode::index_t point_id )
        {
            geode::index_t attr_counter{ 0 };
            for( const auto line_object_position :
                geode::Indices{ split_line } )
            {
                if( line_object_position == header_.xyz_attributes_position[0]
                    || line_object_position
                           == header_.xyz_attributes_position[1]
                    || line_object_position
                           == header_.xyz_attributes_position[2] )
                {
                    continue;
                }
                attributes_[attr_counter]->set_value(
                    point_id, geode::string_to_double(
                                  split_line[line_object_position] ) );
                attr_counter++;
            }
        }

    private:
        std::ifstream file_;
        geode::EdgedCurve3D& curve_;
        std::unique_ptr< geode::EdgedCurveBuilder3D > builder_;
        HeaderData header_;
        std::vector< std::shared_ptr< geode::VariableAttribute< double > > >
            attributes_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        std::unique_ptr< EdgedCurve3D > WellDevInput::read(
            const MeshImpl& impl )
        {
            auto well = EdgedCurve3D::create( impl );
            WellDevInputImpl reader{ this->filename(), *well };
            reader.read_file();
            return well;
        }
    } // namespace detail
} // namespace geode
