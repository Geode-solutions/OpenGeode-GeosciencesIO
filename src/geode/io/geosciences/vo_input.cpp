/*
 * Copyright (c) 2019 - 2022 Geode-solutions
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

#include <absl/strings/str_replace.h>

#include <geode/io/geosciences/private/vo_input.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include <geode/basic/string.h>

#include <geode/geometry/point.h>
#include <geode/basic/attribute_manager.h>

#include <geode/mesh/builder/regular_grid_solid_builder.h>

#include <geode/mesh/core/regular_grid_solid.h>

#include <geode/io/geosciences/private/gocad_common.h>
namespace
{
    class VOInputImpl
    {
    public:
        VOInputImpl( absl::string_view filename, geode::RegularGrid3D& reg_grid )
            : file_{ geode::to_string( filename ) },
              data_file_{ geode::to_string( absl::StrReplaceAll(filename, {{".vo", "__ascii@@"}}) ) },
              reg_grid_( reg_grid ),
              builder_{ geode::RegularGridBuilder3D::create( reg_grid ) }
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename ); 
            OPENGEODE_EXCEPTION(
                data_file_.good(), "Error while opening data file: ", absl::StrReplaceAll(filename, {{".vo", "__ascii@@"}}));
        }

        void read_file()
        {
            std::string line;
            std::getline( file_, line );
            if( !geode::detail::string_starts_with( line, "GOCAD Voxet" ) )
            {
                return;
            }
            const auto header = geode::detail::read_header( file_ );
            builder_->set_name( header.name );
            crs_ = geode::detail::read_CRS( file_ );


            geode::Point3D origin = read_origin();
            std::array< double, 3 > cells_length  = read_cells_length (origin);
            std::array< geode::index_t, 3 > cells_number  = read_cells_number ();
            cells_length  = calculate_cells_length (cells_length, cells_number);
   
            builder_->initialize_grid( origin, cells_number, cells_length );
            create_attributes(); 
            init_points();
        }

    private:
        geode::Point3D read_origin()
        {
            geode::Point3D point;
            auto line = geode::detail::goto_keyword( file_, "AXIS_O" );
            if( geode::detail::string_starts_with( line, "AXIS_O" ) )
            {
                auto translation = read_coord( line, 1 );
                point.set_value( 0, translation.value( 1 ));
                point.set_value( 1, translation.value( 2 ));
                point.set_value( 2, translation.value( 2 ));
            }
            return point;
        }
        
        std::array< double, 3 > calculate_cells_length(std::array< double, 3 > cells_length, std::array< geode::index_t, 3 > cells_number)
        {
            cells_length[0] = cells_length[0] / cells_number[0];
            cells_length[1] = cells_length[1] / cells_number[1];
            cells_length[2] = cells_length[2] / cells_number[2];
            return cells_length;
        }

        std::array< geode::index_t, 3 > read_cells_number()
        {
            std::array< geode::index_t, 3 > cells_number;
            auto line = geode::detail::goto_keyword( file_, "AXIS_N" );
            if(geode::detail::string_starts_with( line, "AXIS_N" ) )
            {
                auto translation = read_coord( line, 1 );
                cells_number[0] = translation.value( 0 );
                cells_number[1] = translation.value( 1 );
                cells_number[2] = translation.value( 2 );
            }
            return cells_number;
        }

        std::array< double, 3 > read_cells_length(geode::Point3D origin)
        {
            std::array< double, 3 > cells_length;
            auto line = geode::detail::goto_keyword( file_, "AXIS_U" );
            if( geode::detail::string_starts_with( line, "AXIS_U" ) )
            {
                auto translation = read_coord( line, 1 );
                cells_length[0] = sqrt(pow((translation.value( 0 ) - origin.value(0)),2)
                    + pow((translation.value( 1 ) - origin.value(1)),2)
                    + pow((translation.value( 2 ) - origin.value(0)),2));
            }
            line = geode::detail::goto_keyword( file_, "AXIS_V" );
            if( geode::detail::string_starts_with( line, "AXIS_V" ) )
            {
                auto translation = read_coord( line, 1 );
                cells_length[1] = sqrt(pow((translation.value( 0 ) - origin.value(0)),2)
                    + pow((translation.value( 1 ) - origin.value(1)),2)
                    + pow((translation.value( 2 ) - origin.value(0)),2));
            }
            line = geode::detail::goto_keyword( file_, "AXIS_W" );
            if( geode::detail::string_starts_with( line, "AXIS_W" ) )
            {
                auto translation = read_coord( line, 1 );
                cells_length[2] = sqrt(pow((translation.value( 0 ) - origin.value(0)),2)
                    + pow((translation.value( 1 ) - origin.value(1)),2)
                    + pow((translation.value( 2 ) - origin.value(0)),2));
            }
            return cells_length;
        }

        geode::Point3D read_coord(
            absl::string_view line, geode::index_t offset ) const
        {
            const auto tokens = geode::string_split( line );
            OPENGEODE_ASSERT( tokens.size() == 3 + offset,
                "[VOInput::read_coord] Wrong number of tokens" );
            return { { geode::string_to_double( tokens[offset] ),
                geode::string_to_double( tokens[1 + offset] ),
                geode::string_to_double( tokens[2 + offset] ) } };
        }

        void create_attributes()
        {
            std::string line;
            std::getline( data_file_, line );
            std::getline( data_file_, line );
            const auto tokens = geode::string_split( line );
            //attributes_.reserve( tokens.size() - 4 );
            int i = 4;
            while((unsigned)i < tokens.size())
            {
                const auto attr_name = tokens[i];
                attributes_.push_back(
                    reg_grid_.cell_attribute_manager()
                        .template find_or_create_attribute<
                            geode::VariableAttribute, double >(
                            attr_name, 0 ) );
                i = i + 1;
            }
        }

        void init_points() {
            std::string line;
            std::getline( data_file_, line );
            while( std::getline( data_file_, line ) )
            {
                const auto split_line = geode::string_split( line );
                OPENGEODE_ASSERT(
                    split_line.size() == attributes_.size() + 3,
                    "[VOInput::create_points] Wrong number of "
                    "split_line: ",
                    split_line.size(), " - ", attributes_.size() );
                const auto point_id = get_point( split_line );
                assign_point_attributes( split_line, point_id );
            }
        }
        
        geode::index_t get_point(
            absl::Span< const absl::string_view > split_line )
        {
            return reg_grid_.cell_index({
                (unsigned int)geode::string_to_double(split_line[0]),
                (unsigned int)geode::string_to_double(split_line[1]),
                (unsigned int)geode::string_to_double(split_line[2])});
        }

        void assign_point_attributes(
            absl::Span< const absl::string_view > split_line,
            geode::index_t point_id )
        {
            geode::index_t attr_counter{ 0 };
            for( const auto line_object_position :
                geode::Indices{ split_line } )
            {
                if( line_object_position == 0
                    || line_object_position
                           == 1
                    || line_object_position
                           == 2 )
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
        std::ifstream data_file_;
        geode::RegularGrid3D& reg_grid_;
        std::unique_ptr< geode::RegularGridBuilder3D > builder_;
        geode::detail::CRSData crs_;
        std::vector< std::shared_ptr< geode::VariableAttribute< double > > >
            attributes_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        std::unique_ptr< RegularGrid3D > VOInput::read( const MeshImpl& impl )
        {
            auto voxet = RegularGrid3D::create( impl );
            VOInputImpl reader{ this->filename(), *voxet };
            reader.read_file();
            return voxet;
        }
    } // namespace detail
} // namespace geode
