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

#include <geode/geosciences_io/mesh/private/vo_input.h>

#include <fstream>
#include <string>

#include <absl/strings/str_replace.h>
#include <absl/types/optional.h>

#include <geode/basic/attribute_manager.h>
#include <geode/basic/file.h>
#include <geode/basic/filename.h>
#include <geode/basic/string.h>

#include <geode/geometry/distance.h>
#include <geode/geometry/vector.h>

#include <geode/mesh/builder/regular_grid_solid_builder.h>
#include <geode/mesh/core/regular_grid_solid.h>
#include <geode/mesh/io/regular_grid_input.h>

#include <geode/geosciences_io/mesh/private/gocad_common.h>

namespace
{

    absl::optional< std::string > get_data_file( std::ifstream& file )
    {
        const auto line =
            geode::goto_keyword_if_it_exists( file, "ASCII_DATA_FILE" );
        if( !line.has_value() )
        {
            return absl::nullopt;
        }
        return absl::StrReplaceAll(
            line.value(), { { "ASCII_DATA_FILE ", "" }, { "\"", "" } } );
    }

    class VOInputImpl
    {
    public:
        VOInputImpl( absl::string_view filename, geode::RegularGrid3D& grid )
            : file_{ geode::to_string( filename ) },
              file_folder_{ geode::filepath_without_filename( filename ) },
              grid_( grid ),
              builder_{ geode::RegularGridBuilder3D::create( grid ) }
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void read_file()
        {
            if( !geode::goto_keyword_if_it_exists( file_, "GOCAD Voxet" ) )
            {
                throw geode::OpenGeodeException{
                    "[VOInput] Cannot find Voxet in the file"
                };
            }
            const auto header = geode::detail::read_header( file_ );
            builder_->set_name( header.name );
            geode::detail::read_CRS( file_ );
            initialize_grid();
            read_data_file();
        }

    private:
        void initialize_grid()
        {
            auto line = geode::goto_keyword( file_, "AXIS_O" );
            auto origin = read_coord( line, 1 );
            const auto grid_size = read_grid_size( origin );
            auto cells_number = read_cells_number();
            auto cells_length = compute_cells_length( grid_size, cells_number );

            builder_->initialize_grid( std::move( origin ),
                std::move( cells_number ), std::move( cells_length ) );
        }

        std::array< double, 3 > read_grid_size( const geode::Point3D& origin )
        {
            std::array< double, 3 > cells_length;
            auto line = geode::goto_keyword( file_, "AXIS_U" );
            cells_length[0] =
                geode::point_point_distance( origin, read_coord( line, 1 ) );
            line = geode::goto_keyword( file_, "AXIS_V" );
            cells_length[1] =
                geode::point_point_distance( origin, read_coord( line, 1 ) );
            line = geode::goto_keyword( file_, "AXIS_W" );
            cells_length[2] =
                geode::point_point_distance( origin, read_coord( line, 1 ) );
            return cells_length;
        }

        std::array< geode::index_t, 3 > read_cells_number()
        {
            auto line = geode::goto_keyword( file_, "AXIS_N" );
            const auto tokens = geode::string_split( line );
            return { geode::string_to_index( tokens[1] ),
                geode::string_to_index( tokens[2] ),
                geode::string_to_index( tokens[3] ) };
        }

        std::array< double, 3 > compute_cells_length(
            const std::array< double, 3 >& grid_size,
            const std::array< geode::index_t, 3 >& nb_cells )
        {
            std::array< double, 3 > result;
            for( const auto axis_id : geode::LRange{ 3 } )
            {
                result[axis_id] = grid_size[axis_id] / nb_cells[axis_id];
            }
            return result;
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

        void read_data_file()
        {
            const auto data_file_name = get_data_file( file_ );
            OPENGEODE_EXCEPTION(
                data_file_name.has_value(), "[VOInput] No data file record" );
            std::ifstream data_file{ absl::StrCat(
                file_folder_, data_file_name.value() ) };
            OPENGEODE_EXCEPTION( data_file.good(),
                "[VOInput] Cannot open data file: ", data_file_name.value() );
            std::string line;
            std::getline( data_file, line );
            std::getline( data_file, line );
            auto tokens = geode::string_split( line );
            absl::FixedArray<
                std::shared_ptr< geode::VariableAttribute< double > > >
                data_attributes( tokens.size() - 4 );
            for( const auto attribute_id : geode::Indices{ data_attributes } )
            {
                data_attributes[attribute_id] =
                    grid_.cell_attribute_manager()
                        .find_or_create_attribute< geode::VariableAttribute,
                            double >( tokens[4 + attribute_id], 0 );
            }
            std::getline( data_file, line );
            while( std::getline( data_file, line ) )
            {
                tokens = geode::string_split( line );
                OPENGEODE_ASSERT( tokens.size() == data_attributes.size() + 3,
                    "[VOInput::read_data_file] Wrong number of tokens in line, "
                    "got",
                    tokens.size(), ", should have ",
                    data_attributes.size() + 3 );
                const auto cell_id =
                    grid_.cell_index( { geode::string_to_index( tokens[0] ),
                        geode::string_to_index( tokens[1] ),
                        geode::string_to_index( tokens[2] ) } );
                for( const auto attribute_id :
                    geode::Indices{ data_attributes } )
                {
                    data_attributes[attribute_id]->set_value( cell_id,
                        geode::string_to_double( tokens[3 + attribute_id] ) );
                }
            }
        }

    private:
        std::ifstream file_;
        std::string file_folder_;
        geode::RegularGrid3D& grid_;
        std::unique_ptr< geode::RegularGridBuilder3D > builder_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        std::unique_ptr< RegularGrid3D > VOInput::read( const MeshImpl& impl )
        {
            auto voxet = RegularGrid3D::create( impl );
            VOInputImpl reader{ filename(), *voxet };
            reader.read_file();
            return voxet;
        }

        RegularGridInput< 3 >::MissingFiles VOInput::check_missing_files() const
        {
            std::ifstream file{ geode::to_string( filename() ) };
            const auto data_file = get_data_file( file );
            file.close();
            if( !data_file.has_value() )
            {
                return {};
            }
            RegularGridInput< 3 >::MissingFiles missing;
            if( !file_exists( data_file.value() ) )
            {
                missing.mandatory_files.push_back( data_file.value() );
            }
            return missing;
        }
    } // namespace detail
} // namespace geode
