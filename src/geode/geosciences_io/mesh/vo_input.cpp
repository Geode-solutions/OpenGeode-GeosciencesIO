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

#include <geode/geosciences_io/mesh/internal/vo_input.hpp>

#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

#include <absl/strings/match.h>
#include <absl/strings/str_replace.h>
#include <absl/strings/strip.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/file.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/string.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/vector.hpp>

#include <geode/mesh/builder/regular_grid_solid_builder.hpp>
#include <geode/mesh/core/regular_grid_solid.hpp>
#include <geode/mesh/io/regular_grid_input.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

namespace
{

    std::string read_data_file_name( std::string_view line )
    {
        return absl::StrReplaceAll( absl::StripAsciiWhitespace( line ),
            { { "ASCII_DATA_FILE ", "" }, { "\"", "" } } );
    }

    std::string data_file_path(
        std::string_view voxet_filename, std::string_view data_file )
    {
        const std::filesystem::path data_path{ geode::to_string( data_file ) };
        if( data_path.is_absolute() )
        {
            return data_path.string();
        }
        return (
            geode::filepath_without_filename( voxet_filename ) / data_path )
            .string();
    }

    std::optional< std::string > get_data_file( std::ifstream& file )
    {
        const auto line =
            geode::goto_keyword_if_it_exists( file, "ASCII_DATA_FILE" );
        if( !line.has_value() )
        {
            return std::nullopt;
        }
        return read_data_file_name( line.value() );
    }

    class VOInputImpl
    {
    public:
        VOInputImpl( std::string_view filename, geode::RegularGrid3D& grid )
            : file_{ geode::to_string( filename ), std::ios::binary },
              filename_{ filename },
              grid_( grid ),
              builder_{ geode::RegularGridBuilder3D::create( grid ) }
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                file_.good(), nullptr, geode::OpenGeodeException::TYPE::data,
                "Error while opening file: ", filename );
        }

        void read_file()
        {
            if( !geode::goto_keyword_if_it_exists( file_, "GOCAD Voxet" ) )
            {
                throw geode::OpenGeodeGeosciencesIOMeshException{ nullptr,
                    geode::OpenGeodeException::TYPE::data,
                    "[VOInput] Cannot find Voxet in the file" };
            }
            const auto header = geode::internal::read_header( file_ );
            if( header.name )
            {
                builder_->set_name( header.name.value() );
            }
            read_voxet_header();
            initialize_grid();
            read_data_file();
        }

    private:
        void read_voxet_header()
        {
            std::string line;
            auto position = file_.tellg();
            while( std::getline( file_, line ) )
            {
                const auto tokens = geode::string_split( line );
                if( tokens.empty() )
                {
                    position = file_.tellg();
                    continue;
                }
                const auto keyword = tokens.front();
                if( keyword == "END" )
                {
                    break;
                }
                if( keyword == "GOCAD_ORIGINAL_COORDINATE_SYSTEM" )
                {
                    file_.seekg( position );
                    crs_ = geode::internal::read_CRS( file_ );
                }
                else if( keyword == "AXIS_O" )
                {
                    origin_ = read_coord( tokens );
                }
                else if( keyword == "AXIS_U" )
                {
                    axes_[0] = geode::Vector3D{ read_coord( tokens ) };
                }
                else if( keyword == "AXIS_V" )
                {
                    axes_[1] = geode::Vector3D{ read_coord( tokens ) };
                }
                else if( keyword == "AXIS_W" )
                {
                    axes_[2] = geode::Vector3D{ read_coord( tokens ) };
                }
                else if( keyword == "AXIS_MIN" )
                {
                    axis_min_ = read_coord( tokens );
                }
                else if( keyword == "AXIS_MAX" )
                {
                    axis_max_ = read_coord( tokens );
                }
                else if( keyword == "AXIS_N" )
                {
                    nb_nodes_ = read_nb_nodes( tokens );
                }
                else if( keyword == "ASCII_DATA_FILE" )
                {
                    data_file_ = read_data_file_name( line );
                }
                else if( tokens.back().back() == '{' )
                {
                    skip_block();
                }
                position = file_.tellg();
            }
            file_.clear();
        }

        void skip_block()
        {
            std::string line;
            while( std::getline( file_, line ) )
            {
                if( absl::StartsWith(
                        absl::StripLeadingAsciiWhitespace( line ), "}" ) )
                {
                    return;
                }
            }
        }

        void initialize_grid()
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                origin_ && axes_[0] && axes_[1] && axes_[2] && nb_nodes_,
                nullptr, geode::OpenGeodeException::TYPE::data,
                "[VOInput] Missing AXIS_O, AXIS_U, AXIS_V, AXIS_W or AXIS_N" );
            const auto z_sign = crs_.z_sign_positive ? 1. : -1.;
            // Voxet values are given at nodes: each node is the center of a
            // grid cell.
            auto origin = origin_.value();
            origin.set_value( 2, z_sign * origin.value( 2 ) );
            std::array< geode::Vector3D, 3 > cell_directions;
            for( const auto axis : geode::LRange{ 3 } )
            {
                auto direction = axes_[axis].value();
                direction.set_value( 2, z_sign * direction.value( 2 ) );
                const auto nb_nodes = nb_nodes_.value()[axis];
                const auto range =
                    axis_max_.value( axis ) - axis_min_.value( axis );
                const auto step =
                    nb_nodes > 1 ? range / ( nb_nodes - 1 ) : range;
                cell_directions[axis] = direction * step;
                geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                    cell_directions[axis].length() > geode::GLOBAL_EPSILON,
                    nullptr, geode::OpenGeodeException::TYPE::data,
                    "[VOInput] Degenerated Voxet axis ", axis );
                origin = origin + direction * axis_min_.value( axis )
                         - cell_directions[axis] / 2.;
            }
            builder_->initialize_grid( std::move( origin ), nb_nodes_.value(),
                std::move( cell_directions ) );
        }

        std::array< geode::index_t, 3 > read_nb_nodes(
            absl::Span< const std::string_view > tokens ) const
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                tokens.size() == 4, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[VOInput] AXIS_N should have 3 values" );
            std::array< geode::index_t, 3 > nb_nodes;
            for( const auto axis : geode::LRange{ 3 } )
            {
                nb_nodes[axis] = geode::string_to_index( tokens[axis + 1] );
                geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                    nb_nodes[axis] > 0, nullptr,
                    geode::OpenGeodeException::TYPE::data,
                    "[VOInput] AXIS_N values should be positive" );
            }
            return nb_nodes;
        }

        geode::Point3D read_coord(
            absl::Span< const std::string_view > tokens ) const
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                tokens.size() == 4, nullptr,
                geode::OpenGeodeException::TYPE::data, "[VOInput] ",
                tokens.front(), " should have 3 values" );
            return geode::Point3D{ { geode::string_to_double( tokens[1] ),
                geode::string_to_double( tokens[2] ),
                geode::string_to_double( tokens[3] ) } };
        }

        void read_data_file()
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                data_file_.has_value(), nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[VOInput] No data file record" );
            const auto data_path =
                data_file_path( filename_, data_file_.value() );
            std::ifstream data_file{ data_path, std::ios::binary };
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                data_file.good(), nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[VOInput] Cannot open data file: ", data_path );
            std::string line;
            std::vector< std::shared_ptr< geode::VariableAttribute< double > > >
                data_attributes;
            auto has_column_header = false;
            geode::index_t nb_outside_cells{ 0 };
            geode::index_t nb_wrong_lines{ 0 };
            while( std::getline( data_file, line ) )
            {
                const auto tokens = geode::string_split( line );
                if( tokens.empty() )
                {
                    continue;
                }
                if( tokens.front().front() == '*'
                    || tokens.front().front() == '#' )
                {
                    if( !has_column_header && tokens.size() >= 4
                        && tokens[1] == "I" )
                    {
                        has_column_header = true;
                        data_attributes = create_data_attributes( tokens );
                    }
                    continue;
                }
                geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                    has_column_header, nullptr,
                    geode::OpenGeodeException::TYPE::data,
                    "[VOInput::read_data_file] Missing \"* I J K ...\" "
                    "column header in data file" );
                if( tokens.size() != data_attributes.size() + 3 )
                {
                    nb_wrong_lines++;
                    continue;
                }
                const std::array< geode::index_t, 3 > cell{
                    geode::string_to_index( tokens[0] ),
                    geode::string_to_index( tokens[1] ),
                    geode::string_to_index( tokens[2] )
                };
                if( !is_cell_inside_grid( cell ) )
                {
                    nb_outside_cells++;
                    continue;
                }
                const auto cell_id = grid_.cell_index( cell );
                for( const auto attribute_id :
                    geode::Indices{ data_attributes } )
                {
                    data_attributes[attribute_id]->set_value( cell_id,
                        geode::string_to_double( tokens[3 + attribute_id] ) );
                }
            }
            if( nb_wrong_lines > 0 )
            {
                geode::Logger::warning( "[VOInput::read_data_file] ",
                    nb_wrong_lines, " data lines do not have ",
                    data_attributes.size() + 3, " values, they are skipped" );
            }
            if( nb_outside_cells > 0 )
            {
                geode::Logger::warning( "[VOInput::read_data_file] ",
                    nb_outside_cells,
                    " data lines reference cells outside of the grid, they "
                    "are skipped" );
            }
        }

        std::vector< std::shared_ptr< geode::VariableAttribute< double > > >
            create_data_attributes(
                absl::Span< const std::string_view > tokens )
        {
            std::vector< std::shared_ptr< geode::VariableAttribute< double > > >
                data_attributes;
            for( const auto token_id : geode::Range{ 4, tokens.size() } )
            {
                geode::AttributeValues< double > attribute_values;
                attribute_values.default_value = 0;
                attribute_values.no_value = 0;
                geode::AttributeProperties attribute_properties;
                attribute_properties.assignable = false;
                attribute_properties.interpolable = false;
                attribute_properties.transferable = true;
                const auto attribute_id =
                    grid_.cell_attribute_manager()
                        .create_attribute< geode::VariableAttribute, double >(
                            tokens[token_id], attribute_values,
                            attribute_properties );
                data_attributes.push_back( grid_.cell_attribute_manager()
                        .find_attribute< geode::VariableAttribute, double >(
                            attribute_id ) );
            }
            return data_attributes;
        }

        bool is_cell_inside_grid(
            const std::array< geode::index_t, 3 >& cell ) const
        {
            for( const auto axis : geode::LRange{ 3 } )
            {
                if( cell[axis] >= grid_.nb_cells_in_direction( axis ) )
                {
                    return false;
                }
            }
            return true;
        }

    private:
        std::ifstream file_;
        std::string filename_;
        geode::RegularGrid3D& grid_;
        std::unique_ptr< geode::RegularGridBuilder3D > builder_;
        geode::internal::CRSData crs_;
        std::optional< geode::Point3D > origin_;
        std::array< std::optional< geode::Vector3D >, 3 > axes_;
        geode::Point3D axis_min_{ { 0, 0, 0 } };
        geode::Point3D axis_max_{ { 1, 1, 1 } };
        std::optional< std::array< geode::index_t, 3 > > nb_nodes_;
        std::optional< std::string > data_file_;
    };
} // namespace

namespace geode::internal
{
    std::unique_ptr< RegularGrid3D > VOInput::read( const MeshImpl& impl )
    {
        auto voxet = RegularGrid3D::create( impl );
        VOInputImpl reader{ filename(), *voxet };
        reader.read_file();
        return voxet;
    }

    auto VOInput::additional_files() const -> AdditionalFiles
    {
        std::ifstream file{ geode::to_string( filename() ), std::ios::binary };
        const auto data_file = get_data_file( file );
        file.close();
        if( !data_file.has_value() )
        {
            return {};
        }
        AdditionalFiles missing;
        missing.mandatory_files.emplace_back( data_file.value(),
            !file_exists( data_file_path( filename(), data_file.value() ) ) );
        return missing;
    }

    Percentage VOInput::is_loadable() const
    {
        std::ifstream file{ to_string( this->filename() ), std::ios::binary };
        if( goto_keyword_if_it_exists( file, "GOCAD Voxet" ) )
        {
            return Percentage{ 1 };
        }
        return Percentage{ 0 };
    }
} // namespace geode::internal
