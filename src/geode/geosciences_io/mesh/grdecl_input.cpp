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

#include <geode/geosciences_io/mesh/internal/grdecl_input.hpp>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>

#include <absl/container/flat_hash_set.h>
#include <absl/strings/ascii.h>
#include <absl/strings/match.h>
#include <absl/strings/numbers.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/file.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/string.hpp>

#include <geode/geometry/nn_search.hpp>

#include <geode/mesh/builder/hybrid_solid_builder.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>

#include <geode/geosciences_io/mesh/common.hpp>
#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

namespace
{
    struct Pillar
    {
        geode::Point3D top;
        geode::Point3D bottom;
    };

    geode::Point3D interpolate_on_pillar(
        const double depth, const Pillar& pillar )
    {
        const auto height = pillar.bottom.value( 2 ) - pillar.top.value( 2 );
        if( std::fabs( height ) <= geode::GLOBAL_EPSILON )
        {
            return geode::Point3D{ { pillar.top.value( 0 ),
                pillar.top.value( 1 ), depth } };
        }
        const auto lambda = ( depth - pillar.top.value( 2 ) ) / height;
        return pillar.bottom * lambda + pillar.top * ( 1 - lambda );
    }

    std::string_view strip_comment( std::string_view line )
    {
        char quote{ 0 };
        for( const auto c : geode::Indices{ line } )
        {
            const auto character = line[c];
            if( quote != 0 )
            {
                if( character == quote )
                {
                    quote = 0;
                }
                continue;
            }
            if( character == '\'' || character == '"' )
            {
                quote = character;
                continue;
            }
            if( character == '-' && c + 1 < line.size() && line[c + 1] == '-' )
            {
                return line.substr( 0, c );
            }
        }
        return line;
    }

    std::uint64_t checked_product(
        std::uint64_t lhs, std::uint64_t rhs, std::string_view what )
    {
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            rhs == 0 || lhs < geode::NO_ID / rhs, nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[GRDECLInput] Grid is too large: number of ", what,
            " exceeds the maximum index" );
        return lhs * rhs;
    }

    /*!
     * Reads the values of a keyword block until the terminating '/'.
     * Repeat counts (N*value, N*) are reported through the callback as
     * (value, count), an empty value meaning the default value.
     */
    template < typename Callback >
    void read_keyword_data( std::istream& file,
        std::string_view keyword,
        std::string first_line,
        Callback&& on_value )
    {
        auto line = std::move( first_line );
        while( true )
        {
            for( auto token : geode::string_split( strip_comment( line ) ) )
            {
                const auto slash = token.find( '/' );
                const auto is_last = slash != std::string_view::npos;
                if( is_last )
                {
                    token = token.substr( 0, slash );
                }
                if( !token.empty() )
                {
                    const auto star = token.find( '*' );
                    if( star == std::string_view::npos )
                    {
                        on_value( token, std::uint64_t{ 1 } );
                    }
                    else
                    {
                        std::uint64_t count{ 0 };
                        geode::OpenGeodeGeosciencesIOMeshException::
                            check_exception(
                                absl::SimpleAtoi(
                                    token.substr( 0, star ), &count ),
                                nullptr, geode::OpenGeodeException::TYPE::data,
                                "[GRDECLInput] Invalid repeat count \"", token,
                                "\" in keyword ", keyword );
                        on_value( token.substr( star + 1 ), count );
                    }
                }
                if( is_last )
                {
                    return;
                }
            }
            if( !std::getline( file, line ) )
            {
                geode::Logger::warning( "[GRDECLInput] Keyword ", keyword,
                    " is not terminated by '/'" );
                return;
            }
        }
    }

    std::vector< double > read_doubles( std::istream& file,
        std::string_view keyword,
        std::string first_line,
        std::uint64_t expected )
    {
        std::vector< double > values;
        bool extra_values{ false };
        read_keyword_data( file, keyword, std::move( first_line ),
            [&]( std::string_view token, std::uint64_t count ) {
                geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                    !token.empty(), nullptr,
                    geode::OpenGeodeException::TYPE::data,
                    "[GRDECLInput] Default values are not supported for "
                    "keyword ",
                    keyword );
                const auto value = geode::string_to_double( token );
                const auto remaining = expected - values.size();
                if( count > remaining )
                {
                    extra_values = true;
                    count = remaining;
                }
                values.insert( values.end(), count, value );
            } );
        if( extra_values )
        {
            geode::Logger::warning( "[GRDECLInput] Keyword ", keyword,
                " has more than the ", expected,
                " expected values, extra values are ignored" );
        }
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            values.size() == expected, nullptr,
            geode::OpenGeodeException::TYPE::data, "[GRDECLInput] Keyword ",
            keyword, " has ", values.size(), " values instead of ", expected );
        return values;
    }

    std::string read_include_path( std::istream& file, std::string first_line )
    {
        auto line = std::move( first_line );
        auto text = absl::StripAsciiWhitespace( strip_comment( line ) );
        while( text.empty() )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                static_cast< bool >( std::getline( file, line ) ), nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[GRDECLInput] Missing file name after INCLUDE keyword" );
            text = absl::StripAsciiWhitespace( strip_comment( line ) );
        }
        std::string path;
        std::string_view remaining;
        if( text.front() == '\'' || text.front() == '"' )
        {
            const auto end = text.find( text.front(), 1 );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                end != std::string_view::npos, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[GRDECLInput] Unterminated quoted INCLUDE file name: ", text );
            path = text.substr( 1, end - 1 );
            remaining = text.substr( end + 1 );
        }
        else
        {
            const auto end = text.find_first_of( " \t/" );
            path = text.substr( 0, end );
            remaining = end == std::string_view::npos ? "" : text.substr( end );
        }
        while( !absl::StrContains( remaining, '/' ) )
        {
            if( !std::getline( file, line ) )
            {
                geode::Logger::warning(
                    "[GRDECLInput] Keyword INCLUDE is not terminated by '/'" );
                break;
            }
            remaining = strip_comment( line );
        }
        return path;
    }

    class GRDECLInputImpl
    {
    public:
        GRDECLInputImpl(
            std::string_view filename, geode::HybridSolid3D& solid )
            : filename_{ filename },
              solid_( solid ),
              builder_{ geode::HybridSolidBuilder< 3 >::create( solid_ ) }
        {
        }

        void read_file()
        {
            parse_file(
                std::filesystem::path{ geode::to_string( filename_ ) } );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                nb_cells_ != 0, nullptr, geode::OpenGeodeException::TYPE::data,
                "[GRDECLInput] Cannot find keyword SPECGRID" );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                !pillars_.empty(), nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[GRDECLInput] Cannot find keyword COORD" );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                !depths_.empty(), nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[GRDECLInput] Cannot find keyword ZCORN" );
            create_cells( pillars_, depths_ );
        }

    private:
        void parse_file( const std::filesystem::path& path )
        {
            std::ifstream file{ path, std::ios::binary };
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                file.is_open(), nullptr, geode::OpenGeodeException::TYPE::data,
                "[GRDECLInput] Cannot open file ", path.string() );
            std::error_code error;
            const auto canonical_path =
                std::filesystem::weakly_canonical( path, error ).string();
            if( !parsed_files_.insert( canonical_path ).second )
            {
                geode::Logger::warning( "[GRDECLInput] File ", path.string(),
                    " is included recursively, skipping it" );
                return;
            }
            std::string line;
            while( std::getline( file, line ) )
            {
                const auto data = strip_comment( line );
                const auto tokens = geode::string_split( data );
                if( tokens.empty() )
                {
                    continue;
                }
                const auto keyword = tokens.front();
                const auto keyword_end =
                    static_cast< std::size_t >( keyword.data() - line.data() )
                    + keyword.size();
                auto rest = line.substr( keyword_end );
                if( keyword == "SPECGRID" || keyword == "DIMENS" )
                {
                    read_dimensions( file, keyword, std::move( rest ) );
                }
                else if( keyword == "COORD" )
                {
                    read_pillars( file, std::move( rest ) );
                }
                else if( keyword == "ZCORN" )
                {
                    read_depths( file, std::move( rest ) );
                }
                else if( keyword == "ACTNUM" )
                {
                    read_active_cells( file, std::move( rest ) );
                }
                else if( keyword == "INCLUDE" )
                {
                    const std::filesystem::path include_path{ read_include_path(
                        file, std::move( rest ) ) };
                    parse_file( include_path.is_absolute()
                                    ? include_path
                                    : path.parent_path() / include_path );
                }
            }
        }

        void read_dimensions(
            std::istream& file, std::string_view keyword, std::string line )
        {
            std::vector< geode::index_t > dimensions;
            read_keyword_data( file, keyword, std::move( line ),
                [&dimensions, keyword](
                    std::string_view token, std::uint64_t count ) {
                    for( std::uint64_t c = 0;
                        c < count && dimensions.size() < 3; c++ )
                    {
                        geode::index_t value{ 0 };
                        geode::OpenGeodeGeosciencesIOMeshException::
                            check_exception(
                                absl::SimpleAtoi( token, &value ) && value > 0,
                                nullptr, geode::OpenGeodeException::TYPE::data,
                                "[GRDECLInput] Invalid grid dimension \"",
                                token, "\" in keyword ", keyword );
                        dimensions.push_back( value );
                    }
                } );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                dimensions.size() == 3, nullptr,
                geode::OpenGeodeException::TYPE::data, "[GRDECLInput] Keyword ",
                keyword, " should define 3 grid dimensions" );
            nx_ = dimensions[0];
            ny_ = dimensions[1];
            nz_ = dimensions[2];
            nb_cells_ = checked_product(
                checked_product( nx_, ny_, "cells" ), nz_, "cells" );
            checked_product( nb_cells_, 8, "ZCORN values" );
            checked_product(
                checked_product( std::uint64_t{ nx_ } + 1, ny_ + 1, "pillars" ),
                6, "COORD values" );
        }

        void check_dimensions_are_read( std::string_view keyword ) const
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                nb_cells_ != 0, nullptr, geode::OpenGeodeException::TYPE::data,
                "[GRDECLInput] Keyword ", keyword,
                " found before grid dimensions (SPECGRID)" );
        }

        void read_pillars( std::istream& file, std::string line )
        {
            check_dimensions_are_read( "COORD" );
            const auto nb_pillars =
                static_cast< geode::index_t >( ( nx_ + 1 ) * ( ny_ + 1 ) );
            const auto values = read_doubles(
                file, "COORD", std::move( line ), 6 * nb_pillars );
            pillars_.resize( nb_pillars );
            for( const auto p : geode::Range{ nb_pillars } )
            {
                for( const auto d : geode::LRange{ 3 } )
                {
                    pillars_[p].top.set_value( d, values[6 * p + d] );
                    pillars_[p].bottom.set_value( d, values[6 * p + 3 + d] );
                }
            }
        }

        void read_active_cells( std::istream& file, std::string line )
        {
            check_dimensions_are_read( "ACTNUM" );
            active_cells_.clear();
            active_cells_.reserve( nb_cells_ );
            bool extra_values{ false };
            read_keyword_data( file, "ACTNUM", std::move( line ),
                [this, &extra_values](
                    std::string_view token, std::uint64_t count ) {
                    int value{ 1 };
                    geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                        token.empty() || absl::SimpleAtoi( token, &value ),
                        nullptr, geode::OpenGeodeException::TYPE::data,
                        "[GRDECLInput] Invalid ACTNUM value \"", token, "\"" );
                    const auto remaining = nb_cells_ - active_cells_.size();
                    if( count > remaining )
                    {
                        extra_values = true;
                        count = remaining;
                    }
                    active_cells_.insert( active_cells_.end(), count,
                        static_cast< std::uint8_t >( value != 0 ) );
                } );
            if( extra_values )
            {
                geode::Logger::warning( "[GRDECLInput] Keyword ACTNUM has "
                                        "more values than cells, extra values "
                                        "are ignored" );
            }
            if( active_cells_.size() != nb_cells_ )
            {
                geode::Logger::warning( "[GRDECLInput] Keyword ACTNUM has ",
                    active_cells_.size(), " values instead of ", nb_cells_,
                    ", missing cells are considered active" );
                active_cells_.resize( nb_cells_, 1 );
            }
        }

        void read_depths( std::istream& file, std::string line )
        {
            check_dimensions_are_read( "ZCORN" );
            depths_ =
                read_doubles( file, "ZCORN", std::move( line ), 8 * nb_cells_ );
        }

        std::array< geode::Point3D, 8 > cell_points(
            const std::array< geode::index_t, 3 >& grid_coordinates,
            absl::Span< const Pillar > pillars,
            absl::Span< const double > depths ) const
        {
            const auto pillars_id = cell_pillars_id( grid_coordinates );
            const auto& bottom_left_pillar = pillars[pillars_id[0]];
            const auto& top_left_pillar = pillars[pillars_id[1]];
            const auto& bottom_right_pillar = pillars[pillars_id[2]];
            const auto& top_right_pillar = pillars[pillars_id[3]];
            std::array< geode::Point3D, 8 > points;
            points[0] = interpolate_on_pillar(
                depths[2 * grid_coordinates[0] + 4 * nx_ * grid_coordinates[1]
                       + 4 * nx_ * ny_ * ( grid_coordinates[2] + 1 ) + 2 * nx_],
                bottom_left_pillar );
            points[1] = interpolate_on_pillar(
                depths[2 * grid_coordinates[0] + 4 * nx_ * grid_coordinates[1]
                       + 4 * nx_ * ny_ * ( grid_coordinates[2] + 1 ) + 2 * nx_
                       + 1],
                bottom_right_pillar );
            points[2] = interpolate_on_pillar(
                depths[2 * grid_coordinates[0] + 4 * nx_ * grid_coordinates[1]
                       + 4 * nx_ * ny_ * ( grid_coordinates[2] + 1 ) + 1],
                top_right_pillar );
            points[3] = interpolate_on_pillar(
                depths[2 * grid_coordinates[0] + 4 * nx_ * grid_coordinates[1]
                       + 4 * nx_ * ny_ * ( grid_coordinates[2] + 1 )],
                top_left_pillar );
            points[4] = interpolate_on_pillar(
                depths[2 * grid_coordinates[0] + 4 * nx_ * grid_coordinates[1]
                       + 2 * nx_ + 4 * nx_ * ny_ * grid_coordinates[2]],
                bottom_left_pillar );
            points[5] = interpolate_on_pillar(
                depths[2 * grid_coordinates[0] + 4 * nx_ * grid_coordinates[1]
                       + 2 * nx_ + 1 + 4 * nx_ * ny_ * grid_coordinates[2]],
                bottom_right_pillar );
            points[6] = interpolate_on_pillar(
                depths[2 * grid_coordinates[0] + 4 * nx_ * grid_coordinates[1]
                       + 1 + 4 * nx_ * ny_ * grid_coordinates[2]],
                top_right_pillar );
            points[7] = interpolate_on_pillar(
                depths[2 * grid_coordinates[0] + 4 * nx_ * grid_coordinates[1]
                       + 4 * nx_ * ny_ * grid_coordinates[2]],
                top_left_pillar );
            return points;
        }

        bool is_cell_active( geode::index_t cell_id ) const
        {
            return active_cells_.empty() || active_cells_[cell_id];
        }

        void create_cells( absl::Span< const Pillar > pillars,
            absl::Span< const double > depths )
        {
            std::vector< geode::Point3D > points;
            points.reserve( 8 * nb_cells_ );
            geode::index_t cell_id{ 0 };
            geode::index_t nb_collapsed_cells{ 0 };
            for( const auto k : geode::Range{ nz_ } )
            {
                for( const auto j : geode::Range{ ny_ } )
                {
                    for( const auto i : geode::Range{ nx_ } )
                    {
                        if( !is_cell_active( cell_id++ ) )
                        {
                            continue;
                        }
                        const auto corners =
                            cell_points( { i, j, 2 * k }, pillars, depths );
                        if( corners[0].inexact_equal( corners[4] )
                            && corners[1].inexact_equal( corners[5] )
                            && corners[2].inexact_equal( corners[6] )
                            && corners[3].inexact_equal( corners[7] ) )
                        {
                            nb_collapsed_cells++;
                            continue;
                        }
                        points.insert(
                            points.end(), corners.begin(), corners.end() );
                    }
                }
            }
            if( nb_collapsed_cells != 0 )
            {
                geode::Logger::warning( "[GRDECLInput] ", nb_collapsed_cells,
                    " active cells have no thickness and are skipped" );
            }
            const auto collocated_mapping =
                geode::NNSearch3D{ points }.colocated_index_mapping(
                    geode::GLOBAL_EPSILON );
            for( const auto& point : collocated_mapping.unique_points )
            {
                builder_->create_point( point );
            }
            const auto& mapping = collocated_mapping.colocated_mapping;
            for( const auto cell : geode::Range{ points.size() / 8 } )
            {
                const auto offset = 8 * cell;
                builder_->create_hexahedron(
                    { mapping[offset], mapping[offset + 1], mapping[offset + 2],
                        mapping[offset + 3], mapping[offset + 4],
                        mapping[offset + 5], mapping[offset + 6],
                        mapping[offset + 7] } );
            }
            builder_->compute_polyhedron_adjacencies();
        }

        std::array< geode::index_t, 4 > cell_pillars_id(
            const std::array< geode::index_t, 3 >& grid_coordinates ) const
        {
            const auto number_of_vertices_on_lines = nx_ + 1;
            const auto bottom_left_pillar_id =
                grid_coordinates[0]
                + number_of_vertices_on_lines * ( grid_coordinates[1] + 1 );
            const auto top_left_pillar_id =
                grid_coordinates[0]
                + number_of_vertices_on_lines * grid_coordinates[1];
            const auto bottom_right_pillar_id =
                grid_coordinates[0]
                + number_of_vertices_on_lines * ( grid_coordinates[1] + 1 ) + 1;
            const auto top_right_pillar_id =
                grid_coordinates[0]
                + number_of_vertices_on_lines * grid_coordinates[1] + 1;
            return { bottom_left_pillar_id, top_left_pillar_id,
                bottom_right_pillar_id, top_right_pillar_id };
        }

    private:
        std::string_view filename_;
        geode::HybridSolid3D& solid_;
        std::unique_ptr< geode::HybridSolidBuilder3D > builder_{ nullptr };
        geode::index_t nx_{ 0 };
        geode::index_t ny_{ 0 };
        geode::index_t nz_{ 0 };
        geode::index_t nb_cells_{ 0 };
        std::vector< Pillar > pillars_;
        std::vector< double > depths_;
        std::vector< std::uint8_t > active_cells_;
        absl::flat_hash_set< std::string > parsed_files_;
    };
} // namespace

namespace geode::internal
{
    std::unique_ptr< HybridSolid3D > GRDECLInput::read( const MeshImpl& impl )
    {
        auto solid = HybridSolid3D::create( impl );
        GRDECLInputImpl reader{ this->filename(), *solid };
        reader.read_file();
        return solid;
    }

    Percentage GRDECLInput::is_loadable() const
    {
        std::ifstream file{ to_string( this->filename() ) };
        if( goto_keyword_if_it_exists( file, "SPECGRID" ) )
        {
            return Percentage{ 1 };
        }
        return Percentage{ 0 };
    }
} // namespace geode::internal