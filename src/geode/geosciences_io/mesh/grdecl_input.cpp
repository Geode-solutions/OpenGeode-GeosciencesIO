/*
 * Copyright (c) 2019 - 2024 Geode-solutions
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

#include <fstream>
#include <optional>
#include <string>

#include <absl/strings/match.h>
#include <absl/strings/str_split.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/file.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/string.hpp>

#include <geode/geometry/nn_search.hpp>

#include <geode/mesh/builder/hybrid_solid_builder.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>

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
        const auto lambda =
            ( depth - pillar.top.value( 2 ) )
            / ( pillar.bottom.value( 2 ) - pillar.top.value( 2 ) );
        return pillar.bottom * lambda + pillar.top * ( 1 - lambda );
    }

    class GRDECLInputImpl
    {
    public:
        GRDECLInputImpl(
            std::string_view filename, geode::HybridSolid3D& solid )
            : file_{ geode::to_string( filename ) },
              filepath_{
                  geode::filepath_without_filename( filename ).string()
              },
              solid_( solid ),
              builder_{ geode::HybridSolidBuilder< 3 >::create( solid ) }
        {
        }

        void read_file()
        {
            read_dimensions();
            get_filenames_and_keywords();
            const auto pillars = keyword_to_filename_map_.contains( "COORD" )
                                     ? read_pillars_with_file()
                                     : read_pillars();
            const auto depths = keyword_to_filename_map_.contains( "ZCORN" )
                                    ? read_depths_with_file()
                                    : read_depths();
            create_cells( pillars, depths );
        }

    private:
        void get_filenames_and_keywords()
        {
            auto line = geode::goto_keyword_if_it_exists( file_, "INCLUDE" );
            while( line != std::nullopt )
            {
                std::getline( file_, line.value() );
                const auto tokens = absl::StrSplit( line.value(), "_" );
                for( const auto& token : tokens )
                {
                    if( !absl::StrContains( token, "." ) )
                    {
                        continue;
                    }
                    const auto second_tokens = absl::StrSplit( token, "." );
                    for( const auto& second_token : second_tokens )
                    {
                        auto string_without_spaces =
                            geode::string_split( line.value() )[0];
                        keyword_to_filename_map_[second_token] =
                            std::string{ string_without_spaces.substr(
                                1, string_without_spaces.size() - 2 ) };
                        break;
                    }
                }

                line = geode::goto_keyword_if_it_exists( file_, "INCLUDE" );
            }
        }

        void read_dimensions()
        {
            auto line = geode::goto_keyword( file_, "SPECGRID" );
            while( !absl::StrContains( line, "F" ) )
            {
                std::getline( file_, line );
            }
            const auto tokens = geode::string_split( line );
            nx_ = geode::string_to_index( tokens[0] );
            ny_ = geode::string_to_index( tokens[1] );
            nz_ = geode::string_to_index( tokens[2] );
        }

        absl::FixedArray< Pillar > read_pillars_from_file( std::ifstream& file )
        {
            absl::FixedArray< Pillar > pillars( ( nx_ + 1 ) * ( ny_ + 1 ) );
            auto line = geode::goto_keyword( file, "COORD" );
            std::getline( file, line );
            geode::index_t pillar_number{ 0 };
            while( line != "/" )
            {
                const auto tokens = geode::string_split( line );
                if( !tokens.empty() )
                {
                    OPENGEODE_ASSERT( tokens.size() == 6,
                        "[GRDECLInput::read_pillars] Wrong "
                        "number of coordinates" );
                    Pillar pillar;
                    pillar.top =
                        geode::Point3D{ { geode::string_to_double( tokens[0] ),
                            geode::string_to_double( tokens[1] ),
                            geode::string_to_double( tokens[2] ) } };
                    pillar.bottom =
                        geode::Point3D{ { geode::string_to_double( tokens[3] ),
                            geode::string_to_double( tokens[4] ),
                            geode::string_to_double( tokens[5] ) } };
                    pillars[pillar_number++] = std::move( pillar );
                }
                std::getline( file, line );
            }
            return pillars;
        }

        absl::FixedArray< Pillar > read_pillars()
        {
            return read_pillars_from_file( file_ );
        }

        absl::FixedArray< Pillar > read_pillars_with_file()
        {
            std::ifstream file{ absl::StrCat(
                filepath_, keyword_to_filename_map_["COORD"] ) };
            return read_pillars_from_file( file );
        }

        absl::FixedArray< double > read_depths_from_file( std::ifstream& file )
        {
            absl::FixedArray< double > depths( 8 * nx_ * ny_ * nz_ );
            auto line = geode::goto_keyword( file, "ZCORN" );
            geode::index_t depths_number{ 0 };
            std::getline( file, line );
            while( line != "/" )
            {
                const auto tokens = geode::string_split( line );
                for( const auto depths_id : geode::Indices{ tokens } )
                {
                    depths[depths_number++] =
                        geode::string_to_double( tokens[depths_id] );
                }
                std::getline( file, line );
            }
            return depths;
        }

        absl::FixedArray< double > read_depths()
        {
            return read_depths_from_file( file_ );
        }

        absl::FixedArray< double > read_depths_with_file()
        {
            std::ifstream file{ absl::StrCat(
                filepath_, keyword_to_filename_map_["ZCORN"] ) };
            return read_depths_from_file( file );
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

        geode::NNSearch3D::ColocatedInfo create_points(
            absl::Span< const Pillar > pillars,
            absl::Span< const double > depths )
        {
            std::vector< geode::Point3D > points;
            points.reserve( 8 * nx_ * ny_ * nz_ );
            for( const auto k : geode::Range{ nz_ } )
            {
                for( const auto j : geode::Range{ ny_ } )
                {
                    for( const auto i : geode::Range{ nx_ } )
                    {
                        const std::array< geode::index_t, 3 > grid_coordinates{
                            i, j, 2 * k
                        };
                        for( const auto point :
                            cell_points( grid_coordinates, pillars, depths ) )
                        {
                            points.push_back( point );
                        }
                    }
                }
            }
            const auto collocated_mapping =
                geode::NNSearch3D{ points }.colocated_index_mapping(
                    geode::GLOBAL_EPSILON );
            for( const auto& point : collocated_mapping.unique_points )
            {
                builder_->create_point( point );
            }
            return collocated_mapping;
        }

        void create_cells( absl::Span< const Pillar > pillars,
            absl::Span< const double > depths )
        {
            const auto collocated_mapping = create_points( pillars, depths );
            for( const auto cell_id : geode::Range{ nx_ * ny_ * nz_ } )
            {
                builder_->create_hexahedron(
                    { collocated_mapping.colocated_mapping[0 + 8 * cell_id],
                        collocated_mapping.colocated_mapping[1 + 8 * cell_id],
                        collocated_mapping.colocated_mapping[2 + 8 * cell_id],
                        collocated_mapping.colocated_mapping[3 + 8 * cell_id],
                        collocated_mapping.colocated_mapping[4 + 8 * cell_id],
                        collocated_mapping.colocated_mapping[5 + 8 * cell_id],
                        collocated_mapping.colocated_mapping[6 + 8 * cell_id],
                        collocated_mapping
                            .colocated_mapping[7 + 8 * cell_id] } );
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
        std::ifstream file_;
        std::string filepath_;
        geode::HybridSolid3D& solid_;
        std::unique_ptr< geode::HybridSolidBuilder3D > builder_{ nullptr };
        geode::index_t nx_{ geode::NO_ID };
        geode::index_t ny_{ geode::NO_ID };
        geode::index_t nz_{ geode::NO_ID };
        absl::flat_hash_map< std::string, std::string >
            keyword_to_filename_map_{};
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::unique_ptr< HybridSolid3D > GRDECLInput::read(
            const MeshImpl& impl )
        {
            auto solid = HybridSolid3D::create( impl );
            GRDECLInputImpl reader{ this->filename(), *solid };
            reader.read_file();
            return solid;
        }

    } // namespace internal
} // namespace geode