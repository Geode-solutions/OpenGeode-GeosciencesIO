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

#include <geode/geosciences_io/mesh/private/grdecl_input.h>

#include <fstream>

#include <geode/basic/attribute_manager.h>
#include <geode/basic/filename.h>
#include <geode/basic/logger.h>
#include <geode/basic/string.h>

#include <geode/mesh/builder/hybrid_solid_builder.h>
#include <geode/mesh/core/hybrid_solid.h>

#include <geode/geosciences_io/mesh/private/gocad_common.h>
#include <geode/geosciences_io/mesh/private/utils.h>

namespace
{
    class GRDECLInputImpl
    {
    public:
        GRDECLInputImpl(
            absl::string_view filename, geode::HybridSolid3D& solid )
            : file_{ geode::to_string( filename ) },
              solid_( solid ),
              builder_{ geode::HybridSolidBuilder< 3 >::create( solid ) }
        {
        }

        void read_file()
        {
            read_dimensions();
            auto pillars = read_pillars();
            auto depths = read_depths();
            for( const auto k : geode::Range{ nz_ } )
            {
                for( const auto j : geode::Range{ ny_ } )
                {
                    for( const auto i : geode::Range{ nx_ } )
                    {
                        auto points =
                            create_cell_points( i, j, 2 * k, depths, pillars );
                        for( const auto point : points )
                        {
                            builder_->create_point( point );
                        }
                        geode::index_t cell_id{ i + nx_ * j + nx_ * ny_ * k };
                        builder_->create_hexahedron( { 0 + 8 * cell_id,
                            1 + 8 * cell_id, 2 + 8 * cell_id, 3 + 8 * cell_id,
                            4 + 8 * cell_id, 5 + 8 * cell_id, 6 + 8 * cell_id,
                            7 + 8 * cell_id } );
                        DEBUG( "CREATE CELL" );
                    }
                }
            }
        }

    private:
        struct Pillar
        {
            geode::Point3D top;
            geode::Point3D bottom;
            Pillar() {}
            Pillar( geode::Point3D top_point, geode::Point3D bottom_point )
            {
                top = top_point;
                bottom = bottom_point;
            }
        };

        void read_dimensions()
        {
            auto line = geode::detail::goto_keyword( file_, "SPECGRID" );
            std::getline( file_, line );
            const auto tokens = geode::string_split( line );
            nx_ = geode::string_to_index( tokens[0] );
            ny_ = geode::string_to_index( tokens[1] );
            nz_ = geode::string_to_index( tokens[2] );
        }

        absl::FixedArray< Pillar > read_pillars()
        {
            absl::FixedArray< Pillar > pillars( 8 * nx_ * ny_ * nz_ );
            auto line = geode::detail::goto_keyword( file_, "COORD" );
            geode::index_t pillar_number{ 0 };
            std::getline( file_, line );
            while( line.compare( "/" ) != 0 )
            {
                // DEBUG( line );
                const auto tokens = geode::string_split( line );

                if( tokens.size() != 0 )
                {
                    geode::Point3D top_point =
                        geode::Point3D( { geode::string_to_double( tokens[0] ),
                            geode::string_to_double( tokens[1] ),
                            geode::string_to_double( tokens[2] ) } );
                    geode::Point3D bottom_point =
                        geode::Point3D( { geode::string_to_double( tokens[3] ),
                            geode::string_to_double( tokens[4] ),
                            geode::string_to_double( tokens[5] ) } );
                    Pillar pillar( top_point, bottom_point );
                    pillars[pillar_number] = pillar;
                    pillar_number += 1;
                    std::getline( file_, line );
                }
                else
                {
                    std::getline( file_, line );
                }
            }

            return pillars;
        }

        absl::FixedArray< double > read_depths()
        {
            absl::FixedArray< double > depths( 8 * nx_ * ny_ * nz_ );
            auto line = geode::detail::goto_keyword( file_, "ZCORN" );
            geode::index_t depths_number{ 0 };
            std::getline( file_, line );
            while( line.compare( "/" ) != 0 )
            {
                const auto tokens = geode::string_split( line );
                for( const auto depths_id : geode::Indices{ tokens } )
                {
                    depths[depths_number] =
                        geode::string_to_double( tokens[depths_id] );
                    depths_number += 1;
                }

                std::getline( file_, line );
            }
            return depths;
        }

        std::array< geode::Point3D, 8 > create_cell_points( geode::index_t i,
            geode::index_t j,
            geode::index_t k,
            absl::FixedArray< double > depths,
            absl::FixedArray< Pillar > pillars )
        {
            const auto bottom_left_pillar =
                pillars[i + ( nx_ + 1 ) * ( j + 1 )];
            const auto top_left_pillar = pillars[i + ( nx_ + 1 ) * ( j )];
            const auto bottom_right_pillar =
                pillars[i + ( nx_ + 1 ) * ( j + 1 ) + 1];
            const auto top_right_pillar = pillars[i + ( nx_ + 1 ) * ( j ) + 1];
            std::array< geode::Point3D, 8 > points;
            points[0] = interpolate_on_pillar(
                geode::Point3D( { bottom_left_pillar.bottom.value( 0 ),
                    bottom_left_pillar.bottom.value( 1 ),
                    depths[2 * i + 4 * nx_ * j + 4 * nx_ * ny_ * ( k + 1 )
                           + 2 * nx_] } ),
                bottom_left_pillar );
            points[1] = interpolate_on_pillar(
                geode::Point3D( { bottom_right_pillar.bottom.value( 0 ),
                    bottom_right_pillar.bottom.value( 1 ),
                    depths[2 * i + 4 * nx_ * j + 4 * nx_ * ny_ * ( k + 1 )
                           + 2 * nx_ + 1] } ),
                bottom_right_pillar );
            points[2] = interpolate_on_pillar(
                geode::Point3D( { top_right_pillar.bottom.value( 0 ),
                    top_right_pillar.bottom.value( 1 ),
                    depths[2 * i + 4 * nx_ * j + 4 * nx_ * ny_ * ( k + 1 )
                           + 1] } ),
                top_right_pillar );
            points[3] = interpolate_on_pillar(
                geode::Point3D( { top_left_pillar.bottom.value( 0 ),
                    top_left_pillar.bottom.value( 1 ),
                    depths[2 * i + 4 * nx_ * j + 4 * nx_ * ny_ * ( k + 1 )] } ),
                top_left_pillar );
            points[4] = interpolate_on_pillar(
                geode::Point3D( { bottom_left_pillar.top.value( 0 ),
                    bottom_left_pillar.top.value( 1 ),
                    depths[2 * i + 4 * nx_ * j + 2 * nx_
                           + 4 * nx_ * ny_ * k] } ),
                bottom_left_pillar );
            points[5] = interpolate_on_pillar(
                geode::Point3D( { bottom_right_pillar.top.value( 0 ),
                    bottom_right_pillar.top.value( 1 ),
                    depths[2 * i + 4 * nx_ * j + 2 * nx_ + 1
                           + 4 * nx_ * ny_ * k] } ),
                bottom_right_pillar );
            points[6] = interpolate_on_pillar(
                geode::Point3D( { top_right_pillar.top.value( 0 ),
                    top_right_pillar.top.value( 1 ),
                    depths[2 * i + 4 * nx_ * j + 1 + 4 * nx_ * ny_ * k] } ),
                top_right_pillar );
            points[7] = interpolate_on_pillar(
                geode::Point3D( { top_left_pillar.top.value( 0 ),
                    top_left_pillar.top.value( 1 ),
                    depths[2 * i + 4 * nx_ * j + 4 * nx_ * ny_ * k] } ),
                top_left_pillar );
            return points;
        }

        geode::Point3D interpolate_on_pillar(
            geode::Point3D point, Pillar pillar )
        {
            const auto lambda =
                ( point.value( 2 ) - pillar.top.value( 2 ) )
                / ( pillar.bottom.value( 2 ) - pillar.top.value( 2 ) );
            return pillar.bottom * lambda + pillar.top * ( 1 - lambda );
        }

    private:
        std::ifstream file_;
        geode::HybridSolid3D& solid_;
        std::unique_ptr< geode::HybridSolidBuilder3D > builder_;
        geode::index_t nx_;
        geode::index_t ny_;
        geode::index_t nz_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        std::unique_ptr< HybridSolid3D > GRDECLInput::read(
            const MeshImpl& impl )
        {
            auto solid = HybridSolid3D::create( impl );
            GRDECLInputImpl reader{ this->filename(), *solid };
            reader.read_file();
            return solid;
        }
    } // namespace detail
} // namespace geode