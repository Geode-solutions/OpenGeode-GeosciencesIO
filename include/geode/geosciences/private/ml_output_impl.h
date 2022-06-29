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

#pragma once

#include <fstream>

#include <absl/strings/str_replace.h>

#include <geode/basic/filename.h>

#include <geode/mesh/core/surface_mesh.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/corner.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/representation/core/brep.h>

#include <geode/geosciences/private/gocad_common.h>

namespace geode
{
    namespace detail
    {
        inline bool check_brep_polygons( const BRep& brep )
        {
            for( const auto& surface : brep.surfaces() )
            {
                const auto& mesh = surface.mesh();
                for( const auto p : Range{ mesh.nb_polygons() } )
                {
                    if( mesh.nb_polygon_vertices( p ) != 3 )
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        inline absl::optional< PolygonEdge > get_one_border_edge(
            const SurfaceMesh3D& mesh )
        {
            for( const auto p : Range{ mesh.nb_polygons() } )
            {
                for( const auto e : LRange{ 3 } )
                {
                    if( mesh.is_edge_on_border( { p, e } ) )
                    {
                        return PolygonEdge{ p, e };
                    }
                }
            }
            return absl::nullopt;
        }

        template < typename Model >
        class MLOutputImpl
        {
        public:
            static constexpr index_t OFFSET_START{ 1 };
            static constexpr char EOL{ '\n' };
            static constexpr char SPACE{ ' ' };

            virtual ~MLOutputImpl() = default;

            void write_file()
            {
                file_ << "GOCAD Model3d 1" << EOL;
                detail::HeaderData header;
                header.name = to_string( model_.name() );
                detail::write_header( file_, header );
                detail::write_CRS( file_, {} );
                write_model_components();
                write_model_surfaces();
            }

        protected:
            MLOutputImpl( absl::string_view filename, const Model& model )
                : file_{ to_string( filename ) },
                  model_( model ),
                  sides_( determine_surface_to_regions_sides( model ) )
            {
                OPENGEODE_EXCEPTION( file_.good(),
                    "[MLOutput] Error while opening file: ", filename );
            }

            index_t& component_id()
            {
                return component_id_;
            }

            absl::flat_hash_map< uuid, index_t >& components()
            {
                return components_;
            }

            std::ofstream& file()
            {
                return file_;
            }

            template < typename Component >
            void write_key_triangle( const Component& component )
            {
                const auto& mesh = component.mesh();
                for( const auto v : LRange{ 3 } )
                {
                    file_ << SPACE << SPACE
                          << mesh.point( mesh.polygon_vertex( { 0, v } ) )
                                 .string()
                          << EOL;
                }
            }

            virtual void write_geological_tsurfs() = 0;

            void write_tsurfs()
            {
                for( const auto& boundary : model_.model_boundaries() )
                {
                    file_ << "TSURF " << component_name( boundary ) << EOL;
                }
                write_geological_tsurfs();
                for( const auto& surface : model_.surfaces() )
                {
                    if( model_.nb_collections( surface.id() ) == 0 )
                    {
                        file_ << "TSURF " << component_name( surface ) << EOL;
                        unclassified_surfaces_.emplace_back( surface.id() );
                    }
                }
            }

            virtual void write_geological_tfaces() = 0;

            void write_tfaces()
            {
                for( const auto& boundary : model_.model_boundaries() )
                {
                    for( const auto& item :
                        model_.model_boundary_items( boundary ) )
                    {
                        if( components_.find( item.id() ) != components_.end() )
                        {
                            Logger::warn( "[MLOutput] A Surface from ",
                                component_name( boundary ),
                                " belongs to several collections. It has been "
                                "exported only once" );
                            continue;
                        }
                        file_ << "TFACE " << component_id_ << SPACE
                              << "boundary" << SPACE
                              << component_name( boundary ) << EOL;
                        write_key_triangle( item );
                        components_.emplace( item.id(), component_id_++ );
                    }
                }
                write_geological_tfaces();
                for( const auto& surface_id : unclassified_surfaces_ )
                {
                    const auto& surface = model_.surface( surface_id );
                    file_ << "TFACE " << component_id_ << SPACE << "boundary"
                          << SPACE << component_name( surface ) << EOL;
                    write_key_triangle( surface );
                    components_.emplace( surface_id, component_id_++ );
                }
            }

            void write_universe()
            {
                file_ << "REGION " << component_id_ << SPACE << SPACE
                      << "Universe " << EOL << SPACE << SPACE;
                index_t counter{ 0 };
                for( const auto& boundary : model_.model_boundaries() )
                {
                    for( const auto& item :
                        model_.model_boundary_items( boundary ) )
                    {
                        char sign{ '-' };
                        if( sides_.universe_surface_sides.at( item.id() ) )
                        {
                            sign = '+';
                        }
                        file_ << sign << components_.at( item.id() ) << SPACE
                              << SPACE;
                        counter++;
                        if( counter % 5 == 0 )
                        {
                            file_ << EOL << SPACE << SPACE;
                        }
                    }
                }
                file_ << 0 << EOL;
                component_id_++;
            }

            void write_regions()
            {
                write_universe();
                for( const auto& region : model_.blocks() )
                {
                    file_ << "REGION " << component_id_ << SPACE
                          << component_name( region ) << EOL << SPACE << SPACE;
                    index_t counter{ 0 };
                    for( const auto& surface : model_.boundaries( region ) )
                    {
                        const auto sign = sides_.regions_surface_sides.at(
                                              { region.id(), surface.id() } )
                                              ? '+'
                                              : '-';
                        file_ << sign << components_.at( surface.id() ) << SPACE
                              << SPACE;
                        counter++;
                        if( counter % 5 == 0 )
                        {
                            file_ << EOL << SPACE << SPACE;
                        }
                    }
                    for( const auto& surface :
                        model_.internal_surfaces( region ) )
                    {
                        file_ << "+" << components_.at( surface.id() ) << SPACE
                              << SPACE;
                        counter++;
                        if( counter % 5 == 0 )
                        {
                            file_ << EOL << SPACE << SPACE;
                        }
                        file_ << "-" << components_.at( surface.id() ) << SPACE
                              << SPACE;
                        counter++;
                        if( counter % 5 == 0 )
                        {
                            file_ << EOL << SPACE << SPACE;
                        }
                    }
                    file_ << 0 << EOL;
                    components_.emplace( region.id(), component_id_++ );
                }
                write_geological_regions();
            }

            virtual void write_geological_regions() = 0;

            void write_model_components()
            {
                write_tsurfs();
                write_tfaces();
                write_regions();
                file_ << "END" << EOL;
            }

            index_t write_surface(
                const Surface3D& surface, const index_t current_offset )
            {
                const auto& mesh = surface.mesh();
                for( const auto v : Range{ mesh.nb_vertices() } )
                {
                    file_ << "VRTX " << current_offset + v << SPACE
                          << mesh.point( v ).string() << EOL;
                }
                for( const auto t : Range{ mesh.nb_polygons() } )
                {
                    file_ << "TRGL "
                          << current_offset + mesh.polygon_vertex( { t, 0 } )
                          << SPACE
                          << current_offset + mesh.polygon_vertex( { t, 1 } )
                          << SPACE
                          << current_offset + mesh.polygon_vertex( { t, 2 } )
                          << EOL;
                }
                return current_offset + mesh.nb_vertices();
            }

            void process_surface_edge( const Surface3D& surface,
                const PolygonEdge& edge,
                const index_t current_offset,
                std::vector< std::array< index_t, 2 > >& line_starts ) const
            {
                const auto& mesh = surface.mesh();
                const auto v0 = mesh.polygon_vertex( edge );
                const auto v1 = mesh.polygon_vertex(
                    { edge.polygon_id, static_cast< local_index_t >(
                                           ( edge.edge_id + 1 ) % 3 ) } );
                const auto uid1 =
                    model_.unique_vertex( { surface.component_id(), v1 } );
                const auto corner_cmvs1 = model_.component_mesh_vertices(
                    uid1, Corner3D::component_type_static() );
                if( !corner_cmvs1.empty() )
                {
                    line_starts.emplace_back( std::array< index_t, 2 >{
                        v1 + current_offset, v0 + current_offset } );
                }
            }

            void add_corners_and_line_starts( const Surface3D& surface,
                const index_t current_offset,
                std::vector< std::array< index_t, 2 > >& line_starts ) const
            {
                // todo several times to process all border edges
                const auto& mesh = surface.mesh();
                const auto first_on_border = get_one_border_edge( mesh );
                if( !first_on_border )
                {
                    return;
                }
                process_surface_edge( surface, first_on_border.value(),
                    current_offset, line_starts );

                auto cur = mesh.previous_on_border( first_on_border.value() );
                while( cur != first_on_border )
                {
                    process_surface_edge(
                        surface, cur, current_offset, line_starts );
                    cur = mesh.previous_on_border( cur );
                }
            }

            void find_boundary_corners_and_line_starts(
                const ModelBoundary3D& surface_collection,
                std::vector< std::array< index_t, 2 > >& line_starts ) const
            {
                index_t current_offset{ OFFSET_START };
                for( const auto& surface :
                    model_.model_boundary_items( surface_collection ) )
                {
                    add_corners_and_line_starts(
                        surface, current_offset, line_starts );
                    current_offset += surface.mesh().nb_vertices();
                }
            }

            template < typename ItemRange >
            void find_corners_and_line_starts( const ItemRange& item_range,
                std::vector< std::array< index_t, 2 > >& line_starts ) const
            {
                index_t current_offset{ OFFSET_START };
                for( const auto& surface : item_range )
                {
                    add_corners_and_line_starts(
                        surface, current_offset, line_starts );
                    current_offset += surface.mesh().nb_vertices();
                }
            }

            std::vector< std::array< index_t, 2 > >
                find_corners_and_line_starts_for_unclassified_surface(
                    const uuid& surface_id ) const
            {
                std::vector< std::array< index_t, 2 > > line_starts;
                const auto& surface = model_.surface( surface_id );
                add_corners_and_line_starts(
                    surface, OFFSET_START, line_starts );
                return line_starts;
            }

            void write_corners(
                absl::Span< const std::array< index_t, 2 > > line_starts )
            {
                for( const auto& line_start : line_starts )
                {
                    file_ << "BSTONE " << line_start[0] << EOL;
                }
            }

            void write_line_starts( index_t current_offset,
                absl::Span< const std::array< index_t, 2 > > line_starts )
            {
                for( const auto& line_start : line_starts )
                {
                    file_ << "BORDER " << current_offset++ << SPACE
                          << line_start[0] << SPACE << line_start[1] << EOL;
                }
            }

            virtual void write_geological_model_surfaces() = 0;

            void write_model_surfaces()
            {
                for( const auto& boundary : model_.model_boundaries() )
                {
                    file_ << "GOCAD TSurf 1" << EOL;
                    detail::HeaderData header;
                    header.name = component_name( boundary );
                    detail::write_header( file_, header );
                    detail::write_CRS( file_, {} );
                    file_ << "GEOLOGICAL_FEATURE " << component_name( boundary )
                          << EOL;
                    file_ << "GEOLOGICAL_TYPE boundary" << EOL;
                    index_t current_offset{ OFFSET_START };
                    for( const auto& item :
                        model_.model_boundary_items( boundary ) )
                    {
                        file_ << "TFACE" << EOL;
                        current_offset = write_surface( item, current_offset );
                    }
                    std::vector< std::array< index_t, 2 > > line_starts;
                    find_boundary_corners_and_line_starts(
                        boundary, line_starts );
                    write_corners( line_starts );
                    write_line_starts( current_offset, line_starts );
                    file_ << "END" << EOL;
                }
                write_geological_model_surfaces();
                for( const auto& surface_id : unclassified_surfaces_ )
                {
                    file_ << "GOCAD TSurf 1" << EOL;
                    const auto& surface = model_.surface( surface_id );
                    detail::HeaderData header;
                    header.name = component_name( surface );
                    detail::write_header( file_, header );
                    detail::write_CRS( file_, {} );
                    file_ << "GEOLOGICAL_FEATURE " << component_name( surface )
                          << EOL;
                    file_ << "GEOLOGICAL_TYPE "
                          << "boundary" << EOL;
                    index_t current_offset{ OFFSET_START };
                    file_ << "TFACE" << EOL;
                    current_offset = write_surface( surface, current_offset );
                    const auto line_starts =
                        find_corners_and_line_starts_for_unclassified_surface(
                            surface_id );
                    write_corners( line_starts );
                    write_line_starts( current_offset, line_starts );
                    file_ << "END" << EOL;
                }
            }

            template < typename Component >
            std::string component_name( const Component& component ) const
            {
                return absl::StrReplaceAll(
                    component.name(), { { " ", "_" } } );
            }

        private:
            std::ofstream file_;
            const Model& model_;
            const RegionSurfaceSide sides_;
            absl::flat_hash_map< uuid, index_t > components_;
            index_t component_id_{ OFFSET_START };
            std::vector< uuid > unclassified_surfaces_;
        };
    } // namespace detail
} // namespace geode
