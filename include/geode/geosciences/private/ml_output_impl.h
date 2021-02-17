/*
 * Copyright (c) 2019 - 2021 Geode-solutions
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

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <queue>

#include <geode/basic/filename.h>

#include <geode/geometry/basic_objects.h>
#include <geode/geometry/bounding_box.h>
#include <geode/geometry/signed_mensuration.h>

#include <geode/mesh/core/polygonal_surface.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/corner.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/mixin/core/vertex_identifier.h>
#include <geode/model/representation/core/brep.h>

#include <geode/geosciences/private/gocad_common.h>

namespace geode
{
    namespace detail{
    bool check_structural_model_polygons(
        const BRep& structural_model )
    {
        for( const auto& surface : structural_model.surfaces() )
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

    absl::optional< PolygonEdge > get_one_border_edge(
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

        void determine_surface_to_regions_signs()
        {
            const auto paired_signs = determine_paired_signs();
            {
                std::vector< uuid > universe_boundaries;
                for( const auto& boundary : model_.model_boundaries() )
                {
                    for( const auto& item :
                        model_.model_boundary_items( boundary ) )
                    {
                        universe_boundaries.push_back( item.id() );
                    }
                }
                const auto relative_signs = determine_relative_signs(
                    universe_boundaries, paired_signs );
                const auto correct =
                    are_correct_sides( universe_boundaries, relative_signs );

                for( const auto& b :
                    Range{ universe_boundaries.size() } )
                {
                    universe_surface_sides_.emplace( universe_boundaries[b],
                        correct ? !relative_signs[b] : relative_signs[b] );
                }
            }
            for( const auto& block : model_.blocks() )
            {
                std::vector< uuid > block_boundaries;
                for( const auto& surface : model_.boundaries( block ) )
                {
                    block_boundaries.push_back( surface.id() );
                }
                const auto relative_signs =
                    determine_relative_signs( block_boundaries, paired_signs );
                const auto correct =
                    are_correct_sides( block_boundaries, relative_signs );

                for( const auto& b : Range{ block_boundaries.size() } )
                {
                    regions_surface_sides_.emplace(
                        std::pair< uuid, uuid >{
                            block.id(), block_boundaries[b] },
                        correct ? relative_signs[b] : !relative_signs[b] );
                }
            }
        }

        void write_file()
        {
            file_ << "GOCAD Model3d 1" << EOL;
            detail::HeaderData header;
            header.name = model_name_;
            detail::write_header( file_, header );
            detail::write_CRS( file_, {} );
            write_model_components();
            write_model_surfaces();
        }

	protected:
        MLOutputImpl( absl::string_view filename, const Model& model )
            : file_( filename.data() ), model_( model )
        {
            OPENGEODE_EXCEPTION( file_.good(),
                "[MLOutput] Error while opening file: ", filename );
            model_name_ = filename_without_extension( filename );
        }

		index_t& component_id()
        {
            return component_id_;
        }

		absl::flat_hash_map< uuid, index_t >& components
        {
            return components_;
        }

        std::ofstream& file()
        {
            return file_;
        }

    private:
        absl::flat_hash_map< std::pair< uuid, uuid >, bool >
            determine_paired_signs() const
        {
            absl::flat_hash_map< std::pair< uuid, uuid >, bool >
                paired_signs;
            for( const auto& line : model_.lines() )
            {
                const auto uid0 =
                    model_.unique_vertex( { line.component_id(), 0 } );
                const auto uid1 =
                    model_.unique_vertex( { line.component_id(), 1 } );
                const auto surface_mcvs0 = model_.mesh_component_vertices(
                    uid0, Surface3D::component_type_static() );
                const auto surface_mcvs1 = model_.mesh_component_vertices(
                    uid1, Surface3D::component_type_static() );
                absl::flat_hash_map< uuid, bool > surface_direct_edges;
                for( const auto& surface_mcv0 : surface_mcvs0 )
                {
                    for( const auto& surface_mcv1 : surface_mcvs1 )
                    {
                        if( surface_mcv1.component_id.id()
                            != surface_mcv0.component_id.id() )
                        {
                            continue;
                        }
                        const auto& surface =
                            model_.surface( surface_mcv0.component_id.id() );
                        const auto& surface_mesh = surface.mesh();
                        const auto v0v1 =
                            surface_mesh.polygon_edge_from_vertices(
                                surface_mcv0.vertex, surface_mcv1.vertex );
                        const auto v1v0 =
                            surface_mesh.polygon_edge_from_vertices(
                                surface_mcv1.vertex, surface_mcv0.vertex );
                        if( v0v1 && !v1v0 )
                        {
                            surface_direct_edges.emplace( surface.id(), true );
                        }
                        else if( v1v0 && !v0v1 )
                        {
                            surface_direct_edges.emplace( surface.id(), false );
                        }
                    }
                }
                const auto nb_adj_surfaces = surface_direct_edges.size();
                if( nb_adj_surfaces < 2 )
                {
                    continue;
                }
                for( const auto& s0 : surface_direct_edges )
                {
                    for( const auto& s1 : surface_direct_edges )
                    {
                        if( s0.first < s1.first )
                        {
                            paired_signs.emplace(
                                std::pair< uuid, uuid >{
                                    s0.first, s1.first },
                                s0.second != s1.second );
                        }
                    }
                }
            }
            return paired_signs;
        }

        absl::FixedArray< bool > determine_relative_signs(
            absl::Span< const uuid > universe_boundaries,
            const absl::flat_hash_map< std::pair< uuid, uuid >,
                bool >& paired_signs ) const
        {
            const auto nb_surfaces = universe_boundaries.size();
            absl::FixedArray< bool > signs( nb_surfaces, true );
            if( nb_surfaces == 1 )
            {
                return signs;
            }
            absl::FixedArray< bool > determined( nb_surfaces, false );
            determined[0] = true;
            std::queue< index_t > to_process;
            to_process.push( 0 );
            while( !to_process.empty() )
            {
                const auto determined_s = to_process.front();
                const auto determined_s_id = universe_boundaries[determined_s];
                to_process.pop();
                for( const auto s : Range{ nb_surfaces } )
                {
                    if( determined[s] )
                    {
                        continue;
                    }
                    const auto s_id = universe_boundaries[s];
                    const auto itr =
                        paired_signs.find( { std::min( determined_s_id, s_id ),
                            std::max( determined_s_id, s_id ) } );
                    if( itr == paired_signs.end() )
                    {
                        continue;
                    }
                    signs[s] = itr->second ? signs[determined_s]
                                           : !signs[determined_s];
                    determined[s] = true;
                    to_process.push( s );
                }
            }
            OPENGEODE_ASSERT( absl::c_count( determined, false ) == 0,
                "All signs should have been found" );
            return signs;
        }

        bool are_correct_sides(
            absl::Span< const uuid > universe_boundaries,
            const absl::FixedArray< bool >& relative_signs ) const
        {
            double signed_volume{ 0 };
            const auto& first_surface_mesh =
                model_.surface( universe_boundaries[0] ).mesh();
            const auto& bbox = first_surface_mesh.bounding_box();
            const auto center = ( bbox.min() + bbox.max() ) * 0.5;
            for( const auto s : Range{ universe_boundaries.size() } )
            {
                const auto sign = relative_signs[s];
                const auto& surface_mesh =
                    model_.surface( universe_boundaries[s] ).mesh();
                for( const auto t : Range{ surface_mesh.nb_polygons() } )
                {
                    std::array< local_index_t, 3 > vertex_order{ 0, 1,
                        2 };
                    if( !sign )
                    {
                        vertex_order[1] = 2;
                        vertex_order[2] = 1;
                    }
                    Tetra tetra{ surface_mesh.point(
                                            surface_mesh.polygon_vertex(
                                                { t, vertex_order[0] } ) ),
                        surface_mesh.point( surface_mesh.polygon_vertex(
                            { t, vertex_order[1] } ) ),
                        surface_mesh.point( surface_mesh.polygon_vertex(
                            { t, vertex_order[2] } ) ),
                        center };
                    signed_volume += tetra_signed_volume( tetra );
                }
            }
            OPENGEODE_ASSERT( std::fabs( signed_volume ) > 0,
                "Null volume block is not valid" );
            return signed_volume > 0;
        }

        virtual void write_geological_tsurfs() = 0;

        void write_tsurfs()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                file_ << "TSURF " << boundary.name() << EOL;
            }
            write_geological_tsurfs();
            bool unclassified_surfaces{ false };
            for( const auto& surface : model_.surfaces() )
            {
                if( model_.nb_collections( surface.id() ) == 0 )
                {
                    unclassified_surfaces = true;
                    break;
                }
            }
            if( unclassified_surfaces )
            {
                file_ << "TSURF " << unclassified_surfaces_name_ << EOL;
            }
        }

template < typename Item >
        void write_key_triangle( const Item& item)
        {
                    const auto& mesh = item.mesh();
                    for( const auto v : LRange{ 3 } )
                    {
                        file_ << SPACE << SPACE << mesh.point( mesh.polygon_vertex( { 0, v } ) ).string()
                              << EOL;
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
                            boundary.name(),
                            " belongs to several collections. It has been "
                            "exported only once" );
                        continue;
                    }
                    file_ << "TFACE " << component_id_ << SPACE << "boundary"
                          << SPACE << boundary.name() << EOL;
                          write_key_triangle(item);
                    components_.emplace( item.id(), component_id_++ );
                }
            }
            write_geological_tfaces();
            if( components_.size() == model_.nb_surfaces() )
            {
                return;
            }
            for( const auto& surface : model_.surfaces() )
            {
                if( components_.find( surface.id() ) != components_.end() )
                {
                    continue;
                }
                file_ << "TFACE " << component_id_ << SPACE << "boundary"
                      << SPACE << unclassified_surfaces_name_ << EOL;
                          write_key_triangle(item);
                components_.emplace( surface.id(), component_id_++ );
                unclassified_surfaces_.emplace_back( surface.id() );
            }
        }

        void write_universe()
        {
            file_ << "REGION " << component_id_ << SPACE << SPACE << "Universe "
                  << EOL << SPACE << SPACE;
            index_t counter{ 0 };
            for( const auto& boundary : model_.model_boundaries() )
            {
                for( const auto& item :
                    model_.model_boundary_items( boundary ) )
                {
                    char sign{ '-' };
                    if( universe_surface_sides_.at( item.id() ) )
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
                file_ << "REGION " << component_id_ << SPACE << region.name()
                      << EOL << SPACE << SPACE;
                index_t counter{ 0 };
                for( const auto& surface : model_.boundaries( region ) )
                {
                    const auto sign = regions_surface_sides_.at(
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
                for( const auto& surface : model_.internal_surfaces( region ) )
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

        index_t write_surface( const Surface3D& surface,
            const index_t current_offset )
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
            const auto corner_mcvs1 = model_.mesh_component_vertices(
                uid1, Corner3D::component_type_static() );
            if( !corner_mcvs1.empty() )
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
            process_surface_edge(
                surface, first_on_border.value(), current_offset, line_starts );

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

        void find_corners_and_line_starts_for_unclassified_surfaces(
            std::vector< std::array< index_t, 2 > >& line_starts ) const
        {
            index_t current_offset{ OFFSET_START };
            for( const auto& surface_id : unclassified_surfaces_ )
            {
                const auto& surface = model_.surface( surface_id );
                add_corners_and_line_starts(
                    surface, current_offset, line_starts );
                current_offset += surface.mesh().nb_vertices();
            }
        }

        void write_corners(
            const std::vector< std::array< index_t, 2 > >& line_starts )
        {
            for( const auto& line_start : line_starts )
            {
                file_ << "BSTONE " << line_start[0] << EOL;
            }
        }

        void write_line_starts( index_t current_offset,
            const std::vector< std::array< index_t, 2 > >& line_starts )
        {
            for( const auto& line_start : line_starts )
            {
                file_ << "BORDER " << current_offset++ << SPACE << line_start[0]
                      << SPACE << line_start[1] << EOL;
            }
        }

        virtual void write_geological_model_surfaces() = 0;

        void write_model_surfaces()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                file_ << "GOCAD TSurf 1" << EOL;
                detail::HeaderData header;
                header.name = boundary.name().data();
                detail::write_header( file_, header );
                detail::write_CRS( file_, {} );
                file_ << "GEOLOGICAL_FEATURE " << boundary.name() << EOL;
                file_ << "GEOLOGICAL_TYPE boundary" << EOL;
                index_t current_offset{ OFFSET_START };
                for( const auto& item :
                    model_.model_boundary_items( boundary ) )
                {
                    file_ << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< index_t, 2 > > line_starts;
                find_boundary_corners_and_line_starts( boundary, line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                file_ << "END" << EOL;
            }
            write_geological_model_surfaces();
            if( unclassified_surfaces_.empty() )
            {
                return;
            }
            file_ << "GOCAD TSurf 1" << EOL;
            detail::HeaderData header;
            header.name = unclassified_surfaces_name_;
            detail::write_header( file_, header );
            detail::write_CRS( file_, {} );
            file_ << "GEOLOGICAL_FEATURE " << unclassified_surfaces_name_
                  << EOL;
            file_ << "GEOLOGICAL_TYPE "
                  << "boundary" << EOL;
            index_t current_offset{ OFFSET_START };
            for( const auto& surface_id : unclassified_surfaces_ )
            {
                file_ << "TFACE" << EOL;
                current_offset = write_surface(
                    model_.surface( surface_id ), current_offset );
            }
            std::vector< std::array< index_t, 2 > > line_starts;
            find_corners_and_line_starts_for_unclassified_surfaces(
                line_starts );
            write_corners( line_starts );
            write_line_starts( current_offset, line_starts );
            file_ << "END" << EOL;
        }

    private:
        std::ofstream file_;
        std::string model_name_;
        const Model& model_;
        absl::flat_hash_map< uuid, index_t > components_;
        absl::flat_hash_map< uuid, bool > universe_surface_sides_;
        absl::flat_hash_map< std::pair< uuid, uuid >, bool >
            regions_surface_sides_;
        index_t component_id_{ OFFSET_START };
        const std::string unclassified_surfaces_name_{
            "unclassified_surfaces"
        };
        std::vector< uuid > unclassified_surfaces_;
    };
    } // namespace detail
} // namespace geode
