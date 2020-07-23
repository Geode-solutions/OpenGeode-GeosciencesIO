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

#include <geode/geosciences/private/ml_output.h>

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

#include <geode/geosciences/private/gocad_common.h>
#include <geode/geosciences/representation/core/structural_model.h>

namespace
{
    bool check_structural_model_polygons(
        const geode::StructuralModel& structural_model )
    {
        for( const auto& surface : structural_model.surfaces() )
        {
            const auto& mesh = surface.mesh();
            for( const auto p : geode::Range{ mesh.nb_polygons() } )
            {
                if( mesh.nb_polygon_vertices( p ) != 3 )
                {
                    return false;
                }
            }
        }
        return true;
    }

    absl::optional< geode::PolygonEdge > get_one_border_edge(
        const geode::SurfaceMesh3D& mesh )
    {
        for( const auto p : geode::Range{ mesh.nb_polygons() } )
        {
            for( const auto e : geode::Range{ 3 } )
            {
                if( mesh.is_edge_on_border( { p, e } ) )
                {
                    return geode::PolygonEdge{ p, e };
                }
            }
        }
        return absl::nullopt;
    }

    class MLOutputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };

        MLOutputImpl(
            absl::string_view filename, const geode::StructuralModel& model )
            : file_( filename.data() ), model_( model )
        {
            file_.precision( 12 );
            OPENGEODE_EXCEPTION( file_.good(),
                "[MLOutput] Error while opening file: ", filename );
            model_name_ = geode::filename_without_extension( filename );
        }

        absl::flat_hash_map< std::pair< geode::uuid, geode::uuid >, bool >
            determine_paired_signs() const
        {
            absl::flat_hash_map< std::pair< geode::uuid, geode::uuid >, bool >
                paired_signs;
            for( const auto& line : model_.lines() )
            {
                const auto uid0 =
                    model_.unique_vertex( { line.component_id(), 0 } );
                const auto uid1 =
                    model_.unique_vertex( { line.component_id(), 1 } );
                const auto surface_mcvs0 = model_.mesh_component_vertices(
                    uid0, geode::Surface3D::component_type_static() );
                const auto surface_mcvs1 = model_.mesh_component_vertices(
                    uid1, geode::Surface3D::component_type_static() );
                absl::flat_hash_map< geode::uuid, bool > surface_direct_edges;
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
                                std::pair< geode::uuid, geode::uuid >{
                                    s0.first, s1.first },
                                s0.second != s1.second );
                        }
                    }
                }
            }
            return paired_signs;
        }

        absl::FixedArray< bool > determine_relative_signs(
            absl::Span< const geode::uuid > universe_boundaries,
            const absl::flat_hash_map< std::pair< geode::uuid, geode::uuid >,
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
            std::queue< geode::index_t > to_process;
            to_process.push( 0 );
            while( !to_process.empty() )
            {
                const auto determined_s = to_process.front();
                const auto determined_s_id = universe_boundaries[determined_s];
                to_process.pop();
                for( const auto s : geode::Range{ nb_surfaces } )
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
            absl::Span< const geode::uuid > universe_boundaries,
            const absl::FixedArray< bool >& relative_signs ) const
        {
            double signed_volume{ 0 };
            const auto& first_surface_mesh =
                model_.surface( universe_boundaries[0] ).mesh();
            const auto& bbox = first_surface_mesh.bounding_box();
            const auto center = ( bbox.min() + bbox.max() ) * 0.5;
            for( const auto s : geode::Range{ universe_boundaries.size() } )
            {
                const auto sign = relative_signs[s];
                const auto& surface_mesh =
                    model_.surface( universe_boundaries[s] ).mesh();
                for( const auto t : geode::Range{ surface_mesh.nb_polygons() } )
                {
                    std::array< geode::index_t, 3 > vertex_order{ 0, 1, 2 };
                    if( !sign )
                    {
                        vertex_order[1] = 2;
                        vertex_order[2] = 1;
                    }
                    geode::Tetra tetra{ surface_mesh.point(
                                            surface_mesh.polygon_vertex(
                                                { t, vertex_order[0] } ) ),
                        surface_mesh.point( surface_mesh.polygon_vertex(
                            { t, vertex_order[1] } ) ),
                        surface_mesh.point( surface_mesh.polygon_vertex(
                            { t, vertex_order[2] } ) ),
                        center };
                    signed_volume += geode::tetra_signed_volume( tetra );
                }
            }
            OPENGEODE_ASSERT( std::fabs( signed_volume ) > 0,
                "Null volume block is not valid" );
            return signed_volume > 0;
        }

        void determine_surface_to_regions_signs()
        {
            const auto paired_signs = determine_paired_signs();
            {
                std::vector< geode::uuid > universe_boundaries;
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
                    geode::Range{ universe_boundaries.size() } )
                {
                    universe_surface_sides_.emplace( universe_boundaries[b],
                        correct ? !relative_signs[b] : relative_signs[b] );
                }
            }
            for( const auto& block : model_.blocks() )
            {
                std::vector< geode::uuid > block_boundaries;
                for( const auto& surface : model_.boundaries( block ) )
                {
                    block_boundaries.push_back( surface.id() );
                }
                const auto relative_signs =
                    determine_relative_signs( block_boundaries, paired_signs );
                const auto correct =
                    are_correct_sides( block_boundaries, relative_signs );

                for( const auto& b : geode::Range{ block_boundaries.size() } )
                {
                    regions_surface_sides_.emplace(
                        std::pair< geode::uuid, geode::uuid >{
                            block.id(), block_boundaries[b] },
                        correct ? relative_signs[b] : !relative_signs[b] );
                }
            }
        }

        void write_file()
        {
            file_ << "GOCAD Model3d 1" << EOL;
            geode::detail::HeaderData header;
            header.name = model_name_;
            geode::detail::write_header( file_, header );
            geode::detail::write_CRS( file_, {} );
            write_model_components();
            write_model_surfaces();
        }

    private:
        void write_tsurfs()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                file_ << "TSURF " << boundary.name() << EOL;
            }
            for( const auto& fault : model_.faults() )
            {
                file_ << "TSURF " << fault.name() << EOL;
            }
            for( const auto& horizon : model_.horizons() )
            {
                file_ << "TSURF " << horizon.name() << EOL;
            }
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

        void write_tfaces()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                for( const auto& item :
                    model_.model_boundary_items( boundary ) )
                {
                    if( components_.find( item.id() ) != components_.end() )
                    {
                        geode::Logger::warn( "[MLOutput] A Surface from ",
                            boundary.name(),
                            " belongs to several collections. It has been "
                            "exported only once" );
                        continue;
                    }
                    file_ << "TFACE " << component_id_ << SPACE << "boundary"
                          << SPACE << boundary.name() << EOL;
                    const auto& mesh = item.mesh();
                    for( const auto v : geode::Range{ 3 } )
                    {
                        const auto& coords =
                            mesh.point( mesh.polygon_vertex( { 0, v } ) );
                        file_ << SPACE << SPACE << coords.value( 0 ) << SPACE
                              << coords.value( 1 ) << SPACE << coords.value( 2 )
                              << EOL;
                    }
                    components_.emplace( item.id(), component_id_++ );
                }
            }
            for( const auto& fault : model_.faults() )
            {
                for( const auto& item : model_.fault_items( fault ) )
                {
                    if( components_.find( item.id() ) != components_.end() )
                    {
                        geode::Logger::warn( "[MLOutput] A Surface from ",
                            fault.name(),
                            " belongs to several collections. It has been "
                            "exported only once" );
                        continue;
                    }
                    file_ << "TFACE " << component_id_ << SPACE
                          << fault_map_.at( fault.type() ) << SPACE
                          << fault.name() << EOL;
                    const auto& mesh = item.mesh();
                    for( const auto v : geode::Range{ 3 } )
                    {
                        const auto& coords =
                            mesh.point( mesh.polygon_vertex( { 0, v } ) );
                        file_ << SPACE << SPACE << coords.value( 0 ) << SPACE
                              << coords.value( 1 ) << SPACE << coords.value( 2 )
                              << EOL;
                    }
                    components_.emplace( item.id(), component_id_++ );
                }
            }
            for( const auto& horizon : model_.horizons() )
            {
                for( const auto& item : model_.horizon_items( horizon ) )
                {
                    if( components_.find( item.id() ) != components_.end() )
                    {
                        geode::Logger::warn( "[MLOutput] A Surface from ",
                            horizon.name(),
                            " belongs to several collections. It has been "
                            "exported only once" );
                        continue;
                    }
                    file_ << "TFACE " << component_id_ << SPACE
                          << horizon_map_.at( horizon.type() ) << SPACE
                          << horizon.name() << EOL;
                    const auto& mesh = item.mesh();
                    for( const auto v : geode::Range{ 3 } )
                    {
                        const auto& coords =
                            mesh.point( mesh.polygon_vertex( { 0, v } ) );
                        file_ << SPACE << SPACE << coords.value( 0 ) << SPACE
                              << coords.value( 1 ) << SPACE << coords.value( 2 )
                              << EOL;
                    }
                    components_.emplace( item.id(), component_id_++ );
                }
            }
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
                const auto& mesh = surface.mesh();
                for( const auto v : geode::Range{ 3 } )
                {
                    const auto& coords =
                        mesh.point( mesh.polygon_vertex( { 0, v } ) );
                    file_ << SPACE << SPACE << coords.value( 0 ) << SPACE
                          << coords.value( 1 ) << SPACE << coords.value( 2 )
                          << EOL;
                }
                components_.emplace( surface.id(), component_id_++ );
                unclassified_surfaces_.emplace_back( surface.id() );
            }
        }

        void write_universe()
        {
            file_ << "REGION " << component_id_ << SPACE << SPACE << "Universe "
                  << EOL << SPACE << SPACE;
            geode::index_t counter{ 0 };
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
                geode::index_t counter{ 0 };
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
        }

        void write_stratigraphic_units()
        {
            for( const auto& stratigraphic_unit : model_.stratigraphic_units() )
            {
                file_ << "LAYER " << stratigraphic_unit.name() << EOL << SPACE
                      << SPACE;
                geode::index_t counter{ 0 };
                for( const auto& item :
                    model_.stratigraphic_unit_items( stratigraphic_unit ) )
                {
                    file_ << components_.at( item.id() ) << SPACE << SPACE;
                    counter++;
                    if( counter % 5 == 0 )
                    {
                        file_ << EOL << SPACE << SPACE;
                    }
                }
                file_ << 0 << EOL;
            }
        }

        void write_fault_blocks()
        {
            for( const auto& fault_block : model_.fault_blocks() )
            {
                file_ << "FAULT_BLOCK " << fault_block.name() << EOL << SPACE
                      << SPACE;
                geode::index_t counter{ 0 };
                for( const auto& item :
                    model_.fault_block_items( fault_block ) )
                {
                    file_ << components_.at( item.id() ) << SPACE << SPACE;
                    counter++;
                    if( counter % 5 == 0 )
                    {
                        file_ << EOL << SPACE << SPACE;
                    }
                }
                file_ << 0 << EOL;
            }
        }

        void write_model_components()
        {
            write_tsurfs();
            write_tfaces();
            write_regions();
            write_stratigraphic_units();
            write_fault_blocks();
            file_ << "END" << EOL;
        }

        geode::index_t write_surface( const geode::Surface3D& surface,
            const geode::index_t current_offset )
        {
            const auto& mesh = surface.mesh();
            for( const auto v : geode::Range{ mesh.nb_vertices() } )
            {
                file_ << "VRTX " << current_offset + v << SPACE
                      << mesh.point( v ).value( 0 ) << SPACE
                      << mesh.point( v ).value( 1 ) << SPACE
                      << mesh.point( v ).value( 2 ) << EOL;
            }
            for( const auto t : geode::Range{ mesh.nb_polygons() } )
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

        void process_surface_edge( const geode::Surface3D& surface,
            const geode::PolygonEdge& edge,
            const geode::index_t current_offset,
            std::vector< std::array< geode::index_t, 2 > >& line_starts ) const
        {
            const auto& mesh = surface.mesh();
            const auto v0 = mesh.polygon_vertex( edge );
            const auto v1 = mesh.polygon_vertex(
                { edge.polygon_id, ( edge.edge_id + 1 ) % 3 } );
            const auto uid1 =
                model_.unique_vertex( { surface.component_id(), v1 } );
            const auto corner_mcvs1 = model_.mesh_component_vertices(
                uid1, geode::Corner3D::component_type_static() );
            if( !corner_mcvs1.empty() )
            {
                line_starts.emplace_back( std::array< geode::index_t, 2 >{
                    v1 + current_offset, v0 + current_offset } );
            }
        }

        void add_corners_and_line_starts( const geode::Surface3D& surface,
            const geode::index_t current_offset,
            std::vector< std::array< geode::index_t, 2 > >& line_starts ) const
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
            const geode::ModelBoundary3D& surface_collection,
            std::vector< std::array< geode::index_t, 2 > >& line_starts ) const
        {
            geode::index_t current_offset{ OFFSET_START };
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
            std::vector< std::array< geode::index_t, 2 > >& line_starts ) const
        {
            geode::index_t current_offset{ OFFSET_START };
            for( const auto& surface : item_range )
            {
                add_corners_and_line_starts(
                    surface, current_offset, line_starts );
                current_offset += surface.mesh().nb_vertices();
            }
        }

        void find_corners_and_line_starts_for_unclassified_surfaces(
            std::vector< std::array< geode::index_t, 2 > >& line_starts ) const
        {
            geode::index_t current_offset{ OFFSET_START };
            for( const auto& surface_id : unclassified_surfaces_ )
            {
                const auto& surface = model_.surface( surface_id );
                add_corners_and_line_starts(
                    surface, current_offset, line_starts );
                current_offset += surface.mesh().nb_vertices();
            }
        }

        void write_corners(
            const std::vector< std::array< geode::index_t, 2 > >& line_starts )
        {
            for( const auto& line_start : line_starts )
            {
                file_ << "BSTONE " << line_start[0] << EOL;
            }
        }

        void write_line_starts( geode::index_t current_offset,
            const std::vector< std::array< geode::index_t, 2 > >& line_starts )
        {
            for( const auto& line_start : line_starts )
            {
                file_ << "BORDER " << current_offset++ << SPACE << line_start[0]
                      << SPACE << line_start[1] << EOL;
            }
        }

        void write_model_surfaces()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                file_ << "GOCAD TSurf 1" << EOL;
                geode::detail::HeaderData header;
                header.name = boundary.name().data();
                geode::detail::write_header( file_, header );
                geode::detail::write_CRS( file_, {} );
                file_ << "GEOLOGICAL_FEATURE " << boundary.name() << EOL;
                file_ << "GEOLOGICAL_TYPE boundary" << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item :
                    model_.model_boundary_items( boundary ) )
                {
                    file_ << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;
                find_boundary_corners_and_line_starts( boundary, line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                file_ << "END" << EOL;
            }
            for( const auto& fault : model_.faults() )
            {
                file_ << "GOCAD TSurf 1" << EOL;
                geode::detail::HeaderData header;
                header.name = fault.name().data();
                geode::detail::write_header( file_, header );
                geode::detail::write_CRS( file_, {} );
                file_ << "GEOLOGICAL_FEATURE " << fault.name() << EOL;
                file_ << "GEOLOGICAL_TYPE " << fault_map_.at( fault.type() )
                      << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.fault_items( fault ) )
                {
                    file_ << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;

                find_corners_and_line_starts(
                    model_.fault_items( fault ), line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                file_ << "END" << EOL;
            }
            for( const auto& horizon : model_.horizons() )
            {
                file_ << "GOCAD TSurf 1" << EOL;
                geode::detail::HeaderData header;
                header.name = horizon.name().data();
                geode::detail::write_header( file_, header );
                geode::detail::write_CRS( file_, {} );
                file_ << "GEOLOGICAL_FEATURE " << horizon.name() << EOL;
                file_ << "GEOLOGICAL_TYPE " << horizon_map_.at( horizon.type() )
                      << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.horizon_items( horizon ) )
                {
                    file_ << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;

                find_corners_and_line_starts(
                    model_.horizon_items( horizon ), line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                file_ << "END" << EOL;
            }
            if( unclassified_surfaces_.empty() )
            {
                return;
            }
            file_ << "GOCAD TSurf 1" << EOL;
            geode::detail::HeaderData header;
            header.name = unclassified_surfaces_name_;
            geode::detail::write_header( file_, header );
            geode::detail::write_CRS( file_, {} );
            file_ << "GEOLOGICAL_FEATURE " << unclassified_surfaces_name_
                  << EOL;
            file_ << "GEOLOGICAL_TYPE "
                  << "boundary" << EOL;
            geode::index_t current_offset{ OFFSET_START };
            for( const auto& surface_id : unclassified_surfaces_ )
            {
                file_ << "TFACE" << EOL;
                current_offset = write_surface(
                    model_.surface( surface_id ), current_offset );
            }
            std::vector< std::array< geode::index_t, 2 > > line_starts;
            find_corners_and_line_starts_for_unclassified_surfaces(
                line_starts );
            write_corners( line_starts );
            write_line_starts( current_offset, line_starts );
            file_ << "END" << EOL;
        }

    private:
        std::ofstream file_;
        std::string model_name_;
        const geode::StructuralModel& model_;
        absl::flat_hash_map< geode::uuid, geode::index_t > components_;
        absl::flat_hash_map< geode::uuid, bool > universe_surface_sides_;
        absl::flat_hash_map< std::pair< geode::uuid, geode::uuid >, bool >
            regions_surface_sides_;
        const absl::flat_hash_map< geode::Fault3D::FAULT_TYPE, std::string >
            fault_map_ = { { geode::Fault3D::FAULT_TYPE::NO_TYPE, "fault" },
                { geode::Fault3D::FAULT_TYPE::NORMAL, "normal_fault" },
                { geode::Fault3D::FAULT_TYPE::REVERSE, "reverse_fault" },
                { geode::Fault3D::FAULT_TYPE::STRIKE_SLIP, "fault" },
                { geode::Fault3D::FAULT_TYPE::LISTRIC, "fault" },
                { geode::Fault3D::FAULT_TYPE::DECOLLEMENT, "fault" } };
        const absl::flat_hash_map< geode::Horizon3D::HORIZON_TYPE, std::string >
            horizon_map_ = { { geode::Horizon3D::HORIZON_TYPE::NO_TYPE,
                                 "none" },
                { geode::Horizon3D::HORIZON_TYPE::CONFORMAL, "top" },
                { geode::Horizon3D::HORIZON_TYPE::TOPOGRAPHY, "topographic" },
                { geode::Horizon3D::HORIZON_TYPE::INTRUSION, "intrusive" },
                { geode::Horizon3D::HORIZON_TYPE::NON_CONFORMAL,
                    "unconformity" } };
        geode::index_t component_id_{ OFFSET_START };
        const std::string unclassified_surfaces_name_{
            "unclassified_surfaces"
        };
        std::vector< geode::uuid > unclassified_surfaces_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        void MLOutput::write() const
        {
            const auto only_triangles =
                check_structural_model_polygons( structural_model() );
            if( !only_triangles )
            {
                geode::Logger::info(
                    "[MLOutput::write] Can not export into .ml a "
                    "StructuralModel with non triangular surface polygons." );
                return;
            }
            MLOutputImpl impl{ filename(), structural_model() };
            impl.determine_surface_to_regions_signs();
            impl.write_file();
        }
    } // namespace detail
} // namespace geode
