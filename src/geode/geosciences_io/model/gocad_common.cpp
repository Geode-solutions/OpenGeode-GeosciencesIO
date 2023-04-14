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

#include <geode/geosciences_io/model/private/gocad_common.h>

#include <queue>

#include <geode/geometry/basic_objects/tetrahedron.h>
#include <geode/geometry/bounding_box.h>
#include <geode/geometry/mensuration.h>

#include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/surface_mesh.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/representation/core/brep.h>

namespace
{
    static constexpr char EOL{ '\n' };
    static constexpr char SPACE{ ' ' };

    using PairedSigns =
        absl::flat_hash_map< std::pair< geode::uuid, geode::uuid >, bool >;

    PairedSigns determine_paired_signs( const geode::BRep& brep )
    {
        PairedSigns paired_signs;
        paired_signs.reserve( brep.nb_lines() * 2 );
        for( const auto& line : brep.lines() )
        {
            const auto& mesh = line.mesh();
            const auto uid0 = brep.unique_vertex(
                { line.component_id(), mesh.edge_vertex( { 0, 0 } ) } );
            const auto uid1 = brep.unique_vertex(
                { line.component_id(), mesh.edge_vertex( { 0, 1 } ) } );
            const auto surface_cmvs0 = brep.component_mesh_vertices(
                uid0, geode::Surface3D::component_type_static() );
            const auto surface_cmvs1 = brep.component_mesh_vertices(
                uid1, geode::Surface3D::component_type_static() );
            absl::flat_hash_map< geode::uuid, bool > surface_direct_edges;
            surface_direct_edges.reserve(
                std::min( surface_cmvs0.size(), surface_cmvs1.size() ) );
            for( const auto& surface_cmv0 : surface_cmvs0 )
            {
                for( const auto& surface_cmv1 : surface_cmvs1 )
                {
                    if( surface_cmv1.component_id.id()
                        != surface_cmv0.component_id.id() )
                    {
                        continue;
                    }
                    const auto& surface =
                        brep.surface( surface_cmv0.component_id.id() );
                    const auto& surface_mesh = surface.mesh();
                    const auto v0v1 = surface_mesh.polygon_edge_from_vertices(
                        surface_cmv0.vertex, surface_cmv1.vertex );
                    const auto v1v0 = surface_mesh.polygon_edge_from_vertices(
                        surface_cmv1.vertex, surface_cmv0.vertex );
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
            if( surface_direct_edges.size() < 2 )
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
                            std::make_pair( s0.first, s1.first ),
                            s0.second != s1.second );
                    }
                }
            }
        }
        return paired_signs;
    }

    absl::FixedArray< bool > determine_relative_signs(
        absl::Span< const geode::uuid > universe_boundaries,
        const PairedSigns& paired_signs )
    {
        const auto nb_surfaces = universe_boundaries.size();
        if( nb_surfaces == 1 )
        {
            return { true };
        }
        absl::FixedArray< bool > signs( nb_surfaces, true );
        absl::FixedArray< bool > determined( nb_surfaces, false );
        determined[0] = true;
        std::queue< geode::index_t > to_process;
        to_process.push( 0 );
        while( !to_process.empty() )
        {
            const auto determined_s = to_process.front();
            const auto& determined_s_id = universe_boundaries[determined_s];
            to_process.pop();
            for( const auto s : geode::Range{ nb_surfaces } )
            {
                if( determined[s] )
                {
                    continue;
                }
                const auto& s_id = universe_boundaries[s];
                const auto itr =
                    paired_signs.find( { std::min( determined_s_id, s_id ),
                        std::max( determined_s_id, s_id ) } );
                if( itr == paired_signs.end() )
                {
                    continue;
                }
                signs[s] =
                    itr->second ? signs[determined_s] : !signs[determined_s];
                determined[s] = true;
                to_process.push( s );
            }
        }
        OPENGEODE_ASSERT( absl::c_count( determined, false ) == 0,
            "All signs should have been found" );
        return signs;
    }

    bool are_correct_sides( const geode::BRep& brep,
        absl::Span< const geode::uuid > universe_boundaries,
        absl::Span< const bool > relative_signs )
    {
        double signed_volume{ 0 };
        const auto& first_surface_mesh =
            brep.surface( universe_boundaries[0] ).mesh();
        const auto& bbox = first_surface_mesh.bounding_box();
        const auto center = ( bbox.min() + bbox.max() ) * 0.5;
        for( const auto s : geode::Indices{ universe_boundaries } )
        {
            const auto sign = relative_signs[s];
            const auto& surface_mesh =
                brep.surface( universe_boundaries[s] ).mesh();
            for( const auto t : geode::Range{ surface_mesh.nb_polygons() } )
            {
                std::array< geode::local_index_t, 3 > vertex_order{ 0, 1, 2 };
                if( !sign )
                {
                    std::swap( vertex_order[1], vertex_order[2] );
                }
                signed_volume += geode::tetrahedron_signed_volume(
                    { surface_mesh.point( surface_mesh.polygon_vertex(
                          { t, vertex_order[0] } ) ),
                        surface_mesh.point( surface_mesh.polygon_vertex(
                            { t, vertex_order[1] } ) ),
                        surface_mesh.point( surface_mesh.polygon_vertex(
                            { t, vertex_order[2] } ) ),
                        center } );
            }
        }
        OPENGEODE_ASSERT(
            std::fabs( signed_volume ) > 0, "Null volume block is not valid" );
        return signed_volume > 0;
    }

    absl::flat_hash_map< geode::uuid, bool > determine_universe_sides(
        const geode::BRep& brep, const PairedSigns& paired_signs )
    {
        std::vector< geode::uuid > universe_boundaries;
        universe_boundaries.reserve( 2 * brep.nb_model_boundaries() );
        for( const auto& boundary : brep.model_boundaries() )
        {
            for( const auto& item : brep.model_boundary_items( boundary ) )
            {
                universe_boundaries.push_back( item.id() );
            }
        }
        const auto relative_signs =
            determine_relative_signs( universe_boundaries, paired_signs );
        const auto correct =
            are_correct_sides( brep, universe_boundaries, relative_signs );

        absl::flat_hash_map< geode::uuid, bool > sides;
        sides.reserve( universe_boundaries.size() );
        for( const auto b : geode::Indices{ universe_boundaries } )
        {
            sides.emplace( universe_boundaries[b],
                correct ? !relative_signs[b] : relative_signs[b] );
        }
        return sides;
    }

    PairedSigns determine_regions_sides(
        const geode::BRep& brep, const PairedSigns& paired_signs )
    {
        PairedSigns sides;
        for( const auto& block : brep.blocks() )
        {
            std::vector< geode::uuid > block_boundaries;
            block_boundaries.reserve( brep.nb_boundaries( block.id() ) );
            for( const auto& surface : brep.boundaries( block ) )
            {
                block_boundaries.push_back( surface.id() );
            }
            const auto relative_signs =
                determine_relative_signs( block_boundaries, paired_signs );
            const auto correct =
                are_correct_sides( brep, block_boundaries, relative_signs );

            for( const auto b : geode::Indices{ block_boundaries } )
            {
                sides.emplace(
                    std::make_pair( block.id(), block_boundaries[b] ),
                    correct ? relative_signs[b] : !relative_signs[b] );
            }
        }
        return sides;
    }
} // namespace

namespace geode
{
    namespace detail
    {
        RegionSurfaceSide determine_surface_to_regions_sides( const BRep& brep )
        {
            RegionSurfaceSide result;
            const auto paired_signs = determine_paired_signs( brep );
            result.universe_surface_sides =
                determine_universe_sides( brep, paired_signs );
            result.regions_surface_sides =
                determine_regions_sides( brep, paired_signs );
            return result;
        }
    } // namespace detail
} // namespace geode
