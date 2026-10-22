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

#include <geode/geosciences_io/model/internal/gocad_common.hpp>

#include <cmath>
#include <optional>
#include <queue>

#include <absl/container/fixed_array.h>

#include <geode/geometry/basic_objects/tetrahedron.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/mensuration.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/line.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/core/brep.hpp>

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
            const auto& cmvs0 = brep.component_mesh_vertices( uid0 );
            const auto& cmvs1 = brep.component_mesh_vertices( uid1 );
            absl::flat_hash_map< geode::uuid, bool > surface_direct_edges;
            surface_direct_edges.reserve(
                std::min( cmvs0.size(), cmvs1.size() ) );
            for( const auto& cmv0 : cmvs0 )
            {
                if( cmv0.component_id.type
                    != geode::Surface3D::component_type_static() )
                {
                    continue;
                }
                for( const auto& cmv1 : cmvs1 )
                {
                    if( cmv1.component_id.type
                        != geode::Surface3D::component_type_static() )
                    {
                        continue;
                    }
                    if( cmv1.component_id.id != cmv0.component_id.id )
                    {
                        continue;
                    }
                    const auto& surface = brep.surface( cmv0.component_id.id );
                    const auto& surface_mesh = surface.mesh();
                    const auto v0v1 = surface_mesh.polygon_edge_from_vertices(
                        cmv0.vertex, cmv1.vertex );
                    const auto v1v0 = surface_mesh.polygon_edge_from_vertices(
                        cmv1.vertex, cmv0.vertex );
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
                            std::pair< geode::uuid, geode::uuid >{
                                s0.first, s1.first },
                            s0.second != s1.second );
                    }
                }
            }
        }
        return paired_signs;
    }

    struct RelativeSigns
    {
        explicit RelativeSigns( geode::index_t nb_surfaces )
            : signs( nb_surfaces, true ),
              components( nb_surfaces, geode::NO_ID )
        {
        }

        absl::FixedArray< bool > signs;
        absl::FixedArray< geode::index_t > components;
        geode::index_t nb_components{ 0 };
    };

    /*!
     * Signs relative to the first Surface of each connected component of
     * Surfaces (connected through shared Lines)
     */
    RelativeSigns determine_relative_signs(
        absl::Span< const geode::uuid > boundaries,
        const PairedSigns& paired_signs )
    {
        const auto nb_surfaces = boundaries.size();
        RelativeSigns result{ static_cast< geode::index_t >( nb_surfaces ) };
        for( const auto start : geode::Range{ nb_surfaces } )
        {
            if( result.components[start] != geode::NO_ID )
            {
                continue;
            }
            const auto component = result.nb_components++;
            result.components[start] = component;
            std::queue< geode::index_t > to_process;
            to_process.push( start );
            while( !to_process.empty() )
            {
                const auto determined_s = to_process.front();
                const auto& determined_s_id = boundaries[determined_s];
                to_process.pop();
                for( const auto s : geode::Range{ nb_surfaces } )
                {
                    if( result.components[s] != geode::NO_ID )
                    {
                        continue;
                    }
                    const auto& s_id = boundaries[s];
                    const auto itr =
                        paired_signs.find( { std::min( determined_s_id, s_id ),
                            std::max( determined_s_id, s_id ) } );
                    if( itr == paired_signs.end() )
                    {
                        continue;
                    }
                    result.signs[s] = itr->second ? result.signs[determined_s]
                                                  : !result.signs[determined_s];
                    result.components[s] = component;
                    to_process.push( s );
                }
            }
        }
        return result;
    }

    std::array< geode::Point3D, 3 > oriented_triangle(
        const geode::SurfaceMesh3D& mesh, geode::index_t polygon, bool sign )
    {
        std::array< geode::local_index_t, 3 > vertex_order{ 0, 1, 2 };
        if( !sign )
        {
            std::swap( vertex_order[1], vertex_order[2] );
        }
        return { mesh.point(
                     mesh.polygon_vertex( { polygon, vertex_order[0] } ) ),
            mesh.point( mesh.polygon_vertex( { polygon, vertex_order[1] } ) ),
            mesh.point( mesh.polygon_vertex( { polygon, vertex_order[2] } ) ) };
    }

    // Signed solid angle of a triangle seen from a point (Van Oosterom and
    // Strackee formula)
    double solid_angle( const std::array< geode::Point3D, 3 >& triangle,
        const geode::Point3D& point )
    {
        const geode::Vector3D a{ point, triangle[0] };
        const geode::Vector3D b{ point, triangle[1] };
        const geode::Vector3D c{ point, triangle[2] };
        const auto la = a.length();
        const auto lb = b.length();
        const auto lc = c.length();
        const auto numerator = a.dot( b.cross( c ) );
        const auto denominator =
            la * lb * lc + a.dot( b ) * lc + a.dot( c ) * lb + b.dot( c ) * la;
        return 2. * std::atan2( numerator, denominator );
    }

    /*!
     * Signs giving, for each connected component (shell), the orientation
     * pointing out of the volume bounded by the Surfaces: outer shells
     * enclose a positive volume, shells nested inside another one
     * (inclusions) enclose a negative volume.
     */
    absl::FixedArray< bool > determine_outward_signs( const geode::BRep& brep,
        absl::Span< const geode::uuid > boundaries,
        const PairedSigns& paired_signs )
    {
        const auto relative =
            determine_relative_signs( boundaries, paired_signs );
        const auto& bbox = brep.surface( boundaries[0] ).mesh().bounding_box();
        const auto center = ( bbox.min() + bbox.max() ) * 0.5;
        absl::FixedArray< double > volumes( relative.nb_components, 0 );
        absl::FixedArray< std::optional< geode::Point3D > > inner_points(
            relative.nb_components );
        for( const auto s : geode::Indices{ boundaries } )
        {
            const auto component = relative.components[s];
            const auto& mesh = brep.surface( boundaries[s] ).mesh();
            for( const auto t : geode::Range{ mesh.nb_polygons() } )
            {
                const auto triangle =
                    oriented_triangle( mesh, t, relative.signs[s] );
                volumes[component] += geode::tetrahedron_signed_volume(
                    { triangle[0], triangle[1], triangle[2], center } );
                if( !inner_points[component] )
                {
                    inner_points[component] =
                        ( triangle[0] + triangle[1] + triangle[2] ) / 3.;
                }
            }
        }
        for( const auto volume : volumes )
        {
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                std::fabs( volume ) > 0, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "Null volume block is not valid" );
        }
        absl::FixedArray< bool > flip( relative.nb_components, false );
        for( const auto component : geode::Range{ relative.nb_components } )
        {
            geode::index_t depth{ 0 };
            if( relative.nb_components > 1 )
            {
                absl::FixedArray< double > windings(
                    relative.nb_components, 0 );
                for( const auto s : geode::Indices{ boundaries } )
                {
                    const auto other = relative.components[s];
                    if( other == component )
                    {
                        continue;
                    }
                    const auto& mesh = brep.surface( boundaries[s] ).mesh();
                    for( const auto t : geode::Range{ mesh.nb_polygons() } )
                    {
                        windings[other] += solid_angle(
                            oriented_triangle( mesh, t, relative.signs[s] ),
                            inner_points[component].value() );
                    }
                }
                for( const auto winding : windings )
                {
                    // Generalized winding number: +-1 inside, 0 outside
                    if( std::fabs( winding ) > 2. * M_PI )
                    {
                        depth++;
                    }
                }
            }
            const auto outer_shell = depth % 2 == 0;
            flip[component] = ( volumes[component] > 0 ) != outer_shell;
        }
        absl::FixedArray< bool > signs( boundaries.size() );
        for( const auto s : geode::Indices{ boundaries } )
        {
            signs[s] = relative.signs[s] != flip[relative.components[s]];
        }
        return signs;
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
        const auto outward_signs =
            determine_outward_signs( brep, universe_boundaries, paired_signs );

        absl::flat_hash_map< geode::uuid, bool > sides;
        sides.reserve( universe_boundaries.size() );
        for( const auto b : geode::Indices{ universe_boundaries } )
        {
            sides.emplace( universe_boundaries[b], !outward_signs[b] );
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
            const auto outward_signs =
                determine_outward_signs( brep, block_boundaries, paired_signs );
            for( const auto b : geode::Indices{ block_boundaries } )
            {
                sides.emplace(
                    std::pair< geode::uuid, geode::uuid >{
                        block.id(), block_boundaries[b] },
                    outward_signs[b] );
            }
        }
        return sides;
    }
} // namespace

namespace geode::internal
{
    RegionSurfaceSide determine_surface_to_regions_sides( const BRep& brep )
    {
        OpenGeodeGeosciencesIOModelException::check_exception(
            brep.nb_model_boundaries() > 0, nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[determine_surface_to_regions_sides] Requires ModelBoundaries "
            "to be defined" );
        RegionSurfaceSide result;
        const auto paired_signs = determine_paired_signs( brep );
        result.universe_surface_sides =
            determine_universe_sides( brep, paired_signs );
        result.regions_surface_sides =
            determine_regions_sides( brep, paired_signs );
        return result;
    }
} // namespace geode::internal
