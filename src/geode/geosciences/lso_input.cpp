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

#include <geode/geosciences/private/lso_input.h>

#include <algorithm>
#include <fstream>
#include <functional>

#include <absl/container/flat_hash_set.h>
#include <absl/strings/numbers.h>
#include <absl/strings/str_split.h>

#include <geode/basic/attribute_manager.h>

#include <geode/geometry/bounding_box.h>
#include <geode/geometry/nn_search.h>
#include <geode/geometry/point.h>
#include <geode/geometry/vector.h>

#include <geode/mesh/builder/edged_curve_builder.h>
#include <geode/mesh/builder/point_set_builder.h>
#include <geode/mesh/builder/tetrahedral_solid_builder.h>
#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/geode_triangulated_surface.h>
#include <geode/mesh/core/mesh_factory.h>
#include <geode/mesh/core/solid_facets.h>
#include <geode/mesh/core/tetrahedral_solid.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/corner.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/mixin/core/vertex_identifier.h>
#include <geode/model/representation/io/brep_output.h>

#include <geode/geosciences/private/gocad_common.h>
#include <geode/geosciences/representation/builder/structural_model_builder.h>
#include <geode/geosciences/representation/core/structural_model.h>

namespace
{
    class LSOInputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };

        LSOInputImpl(
            absl::string_view filename, geode::StructuralModel& model )
            : file_( filename.data() ),
              model_( model ),
              builder_( model ),
              solid_{ geode::TetrahedralSolid3D::create() },
              solid_builder_{ geode::TetrahedralSolidBuilder3D::create(
                  *solid_ ) },
              vertex_id_{
                  solid_->vertex_attribute_manager()
                      .find_or_create_attribute< geode::VariableAttribute,
                          geode::index_t >( "vertex_id", geode::NO_ID )
              }
        {
            solid_->enable_facets();
            OPENGEODE_EXCEPTION( file_.good(),
                "[LSOInput] Error while opening file: ", filename );
        }

        void read_file()
        {
            geode::detail::check_keyword( file_, "GOCAD LightTSolid" );
            geode::detail::read_header( file_ );
            crs_ = geode::detail::read_CRS( file_ );
            read_vertices();
            read_tetra();
            read_surfaces();
            read_blocks();
            build_model_boundaries();
            build_corners();
            build_lines();
            geode::Logger::info(
                "[LSOInput] Properties are not supported yet" );
        }

    private:
        void read_vertices()
        {
            line_ = geode::detail::goto_keywords(
                file_, std::array< absl::string_view, 2 >{ "VRTX", "PVRTX" } );
            geode::index_t nb_unique_vertices{ 0 };
            do
            {
                geode::Point3D point;
                geode::index_t unique_id;
                if( geode::detail::string_starts_with( line_, "SHARED" ) )
                {
                    std::tie( point, unique_id ) = read_shared_point();
                }
                else
                {
                    unique_id = nb_unique_vertices;
                    point = read_point();
                    vertex_mapping_.emplace_back();
                    nb_unique_vertices++;
                }
                const auto id = solid_builder_->create_point( point );
                vertex_id_->set_value( id, unique_id );
                vertex_mapping_[unique_id].push_back( id );
            } while( std::getline( file_, line_ )
                     && absl::StrContains( line_, "VRTX" ) );
            builder_.create_unique_vertices( nb_unique_vertices );
        }

        std::tuple< geode::Point3D, geode::index_t > read_shared_point() const
        {
            const auto tokens = get_tokens();
            geode::index_t value;
            auto ok = absl::SimpleAtoi( tokens[2], &value );
            OPENGEODE_EXCEPTION( ok, "[LSOInput] Error while "
                                     "reading shared point index" );
            const auto unique_id = value - OFFSET_START;
            return std::make_tuple( solid_->point( unique_id ), unique_id );
        }

        geode::Point3D read_point() const
        {
            geode::Point3D point;
            const auto tokens = get_tokens();
            for( const auto i : geode::Range{ 3 } )
            {
                double value;
                auto ok = absl::SimpleAtod( tokens[i + 2], &value );
                OPENGEODE_EXCEPTION( ok, "[LSOInput] Error while "
                                         "reading point coordinates" );
                point.set_value( i, value );
            }
            point.set_value( 2, crs_.z_sign * point.value( 2 ) );
            return point;
        }

        void read_tetra()
        {
            do
            {
                const auto tokens = get_tokens();
                std::array< geode::index_t, 4 > vertices;
                for( const auto i : geode::Range{ 4 } )
                {
                    auto ok = absl::SimpleAtoi( tokens[i + 1], &vertices[i] );
                    OPENGEODE_EXCEPTION(
                        ok, "[LSOInput] Error while reading tetra" );
                    vertices[i] -= OFFSET_START;
                }
                solid_builder_->create_tetrahedron( vertices );
                std::getline( file_, line_ );
            } while( std::getline( file_, line_ )
                     && geode::detail::string_starts_with( line_, "TETRA" ) );
            solid_builder_->compute_polyhedron_adjacencies();
        }

        void read_surfaces()
        {
            if( !geode::detail::string_starts_with( line_, "MODEL" ) )
            {
                geode::detail::goto_keyword( file_, "MODEL" );
            }
            facet_id_ = solid_->facets()
                            .facet_attribute_manager()
                            .find_or_create_attribute< geode::VariableAttribute,
                                geode::uuid >( "facet_id", default_id_ );
            std::getline( file_, line_ );
            while( geode::detail::string_starts_with( line_, "SURFACE" ) )
            {
                std::istringstream iss{ line_ };
                std::string keyword;
                iss >> keyword;
                const auto h_id = builder_.add_horizon();
                builder_.set_horizon_name(
                    h_id, geode::detail::read_name( iss ) );
                const auto& horizon = model_.horizon( h_id );
                read_tfaces( horizon );
            }
        }

        void read_tfaces( const geode::Horizon3D& horizon )
        {
            std::getline( file_, line_ );
            while( geode::detail::string_starts_with( line_, "TFACE" ) )
            {
                const auto id =
                    builder_.add_surface( geode::MeshFactory::default_impl(
                        geode::TriangulatedSurface3D::type_name_static() ) );
                const auto& surface = model_.surface( id );
                builder_.add_surface_in_horizon( surface, horizon );
                surfaces_.push_back( id );
                surface_is_internal_[id] = false;
                read_key_vertices();
                read_triangles( id );
            }
        }

        void read_triangles( const geode::uuid& surface_id )
        {
            absl::flat_hash_map< geode::index_t, geode::index_t >
                vertex_mapping;
            auto builder =
                builder_.surface_mesh_builder< geode::TriangulatedSurface3D >(
                    surface_id );
            const auto component_id =
                model_.surface( surface_id ).component_id();
            while( std::getline( file_, line_ )
                   && geode::detail::string_starts_with( line_, "TRGL" ) )
            {
                const auto tokens = get_tokens();
                std::array< geode::index_t, 3 > facet_vertices;
                std::array< geode::index_t, 3 > vertices;
                for( const auto i : geode::Range{ 3 } )
                {
                    geode::index_t value;
                    auto ok = absl::SimpleAtoi( tokens[i + 1], &value );
                    OPENGEODE_EXCEPTION(
                        ok, "[LSOInput] Error while reading triangles" );
                    value -= OFFSET_START;
                    facet_vertices[i] = value;
                    const auto it = vertex_mapping.find( value );
                    if( it != vertex_mapping.end() )
                    {
                        vertices[i] = it->second;
                    }
                    else
                    {
                        const auto vertex_id =
                            builder->create_point( solid_->point( value ) );
                        vertex_mapping.emplace( value, vertex_id );
                        vertices[i] = vertex_id;
                        builder_.set_unique_vertex( { component_id, vertex_id },
                            vertex_id_->value( value ) );
                    }
                }
                builder->create_triangle( vertices );
                for( const auto& facet :
                    facets_from_key_vertices( facet_vertices ) )
                {
                    facet_id_->set_value( facet.first, surface_id );
                }
            }
            builder->compute_polygon_adjacencies();
        }

        void read_key_vertices()
        {
            std::getline( file_, line_ );
            OPENGEODE_EXCEPTION(
                geode::detail::string_starts_with( line_, "KEYVERTICES" ),
                "[LSOInput] Error while reading key vertices" );
            const auto key_tokens = get_tokens();
            std::array< geode::index_t, 3 > key;
            for( const auto i : geode::Range{ 3 } )
            {
                geode::index_t value;
                auto ok = absl::SimpleAtoi( key_tokens[i + 1], &value );
                OPENGEODE_EXCEPTION(
                    ok, "[LSOInput] Error while reading key vertices" );
                key[i] = value - OFFSET_START;
            }
            surface_keys_.emplace_back( std::move( key ) );
        }

        std::vector<
            std::pair< geode::index_t, std::array< geode::index_t, 3 > > >
            facets_from_key_vertices(
                const std::array< geode::index_t, 3 >& vertices ) const
        {
            std::vector<
                std::pair< geode::index_t, std::array< geode::index_t, 3 > > >
                facets;
            for( const auto v0 : vertex_mapping_[vertices[0]] )
            {
                for( const auto v1 : vertex_mapping_[vertices[1]] )
                {
                    for( const auto v2 : vertex_mapping_[vertices[2]] )
                    {
                        if( const auto facet_id =
                                solid_->facets().facet_from_vertices(
                                    { v0, v1, v2 } ) )
                        {
                            facets.emplace_back( facet_id.value(),
                                std::array< geode::index_t, 3 >{
                                    { v0, v1, v2 } } );
                        }
                    }
                }
            }
            return facets;
        }

        void read_blocks()
        {
            while( geode::detail::string_starts_with( line_, "MODEL_REGION" ) )
            {
                const auto tokens = get_tokens();
                int value;
                auto ok = absl::SimpleAtoi( tokens[2], &value );
                OPENGEODE_EXCEPTION(
                    ok, "[LSOInput] Error while reading tetra" );
                const auto side = value > 0;
                const auto block_id =
                    builder_.add_block( geode::MeshFactory::default_impl(
                        geode::TetrahedralSolid3D::type_name_static() ) );
                builder_.set_block_name( block_id, tokens[1] );
                geode::index_t surface_id = std::abs( value ) - OFFSET_START;
                build_block( block_id, surface_keys_[surface_id], side );
            }
        }

        void build_model_boundaries()
        {
            const auto& id = builder_.add_model_boundary();
            const auto& boundary = model_.model_boundary( id );
            for( const auto& surface : model_.surfaces() )
            {
                if( model_.nb_incidences( surface.id() ) == 1 )
                {
                    builder_.add_surface_in_model_boundary( surface, boundary );
                }
            }
        }

        geode::PolyhedronFacet find_facet(
            const std::array< geode::index_t, 3 >& key_vertices,
            bool side ) const
        {
            for( const auto& facet : facets_from_key_vertices( key_vertices ) )
            {
                for( const auto& polyhedron : solid_->polyhedra_from_facet(
                         solid_->facets().facet_vertices( facet.first ) ) )
                {
                    const auto& key_vertices = facet.second;
                    for( const auto f : geode::Range{
                             solid_->nb_polyhedron_facets( polyhedron ) } )
                    {
                        const geode::PolyhedronFacet polyhedron_facet{
                            polyhedron, f
                        };
                        if( solid_->facets().facet_from_vertices(
                                solid_->polyhedron_facet_vertices(
                                    polyhedron_facet ) )
                            != facet.first )
                        {
                            continue;
                        }
                        if( !facet_matches_key_vertices(
                                key_vertices, polyhedron_facet, side ) )
                        {
                            continue;
                        }
                        return polyhedron_facet;
                    }
                }
            }
            throw geode::OpenGeodeException{
                "[LSOInput] Cannot find starting facet"
            };
            return {};
        }

        std::tuple< absl::flat_hash_set< geode::index_t >,
            absl::flat_hash_map< geode::uuid, geode::index_t > >
            build_tetra( geode::index_t first_tetra )
        {
            absl::flat_hash_map< geode::uuid, geode::index_t >
                surface_relations;
            std::stack< geode::index_t > tetras;
            absl::flat_hash_set< geode::index_t > visited;
            tetras.emplace( first_tetra );
            visited.emplace( first_tetra );
            while( !tetras.empty() )
            {
                const auto tetra = tetras.top();
                tetras.pop();
                for( const auto f : geode::Range{ 4 } )
                {
                    const auto facet_id =
                        solid_->facets()
                            .facet_from_vertices(
                                solid_->polyhedron_facet_vertices(
                                    { tetra, f } ) )
                            .value();
                    const auto& facet_uuid = facet_id_->value( facet_id );
                    if( facet_uuid != default_id_ )
                    {
                        const auto it = surface_relations.find( facet_uuid );
                        if( it != surface_relations.end() )
                        {
                            it->second++;
                        }
                        else
                        {
                            surface_relations.emplace( facet_uuid, 1 );
                        }
                    }
                    else if( const auto adj =
                                 solid_->polyhedron_adjacent( { tetra, f } ) )
                    {
                        if( visited.emplace( adj.value() ).second )
                        {
                            tetras.emplace( adj.value() );
                        }
                    }
                }
            }
            return std::make_tuple( visited, surface_relations );
        }

        void build_solid( const geode::uuid& block_id,
            const absl::flat_hash_set< geode::index_t >& tetras )
        {
            auto builder =
                builder_.block_mesh_builder< geode::TetrahedralSolid3D >(
                    block_id );
            const auto component_id = model_.block( block_id ).component_id();
            absl::flat_hash_map< geode::index_t, geode::index_t >
                vertex_mapping;
            for( const auto tetra : tetras )
            {
                std::array< geode::index_t, 4 > vertices;
                for( const auto i : geode::Range{ 4 } )
                {
                    const auto vertex =
                        solid_->polyhedron_vertex( { tetra, i } );
                    const auto it = vertex_mapping.find( vertex );
                    if( it != vertex_mapping.end() )
                    {
                        vertices[i] = it->second;
                    }
                    else
                    {
                        const auto vertex_id =
                            builder->create_point( solid_->point( vertex ) );
                        vertex_mapping.emplace( vertex, vertex_id );
                        vertices[i] = vertex_id;
                        builder_.set_unique_vertex( { component_id, vertex_id },
                            vertex_id_->value( vertex ) );
                    }
                }
                builder->create_tetrahedron( vertices );
            }
            builder->compute_polyhedron_adjacencies();
        }

        void build_block_relations( const geode::uuid& block_id,
            const absl::flat_hash_map< geode::uuid, geode::index_t >&
                relations )
        {
            const auto& block = model_.block( block_id );
            for( const auto& relation : relations )
            {
                const auto& surface = model_.surface( relation.first );
                const auto nb_relations =
                    relation.second / surface.mesh().nb_polygons();
                if( nb_relations == 1 )
                {
                    builder_.add_surface_block_boundary_relationship(
                        surface, block );
                }
                else
                {
                    OPENGEODE_ASSERT( nb_relations == 2,
                        "[LSOInput] Error in Surface/Block "
                        "relations" );
                    builder_.add_surface_block_internal_relationship(
                        surface, block );
                }
            }
        }

        void build_block( const geode::uuid& block_id,
            const std::array< geode::index_t, 3 >& key_vertices,
            bool side )
        {
            const auto& block = model_.block( block_id );
            const auto component_id = block.component_id();
            const auto polyhedron_facet = find_facet( key_vertices, side );
            const auto info = build_tetra( polyhedron_facet.polyhedron_id );
            build_solid( block_id, std::get< 0 >( info ) );
            build_block_relations( block_id, std::get< 1 >( info ) );
            std::getline( file_, line_ );
        }

        bool facet_matches_key_vertices(
            const std::array< geode::index_t, 3 >& key_vertices,
            const geode::PolyhedronFacet& facet,
            bool side ) const
        {
            const auto pos = static_cast< geode::index_t >(
                std::distance( key_vertices.begin(),
                    absl::c_find( key_vertices,
                        solid_->polyhedron_facet_vertex( { facet, 0 } ) ) ) );
            if( side )
            {
                for( const auto v : geode::Range{ 3 } )
                {
                    const auto v_id =
                        solid_->polyhedron_facet_vertex( { facet, v } );
                    if( v_id != key_vertices[( pos + v ) % 3] )
                    {
                        return false;
                    }
                }
            }
            else
            {
                for( const auto v : geode::Range{ 3 } )
                {
                    const auto v_id =
                        solid_->polyhedron_facet_vertex( { facet, v } );
                    if( v_id != key_vertices[( pos + 3 - v ) % 3] )
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        void build_corners()
        {
            for( const auto& surface : model_.surfaces() )
            {
                const auto& mesh = surface.mesh();
                const auto component_id = surface.component_id();
                for( const auto p : geode::Range{ mesh.nb_polygons() } )
                {
                    for( const auto& border :
                        mesh.polygon_edges_on_border( p ) )
                    {
                        build_corner_from_edge( mesh, component_id, border );
                    }
                }
            }
        }

        void build_corner_from_edge( const geode::SurfaceMesh3D& mesh,
            const geode::ComponentID& component_id,
            const geode::PolygonEdge& border )
        {
            const auto vertex0 = mesh.polygon_edge_vertex( border, 0 );
            const auto unique_id0 =
                model_.unique_vertex( { component_id, vertex0 } );
            const auto surfaces0 = model_.mesh_component_vertices(
                unique_id0, geode::Surface3D::component_type_static() );

            const auto vertex1 = mesh.polygon_edge_vertex( border, 1 );
            const auto unique_id1 =
                model_.unique_vertex( { component_id, vertex1 } );
            const auto surfaces1 = model_.mesh_component_vertices(
                unique_id1, geode::Surface3D::component_type_static() );
            if( surfaces0.size() > surfaces1.size() )
            {
                create_corner( mesh.point( vertex0 ), unique_id0 );
            }
            else if( surfaces0.size() < surfaces1.size() )
            {
                create_corner( mesh.point( vertex1 ), unique_id1 );
            }
            else
            {
                if( !are_equal( surfaces0, surfaces1 ) )
                {
                    create_corner( mesh.point( vertex0 ), unique_id0 );
                    create_corner( mesh.point( vertex1 ), unique_id1 );
                }
            }
        }

        bool are_equal( const std::vector< geode::MeshComponentVertex >& lhs,
            const std::vector< geode::MeshComponentVertex >& rhs ) const
        {
            for( const auto& v0 : lhs )
            {
                bool found{ false };
                for( const auto& v1 : rhs )
                {
                    if( v0.component_id.id() == v1.component_id.id() )
                    {
                        found = true;
                        break;
                    }
                }
                if( !found )
                {
                    return false;
                }
            }
            return true;
        }

        void create_corner( const geode::Point3D& point, geode::index_t id )
        {
            if( !model_.has_mesh_component_vertices(
                    id, geode::Corner3D::component_type_static() ) )
            {
                const auto& corner_id = builder_.add_corner();
                auto builder = builder_.corner_mesh_builder( corner_id );
                builder->create_point( point );
                const auto& corner = model_.corner( corner_id );
                builder_.set_unique_vertex( { corner.component_id(), 0 }, id );
            }
        }

        std::vector< absl::string_view > get_tokens() const
        {
            return absl::StrSplit( line_, ' ' );
        }

        void build_lines()
        {
            for( const auto& surface : model_.surfaces() )
            {
                const auto& mesh = surface.mesh();
                absl::flat_hash_map< geode::uuid, geode::index_t >
                    line_relations;
                for( const auto p : geode::Range{ mesh.nb_polygons() } )
                {
                    for( const auto& border :
                        mesh.polygon_edges_on_border( p ) )
                    {
                        const auto vertex_id0 =
                            mesh.polygon_edge_vertex( border, 0 );
                        const auto unique_id0 = model_.unique_vertex(
                            { surface.component_id(), vertex_id0 } );
                        if( !model_.has_mesh_component_vertices( unique_id0,
                                geode::Corner3D::component_type_static() ) )
                        {
                            continue;
                        }
                        const auto vertex_id1 =
                            mesh.polygon_edge_vertex( border, 1 );
                        const auto unique_id1 = model_.unique_vertex(
                            { surface.component_id(), vertex_id1 } );
                        if( const auto line_id =
                                common_line( unique_id0, unique_id1 ) )
                        {
                            const auto it =
                                line_relations.find( line_id.value() );
                            if( it != line_relations.end() )
                            {
                                it->second++;
                            }
                            else
                            {
                                line_relations.emplace( line_id.value(), 1 );
                            }
                            continue;
                        }
                        build_line( surface, border, line_relations );
                    }
                }
                build_line_relations( surface, line_relations );
            }
        }

        void build_line_relations( const geode::Surface3D& surface,
            absl::flat_hash_map< geode::uuid, geode::index_t >& relations )
        {
            for( const auto& relation : relations )
            {
                const auto& line = model_.line( relation.first );
                if( relation.second == 1 )
                {
                    builder_.add_line_surface_boundary_relationship(
                        line, surface );
                }
                else
                {
                    OPENGEODE_ASSERT( relation.second == 2,
                        "[LSOInput] Error in Line/Surface relations" );
                    builder_.add_line_surface_internal_relationship(
                        line, surface );
                }
            }
        }

        void build_line( const geode::Surface3D& surface,
            const geode::PolygonEdge& border,
            absl::flat_hash_map< geode::uuid, geode::index_t >& line_relations )
        {
            const auto& mesh = surface.mesh();
            const auto& line_id = builder_.add_line();
            const auto& line = model_.line( line_id );
            line_relations.emplace( line_id, 1 );
            auto builder = builder_.line_mesh_builder( line_id );
            const auto vertex_id0 = mesh.polygon_edge_vertex( border, 0 );
            auto v_id = builder->create_point( mesh.point( vertex_id0 ) );
            const auto unique_id0 =
                model_.unique_vertex( { surface.component_id(), vertex_id0 } );
            builder_.set_unique_vertex(
                { line.component_id(), v_id }, unique_id0 );
            auto edge = border;
            auto unique_id = unique_id0;
            const auto corner_type = geode::Corner3D::component_type_static();
            do
            {
                edge = mesh.next_on_border( edge );
                const auto vertex_id = mesh.polygon_edge_vertex( edge, 0 );
                unique_id = model_.unique_vertex(
                    { surface.component_id(), vertex_id } );
                v_id = builder->create_point( mesh.point( vertex_id ) );
                builder_.set_unique_vertex(
                    { line.component_id(), v_id }, unique_id );
                builder->create_edge( v_id - 1, v_id );
            } while(
                !model_.has_mesh_component_vertices( unique_id, corner_type ) );
            const auto corner0 =
                model_.mesh_component_vertices( unique_id0, corner_type )[0];
            const auto corner1 =
                model_.mesh_component_vertices( unique_id, corner_type )[0];
            builder_.add_corner_line_boundary_relationship(
                model_.corner( corner0.component_id.id() ), line );
            builder_.add_corner_line_boundary_relationship(
                model_.corner( corner1.component_id.id() ), line );
        }

        absl::optional< geode::uuid > common_line(
            geode::index_t unique_id0, geode::index_t unique_id1 )
        {
            const auto lines0 = model_.mesh_component_vertices(
                unique_id0, geode::Line3D::component_type_static() );
            const auto lines1 = model_.mesh_component_vertices(
                unique_id1, geode::Line3D::component_type_static() );
            for( const auto& mcv0 : lines0 )
            {
                for( const auto& mcv1 : lines1 )
                {
                    if( mcv0.component_id == mcv1.component_id )
                    {
                        const auto min = std::min( mcv0.vertex, mcv1.vertex );
                        const auto max = std::max( mcv0.vertex, mcv1.vertex );
                        if( max - min == 1 )
                        {
                            return mcv0.component_id.id();
                        }
                    }
                }
            }
            return absl::nullopt;
        }

    private:
        std::ifstream file_;
        std::string line_;
        geode::StructuralModel& model_;
        geode::StructuralModelBuilder builder_;
        geode::detail::CRSData crs_;
        std::unique_ptr< geode::TetrahedralSolid3D > solid_;
        std::unique_ptr< geode::TetrahedralSolidBuilder3D > solid_builder_;
        std::shared_ptr< geode::VariableAttribute< geode::index_t > >
            vertex_id_;
        std::shared_ptr< geode::VariableAttribute< geode::uuid > > facet_id_;
        geode::uuid default_id_;
        std::vector< geode::uuid > surfaces_;
        absl::flat_hash_map< geode::uuid, bool > surface_is_internal_;
        std::vector< std::array< geode::index_t, 3 > > surface_keys_;
        std::vector< absl::InlinedVector< geode::index_t, 1 > > vertex_mapping_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        void LSOInput::read()
        {
            LSOInputImpl impl{ filename(), structural_model() };
            impl.read_file();
        }
    } // namespace detail
} // namespace geode
