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

#include <geode/io/geosciences/private/lso_input.h>

#include <fstream>

#include <absl/container/flat_hash_set.h>
#include <absl/strings/match.h>

#include <geode/basic/attribute_manager.h>
#include <geode/basic/string.h>

#include <geode/geometry/point.h>

#include <geode/mesh/builder/edged_curve_builder.h>
#include <geode/mesh/builder/point_set_builder.h>
#include <geode/mesh/builder/tetrahedral_solid_builder.h>
#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/mesh_factory.h>
#include <geode/mesh/core/solid_facets.h>
#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/mesh/core/triangulated_surface.h>

#include <geode/model/helpers/component_mesh_edges.h>
#include <geode/model/helpers/detail/cut_along_internal_lines.h>
#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/corner.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/surface.h>

#include <geode/geosciences/representation/builder/structural_model_builder.h>
#include <geode/geosciences/representation/core/structural_model.h>
#include <geode/io/geosciences/private/gocad_common.h>

namespace
{
    class LSOInputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr auto block_name_attribute_name =
            "geode_block_name_attribute_name";

        LSOInputImpl(
            absl::string_view filename, geode::StructuralModel& model )
            : file_{ geode::to_string( filename ) },
              model_( model ),
              builder_{ model },
              solid_{ geode::TetrahedralSolid3D::create() },
              solid_builder_{ geode::TetrahedralSolidBuilder3D::create(
                  *solid_ ) },
              vertex_id_{
                  solid_->vertex_attribute_manager()
                      .find_or_create_attribute< geode::VariableAttribute,
                          geode::index_t >( "vertex_id", geode::NO_ID )
              },
              block_name_attribute_{
                  solid_->polyhedron_attribute_manager()
                      .find_or_create_attribute< geode::VariableAttribute,
                          std::string >( block_name_attribute_name, "" )
              }
        {
            solid_->enable_facets();
            OPENGEODE_EXCEPTION( file_.good(),
                "[LSOInput] Error while opening file: ", filename );
        }

        bool read_file()
        {
            if( !geode::detail::goto_keyword_if_it_exists(
                    file_, "GOCAD LightTSolid" ) )
            {
                throw geode::OpenGeodeException{
                    "[LSOInput] Cannot find LightTSolid in the file"
                };
            }
            const auto header = geode::detail::read_header( file_ );
            builder_.set_name( header.name );
            crs_ = geode::detail::read_CRS( file_ );
            vertices_prop_header_ =
                geode::detail::read_prop_header( file_, "" );
            tetrahedra_prop_header_ =
                geode::detail::read_prop_header( file_, "TETRA_" );
            vertices_attributes_.resize( vertices_prop_header_.names.size() );
            tetrahedra_attributes_.resize(
                tetrahedra_prop_header_.names.size() );
            read_vertices();
            read_vertex_region_indicators();
            read_tetrahedra();
            read_tetrahedra_region_indicators();
            read_surfaces();
            read_blocks();
            build_model_boundaries();
            build_corners();
            build_lines();
            cut_on_internal_lines();
            return !inspect_required_;
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
                    read_shared_point_properties();
                }
                else
                {
                    unique_id = nb_unique_vertices;
                    point = read_point();
                    read_unique_point_properties();
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
            for( const auto i : geode::LRange{ 3 } )
            {
                const auto value = geode::string_to_double( tokens[i + 2] );
                point.set_value( i, value );
            }
            point.set_value( 2, crs_.z_sign * point.value( 2 ) );
            return point;
        }

        void read_shared_point_properties()
        {
            read_point_properties( 3 );
        }

        void read_unique_point_properties()
        {
            read_point_properties( 5 );
        }

        void read_point_properties( geode::index_t line_properties_position )
        {
            const auto split_line = get_tokens();
            for( const auto attr_id :
                geode::Indices{ vertices_prop_header_.names } )
            {
                for( const auto item :
                    geode::LRange{ vertices_prop_header_.esizes[attr_id] } )
                {
                    geode_unused( item );
                    OPENGEODE_ASSERT(
                        line_properties_position < split_line.size(),
                        "[LSOInput::read_point_properties] Cannot read "
                        "properties: number of property items is higher than "
                        "number of tokens." );
                    vertices_attributes_[attr_id].push_back(
                        geode::string_to_double(
                            split_line[line_properties_position] ) );
                    line_properties_position++;
                }
            }
        }

        void read_vertex_region_indicators()
        {
            if( geode::detail::string_starts_with(
                    line_, "BEGIN_VERTEX_REGION_INDICATORS" ) )
            {
                geode::detail::goto_keyword(
                    file_, "END_VERTEX_REGION_INDICATORS" );
                std::getline( file_, line_ );
            }
        }

        void read_tetrahedra()
        {
            do
            {
                const auto tokens = get_tokens();
                std::array< geode::index_t, 4 > vertices;
                for( const auto i : geode::LRange{ 4 } )
                {
                    auto ok = absl::SimpleAtoi( tokens[i + 1], &vertices[i] );
                    OPENGEODE_EXCEPTION(
                        ok, "[LSOInput] Error while reading tetra" );
                    vertices[i] -= OFFSET_START;
                }
                const auto tetra_id =
                    solid_builder_->create_tetrahedron( vertices );
                read_tetrahedron_properties();
                std::getline( file_, line_ );
                const auto tokens2 = get_tokens();
                block_name_attribute_->set_value(
                    tetra_id, geode::to_string( tokens2[2] ) );
            } while( std::getline( file_, line_ )
                     && geode::detail::string_starts_with( line_, "TETRA" ) );
            solid_builder_->compute_polyhedron_adjacencies();
        }

        void read_tetrahedron_properties()
        {
            const auto split_line = get_tokens();
            geode::index_t line_properties_position{ 5 };
            for( const auto attr_id :
                geode::Indices{ tetrahedra_prop_header_.names } )
            {
                for( const auto item :
                    geode::LRange{ tetrahedra_prop_header_.esizes[attr_id] } )
                {
                    geode_unused( item );
                    OPENGEODE_ASSERT(
                        line_properties_position < split_line.size(),
                        "[LSOInput::read_tetra_properties] Cannot read "
                        "properties: number of property items is higher than "
                        "number of tokens." );
                    tetrahedra_attributes_[attr_id].push_back(
                        geode::string_to_double(
                            split_line[line_properties_position] ) );
                    line_properties_position++;
                }
            }
        }

        void read_tetrahedra_region_indicators()
        {
            if( geode::detail::string_starts_with(
                    line_, "BEGIN_TETRA_REGION_INDICATORS" ) )
            {
                geode::detail::goto_keyword(
                    file_, "END_TETRA_REGION_INDICATORS" );
                std::getline( file_, line_ );
            }
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
                const auto tokens = geode::string_split( line_ );
                absl::Span< const absl::string_view > remaining_tokens(
                    &tokens[1], tokens.size() - 1 );
                const auto h_id = builder_.add_horizon();
                builder_.set_horizon_name(
                    h_id, geode::detail::read_name( remaining_tokens ) );
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
                builder_.set_surface_name( id, horizon.name() );
                std::getline( file_, line_ );
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
                for( const auto i : geode::LRange{ 3 } )
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
                const auto solid_facets =
                    facets_from_vertices( facet_vertices );
                for( const auto& facet : solid_facets )
                {
                    facet_id_->set_value( facet, surface_id );
                }
                if( solid_facets.empty() )
                {
                    inspect_required_ = true;
                    geode::Logger::warn(
                        "[LSOInput] Surface triangle with vertices [",
                        facet_vertices[0], " ", facet_vertices[1], " ",
                        facet_vertices[2],
                        "] is not conformal to the solid tetrahedra." );
                }
            }
            builder->compute_polygon_adjacencies();
        }

        std::vector< geode::index_t > facets_from_vertices(
            const std::array< geode::index_t, 3 >& vertices ) const
        {
            std::vector< geode::index_t > facets;
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
                            facets.emplace_back( facet_id.value() );
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
                const auto block_id =
                    builder_.add_block( geode::MeshFactory::default_impl(
                        geode::TetrahedralSolid3D::type_name_static() ) );
                builder_.set_block_name( block_id, tokens[1] );
                build_block_mesh( block_id );
                build_block_relations( block_id );
                std::getline( file_, line_ );
            }
        }

        void build_block_mesh( const geode::uuid& block_id )
        {
            auto builder =
                builder_.block_mesh_builder< geode::TetrahedralSolid3D >(
                    block_id );
            const auto component_id = model_.block( block_id ).component_id();
            const auto block_name = model_.block( block_id ).name();
            absl::flat_hash_map< geode::index_t, geode::index_t >
                vertex_mapping;
            std::vector< geode::index_t > inverse_vertex_mapping;
            std::vector< geode::index_t > inverse_tetrahedra_mapping;
            for( const auto tetra : geode::Range{ solid_->nb_polyhedra() } )
            {
                if( block_name_attribute_->value( tetra ) != block_name )
                {
                    continue;
                }
                std::array< geode::index_t, 4 > vertices;
                for( const auto i : geode::LRange{ 4 } )
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
                        inverse_vertex_mapping.push_back( vertex );
                        vertices[i] = vertex_id;
                        builder_.set_unique_vertex( { component_id, vertex_id },
                            vertex_id_->value( vertex ) );
                    }
                }
                builder->create_tetrahedron( vertices );
                inverse_tetrahedra_mapping.push_back( tetra );
            }
            builder->compute_polyhedron_adjacencies();
            create_block_vertices_attributes(
                block_id, inverse_vertex_mapping );
            create_block_tetrahedra_attributes(
                block_id, inverse_tetrahedra_mapping );
        }

        void create_block_vertices_attributes( const geode::uuid& block_id,
            absl::Span< const geode::index_t > inverse_vertex_mapping )
        {
            const auto& block_mesh = model_.block( block_id ).mesh();
            for( const auto attr_id :
                geode::Indices{ vertices_prop_header_.names } )
            {
                const auto nb_attribute_items =
                    vertices_prop_header_.esizes[attr_id];
                if( nb_attribute_items == 1 )
                {
                    auto attribute =
                        block_mesh.vertex_attribute_manager()
                            .find_or_create_attribute< geode::VariableAttribute,
                                double >(
                                vertices_prop_header_.names[attr_id], 0 );
                    for( const auto pt_id :
                        geode::Range{ block_mesh.nb_vertices() } )
                    {
                        attribute->set_value( pt_id,
                            vertices_attributes_
                                [attr_id][inverse_vertex_mapping[pt_id]] );
                    }
                }
                else if( nb_attribute_items == 2 )
                {
                    std::array< double, 2 > container;
                    container.fill( 0 );
                    add_vertices_container_attribute( block_mesh,
                        inverse_vertex_mapping, attr_id, container );
                }
                else if( nb_attribute_items == 3 )
                {
                    std::array< double, 3 > container;
                    container.fill( 0 );
                    add_vertices_container_attribute( block_mesh,
                        inverse_vertex_mapping, attr_id, container );
                }
                else
                {
                    std::vector< double > container( nb_attribute_items, 0 );
                    add_vertices_container_attribute<>( block_mesh,
                        inverse_vertex_mapping, attr_id, container );
                }
            }
        }

        void create_block_tetrahedra_attributes( const geode::uuid& block_id,
            absl::Span< const geode::index_t > inverse_tetrahedra_mapping )
        {
            const auto& block_mesh = model_.block( block_id ).mesh();
            for( const auto attr_id :
                geode::Indices{ tetrahedra_prop_header_.names } )
            {
                const auto nb_attribute_items =
                    tetrahedra_prop_header_.esizes[attr_id];
                if( nb_attribute_items == 1 )
                {
                    auto attribute =
                        block_mesh.polyhedron_attribute_manager()
                            .find_or_create_attribute< geode::VariableAttribute,
                                double >(
                                tetrahedra_prop_header_.names[attr_id], 0 );
                    for( const auto tetra_id :
                        geode::Indices{ inverse_tetrahedra_mapping } )
                    {
                        attribute->set_value( tetra_id,
                            tetrahedra_attributes_
                                [attr_id]
                                [inverse_tetrahedra_mapping[tetra_id]] );
                    }
                }
                else if( nb_attribute_items == 2 )
                {
                    std::array< double, 2 > container;
                    container.fill( 0 );
                    add_tetrahedra_container_attribute( block_mesh,
                        inverse_tetrahedra_mapping, attr_id, container );
                }
                else if( nb_attribute_items == 3 )
                {
                    std::array< double, 3 > container;
                    container.fill( 0 );
                    add_tetrahedra_container_attribute( block_mesh,
                        inverse_tetrahedra_mapping, attr_id, container );
                }
                else
                {
                    std::vector< double > container( nb_attribute_items, 0 );
                    add_tetrahedra_container_attribute( block_mesh,
                        inverse_tetrahedra_mapping, attr_id, container );
                }
            }
        }

        template < typename Container >
        void add_vertices_container_attribute(
            const geode::SolidMesh3D& block_mesh,
            absl::Span< const geode::index_t > inverse_mapping,
            geode::index_t attr_id,
            Container value_array )
        {
            auto attribute =
                block_mesh.vertex_attribute_manager()
                    .template find_or_create_attribute<
                        geode::VariableAttribute, Container >(
                        vertices_prop_header_.names[attr_id], value_array );
            for( const auto pt_id : geode::Range{ block_mesh.nb_vertices() } )
            {
                for( const auto item_id : geode::LRange{ value_array.size() } )
                {
                    value_array[item_id] =
                        vertices_attributes_[attr_id][inverse_mapping[pt_id]
                                                          * value_array.size()
                                                      + item_id];
                }
                attribute->set_value( pt_id, value_array );
            }
        }

        template < typename Container >
        void add_tetrahedra_container_attribute(
            const geode::SolidMesh3D& block_mesh,
            absl::Span< const geode::index_t > inverse_mapping,
            geode::index_t attr_id,
            Container value_array )
        {
            auto attribute =
                block_mesh.polyhedron_attribute_manager()
                    .template find_or_create_attribute<
                        geode::VariableAttribute, Container >(
                        tetrahedra_prop_header_.names[attr_id], value_array );
            for( const auto tetra_id :
                geode::Range{ block_mesh.nb_polyhedra() } )
            {
                for( const auto item_id : geode::LRange{ value_array.size() } )
                {
                    value_array[item_id] =
                        tetrahedra_attributes_[attr_id]
                                              [inverse_mapping[tetra_id]
                                                      * value_array.size()
                                                  + item_id];
                }
                attribute->set_value( tetra_id, value_array );
            }
        }

        absl::flat_hash_map< geode::uuid, geode::index_t > find_block_relations(
            const geode::uuid& block_id )
        {
            absl::flat_hash_map< geode::uuid, geode::index_t >
                surface_relations;
            const auto block_name = model_.block( block_id ).name();
            for( const auto tetra : geode::Range{ solid_->nb_polyhedra() } )
            {
                if( block_name_attribute_->value( tetra ) != block_name )
                {
                    continue;
                }
                for( const auto f : geode::LRange{ 4 } )
                {
                    const auto facet = solid_->facets().facet_from_vertices(
                        solid_->polyhedron_facet_vertices( { tetra, f } ) );
                    const auto& facet_uuid = facet_id_->value( facet.value() );
                    if( facet_uuid == default_id_ )
                    {
                        continue;
                    }
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
            }
            return surface_relations;
        }

        void build_block_relations( const geode::uuid& block_id )
        {
            const auto& block = model_.block( block_id );
            for( const auto& relation : find_block_relations( block_id ) )
            {
                const auto& surface = model_.surface( relation.first );
                const auto nb_relations =
                    relation.second / surface.mesh().nb_polygons();
                if( nb_relations == 1 )
                {
                    builder_.add_surface_block_boundary_relationship(
                        surface, block );
                }
                else if( nb_relations == 2 )
                {
                    builder_.add_surface_block_internal_relationship(
                        surface, block );
                }
                else
                {
                    inspect_required_ = true;
                    geode::Logger::warn( "[LSOInput] Block ", block.name(),
                        " is not conformal to surface ", surface.name(), "." );
                }
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
            const auto surfaces0 = model_.component_mesh_vertices(
                unique_id0, geode::Surface3D::component_type_static() );

            const auto vertex1 = mesh.polygon_edge_vertex( border, 1 );
            const auto unique_id1 =
                model_.unique_vertex( { component_id, vertex1 } );
            const auto surfaces1 = model_.component_mesh_vertices(
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

        bool are_equal( const std::vector< geode::ComponentMeshVertex >& lhs,
            const std::vector< geode::ComponentMeshVertex >& rhs ) const
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
            if( !model_.has_component_mesh_vertices(
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
            return geode::string_split( line_ );
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
                        if( !model_.has_component_mesh_vertices( unique_id0,
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
                !model_.has_component_mesh_vertices( unique_id, corner_type ) );
            const auto corner0 =
                model_.component_mesh_vertices( unique_id0, corner_type )[0];
            const auto corner1 =
                model_.component_mesh_vertices( unique_id, corner_type )[0];
            builder_.add_corner_line_boundary_relationship(
                model_.corner( corner0.component_id.id() ), line );
            builder_.add_corner_line_boundary_relationship(
                model_.corner( corner1.component_id.id() ), line );
        }

        absl::optional< geode::uuid > common_line(
            geode::index_t unique_id0, geode::index_t unique_id1 )
        {
            const auto lines0 = model_.component_mesh_vertices(
                unique_id0, geode::Line3D::component_type_static() );
            const auto lines1 = model_.component_mesh_vertices(
                unique_id1, geode::Line3D::component_type_static() );
            for( const auto& cmv0 : lines0 )
            {
                for( const auto& cmv1 : lines1 )
                {
                    if( cmv0.component_id == cmv1.component_id )
                    {
                        const auto min = std::min( cmv0.vertex, cmv1.vertex );
                        const auto max = std::max( cmv0.vertex, cmv1.vertex );
                        if( max - min == 1 )
                        {
                            return cmv0.component_id.id();
                        }
                    }
                }
            }
            return absl::nullopt;
        }

        void cut_on_internal_lines()
        {
            for( const auto& line : model_.lines() )
            {
                const auto component_edges =
                    geode::component_mesh_edges( model_, line, 0 );
                for( const auto& surface_edges : component_edges.surface_edges )
                {
                    const auto& surface = model_.surface( surface_edges.first );
                    if( model_.is_boundary( line, surface )
                        || model_.is_internal( line, surface ) )
                    {
                        continue;
                    }
                    geode::Logger::warn( "Surface ", surface.name(),
                        " was not cut by one of its internal lines, adding "
                        "the relation and cutting the surface to ensure "
                        "model validity." );
                    builder_.add_line_surface_internal_relationship(
                        line, surface );
                    inspect_required_ = true;
                }
            }

            geode::detail::CutAlongInternalLines< geode::BRep,
                geode::BRepBuilder, 3 >
                cutter{ model_ };
            cutter.cut_all_surfaces();
        }

    private:
        bool inspect_required_{ false };
        std::ifstream file_;
        std::string line_;
        geode::StructuralModel& model_;
        geode::StructuralModelBuilder builder_;
        geode::detail::CRSData crs_;
        geode::detail::PropHeaderData vertices_prop_header_;
        geode::detail::PropHeaderData tetrahedra_prop_header_;
        std::vector< std::vector< double > > vertices_attributes_;
        std::vector< std::vector< double > > tetrahedra_attributes_;
        std::unique_ptr< geode::TetrahedralSolid3D > solid_;
        std::unique_ptr< geode::TetrahedralSolidBuilder3D > solid_builder_;
        std::shared_ptr< geode::VariableAttribute< geode::index_t > >
            vertex_id_;
        std::shared_ptr< geode::VariableAttribute< std::string > >
            block_name_attribute_;
        std::shared_ptr< geode::VariableAttribute< geode::uuid > > facet_id_;
        geode::uuid default_id_;
        std::vector< absl::InlinedVector< geode::index_t, 1 > > vertex_mapping_;
    };
} // namespace

namespace geode
{
    namespace detail
    {
        StructuralModel LSOInput::read()
        {
            StructuralModel structural_model;
            LSOInputImpl impl{ filename(), structural_model };
            const auto file_reading_ok = impl.read_file();
            if( !file_reading_ok )
            {
                this->need_to_inspect_result();
            }
            return structural_model;
        }
    } // namespace detail
} // namespace geode
