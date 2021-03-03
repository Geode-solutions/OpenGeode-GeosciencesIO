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

#include <geode/geosciences/private/ml_input.h>

#include <algorithm>
#include <fstream>
#include <functional>

#include <absl/container/flat_hash_set.h>

#include <geode/geometry/bounding_box.h>
#include <geode/geometry/nn_search.h>
#include <geode/geometry/point.h>
#include <geode/geometry/vector.h>

#include <geode/mesh/builder/edged_curve_builder.h>
#include <geode/mesh/builder/point_set_builder.h>
#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/geode_triangulated_surface.h>

#include <geode/model/mixin/core/corner.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/mixin/core/vertex_identifier.h>

#include <geode/geosciences/private/gocad_common.h>
#include <geode/geosciences/representation/builder/structural_model_builder.h>
#include <geode/geosciences/representation/core/structural_model.h>

namespace
{
    std::tuple< geode::PolygonEdge, bool > find_edge(
        const geode::SurfaceMesh3D& mesh, geode::index_t v0, geode::index_t v1 )
    {
        if( const auto edge = mesh.polygon_edge_from_vertices( v0, v1 ) )
        {
            return std::make_tuple( edge.value(), true );
        }
        if( const auto edge = mesh.polygon_edge_from_vertices( v1, v0 ) )
        {
            return std::make_tuple( edge.value(), false );
        }
        throw geode::OpenGeodeException{
            "[MLInput] Starting edge Line from Corner not found. ",
            "Looking for edge: ", mesh.point( v0 ).string(), " - ",
            mesh.point( v1 ).string()
        };
    }

    class MLInputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };

        MLInputImpl( absl::string_view filename, geode::StructuralModel& model )
            : file_( filename.data() ), model_( model ), builder_( model )
        {
            OPENGEODE_EXCEPTION( file_.good(),
                "[MLInput] Error while opening file: ", filename );
        }

        void read_file()
        {
            geode::detail::check_keyword( file_, "GOCAD Model3d" );
            geode::detail::read_header( file_ );
            geode::detail::read_CRS( file_ );
            read_model_components();
            for( auto& tsurf : tsurfs_ )
            {
                tsurf.data = geode::detail::read_tsurf( file_ ).value();
                build_surfaces( tsurf );
            }
            compute_epsilon();
            build_corners();
            build_lines();
            complete_vertex_identifier();
        }

    private:
        struct TSurfMLData
        {
            TSurfMLData( std::string input_name )
                : name( std::move( input_name ) )
            {
            }

            geode::detail::TSurfData data;
            std::vector< std::reference_wrapper< const geode::uuid > > tfaces;
            std::string feature;
            std::string name = std::string( "unknown" );
        };

        using LinesID =
            absl::InlinedVector< std::reference_wrapper< const geode::uuid >,
                1 >;

        struct LineData
        {
            bool operator!=( const LineData& rhs ) const
            {
                return line != rhs.line;
            }
            bool operator<( const LineData& rhs ) const
            {
                return line < rhs.line;
            }
            geode::uuid corner0, corner1, line, surface;
            std::vector< geode::Point3D > points;
            std::vector< geode::index_t > indices;
        };

        void compute_epsilon()
        {
            const auto box = model_.bounding_box();
            const geode::Vector3D diagonal{ box.min(), box.max() };
            epsilon_ =
                1e-7
                * std::min( diagonal.value( 0 ),
                    std::min( diagonal.value( 1 ), diagonal.value( 2 ) ) );
        }

        void complete_vertex_identifier()
        {
            for( const auto& surface : model_.surfaces() )
            {
                const auto& mesh = surface.mesh();
                const auto surface_id = surface.component_id();
                for( const auto v : geode::Range{ mesh.nb_vertices() } )
                {
                    if( model_.unique_vertex( { surface.component_id(), v } )
                        == geode::NO_ID )
                    {
                        builder_.set_unique_vertex( { surface_id, v },
                            builder_.create_unique_vertex() );
                    }
                }
            }
        }

        void build_corners()
        {
            std::vector< geode::Point3D > corner_points;
            std::vector< geode::MeshComponentVertex > corner_surface_index;
            for( const auto& tsurf : tsurfs_ )
            {
                const auto& data = tsurf.data;
                for( const auto corner : data.bstones )
                {
                    corner_points.push_back( data.points[corner] );
                    const auto tface_id = data.tface_id( corner );
                    const auto& surface =
                        model_.surface( tsurf.tfaces[tface_id] );
                    const auto vertex_id =
                        corner - data.tface_vertices_offset[tface_id];
                    corner_surface_index.emplace_back(
                        surface.component_id(), vertex_id );
                }
            }

            geode::NNSearch3D ann{ std::move( corner_points ) };
            auto colocated_info = ann.colocated_index_mapping( epsilon_ );
            for( auto& point : colocated_info.unique_points )
            {
                const auto& corner_id = builder_.add_corner();
                builder_.corner_mesh_builder( corner_id )
                    ->create_point( point );
                const auto vertex_id = builder_.create_unique_vertex();
                builder_.set_unique_vertex(
                    { model_.corner( corner_id ).component_id(), 0 },
                    vertex_id );
            }
            for( const auto i :
                geode::Indices{ colocated_info.colocated_mapping } )
            {
                if( model_.unique_vertex( corner_surface_index[i] )
                    != geode::NO_ID )
                {
                    geode::Logger::warn(
                        "[MLInput::build_corners] Overriding Corner/Surface "
                        "topological information. Please verify "
                        "StructuralModel validity." );
                }
                builder_.set_unique_vertex(
                    std::move( corner_surface_index[i] ),
                    colocated_info.colocated_mapping[i] );
            }
        }

        void build_lines()
        {
            std::vector< std::pair< geode::MeshComponentVertex,
                geode::MeshComponentVertex > >
                line_surface_index;
            for( const auto& tsurf : tsurfs_ )
            {
                const auto& data = tsurf.data;
                for( const auto border : data.borders )
                {
                    const auto tface_id = data.tface_id( border.corner_id );
                    const auto& surface =
                        model_.surface( tsurf.tfaces[tface_id] );
                    const auto corner_id =
                        border.corner_id - data.tface_vertices_offset[tface_id];
                    const auto next_id =
                        border.next_id - data.tface_vertices_offset[tface_id];
                    line_surface_index.emplace_back(
                        geode::MeshComponentVertex{
                            surface.component_id(), corner_id },
                        geode::MeshComponentVertex{
                            surface.component_id(), next_id } );
                }
            }

            absl::flat_hash_map< geode::uuid, std::vector< LineData > >
                surface2lines;
            for( const auto& line_start : line_surface_index )
            {
                const auto& surface_id = line_start.first.component_id.id();
                const auto& surface = model_.surface( surface_id );
                auto line_data = compute_line( surface, line_start );
                line_data.line = find_or_create_line( line_data );
                surface2lines[surface_id].emplace_back(
                    std::move( line_data ) );
            }

            for( auto& surface2line : surface2lines )
            {
                auto& lines = surface2line.second;
                std::sort( lines.begin(), lines.end() );
                for( size_t l = 1; l < lines.size(); l++ )
                {
                    if( lines[l - 1] != lines[l] )
                    {
                        register_line_data( lines[l - 1] );
                    }
                    else
                    {
                        register_internal_line_data( lines[l - 1] );
                        const auto& surface =
                            model_.surface( lines[l].surface );
                        register_line_surface_vertex_identifier(
                            lines[l], surface.component_id() );
                        l++;
                    }
                }
                const auto size = lines.size() - 1;
                if( lines[size - 1] != lines[size] )
                {
                    register_line_data( lines[size] );
                }
            }
        }

        const geode::uuid& find_or_create_line( LineData& line_data )
        {
            const auto it = corners2line_.find(
                std::make_pair( line_data.corner0, line_data.corner1 ) );
            if( it != corners2line_.end() )
            {
                for( const auto& line_id : it->second )
                {
                    if( are_lines_equal( line_data, line_id ) )
                    {
                        return line_id;
                    }
                }
            }
            const auto it_reverse = corners2line_.find(
                std::make_pair( line_data.corner1, line_data.corner0 ) );
            if( it_reverse != corners2line_.end() )
            {
                for( const auto& line_id : it_reverse->second )
                {
                    if( are_lines_reverse_equal( line_data, line_id ) )
                    {
                        std::reverse( line_data.indices.begin(),
                            line_data.indices.end() );
                        return line_id;
                    }
                }
            }
            return create_line( line_data );
        }

        bool are_lines_equal(
            const LineData& line_data, const geode::uuid& line_id )
        {
            const auto& mesh = model_.line( line_id ).mesh();
            if( mesh.nb_vertices() != line_data.points.size() )
            {
                return false;
            }
            for( const auto v : geode::Range{ mesh.nb_vertices() } )
            {
                if( !mesh.point( v ).inexact_equal(
                        line_data.points[v], epsilon_ ) )
                {
                    return false;
                }
            }
            return true;
        }

        bool are_lines_reverse_equal(
            const LineData& line_data, const geode::uuid& line_id )
        {
            const auto& mesh = model_.line( line_id ).mesh();
            if( mesh.nb_vertices() != line_data.points.size() )
            {
                return false;
            }
            for( const auto v : geode::Range{ mesh.nb_vertices() } )
            {
                if( !mesh.point( v ).inexact_equal(
                        line_data.points[mesh.nb_vertices() - v - 1],
                        epsilon_ ) )
                {
                    return false;
                }
            }
            return true;
        }

        void build_surfaces( const TSurfMLData& tsurf )
        {
            for( const auto i : geode::Range{ tsurf.tfaces.size() } )
            {
                const auto& uuid = tsurf.tfaces[i];
                auto builder =
                    builder_
                        .surface_mesh_builder< geode::TriangulatedSurface3D >(
                            uuid );
                const auto& data = tsurf.data;
                builder->set_name( data.header.name );
                for( const auto p : geode::Range{ data.tface_vertices_offset[i],
                         data.tface_vertices_offset[i + 1] } )
                {
                    builder->create_point( data.points[p] );
                }
                for( const auto t :
                    geode::Range{ data.tface_triangles_offset[i],
                        data.tface_triangles_offset[i + 1] } )
                {
                    builder->create_triangle(
                        { data.triangles[t][0] - data.tface_vertices_offset[i],
                            data.triangles[t][1]
                                - data.tface_vertices_offset[i],
                            data.triangles[t][2]
                                - data.tface_vertices_offset[i] } );
                }

                builder->compute_polygon_adjacencies();
            }
        }

        const geode::uuid& create_line( const LineData& line_data )
        {
            const auto& line_id = builder_.add_line();
            auto it = corners2line_.try_emplace(
                std::make_pair( line_data.corner0, line_data.corner1 ),
                LinesID{ line_id } );
            if( !it.second )
            {
                it.first->second.push_back( line_id );
            }
            create_line_geometry( line_data, line_id );
            create_line_topology( line_data, line_id );
            return line_id;
        }

        void create_line_topology(
            const LineData& line_data, const geode::uuid& line_id )
        {
            const auto& line = model_.line( line_id );
            builder_.add_corner_line_boundary_relationship(
                model_.corner( line_data.corner0 ), line );
            builder_.add_corner_line_boundary_relationship(
                model_.corner( line_data.corner1 ), line );

            const auto line_component_id = line.component_id();
            const auto last_index =
                static_cast< geode::index_t >( line_data.points.size() - 1 );
            builder_.set_unique_vertex( { line_component_id, 0 },
                model_.unique_vertex(
                    { model_.corner( line_data.corner0 ).component_id(),
                        0 } ) );
            builder_.set_unique_vertex( { line_component_id, last_index },
                model_.unique_vertex(
                    { model_.corner( line_data.corner1 ).component_id(),
                        0 } ) );
            for( const auto i : geode::Range{ 1, last_index } )
            {
                builder_.set_unique_vertex(
                    { line_component_id, i }, builder_.create_unique_vertex() );
            }
        }

        void create_line_geometry(
            const LineData& line_data, const geode::uuid& line_id )
        {
            const auto line_builder = builder_.line_mesh_builder( line_id );
            for( const auto& point : line_data.points )
            {
                line_builder->create_point( point );
            }
            for( const auto i : geode::Range( 1, line_data.points.size() ) )
            {
                line_builder->create_edge( i - 1, i );
            }
        }

        void register_line_data( const LineData& line_data )
        {
            const auto& surface = model_.surface( line_data.surface );
            builder_.add_line_surface_boundary_relationship(
                model_.line( line_data.line ), surface );
            register_line_surface_vertex_identifier(
                line_data, surface.component_id() );
        }

        void register_internal_line_data( const LineData& line_data )
        {
            const auto& surface = model_.surface( line_data.surface );
            builder_.add_line_surface_internal_relationship(
                model_.line( line_data.line ), surface );
            register_line_surface_vertex_identifier(
                line_data, surface.component_id() );
        }

        void register_line_surface_vertex_identifier(
            const LineData& line_data, const geode::ComponentID& surface_id )
        {
            for( const auto i : geode::Range{ line_data.indices.size() } )
            {
                const auto unique_id = model_.unique_vertex(
                    { model_.line( line_data.line ).component_id(), i } );
                builder_.set_unique_vertex(
                    { surface_id, line_data.indices[i] }, unique_id );
            }
        }

        LineData compute_line( const geode::Surface3D& surface,
            const std::pair< geode::MeshComponentVertex,
                geode::MeshComponentVertex >& line_start )
        {
            LineData result;
            result.surface = surface.id();
            result.corner0 =
                model_
                    .mesh_component_vertices(
                        model_.unique_vertex( { surface.component_id(),
                            line_start.first.vertex } ),
                        geode::Corner3D::component_type_static() )
                    .front()
                    .component_id.id();
            const auto& mesh = surface.mesh();
            result.indices.push_back( line_start.first.vertex );
            result.points.emplace_back( mesh.point( result.indices.back() ) );
            bool orientation;
            geode::PolygonEdge edge;
            std::tie( edge, orientation ) = find_edge(
                mesh, line_start.second.vertex, line_start.first.vertex );
            if( orientation )
            {
                while( model_.unique_vertex( { surface.component_id(),
                           mesh.polygon_vertex( edge ) } )
                       == geode::NO_ID )
                {
                    result.indices.push_back( mesh.polygon_vertex( edge ) );
                    result.points.emplace_back(
                        mesh.point( result.indices.back() ) );
                    edge = mesh.previous_on_border( edge );
                }
                result.indices.push_back( mesh.polygon_vertex( edge ) );
                result.points.emplace_back(
                    mesh.point( result.indices.back() ) );
            }
            else
            {
                while( model_.unique_vertex( { surface.component_id(),
                           mesh.polygon_edge_vertex( edge, 1 ) } )
                       == geode::NO_ID )
                {
                    result.indices.push_back(
                        mesh.polygon_edge_vertex( edge, 1 ) );
                    result.points.emplace_back(
                        mesh.point( result.indices.back() ) );
                    edge = mesh.next_on_border( edge );
                }
                result.indices.push_back( mesh.polygon_edge_vertex( edge, 1 ) );
                result.points.emplace_back(
                    mesh.point( result.indices.back() ) );
            }

            result.corner1 =
                model_
                    .mesh_component_vertices(
                        model_.unique_vertex(
                            { surface.component_id(), result.indices.back() } ),
                        geode::Corner3D::component_type_static() )
                    .front()
                    .component_id.id();
            return result;
        }

        void read_model_components()
        {
            std::string line;
            while( std::getline( file_, line ) )
            {
                if( geode::detail::string_starts_with( line, "END" ) )
                {
                    create_tsurfs();
                    return;
                }

                std::istringstream iss{ line };
                std::string component_type;
                iss >> component_type;
                if( component_type == "TSURF" )
                {
                    process_TSURF_keyword( iss );
                }
                else if( component_type == "TFACE" )
                {
                    process_TFACE_keyword( iss );
                }
                else if( component_type == "REGION" )
                {
                    process_REGION_keyword( iss );
                }
                else if( component_type == "LAYER" )
                {
                    process_LAYER_keyword( iss );
                }
                else if( component_type == "FAULT_BLOCK" )
                {
                    process_FAULT_BLOCK_keyword( iss );
                }
            }
            throw geode::OpenGeodeException{
                "[MLInput] Cannot find the end of component section"
            };
        }

        void create_tsurfs()
        {
            std::vector< geode::uuid > boundaries;
            for( const auto& tsurf : tsurfs_ )
            {
                if( tsurf.feature.find( "fault" ) != std::string::npos )
                {
                    const auto& fault_uuid =
                        builder_.add_fault( fault_map_.at( tsurf.feature ) );
                    builder_.set_fault_name( fault_uuid, tsurf.name );
                    const auto& fault = model_.fault( fault_uuid );
                    for( const auto& uuid : tsurf.tfaces )
                    {
                        builder_.add_surface_in_fault(
                            model_.surface( uuid ), fault );
                    }
                }
                else if( tsurf.feature == "boundary"
                         || tsurf.feature == "lease" )
                {
                    const auto& model_boundary_uuid =
                        builder_.add_model_boundary();
                    builder_.set_model_boundary_name(
                        model_boundary_uuid, tsurf.name );
                    const auto& model_boundary =
                        model_.model_boundary( model_boundary_uuid );
                    for( const auto& uuid : tsurf.tfaces )
                    {
                        builder_.add_surface_in_model_boundary(
                            model_.surface( uuid ), model_boundary );
                        boundaries.push_back( uuid );
                    }
                }
                else
                {
                    const auto& horizon_uuid = builder_.add_horizon(
                        horizon_map_.at( tsurf.feature ) );
                    builder_.set_horizon_name( horizon_uuid, tsurf.name );
                    const auto& horizon = model_.horizon( horizon_uuid );
                    for( const auto& uuid : tsurf.tfaces )
                    {
                        builder_.add_surface_in_horizon(
                            model_.surface( uuid ), horizon );
                    }
                }
            }
            process_unassigned_model_boundaries( boundaries );
        }

        void process_unassigned_model_boundaries(
            std::vector< geode::uuid >& boundaries )
        {
            std::vector< geode::uuid > diff;
            diff.reserve( boundaries.size() );
            for( auto&& id : geode::EraserRange< geode::uuid >{ boundaries } )
            {
                if( !universe_.contains( id ) )
                {
                    diff.emplace_back( std::move( id ) );
                }
            }
            if( diff.empty() )
            {
                return;
            }
            const auto& model_boundary_uuid = builder_.add_model_boundary();
            builder_.set_model_boundary_name(
                model_boundary_uuid, "undefined boundary" );
            const auto& model_boundary =
                model_.model_boundary( model_boundary_uuid );
            for( const auto& uuid : diff )
            {
                builder_.add_surface_in_model_boundary(
                    model_.surface( uuid ), model_boundary );
            }
        }

        void process_TSURF_keyword( std::istringstream& iss )
        {
            auto name = geode::detail::read_name( iss );
            tsurf_names2index_.emplace( name, tsurfs_.size() );
            tsurfs_.emplace_back( std::move( name ) );
        }

        void process_TFACE_keyword( std::istringstream& iss )
        {
            geode::index_t id;
            std::string feature;
            iss >> id >> feature;
            std::string name = geode::detail::read_name( iss );
            const auto& surface_id = builder_.add_surface(
                geode::OpenGeodeTriangulatedSurface3D::impl_name_static() );
            builder_.set_surface_name( surface_id, name );
            auto& tsurf = tsurfs_[tsurf_names2index_.at( name )];
            tsurf.feature = feature;
            tsurf.tfaces.emplace_back( surface_id );
            surfaces_.emplace_back( surface_id );
        }

        void process_REGION_keyword( std::istringstream& iss )
        {
            geode::index_t id;
            iss >> id;
            auto name = geode::detail::read_name( iss );
            if( name == "Universe" )
            {
                read_universe();
                return;
            }
            const auto& block_id = builder_.add_block();
            builder_.set_block_name( block_id, std::move( name ) );
            create_block_topology( block_id );
            blocks_.emplace_back( block_id );
        }

        void process_LAYER_keyword( std::istringstream& iss )
        {
            auto name = geode::detail::read_name( iss );
            const auto& stratigraphic_unit_id =
                builder_.add_stratigraphic_unit();
            builder_.set_stratigraphic_unit_name(
                stratigraphic_unit_id, std::move( name ) );
            create_stratigraphic_unit_topology( stratigraphic_unit_id );
        }

        void process_FAULT_BLOCK_keyword( std::istringstream& iss )
        {
            auto name = geode::detail::read_name( iss );
            const auto& fault_block_id = builder_.add_fault_block();
            builder_.set_fault_block_name( fault_block_id, std::move( name ) );
            create_fault_block_topology( fault_block_id );
        }

        void read_universe()
        {
            while( true )
            {
                std::string boundary_line;
                std::getline( file_, boundary_line );
                std::istringstream boundary_iss{ boundary_line };
                int surface_id;
                while( boundary_iss >> surface_id )
                {
                    if( surface_id == 0 )
                    {
                        return;
                    }
                    surface_id = std::abs( surface_id ) - OFFSET_START;
                    universe_.emplace( surfaces_[surface_id] );
                }
            }
        }

        void create_block_topology( const geode::uuid& block_id )
        {
            const auto& block = model_.block( block_id );
            std::vector< geode::index_t > surfaces;
            while( true )
            {
                std::string boundary_line;
                std::getline( file_, boundary_line );
                std::istringstream boundary_iss{ boundary_line };
                int surface_id;
                while( boundary_iss >> surface_id )
                {
                    if( surface_id == 0 )
                    {
                        std::sort( surfaces.begin(), surfaces.end() );
                        for( size_t s = 1; s < surfaces.size(); s++ )
                        {
                            if( surfaces[s - 1] != surfaces[s] )
                            {
                                builder_
                                    .add_surface_block_boundary_relationship(
                                        model_.surface(
                                            surfaces_[surfaces[s - 1]] ),
                                        block );
                            }
                            else
                            {
                                builder_
                                    .add_surface_block_internal_relationship(
                                        model_.surface(
                                            surfaces_[surfaces[s]] ),
                                        block );
                                s++;
                            }
                        }
                        const auto size = surfaces.size() - 1;
                        if( surfaces[size - 1] != surfaces[size] )
                        {
                            builder_.add_surface_block_boundary_relationship(
                                model_.surface( surfaces_[surfaces[size]] ),
                                block );
                        }
                        return;
                    }
                    surface_id = std::abs( surface_id ) - OFFSET_START;
                    surfaces.push_back( surface_id );
                }
            }
        }

        void create_stratigraphic_unit_topology(
            const geode::uuid& stratigraphic_unit_id )
        {
            const auto& stratigraphic_unit =
                model_.stratigraphic_unit( stratigraphic_unit_id );
            std::vector< geode::index_t > blocks;
            const auto blocks_offset = OFFSET_START + surfaces_.size();
            while( true )
            {
                std::string boundary_line;
                std::getline( file_, boundary_line );
                std::istringstream boundary_iss{ boundary_line };
                geode::index_t block_id;
                while( boundary_iss >> block_id )
                {
                    if( block_id == 0 )
                    {
                        for( const auto b : blocks )
                        {
                            if( b == 0 )
                            {
                                continue; // Universe
                            }
                            if( b - OFFSET_START < blocks_.size() )
                            {
                                builder_.add_block_in_stratigraphic_unit(
                                    model_.block( blocks_[b - OFFSET_START] ),
                                    stratigraphic_unit );
                            }
                            else
                            {
                                geode::Logger::warn(
                                    "[MLInput] Stated in LAYER ",
                                    stratigraphic_unit.name(), ", Block id ",
                                    b + blocks_offset,
                                    " does not match an existing REGION" );
                            }
                        }
                        return;
                    }
                    block_id = block_id - blocks_offset;
                    blocks.push_back( block_id );
                }
            }
        }

        void create_fault_block_topology( const geode::uuid& fault_block_id )
        {
            const auto& fault_block = model_.fault_block( fault_block_id );
            std::vector< geode::index_t > blocks;
            const auto blocks_offset = OFFSET_START + surfaces_.size();
            while( true )
            {
                std::string boundary_line;
                std::getline( file_, boundary_line );
                std::istringstream boundary_iss{ boundary_line };
                geode::index_t block_id;
                while( boundary_iss >> block_id )
                {
                    if( block_id == 0 )
                    {
                        for( const auto b : blocks )
                        {
                            if( b == 0 )
                            {
                                continue; // Universe
                            }
                            if( b - OFFSET_START < blocks_.size() )
                            {
                                builder_.add_block_in_fault_block(
                                    model_.block( blocks_[b - OFFSET_START] ),
                                    fault_block );
                            }
                            else
                            {
                                geode::Logger::warn(
                                    "[MLInput] Stated in FAULT_BLOCK ",
                                    fault_block.name(), ", Block id ",
                                    b + blocks_offset,
                                    " does not match an existing REGION" );
                            }
                        }
                        return;
                    }
                    block_id = block_id - blocks_offset;
                    blocks.push_back( block_id );
                }
            }
        }

    private:
        std::ifstream file_;
        geode::StructuralModel& model_;
        geode::StructuralModelBuilder builder_;
        absl::flat_hash_map< std::string, geode::index_t > tsurf_names2index_;
        std::vector< TSurfMLData > tsurfs_;
        absl::flat_hash_map< std::pair< geode::uuid, geode::uuid >, LinesID >
            corners2line_;
        std::vector< std::reference_wrapper< const geode::uuid > > surfaces_;
        std::vector< std::reference_wrapper< const geode::uuid > > blocks_;
        absl::flat_hash_set< geode::uuid > universe_;
        double epsilon_;
        const absl::flat_hash_map< std::string, geode::Fault3D::FAULT_TYPE >
            fault_map_ = { { "fault", geode::Fault3D::FAULT_TYPE::NO_TYPE },
                { "reverse_fault", geode::Fault3D::FAULT_TYPE::REVERSE },
                { "normal_fault", geode::Fault3D::FAULT_TYPE::NORMAL } };
        const absl::flat_hash_map< std::string, geode::Horizon3D::HORIZON_TYPE >
            horizon_map_ = { { "top",
                                 geode::Horizon3D::HORIZON_TYPE::CONFORMAL },
                { "none", geode::Horizon3D::HORIZON_TYPE::NO_TYPE },
                { "topographic", geode::Horizon3D::HORIZON_TYPE::TOPOGRAPHY },
                { "intrusive", geode::Horizon3D::HORIZON_TYPE::INTRUSION },
                { "unconformity",
                    geode::Horizon3D::HORIZON_TYPE::NON_CONFORMAL } };
    };
} // namespace

namespace geode
{
    namespace detail
    {
        void MLInput::read()
        {
            MLInputImpl impl{ filename(), structural_model() };
            impl.read_file();
        }
    } // namespace detail
} // namespace geode
