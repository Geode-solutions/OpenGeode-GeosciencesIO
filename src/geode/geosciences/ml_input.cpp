/*
 * Copyright (c) 2019 Geode-solutions
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

#include <geode/geosciences/detail/ml_input.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <unordered_map>

#include <geode/basic/bounding_box.h>
#include <geode/basic/nn_search.h>
#include <geode/basic/point.h>
#include <geode/basic/vector.h>

#include <geode/model/mixin/core/corner.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/mixin/core/vertex_identifier.h>

#include <geode/mesh/builder/edged_curve_builder.h>
#include <geode/mesh/builder/point_set_builder.h>
#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/geode_triangulated_surface.h>

#include <geode/geosciences/representation/builder/structural_model_builder.h>
#include <geode/geosciences/representation/core/structural_model.h>

namespace std
{
    template <>
    struct hash< std::pair< geode::uuid, geode::uuid > >
    {
    public:
        size_t operator()( const pair< geode::uuid, geode::uuid >& x ) const
        {
            return std::hash< geode::uuid >{}( x.first )
                   ^ std::hash< geode::uuid >{}( x.second );
        }
    };
} // namespace std

namespace
{
    bool string_starts_with(
        const std::string& string, const std::string& check )
    {
        return !string.compare( 0, check.length(), check );
    }

    geode::PolygonEdge find_edge( const geode::PolygonalSurface3D& mesh,
        geode::index_t v0,
        geode::index_t v1 )
    {
        bool found;
        geode::PolygonEdge edge;
        std::tie( found, edge ) = mesh.polygon_edge_from_vertices( v0, v1 );
        if( !found )
        {
            throw geode::OpenGeodeException{
                "Starting edge Line from Corner not found.",
                "Looking for edge: ", mesh.point( v0 ), " - ", mesh.point( v1 )
            };
        }
        return edge;
    }

    class MLInputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };

        MLInputImpl(
            const std::string& filename, geode::StructuralModel& model )
            : file_( filename ), model_( model ), builder_( model )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: " + filename );
        }

        void read_file()
        {
            check_keyword( "GOCAD Model3d" );
            read_header();
            read_CRS();
            read_model_components();
            for( auto& tsurf : tsurfs_ )
            {
                read_tsurf( tsurf );
            }
            compute_epsilon();
            build_corners();
            build_lines();
            complete_vertex_identifier();
        }

    private:
        struct TSurfData
        {
            geode::index_t tface_id( geode::index_t vertex_id ) const
            {
                for( auto i : geode::Range{ 1, tface_vertices_offset.size() } )
                {
                    if( vertex_id < tface_vertices_offset[i] )
                    {
                        return i - 1;
                    }
                }
                return tface_vertices_offset.size() - 1;
            }

            std::deque< std::reference_wrapper< const geode::uuid > > tfaces{};
            std::deque< geode::index_t > tface_vertices_offset{ 0 };
            std::string feature;
        };

        struct TopoInfo
        {
            std::vector< geode::Point3D > corner_points;
            std::vector< geode::MeshComponentVertex > corner_surface_index;
            using LineStart = std::pair< geode::MeshComponentVertex,
                geode::MeshComponentVertex >;
            std::vector< LineStart > line_surface_index;
        };

        struct LineData
        {
            geode::uuid corner0, corner1, surface;
            std::vector< geode::Point3D > points;
            std::vector< geode::index_t > indices;
        };

        void compute_epsilon()
        {
            geode::BoundingBox3D box;
            for( const auto& surface : model_.surfaces() )
            {
                const auto& mesh = surface.mesh();
                for( auto v : geode::Range{ mesh.nb_vertices() } )
                {
                    box.add_point( mesh.point( v ) );
                }
            }
            geode::Vector3D diagonal{ box.min(), box.max() };
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
                auto surface_id = surface.component_id();
                for( auto v : geode::Range{ mesh.nb_vertices() } )
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
            geode::NNSearch3D ann{ std::move( info_.corner_points ) };
            auto colocated_info = ann.colocated_index_mapping( epsilon_ );
            for( auto& point : colocated_info.unique_points )
            {
                const auto& corner_id = builder_.add_corner();
                builder_.corner_mesh_builder( corner_id )
                    ->create_point( point );
                auto vertex_id = builder_.create_unique_vertex();
                builder_.set_unique_vertex(
                    { model_.corner( corner_id ).component_id(), 0 },
                    vertex_id );
            }
            for( auto i :
                geode::Range{ colocated_info.colocated_mapping.size() } )
            {
                builder_.set_unique_vertex(
                    std::move( info_.corner_surface_index[i] ),
                    colocated_info.colocated_mapping[i] );
            }
        }

        void build_lines()
        {
            for( auto& line_start : info_.line_surface_index )
            {
                const auto& surface =
                    model_.surface( line_start.first.component_id.id() );
                auto line_data = compute_line( surface, line_start );
                register_line_data(
                    line_data, find_or_create_line( line_data ) );
            }
        }

        const geode::uuid& find_or_create_line( LineData& line_data )
        {
            auto it = corners2line_.find(
                std::make_pair( line_data.corner0, line_data.corner1 ) );
            if( it != corners2line_.end()
                && are_lines_equal( line_data, it->second ) )
            {
                return it->second;
            }
            auto it_reverse = corners2line_.find(
                std::make_pair( line_data.corner1, line_data.corner0 ) );
            if( it_reverse != corners2line_.end()
                && are_lines_reverse_equal( line_data, it_reverse->second ) )
            {
                std::reverse(
                    line_data.indices.begin(), line_data.indices.end() );
                return it_reverse->second;
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
            for( auto v : geode::Range{ mesh.nb_vertices() } )
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
            for( auto v : geode::Range{ mesh.nb_vertices() } )
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

        const geode::uuid& create_line( const LineData& line_data )
        {
            const auto& line_id = builder_.add_line();
            corners2line_.emplace(
                std::make_pair( line_data.corner0, line_data.corner1 ),
                line_id );
            create_line_geometry( line_data, line_id );
            create_line_topology( line_data, line_id );
            return line_id;
        }

        void create_line_topology(
            const LineData& line_data, const geode::uuid& line_id )
        {
            const auto& line = model_.line( line_id );
            builder_.add_corner_line_relationship(
                model_.corner( line_data.corner0 ), line );
            builder_.add_corner_line_relationship(
                model_.corner( line_data.corner1 ), line );

            auto line_component_id = line.component_id();
            auto last_index =
                static_cast< geode::index_t >( line_data.points.size() - 1 );
            builder_.set_unique_vertex( { line_component_id, 0 },
                model_.unique_vertex(
                    { model_.corner( line_data.corner0 ).component_id(),
                        0 } ) );
            builder_.set_unique_vertex( { line_component_id, last_index },
                model_.unique_vertex(
                    { model_.corner( line_data.corner1 ).component_id(),
                        0 } ) );
            for( auto i : geode::Range{ 1, last_index } )
            {
                builder_.set_unique_vertex(
                    { line_component_id, i }, builder_.create_unique_vertex() );
            }
        }

        void create_line_geometry(
            const LineData& line_data, const geode::uuid& line_id )
        {
            auto line_builder = builder_.line_mesh_builder( line_id );
            for( const auto& point : line_data.points )
            {
                line_builder->create_point( point );
            }
            for( auto i : geode::Range( 1, line_data.points.size() ) )
            {
                line_builder->create_edge( i - 1, i );
            }
        }

        void register_line_data(
            const LineData& line_data, const geode::uuid& line_id )
        {
            const auto& surface = model_.surface( line_data.surface );
            builder_.add_line_surface_relationship(
                model_.line( line_id ), surface );
            auto surface_id = surface.component_id();
            for( auto i : geode::Range{ line_data.indices.size() } )
            {
                auto unique_id = model_.unique_vertex(
                    { model_.line( line_id ).component_id(), i } );
                builder_.set_unique_vertex(
                    { surface_id, line_data.indices[i] }, unique_id );
            }
        }

        LineData compute_line( const geode::Surface3D& surface,
            const TopoInfo::LineStart& line_start )
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
            auto edge = find_edge(
                mesh, line_start.second.vertex, line_start.first.vertex );
            while( model_.unique_vertex( { surface.component_id(),
                       mesh.polygon_edge_vertex( edge, 0 ) } )
                   == geode::NO_ID )
            {
                result.indices.push_back( mesh.polygon_edge_vertex( edge, 0 ) );
                result.points.emplace_back(
                    mesh.point( result.indices.back() ) );
                edge = mesh.previous_on_border( edge );
            }

            result.indices.push_back( mesh.polygon_edge_vertex( edge, 0 ) );
            result.points.emplace_back( mesh.point( result.indices.back() ) );
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

        void check_keyword( std::string keyword )
        {
            std::string line;
            std::getline( file_, line );
            OPENGEODE_EXCEPTION( string_starts_with( line, keyword ),
                "Line should starts with \"" + keyword + "\"" );
        }

        void read_header()
        {
            check_keyword( "HEADER" );
            std::string line;
            while( std::getline( file_, line ) )
            {
                if( string_starts_with( line, "}" ) )
                {
                    return;
                }
            }
            throw geode::OpenGeodeException(
                "Cannot find the end of \"HEADER\" section" );
        }

        void read_CRS()
        {
            check_keyword( "GOCAD_ORIGINAL_COORDINATE_SYSTEM" );
            std::string line;
            while( std::getline( file_, line ) )
            {
                if( string_starts_with(
                        line, "END_ORIGINAL_COORDINATE_SYSTEM" ) )
                {
                    return;
                }
            }
            throw geode::OpenGeodeException{
                "Cannot find the end of CRS section"
            };
        }

        void read_model_components()
        {
            std::string line;
            while( std::getline( file_, line ) )
            {
                if( string_starts_with( line, "END" ) )
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
            }
            throw geode::OpenGeodeException{
                "Cannot find the end of component section"
            };
        }

        void create_tsurfs()
        {
            for( const auto& tsurf : tsurfs_ )
            {
                if( tsurf.feature.find( "fault" ) != std::string::npos )
                {
                    const auto& fault_uuid =
                        builder_.add_fault( fault_map_.at( tsurf.feature ) );
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
                    // todo handle ModelBoundary
                }
                else
                {
                    const auto& horizon_uuid = builder_.add_horizon(
                        horizon_map_.at( tsurf.feature ) );
                    const auto& horizon = model_.horizon( horizon_uuid );
                    for( const auto& uuid : tsurf.tfaces )
                    {
                        builder_.add_surface_in_horizon(
                            model_.surface( uuid ), horizon );
                    }
                }
            }
        }

        void process_TSURF_keyword( std::istringstream& iss )
        {
            std::string name;
            iss >> name;
            while( !iss.eof() )
            {
                std::string token;
                iss >> token;
                name += "_" + token;
            }
            tsurf_names2index_.emplace( std::move( name ), tsurfs_.size() );
            tsurfs_.emplace_back();
        }

        void process_TFACE_keyword( std::istringstream& iss )
        {
            geode::index_t id;
            std::string feature, name;
            iss >> id >> feature >> name;
            while( !iss.eof() )
            {
                std::string token;
                iss >> token;
                name += "_" + token;
            }
            const auto& surface_id = builder_.add_surface(
                geode::OpenGeodeTriangulatedSurface3D::type_name_static() );
            auto& tsurf = tsurfs_[tsurf_names2index_.at( name )];
            tsurf.feature = feature;
            tsurf.tfaces.emplace_back( surface_id );
            surfaces_.emplace_back( surface_id );
        }

        void process_REGION_keyword( std::istringstream& iss )
        {
            geode::index_t id;
            std::string name;
            iss >> id >> name;
            if( name == "Universe" )
            {
                return;
            }
            create_block_topology( builder_.add_block() );
        }

        void create_block_topology( const geode::uuid& block_id )
        {
            const auto& block = model_.block( block_id );
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
                    builder_.add_surface_block_relationship(
                        model_.surface( surfaces_[surface_id] ), block );
                }
            }
        }

        void read_tsurf( TSurfData& tsurf )
        {
            check_keyword( "GOCAD TSurf" );
            read_header();
            read_CRS();
            read_surfaces( tsurf );
            read_topological_information( tsurf );
        }

        void read_topological_information( TSurfData& tsurf )
        {
            std::string line;
            while( std::getline( file_, line ) )
            {
                if( string_starts_with( line, "END" ) )
                {
                    return;
                }

                std::istringstream iss{ line };
                std::string keyword;
                iss >> keyword;
                if( keyword == "BSTONE" )
                {
                    process_BSTONE_keyword( iss, tsurf );
                }
                else if( keyword == "BORDER" )
                {
                    process_BORDER_keyword( iss, tsurf );
                }
            }
        }

        void process_BSTONE_keyword( std::istringstream& iss, TSurfData& tsurf )
        {
            geode::index_t corner;
            iss >> corner;
            corner -= OFFSET_START;
            auto tface_id = tsurf.tface_id( corner );
            const auto& surface = model_.surface( tsurf.tfaces[tface_id] );
            auto vertex_id = corner - tsurf.tface_vertices_offset[tface_id];
            info_.corner_points.push_back( surface.mesh().point( vertex_id ) );
            info_.corner_surface_index.emplace_back(
                surface.component_id(), vertex_id );
        }

        void process_BORDER_keyword( std::istringstream& iss, TSurfData& tsurf )
        {
            geode::index_t id, corner, next;
            iss >> id >> corner >> next;
            corner -= OFFSET_START;
            next -= OFFSET_START;
            auto tface_id = tsurf.tface_id( corner );
            const auto& surface = model_.surface( tsurf.tfaces[tface_id] );
            auto corner_id = corner - tsurf.tface_vertices_offset[tface_id];
            auto next_id = next - tsurf.tface_vertices_offset[tface_id];
            info_.line_surface_index.emplace_back(
                geode::MeshComponentVertex{ surface.component_id(), corner_id },
                geode::MeshComponentVertex{ surface.component_id(), next_id } );
        }

        void read_surfaces( TSurfData& tsurf )
        {
            goto_keyword( "TFACE" );
            for( const auto& uuid : tsurf.tfaces )
            {
                std::unique_ptr< geode::TriangulatedSurfaceBuilder3D >
                    mesh_builder{
                        dynamic_cast< geode::TriangulatedSurfaceBuilder3D* >(
                            builder_.surface_mesh_builder( uuid ).release() )
                    };
                std::string line;
                geode::index_t nb_vertices{ 0 };
                auto vertex_offset = tsurf.tface_vertices_offset.back();
                auto previous_file_position = file_.tellg();
                while( std::getline( file_, line ) )
                {
                    std::istringstream iss{ line };
                    std::string keyword;
                    iss >> keyword;
                    if( keyword == "VRTX" || keyword == "PVRTX" )
                    {
                        mesh_builder->create_point(
                            process_VRTX_keyword( iss ) );
                        nb_vertices++;
                    }
                    else if( keyword == "ATOM" || keyword == "PATOM" )
                    {
                        mesh_builder->create_point(
                            process_ATOM_keyword( iss, tsurf ) );
                        nb_vertices++;
                    }
                    else if( keyword == "TRGL" )
                    {
                        geode::index_t v0, v1, v2;
                        iss >> v0 >> v1 >> v2;
                        mesh_builder->create_triangle(
                            { v0 - OFFSET_START - vertex_offset,
                                v1 - OFFSET_START - vertex_offset,
                                v2 - OFFSET_START - vertex_offset } );
                    }
                    else if( keyword == "TFACE" )
                    {
                        tsurf.tface_vertices_offset.push_back(
                            vertex_offset + nb_vertices );
                        break;
                    }
                    else
                    {
                        file_.seekg( previous_file_position );
                        break;
                    }
                    previous_file_position = file_.tellg();
                }
                mesh_builder->compute_polygon_adjacencies();
            }
        }

        geode::Point3D process_VRTX_keyword( std::istringstream& iss )
        {
            geode::index_t dummy;
            double x, y, z;
            iss >> dummy >> x >> y >> z;
            return { { x, y, z } };
        }

        const geode::Point3D& process_ATOM_keyword(
            std::istringstream& iss, TSurfData& tsurf )
        {
            geode::index_t dummy, atom_id;
            iss >> dummy >> atom_id;
            atom_id -= OFFSET_START;
            auto tface_id = tsurf.tface_id( atom_id );
            const auto& surface = model_.surface( tsurf.tfaces[tface_id] );
            auto vertex_id = atom_id - tsurf.tface_vertices_offset[tface_id];
            return surface.mesh().point( vertex_id );
        }

        void goto_keyword( std::string word )
        {
            std::string line;
            while( std::getline( file_, line ) )
            {
                if( string_starts_with( line, word ) )
                {
                    return;
                }
            }
            throw geode::OpenGeodeException(
                "Cannot find the requested word: ", word );
        }

    private:
        std::ifstream file_;
        geode::StructuralModel& model_;
        geode::StructuralModelBuilder builder_;
        std::unordered_map< std::string, geode::index_t > tsurf_names2index_;
        std::deque< TSurfData > tsurfs_;
        TopoInfo info_;
        std::unordered_map< std::pair< geode::uuid, geode::uuid >, geode::uuid >
            corners2line_;
        std::deque< geode::uuid > surfaces_;
        double epsilon_;

        std::unordered_map< std::string, geode::Fault3D::FAULT_TYPE >
            fault_map_ = { { "fault", geode::Fault3D::FAULT_TYPE::NO_TYPE },
                { "reverse_fault", geode::Fault3D::FAULT_TYPE::REVERSE },
                { "normal_fault", geode::Fault3D::FAULT_TYPE::NORMAL } };
        std::unordered_map< std::string, geode::Horizon3D::HORIZON_TYPE >
            horizon_map_ = { { "top", geode::Horizon3D::HORIZON_TYPE::NO_TYPE },
                { "none", geode::Horizon3D::HORIZON_TYPE::NO_TYPE },
                { "topographic", geode::Horizon3D::HORIZON_TYPE::NO_TYPE },
                { "unconformity",
                    geode::Horizon3D::HORIZON_TYPE::NON_CONFORMAL } };
    };
} // namespace

namespace geode
{
    void MLInput::read()
    {
        MLInputImpl impl{ filename(), structural_model() };
        impl.read_file();
    }
} // namespace geode
