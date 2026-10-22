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

#include <geode/geosciences_io/model/internal/lso_output.hpp>

#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>

#include <geode/basic/logger.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/mesh_factory.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/core/triangulated_surface.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/corner.hpp>
#include <geode/model/mixin/core/line.hpp>
#include <geode/model/mixin/core/model_boundary.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/mixin/core/vertex_identifier.hpp>

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>
#include <geode/geosciences_io/model/internal/gocad_common.hpp>

namespace
{
    class LSOOutputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };

        LSOOutputImpl(
            std::string_view file, const geode::StructuralModel& model )
            : file_{ geode::to_string( file ) },
              model_( model ),
              sides_(
                  geode::internal::determine_surface_to_regions_sides( model ) )
        {
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                file_.good(), nullptr, geode::OpenGeodeException::TYPE::data,
                "[LSOOutput] Error while opening file: ", file );
            file_.precision( std::numeric_limits< double >::max_digits10 );
            for( const auto& block : model_.blocks() )
            {
                geode::OpenGeodeGeosciencesIOModelException::check_exception(
                    block.mesh().type_name()
                        == geode::TetrahedralSolid3D::type_name_static(),
                    nullptr, geode::OpenGeodeException::TYPE::data,
                    "[LSOOutput] Only support TetrahedralSolid3D" );
            }
            for( const auto& surface : model_.surfaces() )
            {
                if( surface.mesh().nb_polygons() == 0 )
                {
                    geode::Logger::warning( "[LSOOutput] Surface ",
                        surface.name().value_or( surface.id().string() ),
                        " is empty and will not be exported" );
                    continue;
                }
                geode::OpenGeodeGeosciencesIOModelException::check_exception(
                    surface.mesh().type_name()
                        == geode::TriangulatedSurface3D::type_name_static(),
                    nullptr, geode::OpenGeodeException::TYPE::data,
                    "[LSOOutput] Only support TriangulatedSurface3D" );
            }
        }

        void write_file()
        {
            file_ << "GOCAD LightTSolid 1" << EOL;
            geode::internal::HeaderData header;
            if( const auto name = model_.name() )
            {
                header.name = name.value();
            }
            geode::internal::write_header( file_, header );
            geode::internal::write_CRS( file_, {} );
            write_vertices();
            write_tetrahedron();
            write_model();
            file_ << "END" << EOL;
        }

    private:
        void write_model()
        {
            file_ << "MODEL" << EOL;
            auto nb_tfaces = write_surfaces( model_.horizons(), 1 );
            nb_tfaces = write_surfaces( model_.faults(), nb_tfaces );
            nb_tfaces = write_surfaces( model_.model_boundaries(), nb_tfaces );
            write_remaining_surfaces( nb_tfaces );
            write_regions();
        }

        void write_regions()
        {
            for( const auto& block : model_.blocks() )
            {
                file_ << "MODEL_REGION " << region_name( block ) << " ";
                write_region_boundary( block );
                file_ << EOL;
            }
        }

        void write_region_boundary( const geode::Block3D& block )
        {
            for( const auto& surface : model_.boundaries( block ) )
            {
                const auto tface = exported_surfaces_.find( surface.id() );
                if( tface == exported_surfaces_.end() )
                {
                    continue;
                }
                const auto side = sides_.regions_surface_sides.find(
                    { block.id(), surface.id() } );
                if( side == sides_.regions_surface_sides.end() )
                {
                    continue;
                }
                file_ << ( side->second ? "+" : "-" ) << tface->second;
                return;
            }
            throw geode::OpenGeodeGeosciencesIOModelException{ nullptr,
                geode::OpenGeodeException::TYPE::data, "[LSOOutput] Block ",
                block.name().value_or( block.id().string() ),
                " has no oriented boundary Surface to define its region" };
        }

        static std::string region_name( const geode::Block3D& block )
        {
            auto name = block.name().value_or( block.id().string() );
            if( absl::StrContains( name, ' ' ) )
            {
                return absl::StrCat( "\"", name, "\"" );
            }
            return name;
        }

        template < typename Components >
        geode::index_t write_surfaces(
            const Components& components, geode::index_t nb_tfaces )
        {
            for( const auto& component : components )
            {
                bool all_exported{ true };
                for( const auto& item_id : model_.items( component.id() ) )
                {
                    if( is_exportable( item_id.id ) )
                    {
                        all_exported = false;
                        break;
                    }
                }
                if( all_exported )
                {
                    continue;
                }
                file_ << "SURFACE "
                      << component.name().value_or( component.id().string() )
                      << EOL;
                for( const auto& item_id : model_.items( component.id() ) )
                {
                    if( is_exportable( item_id.id ) )
                    {
                        write_tface( model_.surface( item_id.id ), nb_tfaces );
                    }
                }
            }
            return nb_tfaces;
        }

        void write_remaining_surfaces( geode::index_t nb_tfaces )
        {
            for( const auto& surface : model_.surfaces() )
            {
                if( !is_exportable( surface.id() ) )
                {
                    continue;
                }
                file_ << "SURFACE "
                      << surface.name().value_or( surface.id().string() )
                      << EOL;
                write_tface( surface, nb_tfaces );
            }
        }

        bool is_exportable( const geode::uuid& surface_id ) const
        {
            return !exported_surfaces_.contains( surface_id )
                   && model_.surface( surface_id ).mesh().nb_polygons() > 0;
        }

        void write_tface(
            const geode::Surface3D& surface, geode::index_t& nb_tfaces )
        {
            exported_surfaces_.emplace( surface.id(), nb_tfaces );
            file_ << "TFACE " << nb_tfaces++ << EOL;
            const auto& mesh = surface.mesh();
            const auto& id = surface.component_id();
            file_ << "KEYVERTICES";
            write_triangle( mesh, id, 0 );
            file_ << EOL;
            for( const auto p : geode::Range{ mesh.nb_polygons() } )
            {
                file_ << "TRGL";
                write_triangle( mesh, id, p );
                file_ << EOL;
            }
        }

        void write_triangle( const geode::SurfaceMesh3D& mesh,
            const geode::ComponentID& id,
            geode::index_t p )
        {
            for( const auto v : geode::LRange{ 3 } )
            {
                const auto vertex = mesh.polygon_vertex( { p, v } );
                const auto unique = model_.unique_vertex( { id, vertex } );
                file_ << " " << unique + OFFSET_START;
            }
        }

        void write_tetrahedron()
        {
            for( const auto& block : model_.blocks() )
            {
                const auto name = region_name( block );
                const auto& id = block.component_id();
                const auto& mesh = block.mesh();
                for( const auto p : geode::Range{ mesh.nb_polyhedra() } )
                {
                    file_ << "TETRA";
                    for( const auto v : geode::LRange{ 4 } )
                    {
                        const auto vertex = mesh.polyhedron_vertex( { p, v } );
                        const auto unique =
                            model_.unique_vertex( { id, vertex } );
                        file_ << " "
                              << vertices_.at( unique ).at( { id, vertex } );
                    }
                    file_ << EOL;
                    file_ << "# CTETRA " << name << " none none none none"
                          << EOL;
                }
            }
        }

        void write_vertices()
        {
            geode::index_t count{ 1 };
            auto nb_vertices = model_.nb_unique_vertices();
            vertices_.resize( nb_vertices );
            std::vector< geode::index_t > atoms;
            for( const auto v : geode::Range{ nb_vertices } )
            {
                std::vector< geode::ComponentMeshVertex > block_vertices;
                const auto cmv = model_.component_mesh_vertices( v );
                for( const auto& mv : cmv )
                {
                    if( mv.component_id.type
                        != geode::Block3D::component_type_static() )
                    {
                        continue;
                    }
                    block_vertices.push_back( mv );
                }
                if( block_vertices.empty() )
                {
                    geode::OpenGeodeGeosciencesIOModelException::
                        check_exception( !cmv.empty(), nullptr,
                            geode::OpenGeodeException::TYPE::data,
                            "[LSOOutput] Unique vertex ", v,
                            " is not linked to any component mesh vertex" );
                    file_ << "VRTX " << count++;
                    write_point( component_point( cmv.front() ) );
                    continue;
                }
                auto& cmvs = vertices_[v];
                cmvs.reserve( block_vertices.size() );
                const auto first = first_block( block_vertices );
                const auto& first_vertex = block_vertices[first];
                cmvs.emplace( first_vertex, v + OFFSET_START );
                const auto& block =
                    model_.block( first_vertex.component_id.id );
                const auto& mesh = block.mesh();
                file_ << "VRTX " << count++;
                write_point( mesh.point( first_vertex.vertex ) );
                for( const auto i : geode::Indices{ block_vertices } )
                {
                    if( i == first )
                    {
                        continue;
                    }
                    atoms.push_back( v + OFFSET_START );
                    cmvs.emplace(
                        block_vertices[i], nb_vertices + OFFSET_START );
                    nb_vertices++;
                }
            }
            for( const auto a : atoms )
            {
                file_ << "SHAREDVRTX " << count++ << " " << a << EOL;
            }
        }

        void write_point( const geode::Point3D& point )
        {
            file_ << " " << point.value( 0 ) << " " << point.value( 1 ) << " "
                  << point.value( 2 ) << EOL;
        }

        const geode::Point3D& component_point(
            const geode::ComponentMeshVertex& cmv ) const
        {
            const auto& type = cmv.component_id.type;
            const auto& id = cmv.component_id.id;
            if( type == geode::Corner3D::component_type_static() )
            {
                return model_.corner( id ).mesh().point( cmv.vertex );
            }
            if( type == geode::Line3D::component_type_static() )
            {
                return model_.line( id ).mesh().point( cmv.vertex );
            }
            if( type == geode::Surface3D::component_type_static() )
            {
                return model_.surface( id ).mesh().point( cmv.vertex );
            }
            return model_.block( id ).mesh().point( cmv.vertex );
        }

        geode::index_t first_block(
            const std::vector< geode::ComponentMeshVertex >& block_vertices )
            const
        {
            geode::index_t id{ 0 };
            auto first_uuid = block_vertices.front().component_id.id;
            for( const auto i : geode::Range{ 1, block_vertices.size() } )
            {
                const auto& uuid = block_vertices[i].component_id.id;
                if( uuid < first_uuid )
                {
                    first_uuid = uuid;
                    id = i;
                }
            }
            return id;
        }

    private:
        std::ofstream file_;
        const geode::StructuralModel& model_;
        const geode::internal::RegionSurfaceSide sides_;
        std::vector<
            absl::flat_hash_map< geode::ComponentMeshVertex, geode::index_t > >
            vertices_;
        absl::flat_hash_map< geode::uuid, geode::index_t > exported_surfaces_;
    };
} // namespace

namespace geode::internal
{
    std::vector< std::string > LSOOutput::write(
        const StructuralModel& structural_model ) const
    {
        LSOOutputImpl impl{ filename(), structural_model };
        impl.write_file();
        return { to_string( filename() ) };
    }

    bool LSOOutput::is_saveable( const StructuralModel& structural_model ) const
    {
        for( const auto& surface : structural_model.surfaces() )
        {
            const auto& mesh = surface.mesh();
            if( mesh.nb_polygons() == 0
                || mesh.type_name()
                       != TriangulatedSurface3D::type_name_static() )
            {
                return false;
            }
        }
        for( const auto& block : structural_model.blocks() )
        {
            const auto& mesh = block.mesh();
            if( mesh.nb_polyhedra() == 0
                || mesh.type_name() != TetrahedralSolid3D::type_name_static() )
            {
                return false;
            }
        }
        for( const auto v :
            geode::Range{ structural_model.nb_unique_vertices() } )
        {
            if( structural_model.component_mesh_vertices( v ).empty() )
            {
                return false;
            }
        }
        return true;
    }
} // namespace geode::internal
