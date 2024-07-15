/*
 * Copyright (c) 2019 - 2024 Geode-solutions
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

#include <geode/geosciences_io/model/internal/lso_output.h>

#include <fstream>
#include <string>
#include <vector>

#include <geode/basic/attribute_manager.h>

#include <geode/geometry/point.h>

#include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/mesh_factory.h>
#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/mesh/core/triangulated_surface.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/corner.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/mixin/core/vertex_identifier.h>

#include <geode/geosciences/explicit/representation/core/structural_model.h>
#include <geode/geosciences_io/mesh/internal/gocad_common.h>
#include <geode/geosciences_io/model/internal/gocad_common.h>

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
            OPENGEODE_EXCEPTION(
                file_.good(), "[LSOOutput] Error while opening file: ", file );
            for( const auto& block : model_.blocks() )
            {
                OPENGEODE_EXCEPTION(
                    block.mesh().type_name()
                        == geode::TetrahedralSolid3D::type_name_static(),
                    "[LSOOutput] Only support TetrahedralSolid3D" );
            }
        }

        void write_file()
        {
            file_ << "GOCAD LightTSolid 1" << EOL;
            geode::internal::HeaderData header;
            header.name = geode::to_string( model_.name() );
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
            write_regions();
        }

        void write_regions()
        {
            for( const auto& block : model_.blocks() )
            {
                file_ << "MODEL_REGION " << block.name() << " ";
                const auto& surface = *model_.boundaries( block ).begin();
                if( sides_.regions_surface_sides.at(
                        { block.id(), surface.id() } ) )
                {
                    file_ << "+";
                }
                else
                {
                    file_ << "-";
                }
                file_ << exported_surfaces_.at( surface.id() ) << EOL;
            }
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
                    if( !exported_surfaces_.contains( item_id.id() ) )
                    {
                        all_exported = false;
                        break;
                    }
                }
                if( all_exported )
                {
                    continue;
                }
                file_ << "SURFACE " << component.name() << EOL;
                for( const auto& item_id : model_.items( component.id() ) )
                {
                    if( !exported_surfaces_.emplace( item_id.id(), nb_tfaces )
                             .second )
                    {
                        continue;
                    }
                    file_ << "TFACE " << nb_tfaces++ << EOL;
                    const auto& surface = model_.surface( item_id.id() );
                    const auto& mesh = surface.mesh();
                    file_ << "KEYVERTICES";
                    write_triangle( mesh, item_id, 0 );
                    file_ << EOL;
                    for( const auto p : geode::Range{ mesh.nb_polygons() } )
                    {
                        file_ << "TRGL";
                        write_triangle( mesh, item_id, p );
                        file_ << EOL;
                    }
                }
            }
            return nb_tfaces;
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
                const auto name = block.name();
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
                const auto block_vertices = model_.component_mesh_vertices(
                    v, geode::Block3D::component_type_static() );
                auto& cmvs = vertices_[v];
                cmvs.reserve( block_vertices.size() );
                const auto first = first_block( block_vertices );
                const auto& first_vertex = block_vertices[first];
                cmvs.emplace( first_vertex, v + OFFSET_START );
                const auto& block =
                    model_.block( first_vertex.component_id.id() );
                const auto& mesh = block.mesh();
                file_ << "VRTX " << count++ << " "
                      << mesh.point( first_vertex.vertex ).string() << EOL;
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

        geode::index_t first_block(
            const std::vector< geode::ComponentMeshVertex >& block_vertices )
            const
        {
            geode::index_t id{ 0 };
            auto first_uuid = block_vertices.front().component_id.id();
            for( const auto i : geode::Range{ 1, block_vertices.size() } )
            {
                const auto& uuid = block_vertices[i].component_id.id();
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

namespace geode
{
    namespace internal
    {
        std::vector< std::string > LSOOutput::write(
            const StructuralModel& structural_model ) const
        {
            LSOOutputImpl impl{ filename(), structural_model };
            impl.write_file();
            return { to_string( filename() ) };
        }

        bool LSOOutput::is_saveable(
            const StructuralModel& structural_model ) const
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
                    || mesh.type_name()
                           != TetrahedralSolid3D::type_name_static() )
                {
                    return false;
                }
            }
            return true;
        }
    } // namespace internal
} // namespace geode
