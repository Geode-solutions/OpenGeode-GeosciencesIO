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

#pragma once

#include <fstream>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_replace.h>

#include <geode/basic/filename.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/helpers/detail/build_model_boundaries.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/corner.hpp>
#include <geode/model/mixin/core/line.hpp>
#include <geode/model/mixin/core/model_boundary.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>
#include <geode/geosciences_io/model/internal/gocad_common.hpp>

namespace geode::internal
{
    template < typename Model >
    Model clone_with_model_boundaries( const Model& model )
    {
        Model result;
        typename Model::Builder builder{ result };
        builder.copy( model );
        detail::build_model_boundaries( result, builder );
        return result;
    }

    inline bool check_brep_polygons( const BRep& brep )
    {
        for( const auto& surface : brep.surfaces() )
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

    template < typename Model >
    class MLOutputImpl
    {
    public:
        static constexpr index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };

        virtual ~MLOutputImpl() = default;

        void write_file()
        {
            file_ << "GOCAD Model3d 1" << EOL;
            internal::HeaderData header;
            if( const auto name = model_.name() )
            {
                header.name = name.value();
            }
            internal::write_header( file_, header );
            internal::write_CRS( file_, {} );
            assign_surfaces_to_tsurfs();
            write_model_components();
            write_model_surfaces();
        }

    protected:
        using OwnedSurfaces =
            std::vector< std::reference_wrapper< const Surface3D > >;

        MLOutputImpl( std::string_view filename, const Model& model )
            : file_{ to_string( filename ) },
              model_( model ),
              sides_( determine_surface_to_regions_sides( model ) )
        {
            OpenGeodeGeosciencesIOModelException::check_exception( file_.good(),
                nullptr, OpenGeodeException::TYPE::data,
                "[MLOutput] Error while opening file: ", filename );
        }

        index_t& component_id()
        {
            return component_id_;
        }

        absl::flat_hash_map< uuid, index_t >& components()
        {
            return components_;
        }

        std::ofstream& file()
        {
            return file_;
        }

        /*!
         * In a .ml file, each TFACE belongs to a single TSURF. A Surface
         * shared by several collections (e.g. a topographic Horizon that is
         * also a ModelBoundary) is assigned to the first one, in the export
         * order: ModelBoundaries, then geological collections. A collection
         * without any Surface of its own is not exported.
         */
        void assign_surfaces_to_tsurfs()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                assign_surfaces(
                    boundary, model_.model_boundary_items( boundary ) );
            }
            assign_geological_surfaces();
        }

        virtual void assign_geological_surfaces() {}

        template < typename Collection, typename ItemRange >
        void assign_surfaces(
            const Collection& collection, const ItemRange& items )
        {
            auto& surfaces = owned_surfaces_[collection.id()];
            for( const auto& item : items )
            {
                if( !assigned_surfaces_.insert( item.id() ).second )
                {
                    Logger::warning( "[MLOutput] A Surface from ",
                        component_name( collection ),
                        " belongs to several collections. It has been "
                        "exported only once" );
                    continue;
                }
                surfaces.emplace_back( item );
            }
            if( surfaces.empty() )
            {
                Logger::warning( "[MLOutput] ", component_name( collection ),
                    " has no Surface of its own. It has not been exported" );
            }
        }

        const OwnedSurfaces& owned_surfaces( const uuid& collection_id ) const
        {
            return owned_surfaces_.at( collection_id );
        }

        template < typename Component >
        void write_key_triangle( const Component& component )
        {
            const auto& mesh = component.mesh();
            for( const auto v : LRange{ 3 } )
            {
                file_ << SPACE << SPACE
                      << mesh.point( mesh.polygon_vertex( { 0, v } ) ).string()
                      << EOL;
            }
        }

        virtual void write_geological_tsurfs() = 0;

        /*!
         * Surfaces not owned by any exported collection: each one is
         * exported as its own TSURF
         */
        std::vector< uuid > unclassified_tsurfs() const
        {
            std::vector< uuid > result;
            for( const auto& surface : model_.surfaces() )
            {
                if( !assigned_surfaces_.contains( surface.id() ) )
                {
                    result.push_back( surface.id() );
                }
            }
            return result;
        }

        void write_tsurfs()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                if( !owned_surfaces( boundary.id() ).empty() )
                {
                    file_ << "TSURF " << tsurf_name( boundary ) << EOL;
                }
            }
            write_geological_tsurfs();
            for( const auto& surface_id : unclassified_tsurfs() )
            {
                file_ << "TSURF " << tsurf_name( model_.surface( surface_id ) )
                      << EOL;
                unclassified_surfaces_.emplace_back( surface_id );
            }
        }

        virtual void write_geological_tfaces() = 0;

        void write_collection_tfaces( const OwnedSurfaces& surfaces,
            std::string_view type,
            std::string_view name )
        {
            for( const auto& surface : surfaces )
            {
                file_ << "TFACE " << component_id_ << SPACE << type << SPACE
                      << name << EOL;
                write_key_triangle( surface.get() );
                components_.emplace( surface.get().id(), component_id_++ );
            }
        }

        void write_tfaces()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                write_collection_tfaces( owned_surfaces( boundary.id() ),
                    "boundary", tsurf_name( boundary ) );
            }
            write_geological_tfaces();
            for( const auto& surface_id : unclassified_surfaces_ )
            {
                const auto& surface = model_.surface( surface_id );
                file_ << "TFACE " << component_id_ << SPACE << "boundary"
                      << SPACE << tsurf_name( surface ) << EOL;
                write_key_triangle( surface );
                components_.emplace( surface_id, component_id_++ );
            }
        }

        void write_universe()
        {
            file_ << "REGION " << component_id_ << SPACE << SPACE << "Universe "
                  << EOL << SPACE << SPACE;
            index_t counter{ 0 };
            for( const auto& boundary : model_.model_boundaries() )
            {
                for( const Surface3D& item : owned_surfaces( boundary.id() ) )
                {
                    char sign{ '-' };
                    if( sides_.universe_surface_sides.at( item.id() ) )
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
                file_ << "REGION " << component_id_ << SPACE
                      << component_name( region ) << EOL << SPACE << SPACE;
                index_t counter{ 0 };
                for( const auto& surface : model_.boundaries( region ) )
                {
                    const auto sign = sides_.regions_surface_sides.at(
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

        index_t write_surface(
            const Surface3D& surface, const index_t current_offset )
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
            const auto v0 = mesh.polygon_vertex( PolygonVertex{ edge } );
            const auto v1 = mesh.polygon_vertex( { edge.polygon_id,
                static_cast< local_index_t >( ( edge.edge_id + 1 ) % 3 ) } );
            const auto uid1 =
                model_.unique_vertex( { surface.component_id(), v1 } );
            bool is_corner{ false };
            for( const auto& cmv1 : model_.component_mesh_vertices( uid1 ) )
            {
                if( cmv1.component_id.type
                    == Corner3D::component_type_static() )
                {
                    is_corner = true;
                    break;
                }
            }
            if( is_corner )
            {
                line_starts.emplace_back( std::array< index_t, 2 >{
                    v1 + current_offset, v0 + current_offset } );
            }
        }

        void add_corners_and_line_starts( const Surface3D& surface,
            const index_t current_offset,
            std::vector< std::array< index_t, 2 > >& line_starts ) const
        {
            // Each border loop (outer border, holes) is walked once
            const auto& mesh = surface.mesh();
            std::vector< bool > visited( 3 * mesh.nb_polygons(), false );
            for( const auto p : Range{ mesh.nb_polygons() } )
            {
                for( const auto e : LRange{ 3 } )
                {
                    const PolygonEdge first_on_border{ p, e };
                    if( !mesh.is_edge_on_border( first_on_border )
                        || visited[3 * p + e] )
                    {
                        continue;
                    }
                    auto cur = first_on_border;
                    do
                    {
                        visited[3 * cur.polygon_id + cur.edge_id] = true;
                        process_surface_edge(
                            surface, cur, current_offset, line_starts );
                        cur = mesh.previous_on_border( cur );
                    } while( cur != first_on_border );
                }
            }
        }

        void find_corners_and_line_starts( const OwnedSurfaces& surfaces,
            std::vector< std::array< index_t, 2 > >& line_starts ) const
        {
            index_t current_offset{ OFFSET_START };
            for( const Surface3D& surface : surfaces )
            {
                add_corners_and_line_starts(
                    surface, current_offset, line_starts );
                current_offset += surface.mesh().nb_vertices();
            }
        }

        std::vector< std::array< index_t, 2 > >
            find_corners_and_line_starts_for_unclassified_surface(
                const uuid& surface_id ) const
        {
            std::vector< std::array< index_t, 2 > > line_starts;
            const auto& surface = model_.surface( surface_id );
            add_corners_and_line_starts( surface, OFFSET_START, line_starts );
            return line_starts;
        }

        void write_corners(
            absl::Span< const std::array< index_t, 2 > > line_starts )
        {
            for( const auto& line_start : line_starts )
            {
                file_ << "BSTONE " << line_start[0] << EOL;
            }
        }

        void write_line_starts( index_t current_offset,
            absl::Span< const std::array< index_t, 2 > > line_starts )
        {
            for( const auto& line_start : line_starts )
            {
                file_ << "BORDER " << current_offset++ << SPACE << line_start[0]
                      << SPACE << line_start[1] << EOL;
            }
        }

        virtual void write_geological_model_surfaces() = 0;

        void write_collection_surface( const OwnedSurfaces& surfaces,
            std::string_view type,
            std::string_view name )
        {
            if( surfaces.empty() )
            {
                return;
            }
            file_ << "GOCAD TSurf 1" << EOL;
            internal::HeaderData header;
            header.name = to_string( name );
            internal::write_header( file_, header );
            internal::write_CRS( file_, {} );
            file_ << "GEOLOGICAL_FEATURE " << name << EOL;
            file_ << "GEOLOGICAL_TYPE " << type << EOL;
            index_t current_offset{ OFFSET_START };
            for( const auto& surface : surfaces )
            {
                file_ << "TFACE" << EOL;
                current_offset = write_surface( surface.get(), current_offset );
            }
            std::vector< std::array< index_t, 2 > > line_starts;
            find_corners_and_line_starts( surfaces, line_starts );
            write_corners( line_starts );
            write_line_starts( current_offset, line_starts );
            file_ << "END" << EOL;
        }

        void write_model_surfaces()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                write_collection_surface( owned_surfaces( boundary.id() ),
                    "boundary", tsurf_name( boundary ) );
            }
            write_geological_model_surfaces();
            for( const auto& surface_id : unclassified_surfaces_ )
            {
                file_ << "GOCAD TSurf 1" << EOL;
                const auto& surface = model_.surface( surface_id );
                internal::HeaderData header;
                header.name = tsurf_name( surface );
                internal::write_header( file_, header );
                internal::write_CRS( file_, {} );
                file_ << "GEOLOGICAL_FEATURE " << tsurf_name( surface ) << EOL;
                file_ << "GEOLOGICAL_TYPE " << "boundary" << EOL;
                index_t current_offset{ OFFSET_START };
                file_ << "TFACE" << EOL;
                current_offset = write_surface( surface, current_offset );
                const auto line_starts =
                    find_corners_and_line_starts_for_unclassified_surface(
                        surface_id );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                file_ << "END" << EOL;
            }
        }

        template < typename Component >
        std::string component_name( const Component& component ) const
        {
            return absl::StrReplaceAll(
                component.name().value_or( component.id().string() ),
                { { " ", "_" } } );
        }

        /*!
         * TSURFs are matched by name in a .ml file: the name of a collection
         * (or of an unclassified Surface) is made unique by appending a
         * suffix if another TSURF already uses it.
         */
        template < typename Component >
        const std::string& tsurf_name( const Component& component )
        {
            const auto it = tsurf_names_.find( component.id() );
            if( it != tsurf_names_.end() )
            {
                return it->second;
            }
            auto name = component_name( component );
            if( !used_tsurf_names_.insert( name ).second )
            {
                index_t suffix{ 2 };
                auto unique_name = absl::StrCat( name, "_", suffix );
                while( !used_tsurf_names_.insert( unique_name ).second )
                {
                    unique_name = absl::StrCat( name, "_", ++suffix );
                }
                Logger::warning( "[MLOutput] TSURF name ", name,
                    " is already used. It has been exported as ", unique_name );
                name = std::move( unique_name );
            }
            return tsurf_names_.emplace( component.id(), std::move( name ) )
                .first->second;
        }

    private:
        std::ofstream file_;
        const Model& model_;
        const RegionSurfaceSide sides_;
        absl::flat_hash_map< uuid, index_t > components_;
        index_t component_id_{ OFFSET_START };
        std::vector< uuid > unclassified_surfaces_;
        absl::flat_hash_set< uuid > assigned_surfaces_;
        absl::flat_hash_map< uuid, OwnedSurfaces > owned_surfaces_;
        absl::flat_hash_map< uuid, std::string > tsurf_names_;
        absl::flat_hash_set< std::string > used_tsurf_names_;
    };
} // namespace geode::internal
