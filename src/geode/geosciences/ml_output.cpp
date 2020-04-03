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
#include <unordered_map>

// #include <geode/geometry/bounding_box.h>
// #include <geode/geometry/nn_search.h>
// #include <geode/geometry/point.h>
// #include <geode/geometry/vector.h>

// #include <geode/mesh/builder/edged_curve_builder.h>
// #include <geode/mesh/builder/point_set_builder.h>
// #include <geode/mesh/builder/triangulated_surface_builder.h>
// #include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/polygonal_surface.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/corner.h>
// #include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/mixin/core/vertex_identifier.h>

#include <geode/geosciences/private/gocad_common.h>
// #include
// <geode/geosciences/representation/builder/structural_model_builder.h>
#include <geode/geosciences/mixin/core/fault.h>
#include <geode/geosciences/mixin/core/horizon.h>
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

    template <>
    struct hash< typename geode::Fault3D::FAULT_TYPE >
    {
    public:
        size_t operator()( const geode::Fault3D::FAULT_TYPE& x ) const
        {
            return static_cast< size_t >( x );
        }
    };

    template <>
    struct hash< typename geode::Horizon3D::HORIZON_TYPE >
    {
    public:
        size_t operator()( const geode::Horizon3D::HORIZON_TYPE& x ) const
        {
            return static_cast< size_t >( x );
        }
    };
} // namespace std

namespace
{
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
        }

        void write_file()
        {
            // TODO : need to check that full triangles;
            file_ << "GOCAD Model3d 1" << EOL;
            geode::detail::write_header( file_, {} );
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
            // TODO: create one tsurf for all surfaces not in collection ?
        }

        void write_tfaces()
        {
            for( const auto& boundary : model_.model_boundaries() )
            {
                for( const auto& item : model_.BRep::items( boundary ) )
                {
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
                for( const auto& item : model_.items( fault ) )
                {
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
                for( const auto& item : model_.items( horizon ) )
                {
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
            // TODO: what about all surfaces not in collection ?
            // What about surface that are into several groups : fault and
            // boundary?
        }

        void write_universe()
        {
            file_ << "REGION " << component_id_ << SPACE << "Universe" << EOL
                  << SPACE << SPACE;
            geode::index_t counter{ 0 };
            for( const auto& boundary : model_.model_boundaries() )
            {
                for( const auto& item : model_.BRep::items( boundary ) )
                {
                    file_ << "+" << components_.at( item.id() ) << SPACE
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
                    file_ << "+" << components_.at( surface.id() ) << SPACE
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

        void write_model_components()
        {
            write_tsurfs();
            write_tfaces();
            write_regions();
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

        void toto( const geode::Surface3D& surface,
            const geode::index_t current_offset,
            std::vector< std::array< geode::index_t, 2 > >& line_starts )
        {
            // todo several times to process all border edges
            DEBUG( model_.nb_internal_lines( surface ) );
            const auto& mesh = surface.mesh();
            geode::PolygonEdge first_on_border;
            for( const auto p : geode::Range{ mesh.nb_polygons() } )
            {
                for( const auto e : geode::Range{ 3 } )
                {
                    if( mesh.is_edge_on_border( { p, e } ) )
                    {
                        first_on_border = { p, e };
                        break;
                    }
                }
            }
            if( first_on_border.polygon_id == geode::NO_ID )
            {
                return;
            }
            {
                const auto v0 = mesh.polygon_vertex( first_on_border );
                const auto uid0 =
                    model_.unique_vertex( { surface.component_id(), v0 } );
                const auto v1 =
                    mesh.polygon_vertex( { first_on_border.polygon_id,
                        ( first_on_border.edge_id + 1 ) % 3 } );
                const auto corner_mcvs0 = model_.mesh_component_vertices(
                    uid0, geode::Corner3D::component_type_static() );
                if( !corner_mcvs0.empty() )
                {
                    line_starts.emplace_back( std::array< geode::index_t, 2 >{
                        v0 + current_offset, v1 + current_offset } );
                }
            }

            auto cur = mesh.next_on_border( first_on_border );
            while( cur != first_on_border )
            {
                const auto v0 = mesh.polygon_vertex( cur );
                const auto uid0 =
                    model_.unique_vertex( { surface.component_id(), v0 } );
                const auto v1 = mesh.polygon_vertex(
                    { cur.polygon_id, ( cur.edge_id + 1 ) % 3 } );
                const auto corner_mcvs0 = model_.mesh_component_vertices(
                    uid0, geode::Corner3D::component_type_static() );
                cur = mesh.next_on_border( cur );
                if( corner_mcvs0.empty() )
                {
                    continue;
                }
                line_starts.emplace_back( std::array< geode::index_t, 2 >{
                    v0 + current_offset, v1 + current_offset } );
            }
        }

        void find_boundary_corners_and_line_starts(
            const geode::ModelBoundary3D& surface_collection,
            std::vector< std::array< geode::index_t, 2 > >& line_starts )
        {
            geode::index_t current_offset{ OFFSET_START };
            for( const auto& surface :
                model_.BRep::items( surface_collection ) )
            {
                toto( surface, current_offset, line_starts );
                current_offset += surface.mesh().nb_vertices();
            }
        }

        template < typename CollectionType >
        void find_corners_and_line_starts(
            const CollectionType& surface_collection,
            std::vector< std::array< geode::index_t, 2 > >& line_starts )
        {
            geode::index_t current_offset{ OFFSET_START };
            for( const auto& surface : model_.items( surface_collection ) )
            {
                toto( surface, current_offset, line_starts );
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
                // TODO:: geol feature
                // file_ << "GEOLOGICAL_FEATURE " << boundary.name() << EOL;
                file_ << "GEOLOGICAL_TYPE boundary" << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.BRep::items( boundary ) )
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
                // TODO:: geol feature
                // file_ << "GEOLOGICAL_FEATURE " << boundary.name() << EOL;
                // file_ << "GEOLOGICAL_TYPE boundary" << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.items( fault ) )
                {
                    file_ << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;

                find_corners_and_line_starts( fault, line_starts );
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
                // TODO:: geol feature
                // file_ << "GEOLOGICAL_FEATURE " << horizon.name() << EOL;
                // file_ << "GEOLOGICAL_TYPE boundary" << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.items( horizon ) )
                {
                    file_ << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;

                find_corners_and_line_starts( horizon, line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                file_ << "END" << EOL;
            }
        }

    private:
        std::ofstream file_;
        const geode::StructuralModel& model_;
        std::unordered_map< geode::uuid, geode::index_t > components_;
        const std::unordered_map< geode::Fault3D::FAULT_TYPE, std::string >
            fault_map_ = { { geode::Fault3D::FAULT_TYPE::NO_TYPE, "fault" },
                { geode::Fault3D::FAULT_TYPE::NORMAL, "normal_fault" },
                { geode::Fault3D::FAULT_TYPE::REVERSE, "reverse_fault" },
                { geode::Fault3D::FAULT_TYPE::STRIKE_SLIP, "fault" },
                { geode::Fault3D::FAULT_TYPE::LISTRIC, "fault" },
                { geode::Fault3D::FAULT_TYPE::DECOLLEMENT, "fault" } };
        const std::unordered_map< geode::Horizon3D::HORIZON_TYPE, std::string >
            horizon_map_ = { { geode::Horizon3D::HORIZON_TYPE::NO_TYPE,
                                 "none" },
                { geode::Horizon3D::HORIZON_TYPE::CONFORMAL, "top" },
                { geode::Horizon3D::HORIZON_TYPE::NON_CONFORMAL,
                    "unconformity" } };
        geode::index_t component_id_{ OFFSET_START };
    }; // namespace
} // namespace

namespace geode
{
    namespace detail
    {
        void MLOutput::write() const
        {
            MLOutputImpl impl{ filename(), structural_model() };
            impl.write_file();
        }
    } // namespace detail
} // namespace geode
