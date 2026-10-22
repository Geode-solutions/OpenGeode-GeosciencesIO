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

#include <geode/geosciences_io/mesh/internal/omf_polygonal_surface_input.hpp>

#include <cmath>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>

#include <geode/mesh/builder/polygonal_surface_builder.hpp>
#include <geode/mesh/core/polygonal_surface.hpp>

#include <geode/geosciences_io/mesh/internal/omf_common.hpp>

namespace
{
    constexpr auto OMF_POLYGONAL_SURFACE_ELEMENT_CLASS = "SurfaceElement";
    constexpr auto OMF_POLYGONAL_SURFACE_GEOMETRY_CLASS = "SurfaceGridGeometry";

    std::vector< double > omf_node_coordinates(
        absl::Span< const double > tensor )
    {
        std::vector< double > coordinates( tensor.size() + 1, 0 );
        for( const auto cell : geode::Indices{ tensor } )
        {
            coordinates[cell + 1] = coordinates[cell] + tensor[cell];
        }
        return coordinates;
    }

    /*
     * Grid nodes and cells are ordered with u varying fastest
     */
    class OMFSurfaceGridReader
    {
    public:
        OMFSurfaceGridReader( const geode::internal::OMFFile& file,
            const nlohmann::json& element,
            geode::PolygonalSurface3D& surface,
            geode::PolygonalSurfaceBuilder3D& builder )
            : file_( file ),
              element_( element ),
              geometry_( file.geometry( element ) ),
              surface_( surface ),
              builder_( builder ),
              u_( omf_node_coordinates( geode::internal::OMFFile::read_tensor(
                  geometry_, "tensor_u" ) ) ),
              v_( omf_node_coordinates( geode::internal::OMFFile::read_tensor(
                  geometry_, "tensor_v" ) ) )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                u_.size() * v_.size() < geode::NO_ID, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[OMFInput] Grid surface \"",
                geode::internal::OMFFile::name( element ),
                "\" has too many nodes" );
        }

        void read()
        {
            const auto first_vertex = surface_.nb_vertices();
            create_points();
            const auto nb_u = static_cast< geode::index_t >( u_.size() - 1 );
            const auto nb_v = static_cast< geode::index_t >( v_.size() - 1 );
            const auto first_polygon = surface_.nb_polygons();
            for( const auto v : geode::Range{ nb_v } )
            {
                for( const auto u : geode::Range{ nb_u } )
                {
                    const auto vertex = first_vertex + node( u, v );
                    builder_.create_polygon( { vertex, vertex + 1,
                        vertex + 1 + nb_u + 1, vertex + nb_u + 1 } );
                }
            }
            file_.read_data( element_, "vertices", nb_nodes(),
                surface_.vertex_attribute_manager(),
                [first_vertex]( geode::index_t value ) {
                    return first_vertex + value;
                } );
            file_.read_data( element_, "faces", nb_u * nb_v,
                surface_.polygon_attribute_manager(),
                [first_polygon]( geode::index_t value ) {
                    return first_polygon + value;
                } );
        }

    private:
        geode::index_t nb_nodes() const
        {
            return static_cast< geode::index_t >( u_.size() * v_.size() );
        }

        geode::index_t node( geode::index_t u, geode::index_t v ) const
        {
            return u + v * static_cast< geode::index_t >( u_.size() );
        }

        void create_points()
        {
            const auto axis_u =
                geode::internal::OMFFile::axis( geometry_, "axis_u", 0 );
            const auto axis_v =
                geode::internal::OMFFile::axis( geometry_, "axis_v", 1 );
            const auto axis_w = axis_u.cross( axis_v );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                axis_w.length() > geode::GLOBAL_EPSILON, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[OMFInput] Grid surface \"",
                geode::internal::OMFFile::name( element_ ),
                "\" has collinear axes" );
            const auto origin = file_.origin( geometry_ );
            auto offsets = file_.read_optional_scalars( geometry_, "offset_w" );
            if( !offsets.empty() && offsets.size() != nb_nodes() )
            {
                geode::Logger::warning( "[OMFInput] Grid surface \"",
                    geode::internal::OMFFile::name( element_ ), "\" has ",
                    offsets.size(), " offsets instead of ", nb_nodes(),
                    ": missing offsets are set to 0" );
            }
            offsets.resize( nb_nodes(), 0 );
            for( const auto v : geode::Indices{ v_ } )
            {
                for( const auto u : geode::Indices{ u_ } )
                {
                    const auto offset = offsets[node( u, v )];
                    builder_.create_point( geode::Point3D{
                        origin + axis_u * u_[u] + axis_v * v_[v]
                        + axis_w
                              * ( std::isfinite( offset ) ? offset : 0. ) } );
                }
            }
        }

    private:
        const geode::internal::OMFFile& file_;
        const nlohmann::json& element_;
        const nlohmann::json& geometry_;
        geode::PolygonalSurface3D& surface_;
        geode::PolygonalSurfaceBuilder3D& builder_;
        std::vector< double > u_;
        std::vector< double > v_;
    };
} // namespace

namespace geode::internal
{
    std::unique_ptr< PolygonalSurface3D > OMFPolygonalSurfaceInput::read(
        const MeshImpl& impl )
    {
        auto surface = PolygonalSurface3D::create( impl );
        auto builder = PolygonalSurfaceBuilder3D::create( *surface );
        const OMFFile file{ filename() };
        const auto elements =
            file.elements( OMF_POLYGONAL_SURFACE_ELEMENT_CLASS,
                OMF_POLYGONAL_SURFACE_GEOMETRY_CLASS );
        builder->set_name( file.mesh_name( elements ) );
        for( const auto element_id : Indices{ elements } )
        {
            const auto& element = *elements[element_id];
            const auto first_vertex = surface->nb_vertices();
            const auto first_polygon = surface->nb_polygons();
            OMFSurfaceGridReader{ file, element, *surface, *builder }.read();
            if( elements.size() > 1 )
            {
                set_omf_element_attribute( surface->vertex_attribute_manager(),
                    element_id, first_vertex, surface->nb_vertices() );
                set_omf_element_attribute( surface->polygon_attribute_manager(),
                    element_id, first_polygon, surface->nb_polygons() );
            }
            if( !element.value( "textures", nlohmann::json::array() ).empty() )
            {
                Logger::info( "[OMFInput] Textures of surface \"",
                    OMFFile::name( element ), "\" are not imported" );
            }
        }
        builder->compute_polygon_adjacencies();
        return surface;
    }

    Percentage OMFPolygonalSurfaceInput::is_loadable() const
    {
        return OMFFile::is_loadable( filename(),
            OMF_POLYGONAL_SURFACE_ELEMENT_CLASS,
            OMF_POLYGONAL_SURFACE_GEOMETRY_CLASS );
    }
} // namespace geode::internal
