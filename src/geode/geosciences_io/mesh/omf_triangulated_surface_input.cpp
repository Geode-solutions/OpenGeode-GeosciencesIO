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

#include <geode/geosciences_io/mesh/internal/omf_triangulated_surface_input.hpp>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>

#include <geode/mesh/builder/triangulated_surface_builder.hpp>
#include <geode/mesh/core/triangulated_surface.hpp>

#include <geode/geosciences_io/mesh/internal/omf_common.hpp>

namespace
{
    constexpr auto OMF_TRIANGULATED_SURFACE_ELEMENT_CLASS = "SurfaceElement";
    constexpr auto OMF_TRIANGULATED_SURFACE_GEOMETRY_CLASS = "SurfaceGeometry";
} // namespace

namespace geode::internal
{
    std::unique_ptr< TriangulatedSurface3D > OMFTriangulatedSurfaceInput::read(
        const MeshImpl& impl )
    {
        auto surface = TriangulatedSurface3D::create( impl );
        auto builder = TriangulatedSurfaceBuilder3D::create( *surface );
        const OMFFile file{ filename() };
        const auto elements =
            file.elements( OMF_TRIANGULATED_SURFACE_ELEMENT_CLASS,
                OMF_TRIANGULATED_SURFACE_GEOMETRY_CLASS );
        builder->set_name( file.mesh_name( elements ) );
        for( const auto element_id : Indices{ elements } )
        {
            const auto& element = *elements[element_id];
            const auto& geometry = file.geometry( element );
            const auto points = file.read_points( geometry, "vertices" );
            const auto nb_points = static_cast< index_t >( points.size() );
            const auto triangles =
                file.read_elements< 3 >( geometry, "triangles", nb_points );
            const auto first_vertex = surface->nb_vertices();
            const auto first_polygon = surface->nb_polygons();
            for( const auto& point : points )
            {
                builder->create_point( point );
            }
            std::vector< index_t > face_polygons( triangles.size(), NO_ID );
            for( const auto triangle : Indices{ triangles } )
            {
                auto vertices = triangles[triangle];
                if( vertices[0] == NO_ID )
                {
                    continue;
                }
                for( auto& vertex : vertices )
                {
                    vertex += first_vertex;
                }
                face_polygons[triangle] = builder->create_triangle( vertices );
            }
            file.read_data( element, "vertices", nb_points,
                surface->vertex_attribute_manager(),
                [first_vertex]( index_t value ) {
                    return first_vertex + value;
                } );
            file.read_data( element, "faces",
                static_cast< index_t >( triangles.size() ),
                surface->polygon_attribute_manager(),
                [&face_polygons]( index_t value ) {
                    return face_polygons[value];
                } );
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

    Percentage OMFTriangulatedSurfaceInput::is_loadable() const
    {
        return OMFFile::is_loadable( filename(),
            OMF_TRIANGULATED_SURFACE_ELEMENT_CLASS,
            OMF_TRIANGULATED_SURFACE_GEOMETRY_CLASS );
    }
} // namespace geode::internal
