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

#include <geode/geosciences_io/mesh/internal/omf_edged_curve_input.hpp>

#include <geode/basic/attribute_manager.hpp>

#include <geode/mesh/builder/edged_curve_builder.hpp>
#include <geode/mesh/core/edged_curve.hpp>

#include <geode/geosciences_io/mesh/internal/omf_common.hpp>

namespace
{
    constexpr auto OMF_EDGED_CURVE_ELEMENT_CLASS = "LineSetElement";
} // namespace

namespace geode::internal
{
    std::unique_ptr< EdgedCurve3D > OMFEdgedCurveInput::read(
        const MeshImpl& impl )
    {
        auto curve = EdgedCurve3D::create( impl );
        auto builder = EdgedCurveBuilder3D::create( *curve );
        const OMFFile file{ filename() };
        const auto elements =
            file.elements( OMF_EDGED_CURVE_ELEMENT_CLASS, "" );
        builder->set_name( file.mesh_name( elements ) );
        for( const auto element_id : Indices{ elements } )
        {
            const auto& element = *elements[element_id];
            const auto& geometry = file.geometry( element );
            const auto points = file.read_points( geometry, "vertices" );
            const auto nb_points = static_cast< index_t >( points.size() );
            const auto segments =
                file.read_elements< 2 >( geometry, "segments", nb_points );
            const auto first_vertex = curve->nb_vertices();
            const auto first_edge = curve->nb_edges();
            builder->create_vertices( static_cast< index_t >( points.size() ) );
            for( const auto point : Indices{ points } )
            {
                builder->set_point( first_vertex + point, points[point] );
            }
            std::vector< index_t > segment_edges( segments.size(), NO_ID );
            for( const auto segment : Indices{ segments } )
            {
                const auto& vertices = segments[segment];
                if( vertices[0] != NO_ID )
                {
                    segment_edges[segment] =
                        builder->create_edge( first_vertex + vertices[0],
                            first_vertex + vertices[1] );
                }
            }
            file.read_data( element, "vertices", nb_points,
                curve->vertex_attribute_manager(),
                [first_vertex]( index_t value ) {
                    return first_vertex + value;
                } );
            file.read_data( element, "segments",
                static_cast< index_t >( segments.size() ),
                curve->edge_attribute_manager(),
                [&segment_edges]( index_t value ) {
                    return segment_edges[value];
                } );
            if( elements.size() > 1 )
            {
                set_omf_element_attribute( curve->vertex_attribute_manager(),
                    element_id, first_vertex, curve->nb_vertices() );
                set_omf_element_attribute( curve->edge_attribute_manager(),
                    element_id, first_edge, curve->nb_edges() );
            }
        }
        return curve;
    }

    Percentage OMFEdgedCurveInput::is_loadable() const
    {
        return OMFFile::is_loadable(
            filename(), OMF_EDGED_CURVE_ELEMENT_CLASS, "" );
    }
} // namespace geode::internal
