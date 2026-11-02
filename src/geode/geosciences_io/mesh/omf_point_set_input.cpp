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

#include <geode/geosciences_io/mesh/internal/omf_point_set_input.hpp>

#include <geode/basic/attribute_manager.hpp>

#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/core/point_set.hpp>

#include <geode/geosciences_io/mesh/internal/omf_common.hpp>

namespace
{
    constexpr auto OMF_POINT_SET_ELEMENT_CLASS = "PointSetElement";
} // namespace

namespace geode::internal
{
    std::unique_ptr< PointSet3D > OMFPointSetInput::read( const MeshImpl& impl )
    {
        auto point_set = PointSet3D::create( impl );
        auto builder = PointSetBuilder3D::create( *point_set );
        const OMFFile file{ filename() };
        const auto elements = file.elements( OMF_POINT_SET_ELEMENT_CLASS, "" );
        builder->set_name( file.mesh_name( elements ) );
        auto& manager = point_set->vertex_attribute_manager();
        for( const auto element_id : Indices{ elements } )
        {
            const auto& element = *elements[element_id];
            const auto points =
                file.read_points( file.geometry( element ), "vertices" );
            const auto first = point_set->nb_vertices();
            builder->create_vertices( static_cast< index_t >( points.size() ) );
            for( const auto point : Indices{ points } )
            {
                builder->set_point( first + point, points[point] );
            }
            file.read_data( element, "vertices",
                static_cast< index_t >( points.size() ), manager,
                [first]( index_t value ) {
                    return first + value;
                } );
            if( elements.size() > 1 )
            {
                set_omf_element_attribute(
                    manager, element_id, first, point_set->nb_vertices() );
            }
        }
        return point_set;
    }

    Percentage OMFPointSetInput::is_loadable() const
    {
        return OMFFile::is_loadable(
            filename(), OMF_POINT_SET_ELEMENT_CLASS, "" );
    }
} // namespace geode::internal
