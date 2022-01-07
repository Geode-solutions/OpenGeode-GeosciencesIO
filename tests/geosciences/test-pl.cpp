/*
 * Copyright (c) 2019 - 2022 Geode-solutions
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
#include <geode/basic/assert.h>
#include <geode/basic/attribute_manager.h>
#include <geode/basic/logger.h>
#include <geode/geometry/point.h>
#include <geode/mesh/builder/geode_edged_curve_builder.h>
#include <geode/mesh/core/geode_edged_curve.h>

#include <geode/geosciences/private/pl_output.h>

namespace
{
    void init_edged_curve_geometry( const geode::EdgedCurve3D& edged_curve,
        geode::EdgedCurveBuilder3D& builder )
    {
        builder.create_vertices( 4 );
        builder.set_point( 0, { { 0.5, 0.5, 0.5 } } );
        builder.set_point( 1, { { 1., 1., 1. } } );
        builder.set_point( 2, { { 7.5, 5.2, 6.3 } } );
        builder.set_point( 3, { { 8.7, 1.4, 4.7 } } );
        builder.create_edges( 3 );
        builder.set_edge_vertex( { 0, 0 }, 0 );
        builder.set_edge_vertex( { 0, 1 }, 1 );
        builder.set_edge_vertex( { 1, 0 }, 1 );
        builder.set_edge_vertex( { 1, 1 }, 2 );
        builder.set_edge_vertex( { 2, 0 }, 1 );
        builder.set_edge_vertex( { 2, 1 }, 3 );

        std::shared_ptr< geode::VariableAttribute< double > > attribut =
            edged_curve.vertex_attribute_manager()
                .find_or_create_attribute< geode::VariableAttribute, double >(
                    "attribute_double", -999. );
        attribut->set_value( 0, 0. );
        attribut->set_value( 1, 1. );
        attribut->set_value( 2, 2. );
        attribut->set_value( 3, 3. );
    }

} // namespace

int main()
{
    try
    {
        geode::detail::initialize_geosciences_io();

        auto edged_curve = geode::EdgedCurve3D::create(
            geode::OpenGeodeEdgedCurve3D::impl_name_static() );
        auto builder = geode::EdgedCurveBuilder3D::create( *edged_curve );
        init_edged_curve_geometry( *edged_curve, *builder );

        const auto output_file_native =
            absl::StrCat( "pl3d.", geode::detail::PLOutput::extension() );
        geode::save_edged_curve( *edged_curve, output_file_native );

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}