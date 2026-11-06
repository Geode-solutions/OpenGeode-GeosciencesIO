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

#include <geode/geosciences_io/mesh/internal/pl_input.hpp>

#include <fstream>

#include <absl/algorithm/container.h>

#include <geode/basic/file.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/edged_curve_builder.hpp>
#include <geode/mesh/core/edged_curve.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

namespace
{
    class PLInputImpl
    {
    public:
        PLInputImpl( std::string_view filename, const geode::MeshImpl& impl )
            : file_{ geode::to_string( filename ), std::ios::binary },
              impl_( impl )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                file_.good(), nullptr, geode::OpenGeodeException::TYPE::data,
                "Error while opening file: ", filename );
        }

        std::vector< std::unique_ptr< geode::EdgedCurve3D > > read_file()
        {
            std::vector< std::unique_ptr< geode::EdgedCurve3D > > curves;
            while( const auto ecurve = geode::internal::read_ecurve( file_ ) )
            {
                curves.emplace_back( build_curve( ecurve.value() ) );
            }
            return curves;
        }

    private:
        std::unique_ptr< geode::EdgedCurve3D > build_curve(
            const geode::internal::ECurveData& ecurve ) const
        {
            auto curve = geode::EdgedCurve3D::create( impl_ );
            auto builder = geode::EdgedCurveBuilder3D::create( *curve );
            if( ecurve.header.name )
            {
                builder->set_name( ecurve.header.name.value() );
            }
            builder->create_vertices( ecurve.points.size() );
            for( const auto vertex : geode::Indices{ ecurve.points } )
            {
                builder->set_point( vertex, ecurve.points[vertex] );
            }
            for( const auto& edge : ecurve.edges )
            {
                builder->create_edge( edge[0], edge[1] );
            }
            std::vector< geode::index_t > inverse_vertex_mapping(
                ecurve.points.size() );
            absl::c_iota( inverse_vertex_mapping, 0 );
            geode::internal::create_attributes(
                ecurve.vertices_properties_header,
                ecurve.vertices_attribute_values,
                curve->vertex_attribute_manager(), ecurve.points.size(),
                inverse_vertex_mapping, 0 );
            return curve;
        }

    private:
        std::ifstream file_;
        const geode::MeshImpl& impl_;
    };
} // namespace

namespace geode::internal
{
    std::vector< std::unique_ptr< EdgedCurve3D > > PLInput::read(
        const MeshImpl& impl )
    {
        PLInputImpl reader{ this->filename(), impl };
        return reader.read_file();
    }

    Percentage PLInput::is_loadable() const
    {
        std::ifstream file{ to_string( this->filename() ) };
        if( goto_keyword_if_it_exists( file, "GOCAD PLine" ) )
        {
            return Percentage{ 1 };
        }
        return Percentage{ 0 };
    }
} // namespace geode::internal
