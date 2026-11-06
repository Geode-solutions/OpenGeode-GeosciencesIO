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

#include <geode/geosciences_io/mesh/internal/vs_input.hpp>

#include <fstream>

#include <geode/basic/file.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/core/point_set.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

namespace
{
    class VSInputImpl
    {
    public:
        VSInputImpl( std::string_view filename, const geode::MeshImpl& impl )
            : file_{ geode::to_string( filename ), std::ios::binary },
              impl_( impl )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                file_.good(), nullptr, geode::OpenGeodeException::TYPE::data,
                "Error while opening file: ", filename );
        }

        std::vector< std::unique_ptr< geode::PointSet3D > > read_file()
        {
            std::vector< std::unique_ptr< geode::PointSet3D > > point_sets;
            while( const auto vertex_set =
                       geode::internal::read_vs_points( file_ ) )
            {
                point_sets.emplace_back(
                    build_point_set( vertex_set.value() ) );
            }
            return point_sets;
        }

    private:
        std::unique_ptr< geode::PointSet3D > build_point_set(
            const geode::internal::VSetData& vertex_set ) const
        {
            auto point_set = geode::PointSet3D::create( impl_ );
            auto builder = geode::PointSetBuilder3D::create( *point_set );
            if( vertex_set.header.name )
            {
                builder->set_name( vertex_set.header.name.value() );
            }
            builder->create_vertices( vertex_set.points.size() );
            for( const auto vertex : geode::Indices{ vertex_set.points } )
            {
                builder->set_point( vertex, vertex_set.points[vertex] );
            }
            std::vector< geode::index_t > inverse_vertex_mapping(
                vertex_set.points.size() );
            absl::c_iota( inverse_vertex_mapping, 0 );
            geode::internal::create_attributes(
                vertex_set.vertices_properties_header,
                vertex_set.vertices_attribute_values,
                point_set->vertex_attribute_manager(), vertex_set.points.size(),
                inverse_vertex_mapping, 0 );
            return point_set;
        }

    private:
        std::ifstream file_;
        const geode::MeshImpl& impl_;
    };
} // namespace

namespace geode::internal
{
    std::vector< std::unique_ptr< PointSet3D > > VSInput::read(
        const MeshImpl& impl )
    {
        VSInputImpl reader{ this->filename(), impl };
        return reader.read_file();
    }

    Percentage VSInput::is_loadable() const
    {
        std::ifstream file{ to_string( this->filename() ) };
        if( goto_keyword_if_it_exists( file, "GOCAD VSet" ) )
        {
            return Percentage{ 1 };
        }
        return Percentage{ 0 };
    }
} // namespace geode::internal
