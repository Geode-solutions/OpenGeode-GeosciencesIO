/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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
        VSInputImpl( std::string_view filename, geode::PointSet3D& point_set )
            : file_{ geode::to_string( filename ), std::ios::binary },
              point_set_( point_set ),
              builder_( geode::PointSetBuilder< 3 >::create( point_set ) )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void read_file()
        {
            while( const auto vertex_set =
                       geode::internal::read_vs_points( file_ ) )
            {
                build_point_set( vertex_set.value() );
            }
        }

    private:
        void build_point_set( const geode::internal::VSetData& vertex_set )
        {
            if( vertex_set.header.name )
            {
                builder_->set_name( vertex_set.header.name.value() );
            }
            for( const auto& point : vertex_set.points )
            {
                builder_->create_point( point );
            }
            std::vector< geode::index_t > inverse_vertex_mapping(
                vertex_set.points.size() );
            absl::c_iota( inverse_vertex_mapping, 0 );
            geode::internal::create_attributes(
                vertex_set.vertices_properties_header,
                vertex_set.vertices_attribute_values,
                point_set_.vertex_attribute_manager(), vertex_set.points.size(),
                inverse_vertex_mapping );
        }

    private:
        std::ifstream file_;
        geode::PointSet3D& point_set_;
        std::unique_ptr< geode::PointSetBuilder3D > builder_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::unique_ptr< PointSet3D > VSInput::read( const MeshImpl& impl )
        {
            auto surface = PointSet3D::create( impl );
            VSInputImpl reader{ this->filename(), *surface };
            reader.read_file();
            return surface;
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
    } // namespace internal
} // namespace geode
