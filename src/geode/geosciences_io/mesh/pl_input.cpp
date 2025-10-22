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

#include <geode/geosciences_io/mesh/internal/pl_input.hpp>

#include <fstream>

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
        PLInputImpl( std::string_view filename, geode::EdgedCurve3D& curve )
            : file_{ geode::to_string( filename ),std::ios::binary },
              curve_( curve ),
              builder_( geode::EdgedCurveBuilder< 3 >::create( curve ) )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void read_file()
        {
            while( const auto ecurve = geode::internal::read_ecurve( file_ ) )
            {
                build_curve( ecurve.value() );
            }
        }

    private:
        void build_curve( const geode::internal::ECurveData& ecurve )
        {
            const auto offset = curve_.nb_vertices();
            if( ecurve.header.name )
            {
                builder_->set_name( ecurve.header.name.value() );
            }
            for( const auto& point : ecurve.points )
            {
                builder_->create_point( point );
            }
            for( auto edge : ecurve.edges )
            {
                for( auto& vertex : edge )
                {
                    vertex += offset;
                }
                builder_->create_edge( edge[0], edge[1] );
            }
        }

    private:
        std::ifstream file_;
        geode::EdgedCurve3D& curve_;
        std::unique_ptr< geode::EdgedCurveBuilder3D > builder_;
    };

} // namespace
namespace geode
{
    namespace internal
    {
        std::unique_ptr< EdgedCurve3D > PLInput::read( const MeshImpl& impl )
        {
            auto curve = EdgedCurve3D::create( impl );
            PLInputImpl reader{ this->filename(), *curve };
            reader.read_file();
            return curve;
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
    } // namespace internal
} // namespace geode
