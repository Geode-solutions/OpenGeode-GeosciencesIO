/*
 * Copyright (c) 2019 - 2024 Geode-solutions
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

#include <geode/geosciences_io/mesh/internal/wl_input.hpp>

#include <fstream>

#include <geode/basic/file.hpp>
#include <geode/basic/string.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/edged_curve_builder.hpp>
#include <geode/mesh/core/edged_curve.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

namespace
{
    class WLInputImpl
    {
    public:
        WLInputImpl( std::string_view filename, geode::EdgedCurve3D& curve )
            : file_{ geode::to_string( filename ) },
              curve_( curve ),
              builder_( geode::EdgedCurveBuilder3D::create( curve ) )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void read_file()
        {
            if( !geode::goto_keyword_if_it_exists( file_, "GOCAD Well" ) )
            {
                throw geode::OpenGeodeException{
                    "[WLInput] Cannot find Well in the file"
                };
            }
            const auto header = geode::internal::read_header( file_ );
            builder_->set_name( header.name );
            crs_ = geode::internal::read_CRS( file_ );
            const auto ref = read_ref();
            builder_->create_point( ref );
            read_paths( ref );
        }

    private:
        geode::Point3D read_ref()
        {
            const auto line = geode::goto_keyword( file_, "WREF" );
            auto ref = read_coord( line, 1 );
            ref.set_value(
                2, ref.value( 2 ) * ( crs_.z_sign_positive ? 1. : -1. ) );
            return ref;
        }

        geode::Point3D read_coord(
            std::string_view line, geode::index_t offset ) const
        {
            const auto tokens = geode::string_split( line );
            OPENGEODE_ASSERT( tokens.size() == 3 + offset,
                "[WLInput::read_coord] Wrong number of tokens" );
            return geode::Point3D{ { geode::string_to_double( tokens[offset] ),
                geode::string_to_double( tokens[1 + offset] ),
                geode::string_to_double( tokens[2 + offset] ) } };
        }

        void read_paths( const geode::Point3D& ref )
        {
            auto line = geode::goto_keyword( file_, "PATH" );
            while( std::getline( file_, line ) )
            {
                if( !geode::string_starts_with( line, "PATH" ) )
                {
                    return;
                }
                auto translation = read_coord( line, 2 );
                geode::Point3D point;
                point.set_value( 0, translation.value( 1 ) + ref.value( 0 ) );
                point.set_value( 1, translation.value( 2 ) + ref.value( 1 ) );
                point.set_value( 2, translation.value( 0 )
                                        * ( crs_.z_sign_positive ? 1. : -1. ) );
                const auto id = builder_->create_point( point );
                builder_->create_edge( id - 1, id );
            }
        }

    private:
        std::ifstream file_;
        geode::EdgedCurve3D& curve_;
        std::unique_ptr< geode::EdgedCurveBuilder3D > builder_;
        geode::internal::CRSData crs_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::unique_ptr< EdgedCurve3D > WLInput::read( const MeshImpl& impl )
        {
            auto well = EdgedCurve3D::create( impl );
            WLInputImpl reader{ this->filename(), *well };
            reader.read_file();
            return well;
        }
    } // namespace internal
} // namespace geode
