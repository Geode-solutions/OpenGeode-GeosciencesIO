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

#include <fstream>

#include <absl/strings/str_split.h>

#include <geode/basic/filename.h>

#include <geode/geometry/point.h>

#include <geode/mesh/builder/edged_curve_builder.h>
#include <geode/mesh/core/edged_curve.h>

namespace geode
{
    namespace detail
    {
        class WellInputImpl
        {
        public:
            WellInputImpl(
                absl::string_view filename, geode::EdgedCurve3D& curve )
                : file_{ geode::to_string( filename ) },
                  curve_( curve ),
                  builder_( geode::EdgedCurveBuilder3D::create( curve ) )
            {
                OPENGEODE_EXCEPTION(
                    file_.good(), "Error while opening file: ", filename );
                builder_->set_name(
                    geode::filename_without_extension( filename ) );
            }

            void read_file()
            {
                std::string line;
                while( std::getline( file_, line ) )
                {
                    builder_->create_point( read_coord( line ) );
                }
                for( const auto pt_id : Range{ curve_.nb_vertices() - 1 } )
                {
                    builder_->create_edge( pt_id, pt_id + 1 );
                }
            }

        private:
            geode::Point3D read_coord( absl::string_view line ) const
            {
                geode::Point3D coord;
                std::vector< absl::string_view > tokens =
                    absl::StrSplit( absl::StripAsciiWhitespace( line ),
                        absl::ByAnyChar( " 	" ), absl::SkipEmpty() );
                OPENGEODE_ASSERT( tokens.size() == 3,
                    "[WellInput::read_coord] Wrong number of tokens" );
                for( const auto i : geode::Range{ 3 } )
                {
                    double value;
                    const auto ok = absl::SimpleAtod( tokens[i], &value );
                    OPENGEODE_EXCEPTION( ok,
                        "[WellInput::read_coord] Error while reading "
                        "point coordinates, on axis ",
                        i, " with value '", tokens[i], "'" );
                    coord.set_value( i, value );
                }
                return coord;
            }

        private:
            std::ifstream file_;
            geode::EdgedCurve3D& curve_;
            std::unique_ptr< geode::EdgedCurveBuilder3D > builder_;
        };
    } // namespace detail
} // namespace geode
