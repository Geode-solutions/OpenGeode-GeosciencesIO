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

#pragma once

#include <fstream>

#include <geode/basic/filename.hpp>
#include <geode/basic/string.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/edged_curve_builder.hpp>
#include <geode/mesh/core/edged_curve.hpp>

namespace geode
{
    namespace internal
    {
        class WellInputImpl
        {
        public:
            WellInputImpl( std::string_view filename, EdgedCurve3D& curve )
                : file_{ to_string( filename ) },
                  curve_( curve ),
                  builder_( EdgedCurveBuilder3D::create( curve ) )
            {
                OPENGEODE_EXCEPTION(
                    file_.good(), "Error while opening file: ", filename );
                builder_->set_name(
                    filename_without_extension( filename ).string() );
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
            geode::Point3D read_coord( std::string_view line ) const
            {
                const auto tokens = string_split( line );
                OPENGEODE_ASSERT( tokens.size() == 3,
                    "[WellInput::read_coord] Wrong number of tokens" );
                return Point3D{ { string_to_double( tokens[0] ),
                    string_to_double( tokens[1] ),
                    string_to_double( tokens[2] ) } };
            }

        private:
            std::ifstream file_;
            EdgedCurve3D& curve_;
            std::unique_ptr< EdgedCurveBuilder3D > builder_;
        };
    } // namespace internal
} // namespace geode
