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

#include <geode/geosciences_io/mesh/internal/vs_output.hpp>

#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/types.hpp>

#include <geode/basic/filename.hpp>
#include <geode/mesh/core/point_set.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

namespace
{
    class VSOutputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };
        VSOutputImpl(
            std::string_view filename, const geode::PointSet3D& pointset )
            : file_{ geode::to_string( filename ) },
              pointset_( pointset ),
              properties_( pointset.vertex_attribute_manager() )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                file_.good(), nullptr, geode::OpenGeodeException::TYPE::data,
                "Error while opening file: ", filename );
            file_.precision( std::numeric_limits< double >::max_digits10 );
        }

        void write_file()
        {
            geode::Logger::info( "[VSOutput::write] Writing vs file." );
            file_ << "GOCAD VSet 1" << EOL;
            geode::internal::HeaderData header;
            if( const auto name = pointset_.name() )
            {
                header.name = name.value();
            }
            geode::internal::write_header( file_, header );
            geode::internal::CRSData crs_data;
            geode::internal::write_CRS( file_, crs_data );
            write_prop_header();
            write_vset();
            file_ << "END" << EOL;
        }

        void write_prop_header()
        {
            if( !properties_.empty() )
            {
                VRTX_KEYWORD = "PVRTX";
                geode::internal::write_prop_header(
                    file_, properties_.header() );
            }
            write_xyz_prop_class_header();
            for( const auto& class_header : properties_.class_headers() )
            {
                write_property_class_header( file_, class_header );
            }
        }

        void write_xyz_prop_class_header()
        {
            geode::internal::PropClassHeaderData x_prop_header;
            x_prop_header.name = "X";
            x_prop_header.kind = "X";
            x_prop_header.unit = "m";
            x_prop_header.is_z = false;
            write_property_class_header( file_, x_prop_header );

            geode::internal::PropClassHeaderData y_prop_header;
            y_prop_header.name = "Y";
            y_prop_header.kind = "Y";
            y_prop_header.unit = "m";
            y_prop_header.is_z = false;
            write_property_class_header( file_, y_prop_header );

            geode::internal::PropClassHeaderData z_prop_header;
            z_prop_header.name = "Z";
            z_prop_header.kind = "Z";
            z_prop_header.unit = "m";
            z_prop_header.is_z = true;
            write_property_class_header( file_, z_prop_header );
        }

        void write_vset()
        {
            for( const auto v : geode::Range{ pointset_.nb_vertices() } )
            {
                write_vrtx( v );
            }
        }

        void write_vrtx( const geode::index_t vertex_id )
        {
            file_ << VRTX_KEYWORD << SPACE << vertex_id << SPACE
                  << pointset_.point( vertex_id ).string();
            properties_.write_values( file_, vertex_id );
            file_ << EOL;
        }

    private:
        std::ofstream file_;
        const geode::PointSet3D& pointset_;
        geode::internal::GocadVertexProperties properties_;
        std::string VRTX_KEYWORD{ "VRTX" };
    };
} // namespace

namespace geode::internal
{
    std::vector< std::string > VSOutput::write(
        const PointSet3D& pointset ) const
    {
        VSOutputImpl impl{ filename(), pointset };
        impl.write_file();
        return { to_string( filename() ) };
    }
} // namespace geode::internal
