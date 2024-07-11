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

#include <geode/geosciences_io/mesh/internal/ts_output.h>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <geode/basic/attribute_manager.h>
#include <geode/basic/types.h>

#include <geode/basic/filename.h>
#include <geode/mesh/core/triangulated_surface.h>

#include <geode/geosciences_io/mesh/internal/gocad_common.h>

namespace
{
    class TSOutputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };
        TSOutputImpl( std::string_view filename,
            const geode::TriangulatedSurface3D& surface )
            : file_{ geode::to_string( filename ) }, surface_( surface )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void write_file()
        {
            geode::Logger::info( "[TSOutput::write] Writing ts file." );
            file_ << "GOCAD TSurf 1" << EOL;
            geode::internal::HeaderData header;
            header.name = geode::to_string( surface_.name() );
            geode::internal::write_header( file_, header );
            geode::internal::write_CRS( file_, {} );
            write_prop_header();
            write_tface();
            file_ << "END" << EOL;
        }

        void write_prop_header()
        {
            const auto names =
                surface_.vertex_attribute_manager().attribute_names();
            geode::internal::PropHeaderData prop_header;
            std::vector< geode::internal::PropClassHeaderData >
                header_properties_data;
            header_properties_data.reserve( names.size() );
            generic_att_.reserve( names.size() );

            for( const auto& name : names )
            {
                VRTX_KEYWORD = "PVRTX";
                const auto attribute =
                    surface_.vertex_attribute_manager().find_generic_attribute(
                        name );
                if( !attribute || !attribute->is_genericable()
                    || name == "points" )
                {
                    continue;
                }
                generic_att_.push_back( attribute );
                prop_header.names.emplace_back( geode::to_string( name ) );
                prop_header.prop_legal_ranges.push_back(
                    std::make_pair( "**none**", "**none**" ) );
                prop_header.no_data_values.push_back( -99999. );
                prop_header.property_classes.emplace_back(
                    geode::to_string( name ) );
                prop_header.kinds.push_back( "Real Number" );
                prop_header.property_subclass.push_back(
                    std::make_pair( "QUANTITY", "Float" ) );
                prop_header.esizes.push_back( 1 );
                prop_header.units.push_back( "unitless" );

                geode::internal::PropClassHeaderData tmp_prop_class_header;
                tmp_prop_class_header.name = geode::to_string( name );

                header_properties_data.push_back( tmp_prop_class_header );
            }
            if( !prop_header.empty() )
            {
                geode::internal::write_prop_header( file_, prop_header );
            }
            write_xyz_prop_class_header();
            for( const auto& property_data : header_properties_data )
            {
                write_property_class_header( file_, property_data );
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

        void write_tface()
        {
            file_ << "TFACE" << EOL;
            for( const auto v : geode::Range{ surface_.nb_vertices() } )
            {
                write_vrtx( v );
            }
            for( const auto triangle_id :
                geode::Range{ surface_.nb_polygons() } )
            {
                write_triangle( triangle_id );
            }
        }

        void write_vrtx( const geode::index_t vertex_id )
        {
            file_ << VRTX_KEYWORD << SPACE << vertex_id << SPACE
                  << surface_.point( vertex_id ).string();
            for( const auto& att : generic_att_ )
            {
                file_ << SPACE << att->generic_value( vertex_id );
            }
            file_ << EOL;
        }

        void write_triangle( const geode::index_t triangle_id )
        {
            const auto& vertices = surface_.polygon_vertices( triangle_id );
            file_ << "TRGL" << SPACE << vertices[0] << SPACE << vertices[1]
                  << SPACE << vertices[2];
            file_ << EOL;
        }

    private:
        std::ofstream file_;
        const geode::TriangulatedSurface3D& surface_;
        std::vector< std::shared_ptr< geode::AttributeBase > > generic_att_{};
        std::string VRTX_KEYWORD{ "VRTX" };
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::vector< std::string > TSOutput::write(
            const TriangulatedSurface3D& surface ) const
        {
            TSOutputImpl impl{ filename(), surface };
            impl.write_file();
            return { to_string( filename() ) };
        }
    } // namespace internal
} // namespace geode
