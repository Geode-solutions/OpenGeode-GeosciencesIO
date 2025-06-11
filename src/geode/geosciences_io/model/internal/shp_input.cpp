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

#include <geode/geosciences_io/model/internal/shp_input.hpp>

#include <ogrsf_frmts.h>

#include <geode/basic/file.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/basic_objects/polygon.hpp>
#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/edged_curve_builder.hpp>
#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/builder/surface_mesh_builder.hpp>
#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/corner.hpp>
#include <geode/model/mixin/core/line.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/builder/section_builder.hpp>
#include <geode/model/representation/core/section.hpp>

namespace
{
    class SHPInputImpl
    {
    public:
        SHPInputImpl( geode::Section& section, std::string_view filename )
            : section_( section ),
              builder_{ section },
              gdal_data_{ GDALDataset::Open(
                  geode::to_string( filename ).c_str(), GDAL_OF_VECTOR ) }
        {
            OPENGEODE_EXCEPTION(
                gdal_data_, "[SHPInput] Failed to open file ", filename );
        }

        void read_file()
        {
            for( auto* layer : gdal_data_->GetLayers() )
            {
                if( layer->GetGeomType() == OGRwkbGeometryType::wkbPoint )
                {
                    create_corner( *layer );
                }
                else if( layer->GetGeomType()
                         == OGRwkbGeometryType::wkbLineString )
                {
                    create_line( *layer );
                }
                else if( layer->GetGeomType()
                         == OGRwkbGeometryType::wkbPolygon )
                {
                    create_surface( *layer );
                }
                else
                {
                    geode::Logger::warn( "[SHPInput] Unknown Layer type: ",
                        layer->GetGeomType() );
                }
            }
        }

    private:
        void create_corner( OGRLayer& layer )
        {
            const auto& id = builder_.add_corner();
            builder_.set_corner_name( id, layer.GetName() );
            const auto& corner = section_.corner( id );
            auto corner_builder = builder_.corner_mesh_builder( id );
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                OPENGEODE_EXCEPTION( geometry,
                    "[SHPInput::create_corner] Failed to retrieve geometry "
                    "data" );
                if( geometry->getGeometryType()
                    == OGRwkbGeometryType::wkbPoint )
                {
                    const auto* point = geometry->toPoint();
                    corner_builder->create_point(
                        geode::Point2D{ { point->getX(), point->getY() } } );
                    builder_.set_unique_vertex( { corner.component_id(), 0 },
                        builder_.create_unique_vertex() );
                }
                else
                {
                    geode::Logger::warn(
                        "[SHPInput::create_corner] Unknown geometry type: ",
                        geometry->getGeometryType() );
                }
            }
        }

        void create_line( OGRLayer& layer )
        {
            const auto& id = builder_.add_line();
            builder_.set_line_name( id, layer.GetName() );
            const auto& line = section_.line( id );
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                OPENGEODE_EXCEPTION( geometry,
                    "[SHPInput::create_line] Failed to retrieve geometry "
                    "data" );
                if( geometry->getGeometryType()
                    == OGRwkbGeometryType::wkbLineString )
                {
                    read_line( *geometry->toLineString(), line );
                }
                else if( geometry->getGeometryType()
                         == OGRwkbGeometryType::wkbMultiLineString )
                {
                    for( const auto* line_string :
                        *geometry->toMultiLineString() )
                    {
                        read_line( *line_string, line );
                    }
                }
                else
                {
                    geode::Logger::warn( "[SHPInput::create_line] "
                                         "Unknown geometry type: ",
                        geometry->getGeometryType() );
                }
            }
        }

        void read_line(
            const OGRLineString& line_string, const geode::Line2D& line )
        {
            const auto& curve = line.mesh();
            auto curve_builder = builder_.line_mesh_builder( line.id() );
            const auto start = read_points( line_string, line, *curve_builder );
            for( const auto p : geode::Range{ start, curve.nb_vertices() - 1 } )
            {
                curve_builder->create_edge( p, p + 1 );
            }
            if( line_string.get_IsClosed() )
            {
                curve_builder->create_edge( start, curve.nb_vertices() - 1 );
            }
        }

        template < typename Builder >
        geode::index_t read_points( const OGRLineString& line_string,
            const geode::Component2D& component,
            Builder& mesh_builder )
        {
            const auto nb_points = line_string.get_IsClosed()
                                       ? line_string.getNumPoints() - 1
                                       : line_string.getNumPoints();
            const auto start = mesh_builder.create_vertices( nb_points );
            const auto unique_start =
                builder_.create_unique_vertices( nb_points );
            for( const auto p : geode::Range{ nb_points } )
            {
                OGRPoint point;
                line_string.getPoint( p, &point );
                const auto vertex = start + p;
                mesh_builder.set_point(
                    vertex, geode::Point2D{ { point.getX(), point.getY() } } );
                const auto unique_vertex = unique_start + p;
                builder_.set_unique_vertex(
                    { component.component_id(), vertex }, unique_vertex );
            }
            return start;
        }

        void create_surface( OGRLayer& layer )
        {
            const auto& id = builder_.add_surface();
            builder_.set_surface_name( id, layer.GetName() );
            const auto& surface = section_.surface( id );
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                OPENGEODE_EXCEPTION( geometry,
                    "[SHPInput::create_surface] Failed to retrieve geometry "
                    "data" );
                if( geometry->getGeometryType()
                    == OGRwkbGeometryType::wkbPolygon )
                {
                    read_polygon( *geometry->toPolygon(), surface );
                }
                else if( geometry->getGeometryType()
                         == OGRwkbGeometryType::wkbMultiPolygon )
                {
                    for( const auto* polygon : *geometry->toMultiPolygon() )
                    {
                        read_polygon( *polygon, surface );
                    }
                }
                else
                {
                    geode::Logger::warn( "[SHPInput::create_surface] "
                                         "Unknown geometry type: ",
                        geometry->getGeometryType() );
                }
            }
        }

        void read_polygon(
            const OGRPolygon& polygon_string, const geode::Surface2D& surface )
        {
            const auto& mesh = surface.mesh();
            auto surface_builder =
                builder_.surface_mesh_builder( surface.id() );
            for( const auto* line_string : polygon_string )
            {
                const auto start =
                    read_points( *line_string, surface, *surface_builder );
                const auto nb_vertices = line_string->getNumPoints() - 1;
                if( nb_vertices < geode::NO_LID )
                {
                    absl::FixedArray< geode::index_t > vertices( nb_vertices );
                    absl::c_iota( vertices, start );
                    surface_builder->create_polygon( vertices );
                    continue;
                }
                std::vector< std::reference_wrapper< const geode::Point2D > >
                    points;
                points.reserve( nb_vertices );
                for( const auto p : geode::Range{ nb_vertices } )
                {
                    points.emplace_back( mesh.point( start + p ) );
                }
                geode::Polygon2D polygon{ std::move( points ) };
                for( auto& triangle : polygon.triangulate() )
                {
                    for( const auto v : geode::LRange{ 3 } )
                    {
                        triangle[v] = start + triangle[v];
                    }
                    surface_builder->create_polygon( triangle );
                }
            }
        }

    private:
        geode::Section& section_;
        geode::SectionBuilder builder_;
        GDALDatasetUniquePtr gdal_data_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        Section SHPInput::read()
        {
            Section section;
            SHPInputImpl impl{ section, filename() };
            impl.read_file();
            return section;
        }

        SectionInput::MissingFiles SHPInput::check_missing_files() const
        {
            const auto file_path = filepath_without_extension( filename() );
            const auto shx_file = absl::StrCat( file_path.string(), ".shx" );
            SectionInput::MissingFiles missing;
            if( !file_exists( shx_file ) )
            {
                missing.mandatory_files.push_back( shx_file );
            }
            return missing;
        }

    } // namespace internal
} // namespace geode