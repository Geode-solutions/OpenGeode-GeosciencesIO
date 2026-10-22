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

#include <geode/io/image/detail/gdal_file.hpp>

namespace
{
    class SHPInputImpl : public geode::detail::GDALFile
    {
    public:
        SHPInputImpl( geode::Section& section, std::string_view filename )
            : geode::detail::GDALFile{ filename },
              section_( section ),
              builder_{ section }
        {
        }

        void read_file()
        {
            for( auto* layer : dataset().GetLayers() )
            {
                const auto type = wkbFlatten( layer->GetGeomType() );
                if( type == OGRwkbGeometryType::wkbPoint
                    || type == OGRwkbGeometryType::wkbMultiPoint )
                {
                    create_corner( *layer );
                }
                else if( type == OGRwkbGeometryType::wkbLineString
                         || type == OGRwkbGeometryType::wkbMultiLineString )
                {
                    create_line( *layer );
                }
                else if( type == OGRwkbGeometryType::wkbPolygon
                         || type == OGRwkbGeometryType::wkbMultiPolygon )
                {
                    create_surface( *layer );
                }
                else
                {
                    geode::Logger::warning( "[SHPInput] Unknown Layer type: ",
                        layer->GetGeomType() );
                }
            }
        }

    private:
        void create_corner( OGRLayer& layer )
        {
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                if( geometry == nullptr )
                {
                    geode::Logger::warning(
                        "[SHPInput::create_corner] Skipping feature "
                        "without geometry" );
                    continue;
                }
                const auto type = wkbFlatten( geometry->getGeometryType() );
                if( type == OGRwkbGeometryType::wkbPoint )
                {
                    read_point( *geometry->toPoint(), layer.GetName() );
                }
                else if( type == OGRwkbGeometryType::wkbMultiPoint )
                {
                    for( const auto* point : *geometry->toMultiPoint() )
                    {
                        read_point( *point, layer.GetName() );
                    }
                }
                else
                {
                    geode::Logger::warning(
                        "[SHPInput::create_corner] Unknown geometry type: ",
                        geometry->getGeometryType() );
                }
            }
        }

        void read_point( const OGRPoint& point, std::string_view name )
        {
            const auto& id = builder_.add_corner();
            const auto& corner = section_.corner( id );
            builder_.set_corner_name( corner, name );
            builder_.corner_mesh_builder( corner )->create_point(
                geode::Point2D{ { point.getX(), point.getY() } } );
            builder_.set_unique_vertex(
                { corner.component_id(), 0 }, builder_.create_unique_vertex() );
        }

        void create_line( OGRLayer& layer )
        {
            const auto& id = builder_.add_line();
            const auto& line = section_.line( id );
            builder_.set_line_name( line, layer.GetName() );
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                if( geometry == nullptr )
                {
                    geode::Logger::warning(
                        "[SHPInput::create_line] Skipping feature "
                        "without geometry" );
                    continue;
                }
                const auto type = wkbFlatten( geometry->getGeometryType() );
                if( type == OGRwkbGeometryType::wkbLineString )
                {
                    read_line( *geometry->toLineString(), line );
                }
                else if( type == OGRwkbGeometryType::wkbMultiLineString )
                {
                    for( const auto* line_string :
                        *geometry->toMultiLineString() )
                    {
                        read_line( *line_string, line );
                    }
                }
                else
                {
                    geode::Logger::warning( "[SHPInput::create_line] "
                                            "Unknown geometry type: ",
                        geometry->getGeometryType() );
                }
            }
        }

        void read_line(
            const OGRLineString& line_string, const geode::Line2D& line )
        {
            const auto nb_points = nb_line_points( line_string );
            if( nb_points < 2 )
            {
                geode::Logger::warning( "[SHPInput::read_line] Skipping a line "
                                        "part with less than 2 points" );
                return;
            }
            auto curve_builder = builder_.line_mesh_builder( line );
            const auto start = read_points( line_string, line, *curve_builder );
            for( const auto p : geode::Range{ nb_points - 1 } )
            {
                curve_builder->create_edge( start + p, start + p + 1 );
            }
            if( line_string.get_IsClosed() && nb_points > 2 )
            {
                curve_builder->create_edge( start, start + nb_points - 1 );
            }
        }

        static geode::index_t nb_line_points( const OGRLineString& line_string )
        {
            const auto nb_points =
                static_cast< geode::index_t >( line_string.getNumPoints() );
            if( nb_points > 1 && line_string.get_IsClosed() )
            {
                return nb_points - 1;
            }
            return nb_points;
        }

        template < typename Builder >
        geode::index_t read_points( const OGRLineString& line_string,
            const geode::Component2D& component,
            Builder& mesh_builder )
        {
            const auto nb_points = nb_line_points( line_string );
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
            const auto& surface = section_.surface( id );
            builder_.set_surface_name( surface, layer.GetName() );
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                if( geometry == nullptr )
                {
                    geode::Logger::warning(
                        "[SHPInput::create_surface] Skipping feature "
                        "without geometry" );
                    continue;
                }
                const auto type = wkbFlatten( geometry->getGeometryType() );
                if( type == OGRwkbGeometryType::wkbPolygon )
                {
                    read_polygon( *geometry->toPolygon(), surface );
                }
                else if( type == OGRwkbGeometryType::wkbMultiPolygon )
                {
                    for( const auto* polygon : *geometry->toMultiPolygon() )
                    {
                        read_polygon( *polygon, surface );
                    }
                }
                else
                {
                    geode::Logger::warning( "[SHPInput::create_surface] "
                                            "Unknown geometry type: ",
                        geometry->getGeometryType() );
                }
            }
            if( nb_ignored_holes_ > 0 )
            {
                geode::Logger::warning( "[SHPInput::create_surface] Polygon "
                                        "holes are not supported: ",
                    nb_ignored_holes_, " interior ring(s) of layer ",
                    layer.GetName(), " ignored, holes are filled" );
                nb_ignored_holes_ = 0;
            }
        }

        void read_polygon(
            const OGRPolygon& polygon_string, const geode::Surface2D& surface )
        {
            const auto* exterior = polygon_string.getExteriorRing();
            if( exterior == nullptr )
            {
                return;
            }
            nb_ignored_holes_ += polygon_string.getNumInteriorRings();
            const auto nb_vertices = nb_line_points( *exterior );
            if( nb_vertices < 3 )
            {
                geode::Logger::warning(
                    "[SHPInput::read_polygon] Skipping a "
                    "polygon ring with less than 3 points" );
                return;
            }
            const auto& mesh = surface.mesh();
            auto surface_builder = builder_.surface_mesh_builder( surface );
            const auto start =
                read_points( *exterior, surface, *surface_builder );
            if( nb_vertices < geode::NO_LID )
            {
                absl::FixedArray< geode::index_t > vertices( nb_vertices );
                absl::c_iota( vertices, start );
                surface_builder->create_polygon( vertices );
                return;
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

    private:
        geode::Section& section_;
        geode::SectionBuilder builder_;
        geode::index_t nb_ignored_holes_{ 0 };
    };
} // namespace

namespace geode::internal
{
    Section SHPInput::read()
    {
        Section section;
        SHPInputImpl impl{ section, filename() };
        impl.read_file();
        return section;
    }

    auto SHPInput::additional_files() const -> AdditionalFiles
    {
        const auto file_path = filepath_without_extension( filename() );
        const auto shx_file = absl::StrCat( file_path.string(), ".shx" );
        AdditionalFiles missing;
        missing.mandatory_files.emplace_back(
            shx_file, file_exists( shx_file ) );
        return missing;
    }

    Percentage SHPInput::is_loadable() const
    {
        detail::GDALFile reader{ filename() };
        if( reader.dataset().GetLayerCount() == 0 )
        {
            return Percentage{ 0 };
        }
        return Percentage{ 1 };
    }
} // namespace geode::internal