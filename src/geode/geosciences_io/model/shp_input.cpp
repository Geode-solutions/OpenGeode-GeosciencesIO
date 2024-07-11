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

#include <geode/geosciences_io/model/internal/shp_input.h>

#include <ogrsf_frmts.h>

#include <geode/basic/file.h>
#include <geode/basic/filename.h>
#include <geode/basic/logger.h>

#include <geode/geometry/point.h>

#include <geode/mesh/builder/edged_curve_builder.h>
#include <geode/mesh/builder/point_set_builder.h>
#include <geode/mesh/builder/surface_mesh_builder.h>
#include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/surface_mesh.h>

#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/representation/builder/section_builder.h>
#include <geode/model/representation/core/section.h>

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
                switch( layer->GetGeomType() )
                {
                    case OGRwkbGeometryType::wkbPoint: {
                        create_corner( *layer );
                        break;
                    }
                    case OGRwkbGeometryType::wkbLineString: {
                        create_line( *layer );
                        break;
                    }
                    case OGRwkbGeometryType::wkbPolygon: {
                        create_surface( *layer );
                        break;
                    }
                    default: {
                        geode::Logger::warn( "[SHPInput] Unknown Layer type: ",
                            layer->GetGeomType() );
                        break;
                    }
                }
            }
        }

    private:
        void create_corner( OGRLayer& layer )
        {
            const auto& id = builder_.add_corner();
            builder_.set_corner_name( id, layer.GetName() );
            auto corner_builder = builder_.corner_mesh_builder( id );
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                OPENGEODE_EXCEPTION( geometry,
                    "[SHPInput::create_corner] Failed to retrieve geometry "
                    "data" );
                switch( geometry->getGeometryType() )
                {
                    case OGRwkbGeometryType::wkbPoint: {
                        const auto* point = geometry->toPoint();
                        corner_builder->create_point( geode::Point2D{
                            { point->getX(), point->getY() } } );
                        break;
                    }
                    default: {
                        geode::Logger::warn(
                            "[SHPInput::create_corner] Unknown geometry type: ",
                            geometry->getGeometryType() );
                        break;
                    }
                }
            }
        }

        void create_line( OGRLayer& layer )
        {
            const auto& id = builder_.add_line();
            builder_.set_line_name( id, layer.GetName() );
            auto curve_builder = builder_.line_mesh_builder( id );
            const auto& curve = section_.line( id ).mesh();
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                OPENGEODE_EXCEPTION( geometry,
                    "[SHPInput::create_line] Failed to retrieve geometry "
                    "data" );
                switch( geometry->getGeometryType() )
                {
                    case OGRwkbGeometryType::wkbLineString: {
                        read_line(
                            *geometry->toLineString(), curve, *curve_builder );
                        break;
                    }
                    case OGRwkbGeometryType::wkbMultiLineString: {
                        for( const auto* line : *geometry->toMultiLineString() )
                        {
                            read_line( *line, curve, *curve_builder );
                        }
                        break;
                    }
                    default: {
                        geode::Logger::warn( "[SHPInput::create_line] "
                                             "Unknown geometry type: ",
                            geometry->getGeometryType() );
                        break;
                    }
                }
            }
        }

        void read_line( const OGRLineString& line,
            const geode::EdgedCurve2D& curve,
            geode::EdgedCurveBuilder2D& curve_builder )
        {
            const auto start = read_points( line, curve_builder );
            for( const auto p : geode::Range{ start, curve.nb_vertices() - 1 } )
            {
                curve_builder.create_edge( p, p + 1 );
            }
            if( line.get_IsClosed() )
            {
                curve_builder.create_edge( start, curve.nb_vertices() - 1 );
            }
        }

        template < typename Builder >
        geode::index_t read_points(
            const OGRLineString& line, Builder& mesh_builder )
        {
            const auto nb_points = line.get_IsClosed() ? line.getNumPoints() - 1
                                                       : line.getNumPoints();
            const auto start = mesh_builder.create_vertices( nb_points );
            for( const auto p : geode::Range{ nb_points } )
            {
                OGRPoint point;
                line.getPoint( p, &point );
                mesh_builder.set_point( start + p,
                    geode::Point2D{ { point.getX(), point.getY() } } );
            }
            return start;
        }

        void create_surface( OGRLayer& layer )
        {
            const auto& id = builder_.add_surface();
            builder_.set_surface_name( id, layer.GetName() );
            auto surface_builder = builder_.surface_mesh_builder( id );
            for( const auto& feature : layer )
            {
                const auto* geometry = feature->GetGeometryRef();
                OPENGEODE_EXCEPTION( geometry,
                    "[SHPInput::create_surface] Failed to retrieve geometry "
                    "data" );
                switch( geometry->getGeometryType() )
                {
                    case OGRwkbGeometryType::wkbPolygon: {
                        read_polygon(
                            *geometry->toPolygon(), *surface_builder );
                        break;
                    }
                    case OGRwkbGeometryType::wkbMultiPolygon: {
                        for( const auto* polygon : *geometry->toMultiPolygon() )
                        {
                            read_polygon( *polygon, *surface_builder );
                        }
                        break;
                    }
                    default: {
                        geode::Logger::warn( "[SHPInput::create_surface] "
                                             "Unknown geometry type: ",
                            geometry->getGeometryType() );
                        break;
                    }
                }
            }
        }

        void read_polygon( const OGRPolygon& polygon,
            geode::SurfaceMeshBuilder2D& surface_builder )
        {
            for( const auto* line : polygon )
            {
                const auto start = read_points( *line, surface_builder );
                absl::FixedArray< geode::index_t > vertices(
                    line->getNumPoints() - 1 );
                absl::c_iota( vertices, start );
                surface_builder.create_polygon( vertices );
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
            const auto shx_file = absl::StrCat( file_path, ".shx" );
            SectionInput::MissingFiles missing;
            if( !file_exists( shx_file ) )
            {
                missing.mandatory_files.push_back( shx_file );
            }
            return missing;
        }

    } // namespace internal
} // namespace geode