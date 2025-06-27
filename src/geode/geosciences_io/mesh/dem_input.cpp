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

#include <geode/geosciences_io/mesh/internal/dem_input.hpp>

#include <gdal.h>
#include <gdal_priv.h>
#include <gdal_utils.h>

#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/builder/polygonal_surface_builder.hpp>
#include <geode/mesh/core/polygonal_surface.hpp>

namespace
{
    class DEMInputImpl
    {
    public:
        DEMInputImpl(
            geode::PolygonalSurface3D& surface, std::string_view filename )
            : surface_{ surface },
              builder_{ geode::PolygonalSurfaceBuilder3D::create( surface ) },
              gdal_data_{ GDALDataset::Open(
                  geode::to_string( filename ).c_str(), GDAL_OF_READONLY ) }
        {
            OPENGEODE_EXCEPTION(
                gdal_data_, "[DEMInput] Failed to open file ", filename );
        }

        void read_file()
        {
            read_metadata();
            const auto vertices = read_vertices();
            create_polygons( vertices );
        }

    private:
        void create_polygons( absl::Span< const geode::index_t > vertices )
        {
            for( const auto i : geode::Range{ height_ - 1 } )
            {
                const auto i_offset = i * width_;
                for( const auto j : geode::Range{ width_ - 1 } )
                {
                    const auto j_offset = j + width_;
                    std::array polygon{ i_offset + j, i_offset + j + 1,
                        i_offset + j_offset, i_offset + j_offset + 1 };
                    bool found{ true };
                    for( auto& vertex : polygon )
                    {
                        const auto vertex_id = vertices[vertex];
                        if( vertex_id == geode::NO_ID )
                        {
                            found = false;
                            break;
                        }
                        vertex = vertex_id;
                    }
                    if( !found )
                    {
                        continue;
                    }
                    builder_->create_polygon( polygon );
                }
            }
        }

        absl::FixedArray< geode::index_t > read_vertices()
        {
            const auto nb_vertices = width_ * height_;
            const auto nb_bands = gdal_data_->GetRasterCount();
            OPENGEODE_EXCEPTION( nb_bands > 0, "[DEMInput] No bands found" );
            absl::FixedArray< float > elevation( nb_vertices );
            absl::FixedArray< geode::index_t > vertices(
                nb_vertices, geode::NO_ID );
            const auto band = gdal_data_->GetRasterBand( 1 );
            const auto status = band->RasterIO( GF_Read, 0, 0, width_, height_,
                elevation.data(), width_, height_, GDT_Float32, 0, 0 );
            const auto no_data_value = band->GetNoDataValue();
            OPENGEODE_EXCEPTION(
                status == CE_None, "[DEMInput] Failed to read elevation" );
            for( const auto i : geode::Range{ height_ } )
            {
                const auto i_contribution = y_direction_ * i;
                for( const auto j : geode::Range{ width_ } )
                {
                    const auto vertex = i * width_ + j;
                    const auto current_elevation = elevation[vertex];
                    if( current_elevation == no_data_value )
                    {
                        continue;
                    }
                    const auto j_contribution = x_direction_ * j;
                    const auto point =
                        origin_ + i_contribution + j_contribution;
                    vertices[vertex] = builder_->create_point(
                        geode::Point3D{ { point.value( 0 ), point.value( 1 ),
                            current_elevation } } );
                }
            }
            return vertices;
        }

        void read_metadata()
        {
            std::array< double, 6 > geo_transform;
            const auto status =
                gdal_data_->GetGeoTransform( geo_transform.data() );
            OPENGEODE_EXCEPTION(
                status == CE_None, "[DEMInput] Failed to read geotransform" );
            for( const auto axis : geode::LRange{ 2 } )
            {
                origin_.set_value( axis, geo_transform[3 * axis] );
                x_direction_.set_value( axis, geo_transform[3 * axis + 1] );
                y_direction_.set_value( axis, geo_transform[3 * axis + 2] );
            }
            width_ = gdal_data_->GetRasterXSize();
            height_ = gdal_data_->GetRasterYSize();
        }

    private:
        const geode::PolygonalSurface3D& surface_;
        std::unique_ptr< geode::PolygonalSurfaceBuilder3D > builder_;
        GDALDatasetUniquePtr gdal_data_;
        geode::Point2D origin_{ { 0., 0. } };
        geode::Vector2D x_direction_{ { 1., 0. } };
        geode::Vector2D y_direction_{ { 0., 1. } };
        geode::index_t width_{ 0 };
        geode::index_t height_{ 0 };
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::unique_ptr< PolygonalSurface3D > DEMInput::read(
            const MeshImpl& impl )
        {
            auto surface = PolygonalSurface3D::create( impl );
            DEMInputImpl reader{ *surface, this->filename() };
            reader.read_file();
            return surface;
        }
    } // namespace internal
} // namespace geode