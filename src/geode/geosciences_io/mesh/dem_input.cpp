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

#include <geode/geometry/coordinate_system.hpp>
#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/builder/polygonal_surface_builder.hpp>
#include <geode/mesh/core/polygonal_surface.hpp>

#include <geode/io/image/detail/gdal_file.hpp>

namespace
{
    class DEMInputImpl : public geode::detail::GDALFile
    {
    public:
        DEMInputImpl(
            geode::PolygonalSurface3D& surface, std::string_view filename )
            : geode::detail::GDALFile{ filename },
              builder_{ geode::PolygonalSurfaceBuilder3D::create( surface ) }
        {
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
            const auto nb_bands = dataset().GetRasterCount();
            OPENGEODE_EXCEPTION( nb_bands > 0, "[DEMInput] No bands found" );
            absl::FixedArray< float > elevation( nb_vertices );
            absl::FixedArray< geode::index_t > vertices(
                nb_vertices, geode::NO_ID );
            const auto band = dataset().GetRasterBand( 1 );
            const auto status = band->RasterIO( GF_Read, 0, 0, width_, height_,
                elevation.data(), width_, height_, GDT_Float32, 0, 0 );
            const auto no_data_value = band->GetNoDataValue();
            OPENGEODE_EXCEPTION(
                status == CE_None, "[DEMInput] Failed to read elevation" );
            for( const auto i : geode::Range{ height_ } )
            {
                const auto i_contribution =
                    coordinate_system_.direction( 1 ) * i;
                for( const auto j : geode::Range{ width_ } )
                {
                    const auto vertex = i * width_ + j;
                    const auto current_elevation = elevation[vertex];
                    if( current_elevation == no_data_value )
                    {
                        continue;
                    }
                    const auto j_contribution =
                        coordinate_system_.direction( 0 ) * j;
                    const auto point = coordinate_system_.origin()
                                       + i_contribution + j_contribution;
                    vertices[vertex] = builder_->create_point(
                        geode::Point3D{ { point.value( 0 ), point.value( 1 ),
                            current_elevation } } );
                }
            }
            return vertices;
        }

        void read_metadata()
        {
            coordinate_system_ = read_coordinate_system();
            width_ = dataset().GetRasterXSize();
            height_ = dataset().GetRasterYSize();
        }

    private:
        std::unique_ptr< geode::PolygonalSurfaceBuilder3D > builder_;
        geode::CoordinateSystem2D coordinate_system_;
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

        auto DEMInput::additional_files() const -> AdditionalFiles
        {
            detail::GDALFile reader{ this->filename() };
            return reader.additional_files< AdditionalFiles >();
        }
    } // namespace internal
} // namespace geode