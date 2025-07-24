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

#include <geode/geosciences_io/mesh/internal/polytiff_input.hpp>

#include <gdal_priv.h>

#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/polygonal_surface_builder.hpp>
#include <geode/mesh/core/light_regular_grid.hpp>
#include <geode/mesh/core/polygonal_surface.hpp>
#include <geode/mesh/helpers/convert_surface_mesh.hpp>
#include <geode/mesh/io/light_regular_grid_input.hpp>

namespace
{
    class PolyTIFFInputImpl
    {
    public:
        PolyTIFFInputImpl( std::string_view filename )
            : filename_{ filename },
              gdal_data_{ GDALDataset::Open(
                  geode::to_string( filename ).c_str(), GDAL_OF_READONLY ) }
        {
            OPENGEODE_EXCEPTION( gdal_data_,
                "[PolyTIFFInputImpl] Failed to open file ", filename );
        }

        std::unique_ptr< geode::PolygonalSurface3D > read_file()
        {
            const auto level = geode::Logger::level();
            geode::Logger::set_level( geode::Logger::LEVEL::critical );
            const auto grid = geode::load_light_regular_grid< 2 >( filename_ );
            geode::Logger::set_level( level );
            auto surface2d = geode::convert_grid_into_polygonal_surface( grid );
            const auto nb_bands = gdal_data_->GetRasterCount();
            OPENGEODE_EXCEPTION(
                nb_bands > 0, "[PolyTIFFInput] No bands found" );
            absl::FixedArray< float > elevation( surface2d->nb_vertices() );
            const auto band = gdal_data_->GetRasterBand( 1 );
            const auto width = gdal_data_->GetRasterXSize();
            const auto height = gdal_data_->GetRasterYSize();
            const auto status = band->RasterIO( GF_Read, 0, 0, width, height,
                elevation.data(), width, height, GDT_Float32, 0, 0 );
            OPENGEODE_EXCEPTION(
                status == CE_None, "[PolyTIFFInput] Failed to read elevation" );
            const auto no_data_value = band->GetNoDataValue();
            auto surface3d = geode::convert_polygonal_surface2d_into_3d(
                *surface2d, 2, no_data_value );
            auto builder =
                geode::PolygonalSurfaceBuilder3D::create( *surface3d );
            std::vector< geode::index_t > vertices_to_delete;
            for( const auto vertex : geode::Range{ surface3d->nb_vertices() } )
            {
                const auto current_elevation = elevation[vertex];
                if( current_elevation == no_data_value )
                {
                    vertices_to_delete.push_back( vertex );
                    continue;
                }
                auto point = surface3d->point( vertex );
                point.set_value( 2, current_elevation );
                builder->set_point( vertex, point );
            }
            std::vector< bool > polygon_to_delete(
                surface3d->nb_polygons(), false );
            for( const auto vertex : vertices_to_delete )
            {
                for( const auto [polygon, _] :
                    surface3d->polygons_around_vertex( vertex ) )
                {
                    polygon_to_delete[polygon] = true;
                }
            }
            builder->delete_polygons( polygon_to_delete );
            builder->delete_isolated_vertices();
            return surface3d;
        }

    private:
        std::string_view filename_;
        GDALDatasetUniquePtr gdal_data_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::unique_ptr< PolygonalSurface3D > PolyTIFFInput::read(
            const MeshImpl& /*impl*/ )
        {
            PolyTIFFInputImpl geo_reader( filename() );
            return geo_reader.read_file();
        }
    } // namespace internal
} // namespace geode