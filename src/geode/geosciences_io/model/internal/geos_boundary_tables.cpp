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

#include <geode/geosciences_io/model/internal/geos_boundary_tables.hpp>

#include <algorithm>
#include <fstream>
#include <typeinfo>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_format.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/range.hpp>

#include <geode/geometry/nn_search.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/core/brep.hpp>

namespace
{
    using geode::internal::TableCoordinates;

    // Sorted coordinates of the given polygon centers along an axis, the ones
    // closer than GLOBAL_EPSILON being merged
    std::vector< double > axis_coordinates(
        absl::Span< const geode::Point3D > centers,
        absl::Span< const geode::index_t > polygons,
        geode::local_index_t axis )
    {
        std::vector< double > values;
        values.reserve( polygons.size() );
        for( const auto polygon_id : polygons )
        {
            values.push_back( centers[polygon_id].value( axis ) );
        }
        absl::c_sort( values );
        values.erase( std::unique( values.begin(), values.end(),
                          []( double kept, double value ) {
                              return value - kept <= geode::GLOBAL_EPSILON;
                          } ),
            values.end() );
        return values;
    }

    // Replicates GEOS TableFunction nearest lookup along one axis
    geode::index_t geos_nearest_index(
        absl::Span< const double > coordinates, double input )
    {
        if( input <= coordinates.front() )
        {
            return 0;
        }
        const auto last =
            static_cast< geode::index_t >( coordinates.size() - 1 );
        if( input >= coordinates.back() )
        {
            return last;
        }
        const auto upper = static_cast< geode::index_t >(
            absl::c_lower_bound( coordinates, input ) - coordinates.begin() );
        if( input - coordinates[upper - 1] <= coordinates[upper] - input )
        {
            return upper - 1;
        }
        return upper;
    }

    geode::index_t nb_table_nodes( const TableCoordinates& coordinates )
    {
        std::size_t nb_nodes{ 1 };
        for( const auto& axis_coordinates : coordinates )
        {
            nb_nodes *= axis_coordinates.size();
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                nb_nodes < geode::NO_ID, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[GeosExporter] Boundary condition table is too large: ",
                coordinates[0].size(), " x ", coordinates[1].size(), " x ",
                coordinates[2].size(), " nodes." );
        }
        return static_cast< geode::index_t >( nb_nodes );
    }

    // Polygon centers flattened on the table axes holding a single coordinate
    std::vector< geode::Point3D > table_projected_centers(
        const TableCoordinates& coordinates,
        absl::Span< const geode::Point3D > centers )
    {
        std::vector< geode::Point3D > projected_centers{ centers.begin(),
            centers.end() };
        for( const auto axis : geode::LRange{ 3 } )
        {
            if( coordinates[axis].size() != 1 )
            {
                continue;
            }
            for( auto& center : projected_centers )
            {
                center.set_value( axis, coordinates[axis].front() );
            }
        }
        return projected_centers;
    }

    std::vector< geode::index_t > polygons_sharing_node(
        absl::Span< const geode::index_t > polygon_nodes )
    {
        absl::flat_hash_map< geode::index_t, geode::index_t > node_counts;
        for( const auto node : polygon_nodes )
        {
            node_counts[node]++;
        }
        std::vector< geode::index_t > polygons;
        for( const auto polygon_id : geode::Indices{ polygon_nodes } )
        {
            if( node_counts.at( polygon_nodes[polygon_id] ) > 1 )
            {
                polygons.push_back( polygon_id );
            }
        }
        return polygons;
    }
} // namespace

namespace geode::internal
{
    geode::Point3D face_center(
        const geode::SurfaceMesh3D& mesh, geode::index_t polygon_id )
    {
        const auto vertices = mesh.polygon_vertices( polygon_id );
        if( vertices.size() == 3 )
        {
            return mesh.polygon_barycenter( polygon_id );
        }
        const auto& first = mesh.point( vertices[0] );
        geode::Point3D weighted_centers;
        double total_area{ 0 };
        for( const auto vertex : geode::LRange{ 1, vertices.size() - 1 } )
        {
            const auto& second = mesh.point( vertices[vertex] );
            const auto& third = mesh.point( vertices[vertex + 1] );
            const auto area = geode::Vector3D{ first, second }
                                  .cross( geode::Vector3D{ first, third } )
                                  .length();
            weighted_centers += ( first + second + third ) * ( area / 3. );
            total_area += area;
        }
        if( total_area == 0 )
        {
            return mesh.polygon_barycenter( polygon_id );
        }
        return weighted_centers / total_area;
    }

    std::vector< geode::index_t > table_node_polygons(
        const TableCoordinates& coordinates,
        absl::Span< const geode::index_t > polygon_nodes,
        absl::Span< const geode::Point3D > centers )
    {
        const auto nb_nodes = nb_table_nodes( coordinates );
        std::vector< geode::index_t > node_polygons( nb_nodes, geode::NO_ID );
        for( const auto polygon_id : geode::Indices{ polygon_nodes } )
        {
            node_polygons[polygon_nodes[polygon_id]] = polygon_id;
        }
        if( polygon_nodes.size() == nb_nodes )
        {
            return node_polygons;
        }
        const geode::NNSearch3D search{ table_projected_centers(
            coordinates, centers ) };
        const auto nb_x = coordinates[0].size();
        const auto nb_xy = nb_x * coordinates[1].size();
        for( const auto node : geode::Range{ nb_nodes } )
        {
            if( node_polygons[node] != geode::NO_ID )
            {
                continue;
            }
            node_polygons[node] = search.closest_neighbor(
                geode::Point3D{ { coordinates[0][node % nb_x],
                    coordinates[1][( node / nb_x ) % coordinates[1].size()],
                    coordinates[2][node / nb_xy] } } );
        }
        return node_polygons;
    }

    std::vector< geode::index_t > table_polygon_nodes(
        const TableCoordinates& coordinates,
        absl::Span< const geode::Point3D > centers )
    {
        const auto nb_x = coordinates[0].size();
        const auto nb_xy = nb_x * coordinates[1].size();
        std::vector< geode::index_t > polygon_nodes;
        polygon_nodes.reserve( centers.size() );
        for( const auto& center : centers )
        {
            polygon_nodes.push_back(
                geos_nearest_index( coordinates[0], center.value( 0 ) )
                + geos_nearest_index( coordinates[1], center.value( 1 ) ) * nb_x
                + geos_nearest_index( coordinates[2], center.value( 2 ) )
                      * nb_xy );
        }
        return polygon_nodes;
    }

    std::optional< TableCoordinates > separating_table_coordinates(
        absl::Span< const geode::Point3D > centers )
    {
        TableCoordinates coordinates;
        for( const auto axis : geode::LRange{ 3 } )
        {
            coordinates[axis] = { centers.front().value( axis ) };
        }
        std::vector< geode::index_t > colliding_polygons( centers.size() );
        absl::c_iota( colliding_polygons, 0 );
        for( const auto axis : geode::LRange{ 3 } )
        {
            coordinates[axis] =
                axis_coordinates( centers, colliding_polygons, axis );
            colliding_polygons = polygons_sharing_node(
                table_polygon_nodes( coordinates, centers ) );
            if( colliding_polygons.empty() )
            {
                return coordinates;
            }
        }
        return std::nullopt;
    }

    void write_table_values(
        const std::string& filename, absl::Span< const double > values )
    {
        std::ofstream file{ filename };
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            file.good(), nullptr, geode::OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot open file ", filename );
        std::string content;
        for( const auto value : values )
        {
            absl::StrAppendFormat( &content, "%.17g\n", value );
        }
        file << content;
        file.close();
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            !file.fail(), nullptr, geode::OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot write file ", filename );
    }

} // namespace geode::internal
