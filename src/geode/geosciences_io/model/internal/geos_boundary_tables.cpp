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

    template < typename T >
    std::string_view value_type_name();

    template <>
    std::string_view value_type_name< double >()
    {
        return "double";
    }

    template <>
    std::string_view value_type_name< std::array< double, 3 > >()
    {
        return "std::array< double, 3 >";
    }

    template < typename T >
    bool is_uniform( absl::Span< const T > values )
    {
        return absl::c_all_of( values, [&values]( const T& value ) {
            return value == values.front();
        } );
    }

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
    template < typename T >
    std::optional< std::vector< T > > surface_values(
        const geode::Surface3D& surface, const geode::uuid& attribute_id )
    {
        const auto& mesh = surface.mesh();
        const auto& manager = mesh.polygon_attribute_manager();
        if( !manager.attribute_exists( attribute_id )
            || mesh.nb_polygons() == 0 )
        {
            return std::nullopt;
        }
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            manager.attribute_type( attribute_id ) == typeid( T ).name(),
            nullptr, geode::OpenGeodeException::TYPE::data,
            "[GeosExporter] Boundary condition attribute ",
            attribute_id.string(), " on Surface ", surface.id().string(),
            " must store ", value_type_name< T >(), " values." );
        const auto attribute =
            manager.find_read_only_attribute< T >( attribute_id );
        if( const auto* constant =
                dynamic_cast< const geode::ConstantAttribute< T >* >(
                    attribute.get() ) )
        {
            return std::vector< T >{ constant->value() };
        }
        std::vector< T > values;
        values.reserve( mesh.nb_polygons() );
        for( const auto polygon_id : geode::Range{ mesh.nb_polygons() } )
        {
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                attribute->has_value( polygon_id ), nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[GeosExporter] Boundary condition attribute ",
                attribute_id.string(), " has no value on polygon ", polygon_id,
                " of Surface ", surface.id().string(),
                ". Boundary conditions must be defined on every polygon of "
                "a Surface." );
            values.push_back( attribute->value( polygon_id ) );
        }
        return values;
    }

    template < typename T >
    SurfaceValueGroups< T > group_surface_values( const BRep& model,
        const absl::flat_hash_map< geode::uuid, geode::index_t >&
            surface_regions,
        const geode::uuid& attribute_id )
    {
        SurfaceValueGroups< T > groups;
        for( const auto& surface : model.surfaces() )
        {
            auto values = surface_values< T >( surface, attribute_id );
            if( !values )
            {
                continue;
            }
            const auto region_id = surface_regions.at( surface.id() );
            if( is_uniform< T >( values.value() ) )
            {
                groups.uniform_regions[values->front()].push_back( region_id );
                continue;
            }
            groups.variable_surfaces.push_back(
                { &surface, region_id, std::move( values.value() ) } );
        }
        return groups;
    }

    std::vector< double > axis_values(
        absl::Span< const std::array< double, 3 > > values,
        geode::local_index_t axis )
    {
        std::vector< double > component_values;
        component_values.reserve( values.size() );
        for( const auto& value : values )
        {
            component_values.push_back( value[axis] );
        }
        return component_values;
    }

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

    template std::optional< std::vector< double > > surface_values(
        const Surface3D&, const uuid& );
    template std::optional< std::vector< std::array< double, 3 > > >
        surface_values( const Surface3D&, const uuid& );
    template SurfaceValueGroups< double > group_surface_values(
        const BRep&, const absl::flat_hash_map< uuid, index_t >&, const uuid& );
    template SurfaceValueGroups< std::array< double, 3 > > group_surface_values(
        const BRep&, const absl::flat_hash_map< uuid, index_t >&, const uuid& );
} // namespace geode::internal
