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

#include <geode/geosciences_io/model/internal/geos_boundary_conditions.hpp>

#include <algorithm>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_map.h>

#include <geode/basic/range.hpp>

#include <geode/geometry/nn_search.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>

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
        return static_cast< geode::index_t >( coordinates[0].size()
                                              * coordinates[1].size()
                                              * coordinates[2].size() );
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
    const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME, GeosBoundaryField >&
        geos_boundary_condition_fields()
    {
        static const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME,
            GeosBoundaryField >
            fields{
                { PHYSICAL_PROPERTY_NAME::boundary_pressure,
                    { "pressure", "boundaryPressure" } },
                { PHYSICAL_PROPERTY_NAME::boundary_temperature,
                    { "temperature", "boundaryTemperature" } },
            };
        return fields;
    }

    const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME, GeosFluidFraction >&
        geos_fluid_fractions()
    {
        static const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME,
            GeosFluidFraction >
            fractions{
                { PHYSICAL_PROPERTY_NAME::boundary_oil_fraction,
                    { "oil", PHASE::oil } },
                { PHYSICAL_PROPERTY_NAME::boundary_gas_fraction,
                    { "gas", PHASE::gas } },
                { PHYSICAL_PROPERTY_NAME::boundary_water_fraction,
                    { "water", PHASE::water } },
            };
        return fractions;
    }

    std::optional< local_index_t > geos_boundary_fraction_component(
        const PhysicalProblem& problem, const GeosFluidFraction& fraction )
    {
        if( !problem.has_flow() || !problem.flow().has_fluid() )
        {
            return std::nullopt;
        }
        const auto component =
            phase_component( problem.fluid(), fraction.phase );
        if( !component )
        {
            return std::nullopt;
        }
        return static_cast< local_index_t >( component.value() );
    }

    std::vector< PHYSICAL_PROPERTY_NAME > geos_boundary_properties(
        const PhysicalProblem& problem )
    {
        std::vector< PHYSICAL_PROPERTY_NAME > properties{
            PHYSICAL_PROPERTY_NAME::boundary_pressure
        };
        if( !problem.has_flow() )
        {
            return properties;
        }
        const auto is_multiphase =
            dynamic_cast< const MultiphaseFlow* >( &problem.flow() ) != nullptr;
        if( is_multiphase || problem.thermal() )
        {
            properties.push_back(
                PHYSICAL_PROPERTY_NAME::boundary_temperature );
        }
        if( !is_multiphase || !problem.flow().has_fluid() )
        {
            return properties;
        }
        for( const auto& [property, fraction] : geos_fluid_fractions() )
        {
            if( absl::c_linear_search(
                    problem.fluid().phases(), fraction.phase ) )
            {
                properties.push_back( property );
            }
        }
        return properties;
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
} // namespace geode::internal
