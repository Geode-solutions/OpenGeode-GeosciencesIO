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

#include <geode/geosciences_io/model/internal/simulation_boundary_values.hpp>

#include <array>
#include <optional>
#include <string_view>
#include <typeinfo>
#include <utility>
#include <vector>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_map.h>
#include <absl/types/span.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/range.hpp>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/simulation/core/fluid.hpp>

namespace
{
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

    constexpr std::
        array< std::pair< geode::PHYSICAL_PROPERTY_NAME, geode::PHASE >, 3 >
            BOUNDARY_FRACTION_PHASES{ {
                { geode::PHYSICAL_PROPERTY_NAME::boundary_oil_fraction,
                    geode::PHASE::oil },
                { geode::PHYSICAL_PROPERTY_NAME::boundary_gas_fraction,
                    geode::PHASE::gas },
                { geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction,
                    geode::PHASE::water },
            } };
} // namespace

namespace geode::internal
{
    template < typename T >
    std::optional< std::vector< T > > surface_values(
        const geode::Surface3D& surface,
        const geode::uuid& attribute_id,
        std::string_view context )
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
            nullptr, geode::OpenGeodeException::TYPE::data, context,
            " Boundary condition attribute ", attribute_id.string(),
            " on Surface ", surface.id().string(), " must store ",
            value_type_name< T >(), " values." );
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
                geode::OpenGeodeException::TYPE::data, context,
                " Boundary condition attribute ", attribute_id.string(),
                " has no value on polygon ", polygon_id, " of Surface ",
                surface.id().string(),
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
        const geode::uuid& attribute_id,
        std::string_view context )
    {
        SurfaceValueGroups< T > groups;
        for( const auto& surface : model.surfaces() )
        {
            auto values = surface_values< T >( surface, attribute_id, context );
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

    ModelBoundaryConditions model_boundary_conditions( const BRep& model )
    {
        ModelBoundaryConditions conditions;
        conditions.has_pressure = model.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_pressure );
        conditions.has_temperature = model.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_temperature );
        for( const auto& [property, phase] : BOUNDARY_FRACTION_PHASES )
        {
            if( model.has_physical_property( property ) )
            {
                conditions.fraction_phases.push_back( phase );
            }
        }
        return conditions;
    }

    template std::optional< std::vector< double > > surface_values(
        const Surface3D&, const uuid&, std::string_view );
    template std::optional< std::vector< std::array< double, 3 > > >
        surface_values( const Surface3D&, const uuid&, std::string_view );
    template SurfaceValueGroups< double > group_surface_values( const BRep&,
        const absl::flat_hash_map< uuid, index_t >&,
        const uuid&,
        std::string_view );
    template SurfaceValueGroups< std::array< double, 3 > > group_surface_values(
        const BRep&,
        const absl::flat_hash_map< uuid, index_t >&,
        const uuid&,
        std::string_view );
} // namespace geode::internal
