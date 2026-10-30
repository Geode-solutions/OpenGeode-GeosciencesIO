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

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <typeinfo>
#include <utility>

#include <pugixml.hpp>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_map.h>
#include <absl/container/linked_hash_map.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_join.h>
#include <absl/types/span.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/filename.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/mapping.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/aabb.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/distance.hpp>
#include <geode/geometry/nn_search.hpp>
#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/geosciences_io/model/common.hpp>

#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/polyhedral_solid.hpp>
#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/surface_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/helpers/aabb_solid_helpers.hpp>

#include <geode/mesh/io/point_set_output.hpp>

#include <geode/model/helpers/component_mesh_polygons.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/block_collection.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>
#include <geode/model/mixin/core/surface.hpp>

#include <geode/io/mesh/detail/vtu_solid_output_impl.hpp>

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/solid_mechanics.hpp>
#include <geode/geosciences/simulation/helpers/physical_problem_validation.hpp>
#include <geode/geosciences_io/model/internal/geos_export.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_validation.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_writer.hpp>
#include <geode/geosciences_io/model/internal/geos_xml.hpp>
#include <geode/model/representation/core/brep.hpp>

namespace
{
    // Writes the solid polyhedra followed by the surface polygons in a single
    // UnstructuredGrid, all the cells sharing the solid vertices
    template < typename SolidOutputImpl >
    class GeosVTUOutputImpl : public SolidOutputImpl
    {
    public:
        template < typename Solid >
        GeosVTUOutputImpl( std::string_view filename,
            const Solid& solid,
            absl::Span< const geode::PolygonVertices > surface_cells,
            const geode::AttributeManager& cell_attributes )
            : SolidOutputImpl{ filename, solid },
              surface_cells_{ surface_cells },
              cell_attributes_( cell_attributes )
        {
        }

    private:
        [[nodiscard]] geode::index_t nb_additional_polygons() const override
        {
            return surface_cells_.size();
        }

        [[nodiscard]] absl::Span< const geode::index_t >
            additional_polygon_vertices(
                geode::index_t polygon_id ) const override
        {
            return surface_cells_[polygon_id];
        }

        [[nodiscard]] const geode::AttributeManager&
            cell_attribute_manager() const override
        {
            return cell_attributes_;
        }

        absl::Span< const geode::PolygonVertices > surface_cells_;
        const geode::AttributeManager& cell_attributes_;
    };

    template < typename SolidOutputImpl, typename Solid >
    void save_geos_vtu( std::string_view filename,
        const Solid& solid,
        absl::Span< const geode::PolygonVertices > surface_cells,
        const geode::AttributeManager& cell_attributes )
    {
        GeosVTUOutputImpl< SolidOutputImpl > writer{ filename, solid,
            surface_cells, cell_attributes };
        writer.write_file();
    }

    struct GeosBoundaryField
    {
        std::string_view field_name;
        std::string_view name_prefix;
    };

    const absl::linked_hash_map< geode::PHYSICAL_PROPERTY_NAME,
        GeosBoundaryField >
        BOUNDARY_CONDITION_GEOS_FIELDS{
            { geode::PHYSICAL_PROPERTY_NAME::boundary_pressure,
                { "pressure", "boundaryPressure" } },
            { geode::PHYSICAL_PROPERTY_NAME::boundary_temperature,
                { "temperature", "boundaryTemperature" } },
        };

    struct GeosFluidFraction
    {
        std::string_view fluid;
        geode::PHASE phase;
        // Component index used without PhysicalProblem
        geode::local_index_t legacy_component;
    };

    const absl::linked_hash_map< geode::PHYSICAL_PROPERTY_NAME,
        GeosFluidFraction >
        FLUID_FRACTION_GEOS_NAMES{
            { geode::PHYSICAL_PROPERTY_NAME::boundary_oil_fraction,
                { "oil", geode::PHASE::oil, 0 } },
            { geode::PHYSICAL_PROPERTY_NAME::boundary_gas_fraction,
                { "gas", geode::PHASE::gas, 1 } },
            { geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction,
                { "water", geode::PHASE::water, 2 } },
        };

    // GEOS component of a boundary phase fraction: the component carried by
    // the phase in the problem fluid, or a fixed oil=0, gas=1, water=2 index
    // without PhysicalProblem
    std::optional< geode::local_index_t > boundary_fraction_component(
        const geode::PhysicalProblem* problem,
        const GeosFluidFraction& fraction )
    {
        if( problem == nullptr )
        {
            return fraction.legacy_component;
        }
        if( !problem->has_flow() || !problem->flow().has_fluid() )
        {
            return std::nullopt;
        }
        const auto component =
            geode::phase_component( problem->fluid(), fraction.phase );
        if( !component )
        {
            return std::nullopt;
        }
        return static_cast< geode::local_index_t >( component.value() );
    }

    void check_no_issue( absl::Span< const std::string > issues )
    {
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            issues.empty(), nullptr, geode::OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot write the GEOS simulation deck:\n - ",
            absl::StrJoin( issues, "\n - " ) );
    }

    constexpr auto BOUNDARY_TABLES_DIRECTORY = "boundary_tables";
    constexpr std::array< std::string_view, 3 > AXIS_NAMES{ "x", "y", "z" };

    using TableCoordinates = std::array< std::vector< double >, 3 >;

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

    // Value of the attribute on each Surface polygon, reduced to a single value
    // for a ConstantAttribute
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

    bool is_uniform( absl::Span< const double > values )
    {
        return absl::c_all_of( values, [&values]( double value ) {
            return value == values.front();
        } );
    }

    // Free displacement components are NaN, which equals no value
    bool same_components( const std::array< double, 3 >& first,
        const std::array< double, 3 >& second )
    {
        for( const auto axis : geode::LRange{ 3 } )
        {
            if( first[axis] != second[axis]
                && !(
                    std::isnan( first[axis] ) && std::isnan( second[axis] ) ) )
            {
                return false;
            }
        }
        return true;
    }

    bool is_uniform( absl::Span< const std::array< double, 3 > > values )
    {
        return absl::c_all_of(
            values, [&values]( const std::array< double, 3 >& value ) {
                return same_components( value, values.front() );
            } );
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

    // Polygon whose value is stored at each table node, in Fortran order.
    // Nodes looked up by no polygon center take the value of the closest one.
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

    // Adds the x, y and z axes until every polygon center is looked up by
    // GEOS on its own node. Each added axis only gets the coordinates of the
    // polygons still sharing a node, which keeps the table small.
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
        for( const auto value : values )
        {
            file << absl::StrFormat( "%.17g", value ) << "\n";
        }
    }

    // Mean of one item of a numeric attribute over the first elements
    double mean_attribute_item( const geode::AttributeBase& attribute,
        geode::index_t nb_elements,
        geode::local_index_t item )
    {
        double sum{ 0 };
        for( const auto element : geode::Range{ nb_elements } )
        {
            sum += attribute.generic_item_value( element, item );
        }
        return nb_elements == 0 ? 0. : sum / nb_elements;
    }

    // Mean and minimum of a scalar cell attribute over the first elements,
    // read in double precision when the attribute stores doubles
    std::pair< double, double > attribute_mean_and_min(
        const geode::AttributeManager& manager,
        const geode::uuid& attribute_id,
        geode::index_t nb_elements )
    {
        const auto is_double =
            manager.attribute_type( attribute_id ) == typeid( double ).name();
        const auto generic = manager.find_generic_attribute( attribute_id );
        const auto typed =
            is_double
                ? manager.find_read_only_attribute< double >( attribute_id )
                : nullptr;
        double sum{ 0 };
        double min{ std::numeric_limits< double >::max() };
        for( const auto element : geode::Range{ nb_elements } )
        {
            const auto value = is_double
                                   ? typed->value( element )
                                   : static_cast< double >(
                                       generic->generic_value( element ) );
            sum += value;
            min = std::min( min, value );
        }
        return { nb_elements == 0 ? 0. : sum / nb_elements, min };
    }

    template < typename Model >
    const geode::uuid& boundary_condition_attribute_id(
        const Model& model, geode::PHYSICAL_PROPERTY_NAME property )
    {
        const auto& property_info = model.physical_property_info( property );
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            property_info.component_type
                == geode::Surface3D::component_type_static(),
            nullptr, geode::OpenGeodeException::TYPE::data,
            "[GeosExporter] Physical property attribute ",
            property_info.attribute_id.string(),
            " must be defined on Surfaces to be exported as a boundary "
            "condition." );
        return property_info.attribute_id;
    }
} // namespace

namespace geode::internal
{
    constexpr auto REGION_ID_ATTRIBUTE_NAME = "attribute";

    const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME, std::string_view >
        PHYSICAL_PROPERTY_GEOS_NAMES{
            { PHYSICAL_PROPERTY_NAME::permeability, "rockPerm_permeability" },
            { PHYSICAL_PROPERTY_NAME::porosity,
                "rockPorosity_referencePorosity" },
        };

    constexpr std::array< PHYSICAL_PROPERTY_NAME, 3 > MECHANICS_PROPERTIES{
        PHYSICAL_PROPERTY_NAME::bulk_modulus,
        PHYSICAL_PROPERTY_NAME::shear_modulus, PHYSICAL_PROPERTY_NAME::density
    };

    constexpr auto INITIAL_STRESS_TENSOR = "initial_stress_tensor";

    using StressTensor = std::array< double, 6 >;

    // Block values are read one by one: importing one ConstantAttribute per
    // Block onto the solid would keep only the last Block value.
    template < typename Model >
    void transfer_stress_components( const Model& model,
        const ModelToMeshMappings& model2solid,
        const PhysicalProperties::Info& property_info,
        std::string_view name,
        local_index_t first_component,
        VariableAttribute< StressTensor >& tensor )
    {
        const auto& attribute_id = property_info.attribute_id;
        OpenGeodeGeosciencesIOModelException::check_exception(
            property_info.component_type == Block3D::component_type_static(),
            nullptr, OpenGeodeException::TYPE::data, "[GeosExporter] ", name,
            " physical property attribute ", attribute_id.string(),
            " must be defined on Blocks to be exported." );
        const auto& polyhedra_mapping = model2solid.solid_polyhedra_mapping;
        for( const auto& block : model.blocks() )
        {
            const auto& block_mesh = block.mesh();
            const auto& block_manager =
                block_mesh.polyhedron_attribute_manager();
            OpenGeodeGeosciencesIOModelException::check_exception(
                block_manager.attribute_exists( attribute_id ), nullptr,
                OpenGeodeException::TYPE::data, "[GeosExporter] ", name,
                " physical property attribute ", attribute_id.string(),
                " is missing on Block ", block.id().string(), "." );
            OpenGeodeGeosciencesIOModelException::check_exception(
                block_manager.attribute_type( attribute_id )
                    == typeid( std::array< double, 3 > ).name(),
                nullptr, OpenGeodeException::TYPE::data, "[GeosExporter] ",
                name, " physical property attribute ", attribute_id.string(),
                " must store std::array< double, 3 > values." );
            const auto stress = block_manager.template find_read_only_attribute<
                std::array< double, 3 > >( attribute_id );
            for( const auto polyhedron_id : Range{ block_mesh.nb_polyhedra() } )
            {
                const MeshElement block_polyhedron{ block.id(), polyhedron_id };
                if( !polyhedra_mapping.has_mapping_input( block_polyhedron ) )
                {
                    continue;
                }
                const auto& value = stress->value( polyhedron_id );
                tensor.modify_value(
                    polyhedra_mapping.in2out( block_polyhedron ),
                    [&value, first_component]( StressTensor& tensor_value ) {
                        // GEOS counts tractions positively
                        for( const auto axis : LRange{ 3 } )
                        {
                            tensor_value[first_component + axis] = -value[axis];
                        }
                    } );
            }
        }
    }

    template < typename Model >
    std::string transfer_block_attribute( const Model& model,
        SolidMesh3D& solid,
        const ModelToMeshMappings& model2solid,
        const PhysicalProperties::Info& property_info )
    {
        const auto& attribute_id = property_info.attribute_id;
        OpenGeodeGeosciencesIOModelException::check_exception(
            property_info.component_type == Block3D::component_type_static(),
            nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Physical property attribute ",
            attribute_id.string(),
            " must be defined on Blocks to be exported." );
        auto& solid_manager = solid.polyhedron_attribute_manager();
        const auto& polyhedra_mapping = model2solid.solid_polyhedra_mapping;
        for( const auto& block : model.blocks() )
        {
            const auto& block_mesh = block.mesh();
            const auto& block_manager =
                block_mesh.polyhedron_attribute_manager();
            OpenGeodeGeosciencesIOModelException::check_exception(
                block_manager.attribute_exists( attribute_id ), nullptr,
                OpenGeodeException::TYPE::data,
                "[GeosExporter] Physical property attribute ",
                attribute_id.string(), " is missing on Block ",
                block.id().string(), "." );
            GenericMapping< index_t > block2solid;
            for( const auto polyhedron_id : Range( block_mesh.nb_polyhedra() ) )
            {
                const MeshElement block_polyhedron{ block.id(), polyhedron_id };
                if( !polyhedra_mapping.has_mapping_input( block_polyhedron ) )
                {
                    continue;
                }
                block2solid.map( polyhedron_id,
                    polyhedra_mapping.in2out( block_polyhedron ) );
            }
            solid_manager.import( block_manager, block2solid, attribute_id );
        }
        const auto attribute =
            solid_manager.find_generic_attribute( attribute_id );
        OpenGeodeGeosciencesIOModelException::check_exception(
            attribute && attribute->name().has_value(), nullptr,
            OpenGeodeException::TYPE::data,
            "[GeosExporter] Physical property attribute ",
            attribute_id.string(), " has no name." );
        return attribute->name().value();
    }

    template < typename Model >
    GeosExporterImpl< Model >::GeosExporterImpl(
        std::string_view files_directory, const Model& model )
        : model_( model ),
          files_directory_{
              std::filesystem::path{ to_string( files_directory ) }.string()
          },
          prefix_{ filename_without_extension( files_directory ).string() }
    {
        if( std::filesystem::path{ to_string( files_directory ) }
                .is_relative() )
        {
            std::filesystem::create_directory(
                std::filesystem::current_path() / files_directory_ );
        }
        else
        {
            std::filesystem::create_directory( files_directory_ );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::convert_model()
    {
        auto curve_conversion = convert_brep_into_curve( model_ );
        model_curve_ = std::move( std::get< 0 >( curve_conversion ) );
        std::tie( model_solid_, model2solid_ ) =
            convert_brep_into_solid( model_ );
        AttributeValues< index_t > region_attribute_values;
        region_attribute_values.default_value = NO_ID;
        region_attribute_values.no_value = NO_ID;
        AttributeProperties region_attribute_properties;
        region_attribute_properties.assignable = false;
        region_attribute_properties.interpolable = false;
        region_attribute_properties.transferable = true;
        region_attribute_id_ =
            model_solid_->polyhedron_attribute_manager()
                .template create_attribute< VariableAttribute, index_t >(
                    REGION_ID_ATTRIBUTE_NAME, region_attribute_values,
                    region_attribute_properties );
        region_attribute_ = model_solid_->polyhedron_attribute_manager()
                                .find_attribute< VariableAttribute, index_t >(
                                    region_attribute_id_ );
    }

    template < typename Model >
    GeosExporterImpl< Model >::~GeosExporterImpl() = default;

    template < typename Model >
    void GeosExporterImpl< Model >::set_physical_problem(
        const PhysicalProblem& problem )
    {
        problem_ = &problem;
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_files() const
    {
        pugi::xml_document doc_xml;
        auto pb_node = doc_xml.append_child( "Problem" );

        auto mesh_node = pb_node.append_child( "Mesh" );
        write_mesh_files( mesh_node );

        if( !well_perforations_.empty() )
        {
            auto geometry_node = pb_node.append_child( "Geometry" );
            write_well_perforations_boxes( geometry_node );
            write_well_perforation_file();
        }

        if( writer_ )
        {
            write_simulation_sections( pb_node );
        }
        else
        {
            if( !boundary_conditions_.empty() )
            {
                auto field_specifications_node =
                    pb_node.append_child( "FieldSpecifications" );
                write_boundary_conditions( field_specifications_node );
            }

            if( !boundary_tables_.empty() )
            {
                auto functions_node = pb_node.append_child( "Functions" );
                write_boundary_functions( functions_node );
            }
        }
        if( !boundary_tables_.empty() )
        {
            write_boundary_table_files();
        }

        const auto filename_xml =
            absl::StrCat( files_directory(), "/", prefix(), "_simulation.xml" );
        doc_xml.save_file( filename_xml.c_str(), PUGIXML_TEXT( "    " ),
            pugi::format_indent_attributes );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::add_well_perforations(
        const PointSet3D& perforations, std::string_view name )
    {
        OpenGeodeGeosciencesIOModelException::check_exception(
            absl::c_find( well_names_, name ) == well_names_.end(), nullptr,
            OpenGeodeException::TYPE::data,
            "[GeosExporter] Well perforations named ", name,
            " already added." );
        well_perforations_.push_back( perforations.clone() );
        well_names_.emplace_back( to_string( name ) );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::prepare_export()
    {
        if( problem_ != nullptr )
        {
            validate_problem();
        }
        convert_model();
        const auto nb_solid_regions = initialize_solid_region_attribute();
        initialize_surface_cells( nb_solid_regions );
        transfer_physical_properties();
        transfer_boundary_conditions();
        if( problem_ != nullptr && problem_->solid_mechanics() )
        {
            transfer_mechanics_properties();
            transfer_mechanics_boundary_conditions();
        }
        if( problem_ != nullptr )
        {
            prepare_simulation();
        }
    }

    template < typename Model >
    std::string_view GeosExporterImpl< Model >::files_directory() const
    {
        return files_directory_;
    }
    template < typename Model >
    std::string_view GeosExporterImpl< Model >::prefix() const
    {
        return prefix_;
    }

    template < typename Model >
    index_t GeosExporterImpl< Model >::initialize_solid_region_attribute()
    {
        region_definitions_ = create_region_definitions( model_ );
        add_block_collection_regions();
        const auto& block_regions = region_definitions_.block_regions;
        for( const auto polyhedron_id : Range( model_solid_->nb_polyhedra() ) )
        {
            const auto& block_id =
                model2solid_.solid_polyhedra_mapping.out2in( polyhedron_id )
                    .mesh_id;
            const auto region = block_regions.find( block_id );
            OpenGeodeGeosciencesIOModelException::check_exception(
                region != block_regions.end(), nullptr,
                OpenGeodeException::TYPE::internal, "[GeosExporter] Block ",
                block_id.string(), " belongs to no region." );
            region_attribute_->set_value( polyhedron_id, region->second );
        }
        return static_cast< index_t >( region_definitions_.regions.size() );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::add_block_collection_regions()
    {
        for( const auto& collection : model_.block_collections() )
        {
            std::vector< index_t > collection_regions;
            bool has_merged_block{ false };
            for( const auto& block :
                model_.block_collection_items( collection ) )
            {
                if( region_definitions_.merged_blocks.contains( block.id() ) )
                {
                    has_merged_block = true;
                    break;
                }
                collection_regions.push_back(
                    region_definitions_.block_regions.at( block.id() ) );
            }
            if( !has_merged_block )
            {
                region_definitions_.regions_by_uuid.emplace(
                    collection.id(), std::move( collection_regions ) );
            }
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::initialize_surface_cells(
        index_t first_surface_region_id )
    {
        auto region_id = first_surface_region_id;
        for( const auto& surface : model_.surfaces() )
        {
            for( const auto polygon_id : Range{ surface.mesh().nb_polygons() } )
            {
                auto vertices =
                    polygon_unique_vertices( model_, surface, polygon_id );
                for( auto& vertex : vertices )
                {
                    vertex =
                        model2solid_.unique_vertices_mapping.in2out( vertex );
                }
                surface_cells_.emplace_back( std::move( vertices ) );
                surface_cells_region_.push_back( region_id );
            }
            surface_regions_.emplace( surface.id(), region_id );
            region_id++;
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_well_perforations_boxes(
        pugi::xml_node& root ) const
    {
        auto aabb = create_aabb_tree( *model_solid_ );
        for( const auto well_id : Indices{ well_perforations_ } )
        {
            const auto& well = well_perforations_[well_id];
            BoundingBox3D perf_box;
            for( const auto point : Range( well->nb_vertices() ) )
            {
                const auto neigh_cells =
                    aabb.containing_boxes( well->point( point ) );
                if( neigh_cells.empty() )
                {
                    continue;
                }
                double distance_to_nearest_cell_center{
                    std::numeric_limits< double >::max()
                };
                index_t selected_cell_id{ NO_ID };
                for( const auto cell_id : neigh_cells )
                {
                    const auto tmp_dist =
                        point_point_distance< 3 >( well->point( point ),
                            model_solid_->polyhedron_barycenter( cell_id ) );
                    if( distance_to_nearest_cell_center < tmp_dist )
                    {
                        continue;
                    }
                    distance_to_nearest_cell_center = tmp_dist;
                    selected_cell_id = cell_id;
                }
                for( auto& vertex_id :
                    model_solid_->polyhedron_vertices( selected_cell_id ) )
                {
                    perf_box.add_point( model_solid_->point( vertex_id ) );
                }
            }
            auto box_node = root.append_child( "Box" );
            box_node.append_attribute( "name" ).set_value(
                well_names_[well_id].c_str() );
            static constexpr auto SAFETY_OFFSET = 100. * GLOBAL_EPSILON;
            box_node.append_attribute( "xMin" ).set_value(
                absl::StrCat( "{", perf_box.min().value( 0 ) - SAFETY_OFFSET,
                    ", ", perf_box.min().value( 1 ) - SAFETY_OFFSET, ", ",
                    perf_box.min().value( 2 ) - SAFETY_OFFSET, "}" )
                    .c_str() );
            box_node.append_attribute( "xMax" ).set_value(
                absl::StrCat( "{", perf_box.max().value( 0 ) + SAFETY_OFFSET,
                    ", ", perf_box.max().value( 1 ) + SAFETY_OFFSET, ", ",
                    perf_box.max().value( 2 ) + SAFETY_OFFSET, "}" )
                    .c_str() );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_mesh_files(
        pugi::xml_node& root ) const
    {
        const auto file_vtu = write_solid_file();

        auto vtk_mesh_node = root.append_child( "VTKMesh" );
        vtk_mesh_node.append_attribute( "name" ).set_value(
            to_string( prefix() ).c_str() );
        vtk_mesh_node.append_attribute( "file" ).set_value(
            absl::StrCat( "./", file_vtu ).c_str() );
        if( writer_ )
        {
            std::vector< std::string > vtu_names;
            std::vector< std::string > geos_names_list;
            for( const auto& field : writer_->imported_fields() )
            {
                vtu_names.push_back( field.vtu_name );
                geos_names_list.push_back( field.geos_name );
            }
            if( !vtu_names.empty() )
            {
                set_attribute(
                    vtk_mesh_node, "fieldsToImport", geos_names( vtu_names ) );
                set_attribute( vtk_mesh_node, "fieldNamesInGEOS",
                    geos_names( geos_names_list ) );
            }
            return;
        }
        if( imported_fields_.empty() )
        {
            return;
        }
        const auto join_fields = [this]( bool vtu_names ) {
            return absl::StrCat( "{ ",
                absl::StrJoin( imported_fields_, ", ",
                    [vtu_names]( std::string* out, const auto& field ) {
                        absl::StrAppend(
                            out, vtu_names ? field.first : field.second );
                    } ),
                " }" );
        };
        vtk_mesh_node.append_attribute( "fieldsToImport" )
            .set_value( join_fields( true ).c_str() );
        vtk_mesh_node.append_attribute( "fieldNamesInGEOS" )
            .set_value( join_fields( false ).c_str() );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_boundary_conditions(
        pugi::xml_node& root ) const
    {
        for( const auto& condition : boundary_conditions_ )
        {
            const auto is_traction = condition.kind != BOUNDARY_KIND::field;
            auto specification_node = root.append_child(
                is_traction ? "Traction" : "FieldSpecification" );
            specification_node.append_attribute( "name" ).set_value(
                condition.name.c_str() );
            specification_node.append_attribute( "objectPath" )
                .set_value( condition.object_path.c_str() );
            if( is_traction )
            {
                specification_node.append_attribute( "tractionType" )
                    .set_value( condition.kind == BOUNDARY_KIND::vector_traction
                                    ? "vector"
                                    : "normal" );
            }
            else
            {
                specification_node.append_attribute( "fieldName" )
                    .set_value( condition.field_name.c_str() );
            }
            if( condition.direction )
            {
                const auto& direction = condition.direction.value();
                specification_node.append_attribute( "direction" )
                    .set_value( absl::StrFormat( "{ %.15g, %.15g, %.15g }",
                        direction[0], direction[1], direction[2] )
                                    .c_str() );
            }
            if( condition.component )
            {
                specification_node.append_attribute( "component" )
                    .set_value(
                        static_cast< unsigned int >( *condition.component ) );
            }
            if( condition.function_name )
            {
                specification_node.append_attribute( "functionName" )
                    .set_value( condition.function_name->c_str() );
            }
            specification_node.append_attribute( "scale" ).set_value(
                absl::StrFormat( "%.15g", condition.scale ).c_str() );
            specification_node.append_attribute( "setNames" )
                .set_value( absl::StrCat(
                    "{ ", absl::StrJoin( condition.region_ids, ", " ), " }" )
                                .c_str() );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_boundary_functions(
        pugi::xml_node& root ) const
    {
        for( const auto& table : boundary_tables_ )
        {
            const auto& grid = boundary_grids_.at( table.surface_id );
            auto function_node = root.append_child( "TableFunction" );
            function_node.append_attribute( "name" ).set_value(
                table.function_name.c_str() );
            function_node.append_attribute( "inputVarNames" )
                .set_value( "{ faceCenter }" );
            function_node.append_attribute( "coordinateFiles" )
                .set_value( absl::StrCat( "{ ",
                    absl::StrJoin( grid.files, ", ",
                        []( std::string* out, const std::string& file ) {
                            absl::StrAppend( out, "./", file );
                        } ),
                    " }" )
                                .c_str() );
            function_node.append_attribute( "voxelFile" )
                .set_value( absl::StrCat( "./", table.voxel_file ).c_str() );
            function_node.append_attribute( "interpolation" )
                .set_value( "nearest" );
        }
    }

    template < typename Model >
    bool GeosExporterImpl< Model >::check_property_name(
        std::string_view property_name ) const
    {
        for( const auto& block : model_.blocks() )
        {
            if( !block.mesh()
                     .polyhedron_attribute_manager()
                     .attribute_ids_matching_name( property_name )
                     .has_value() )
            {
                Logger::info( "The property ", property_name,
                    " will not be exported because it is not defined on "
                    "every block of the model." );
                return false;
            }
        }
        return true;
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_physical_properties()
    {
        for( const auto& [property, geos_name] : PHYSICAL_PROPERTY_GEOS_NAMES )
        {
            if( !model_.has_physical_property( property ) )
            {
                continue;
            }
            const auto& property_info =
                model_.physical_property_info( property );
            auto attribute_name = transfer_block_attribute(
                model_, *model_solid_, model2solid_, property_info );
            imported_fields_.emplace_back(
                std::move( attribute_name ), to_string( geos_name ) );
            if( property == PHYSICAL_PROPERTY_NAME::porosity )
            {
                porosity_attribute_id_ = property_info.attribute_id;
            }
            else
            {
                permeability_attribute_id_ = property_info.attribute_id;
            }
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_mechanics_properties()
    {
        for( const auto property : MECHANICS_PROPERTIES )
        {
            if( !model_.has_physical_property( property ) )
            {
                continue;
            }
            const auto& property_info =
                model_.physical_property_info( property );
            transfer_block_attribute(
                model_, *model_solid_, model2solid_, property_info );
            mechanics_attribute_ids_.emplace(
                property, property_info.attribute_id );
        }
        transfer_initial_stress();
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_initial_stress()
    {
        auto& manager = model_solid_->polyhedron_attribute_manager();
        // Index of the first tensor component given by each property
        for( const auto& [property, name, first_component] :
            { std::tuple{ PHYSICAL_PROPERTY_NAME::initial_stress,
                  "initial_stress", local_index_t{ 0 } },
                std::tuple{ PHYSICAL_PROPERTY_NAME::initial_shear_stress,
                    "initial_shear_stress", local_index_t{ 3 } } } )
        {
            if( !model_.has_physical_property( property ) )
            {
                continue;
            }
            if( !initial_stress_attribute_id_ )
            {
                AttributeValues< StressTensor > values;
                values.default_value = { 0, 0, 0, 0, 0, 0 };
                initial_stress_attribute_id_ =
                    manager.template create_attribute< VariableAttribute,
                        StressTensor >(
                        INITIAL_STRESS_TENSOR, values, AttributeProperties{} );
            }
            auto tensor = manager.template find_attribute< VariableAttribute,
                StressTensor >( initial_stress_attribute_id_.value() );
            transfer_stress_components( model_, model2solid_,
                model_.physical_property_info( property ), name,
                first_component, *tensor );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_boundary_conditions()
    {
        for( const auto& [property, field] : BOUNDARY_CONDITION_GEOS_FIELDS )
        {
            if( !model_.has_physical_property( property ) )
            {
                continue;
            }
            BoundaryCondition prototype;
            prototype.field_name = to_string( field.field_name );
            transfer_boundary_condition(
                boundary_condition_attribute_id( model_, property ),
                field.name_prefix, prototype, 1. );
        }
        for( const auto& [property, fraction] : FLUID_FRACTION_GEOS_NAMES )
        {
            if( !model_.has_physical_property( property ) )
            {
                continue;
            }
            boundary_fraction_phases_.push_back( fraction.phase );
            const auto component =
                boundary_fraction_component( problem_, fraction );
            if( !component )
            {
                // Reported by validate_geos_deck
                continue;
            }
            BoundaryCondition prototype;
            prototype.field_name = "globalCompFraction";
            prototype.component = component;
            transfer_boundary_condition(
                boundary_condition_attribute_id( model_, property ),
                absl::StrCat( "boundaryComposition_", fraction.fluid ),
                prototype, 1. );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_mechanics_boundary_conditions()
    {
        if( model_.has_physical_property(
                PHYSICAL_PROPERTY_NAME::boundary_displacement ) )
        {
            transfer_displacements( boundary_condition_attribute_id(
                model_, PHYSICAL_PROPERTY_NAME::boundary_displacement ) );
        }
        if( model_.has_physical_property(
                PHYSICAL_PROPERTY_NAME::boundary_traction ) )
        {
            transfer_tractions( boundary_condition_attribute_id(
                model_, PHYSICAL_PROPERTY_NAME::boundary_traction ) );
        }
        if( model_.has_physical_property(
                PHYSICAL_PROPERTY_NAME::boundary_normal_traction ) )
        {
            BoundaryCondition prototype;
            prototype.kind = BOUNDARY_KIND::normal_traction;
            // GEOS normal tractions pull the faces outwards
            transfer_boundary_condition(
                boundary_condition_attribute_id(
                    model_, PHYSICAL_PROPERTY_NAME::boundary_normal_traction ),
                "normalTraction", prototype, -1. );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_displacements(
        const uuid& attribute_id )
    {
        absl::linked_hash_map< std::pair< local_index_t, double >,
            std::vector< index_t > >
            value_regions;
        for( const auto& surface : model_.surfaces() )
        {
            const auto values = surface_values< std::array< double, 3 > >(
                surface, attribute_id );
            if( !values )
            {
                continue;
            }
            OpenGeodeGeosciencesIOModelException::check_exception(
                is_uniform( values.value() ), nullptr,
                OpenGeodeException::TYPE::data,
                "[GeosExporter] Displacement boundary condition on Surface ",
                surface.id().string(), " must be the same on every polygon." );
            const auto& displacement = values->front();
            const auto region_id = surface_regions_.at( surface.id() );
            for( const auto axis : LRange{ 3 } )
            {
                if( std::isnan( displacement[axis] ) )
                {
                    continue;
                }
                value_regions[{ axis, displacement[axis] }].push_back(
                    region_id );
                displacement_axes_[axis] = true;
            }
        }
        for( auto& [axis_value, region_ids] : value_regions )
        {
            BoundaryCondition condition;
            condition.name = absl::StrCat( "boundaryDisplacement_",
                AXIS_NAMES[axis_value.first], "_", region_ids.front() );
            condition.field_name = "totalDisplacement";
            condition.component = axis_value.first;
            condition.region_ids = std::move( region_ids );
            condition.scale = axis_value.second;
            condition.object_path = "nodeManager";
            boundary_conditions_.push_back( std::move( condition ) );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_tractions(
        const uuid& attribute_id )
    {
        absl::linked_hash_map< std::array< double, 3 >, std::vector< index_t > >
            vector_regions;
        for( const auto& surface : model_.surfaces() )
        {
            const auto values = surface_values< std::array< double, 3 > >(
                surface, attribute_id );
            if( !values )
            {
                continue;
            }
            const auto region_id = surface_regions_.at( surface.id() );
            if( is_uniform( values.value() ) )
            {
                vector_regions[values->front()].push_back( region_id );
                continue;
            }
            for( const auto axis : LRange{ 3 } )
            {
                const auto name = absl::StrCat(
                    "traction_", AXIS_NAMES[axis], "_surface", region_id );
                BoundaryCondition condition;
                condition.kind = BOUNDARY_KIND::vector_traction;
                condition.direction = std::array< double, 3 >{ 0, 0, 0 };
                condition.direction->at( axis ) = 1;
                condition.function_name = add_boundary_table(
                    surface, axis_values( values.value(), axis ), name );
                condition.name = name;
                condition.region_ids = { region_id };
                // The table stores the actual values, hence a unit scale
                condition.scale = 1.;
                boundary_conditions_.push_back( std::move( condition ) );
            }
        }
        index_t condition_id{ 0 };
        for( auto& [traction, region_ids] : vector_regions )
        {
            const Vector3D vector{ traction };
            const auto norm = vector.length();
            if( norm == 0 )
            {
                continue;
            }
            BoundaryCondition condition;
            condition.name = absl::StrCat( "traction_", condition_id++ );
            condition.kind = BOUNDARY_KIND::vector_traction;
            const auto direction = vector / norm;
            condition.direction = std::array< double, 3 >{ direction.value( 0 ),
                direction.value( 1 ), direction.value( 2 ) };
            condition.region_ids = std::move( region_ids );
            condition.scale = norm;
            boundary_conditions_.push_back( std::move( condition ) );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_boundary_condition(
        const uuid& attribute_id,
        std::string_view name_prefix,
        const BoundaryCondition& prototype,
        double factor )
    {
        absl::linked_hash_map< double, std::vector< index_t > > value_regions;
        for( const auto& surface : model_.surfaces() )
        {
            auto values = surface_values< double >( surface, attribute_id );
            if( !values )
            {
                continue;
            }
            for( auto& value : values.value() )
            {
                value *= factor;
            }
            const auto region_id = surface_regions_.at( surface.id() );
            if( is_uniform( values.value() ) )
            {
                value_regions[values->front()].push_back( region_id );
                continue;
            }
            auto condition = prototype;
            condition.name = absl::StrCat( name_prefix, "_surface", region_id );
            condition.function_name =
                add_boundary_table( surface, values.value(), condition.name );
            condition.region_ids = { region_id };
            // The table stores the actual values, hence a unit scale
            condition.scale = 1.;
            boundary_conditions_.push_back( std::move( condition ) );
        }
        index_t condition_id{ 0 };
        for( auto& [value, region_ids] : value_regions )
        {
            auto condition = prototype;
            condition.name = absl::StrCat( name_prefix, "_", condition_id++ );
            condition.region_ids = std::move( region_ids );
            condition.scale = value;
            boundary_conditions_.push_back( std::move( condition ) );
        }
    }

    template < typename Model >
    auto GeosExporterImpl< Model >::boundary_grid( const Surface3D& surface )
        -> const BoundaryGrid&
    {
        if( const auto grid = boundary_grids_.find( surface.id() );
            grid != boundary_grids_.end() )
        {
            return grid->second;
        }
        const auto& mesh = surface.mesh();
        std::vector< Point3D > centers;
        centers.reserve( mesh.nb_polygons() );
        for( const auto polygon_id : Range{ mesh.nb_polygons() } )
        {
            centers.push_back( mesh.polygon_barycenter( polygon_id ) );
        }
        auto coordinates = separating_table_coordinates( centers );
        OpenGeodeGeosciencesIOModelException::check_exception(
            coordinates.has_value(), nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot build an exact boundary condition table "
            "for Surface ",
            surface.id().string(), ": some polygons share the same center." );
        const auto region_id = surface_regions_.at( surface.id() );
        BoundaryGrid grid;
        for( const auto axis : LRange{ 3 } )
        {
            grid.coordinates[axis] = std::move( ( *coordinates )[axis] );
            grid.files[axis] = absl::StrCat( BOUNDARY_TABLES_DIRECTORY,
                "/surface", region_id, "_", AXIS_NAMES[axis], "lin.geos" );
        }
        grid.node_polygons = table_node_polygons( grid.coordinates,
            table_polygon_nodes( grid.coordinates, centers ), centers );
        return boundary_grids_.emplace( surface.id(), std::move( grid ) )
            .first->second;
    }

    template < typename Model >
    std::string GeosExporterImpl< Model >::add_boundary_table(
        const Surface3D& surface,
        absl::Span< const double > values,
        std::string_view name )
    {
        const auto& grid = boundary_grid( surface );
        BoundaryTable table;
        table.function_name = absl::StrCat( name, "_function" );
        table.surface_id = surface.id();
        table.voxel_file =
            absl::StrCat( BOUNDARY_TABLES_DIRECTORY, "/", name, ".geos" );
        table.values.reserve( grid.node_polygons.size() );
        for( const auto polygon_id : grid.node_polygons )
        {
            table.values.push_back( values[polygon_id] );
        }
        auto function_name = table.function_name;
        boundary_tables_.push_back( std::move( table ) );
        return function_name;
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_boundary_table_files() const
    {
        std::filesystem::create_directories(
            absl::StrCat( files_directory(), "/", BOUNDARY_TABLES_DIRECTORY ) );
        for( const auto& [surface_id, grid] : boundary_grids_ )
        {
            for( const auto axis : LRange{ 3 } )
            {
                write_table_values(
                    absl::StrCat( files_directory(), "/", grid.files[axis] ),
                    grid.coordinates[axis] );
            }
        }
        for( const auto& table : boundary_tables_ )
        {
            write_table_values(
                absl::StrCat( files_directory(), "/", table.voxel_file ),
                table.values );
        }
    }

    template < typename Model >
    std::string GeosExporterImpl< Model >::write_solid_file() const
    {
        const auto filename = absl::StrCat( prefix(), ".vtu" );
        const auto file_vtu = absl::StrCat( files_directory(), "/", filename );
        const auto nb_polyhedra = model_solid_->nb_polyhedra();
        AttributeManager cell_attributes;
        cell_attributes.copy( model_solid_->polyhedron_attribute_manager() );
        cell_attributes.resize( nb_polyhedra + surface_cells_.size() );
        auto region_attribute =
            cell_attributes.find_attribute< VariableAttribute, index_t >(
                region_attribute_id_ );
        for( const auto surface_cell : Indices{ surface_cells_region_ } )
        {
            region_attribute->set_value( nb_polyhedra + surface_cell,
                surface_cells_region_[surface_cell] );
        }
        duplicate_imported_fields( cell_attributes );
        if( const auto* tetra = dynamic_cast< const TetrahedralSolid3D* >(
                model_solid_.get() ) )
        {
            save_geos_vtu< detail::VTUTetrahedralOutputImpl >(
                file_vtu, *tetra, surface_cells_, cell_attributes );
        }
        else if( const auto* hybrid = dynamic_cast< const HybridSolid3D* >(
                     model_solid_.get() ) )
        {
            save_geos_vtu< detail::VTUHybridOutputImpl >(
                file_vtu, *hybrid, surface_cells_, cell_attributes );
        }
        else if( const auto* poly = dynamic_cast< const PolyhedralSolid3D* >(
                     model_solid_.get() ) )
        {
            save_geos_vtu< detail::VTUPolyhedralOutputImpl >(
                file_vtu, *poly, surface_cells_, cell_attributes );
        }
        else
        {
            throw OpenGeodeGeosciencesIOModelException{ nullptr,
                OpenGeodeException::TYPE::data,
                "[Blocks::save_geos] Cannot find the explicit "
                "SolidMesh type" };
        }
        return filename;
    }
    template < typename Model >
    void GeosExporterImpl< Model >::write_well_perforation_file() const
    {
        index_t well_id{ 0 };
        for( const auto& well : well_perforations_ )
        {
            const auto file = absl::StrCat(
                files_directory(), "/", prefix(), "_well", well_id, ".vtp" );
            save_point_set( *well, file );
            well_id++;
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::duplicate_imported_fields(
        AttributeManager& cell_attributes ) const
    {
        if( !writer_ )
        {
            return;
        }
        for( const auto& field : writer_->imported_fields() )
        {
            if( !field.is_copy() )
            {
                continue;
            }
            const auto& source_id = field.source.attribute_id;
            OpenGeodeGeosciencesIOModelException::check_exception(
                cell_attributes.attribute_type( source_id )
                    == typeid( double ).name(),
                nullptr, OpenGeodeException::TYPE::data, "[GeosExporter] ",
                field.source.name,
                " attribute must store double values to be imported "
                "several times." );
            const auto source =
                cell_attributes.find_read_only_attribute< double >( source_id );
            AttributeValues< double > values;
            values.default_value = 0.;
            const auto copy_id =
                cell_attributes
                    .template create_attribute< VariableAttribute, double >(
                        field.vtu_name, values, AttributeProperties{} );
            auto copy =
                cell_attributes.find_attribute< VariableAttribute, double >(
                    copy_id );
            for( const auto cell : Range{ cell_attributes.nb_elements() } )
            {
                copy->set_value( cell, source->value( cell ) );
            }
        }
    }

    template < typename Model >
    GeosDeckContext GeosExporterImpl< Model >::build_deck_context() const
    {
        GeosDeckContext context;
        context.solid_regions = region_definitions_.regions;
        context.regions_by_uuid = region_definitions_.regions_by_uuid;
        context.merged_blocks = region_definitions_.merged_blocks;
        const auto& manager = model_solid_->polyhedron_attribute_manager();
        const auto nb_polyhedra = model_solid_->nb_polyhedra();
        if( porosity_attribute_id_ )
        {
            const auto attribute = manager.find_generic_attribute(
                porosity_attribute_id_.value() );
            context.porosity_array =
                GeosCellArray{ porosity_attribute_id_.value(),
                    attribute->name().value() };
            context.mean_porosity =
                mean_attribute_item( *attribute, nb_polyhedra, 0 );
        }
        if( permeability_attribute_id_ )
        {
            const auto attribute = manager.find_generic_attribute(
                permeability_attribute_id_.value() );
            context.permeability_array =
                GeosCellArray{ permeability_attribute_id_.value(),
                    attribute->name().value() };
            std::array< double, 3 > mean{ 0, 0, 0 };
            const auto is_isotropic = attribute->nb_items() < 3;
            for( const auto axis : LRange{ 3 } )
            {
                mean[axis] = mean_attribute_item(
                    *attribute, nb_polyhedra, is_isotropic ? 0 : axis );
            }
            context.mean_permeability = mean;
        }
        for( const auto& [property, attribute_id] : mechanics_attribute_ids_ )
        {
            const auto [mean, min] =
                attribute_mean_and_min( manager, attribute_id, nb_polyhedra );
            GeosCellProperty cell_property{
                GeosCellArray{ attribute_id,
                    manager.find_generic_attribute( attribute_id )
                        ->name()
                        .value() },
                mean, min
            };
            if( property == PHYSICAL_PROPERTY_NAME::bulk_modulus )
            {
                context.bulk_modulus = std::move( cell_property );
            }
            else if( property == PHYSICAL_PROPERTY_NAME::shear_modulus )
            {
                context.shear_modulus = std::move( cell_property );
            }
            else
            {
                context.density = std::move( cell_property );
            }
        }
        if( initial_stress_attribute_id_ )
        {
            context.initial_stress =
                GeosCellArray{ initial_stress_attribute_id_.value(),
                    INITIAL_STRESS_TENSOR };
        }
        context.displacement_axes = displacement_axes_;
        context.boundary_fraction_phases = boundary_fraction_phases_;
        context.has_boundary_pressure = model_.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_pressure );
        context.has_boundary_temperature = model_.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_temperature );
        context.has_flow_boundary = context.has_boundary_pressure
                                    || context.has_boundary_temperature
                                    || !boundary_fraction_phases_.empty();
        for( const auto property :
            { PHYSICAL_PROPERTY_NAME::boundary_displacement,
                PHYSICAL_PROPERTY_NAME::boundary_traction,
                PHYSICAL_PROPERTY_NAME::boundary_normal_traction } )
        {
            context.has_mechanics_boundary =
                context.has_mechanics_boundary
                || model_.has_physical_property( property );
        }
        context.well_names = well_names_;
        return context;
    }

    template < typename Model >
    void GeosExporterImpl< Model >::validate_problem() const
    {
        auto issues = validate_physical_problem( *problem_ );
        for( auto& issue :
            validate_physical_problem_regions( *problem_, model_ ) )
        {
            issues.push_back( std::move( issue ) );
        }
        for( auto& issue : validate_geos_problem( *problem_, well_names_ ) )
        {
            issues.push_back( std::move( issue ) );
        }
        check_no_issue( issues );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::prepare_simulation()
    {
        context_ = build_deck_context();
        auto issues = validate_geos_deck( *problem_, context_ );
        // Solid mechanics data are cell properties: no material to assign
        GeosRegionAssignment materials;
        GeosRegionAssignment initial_conditions;
        if( problem_->has_flow() )
        {
            materials = assign_geos_materials( *problem_, context_ );
            initial_conditions =
                assign_geos_initial_conditions( *problem_, context_ );
        }
        for( const auto* assignment : { &materials, &initial_conditions } )
        {
            issues.insert( issues.end(), assignment->issues.begin(),
                assignment->issues.end() );
        }
        check_no_issue( issues );
        if( problem_->has_flow() && !context_.has_boundary_pressure )
        {
            Logger::warning(
                "[GeosExporter] The model has no boundary pressure: "
                "nothing drives the flow in the GEOS deck." );
        }
        if( const auto temperature =
                geos_default_boundary_temperature( *problem_, context_ ) )
        {
            add_isothermal_boundary_temperature( temperature.value() );
        }
        writer_ = std::make_unique< GeosSimulationWriter >(
            *problem_, context_, materials, initial_conditions );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::add_isothermal_boundary_temperature(
        double temperature )
    {
        std::vector< index_t > region_ids;
        for( const auto& condition : boundary_conditions_ )
        {
            if( condition.field_name == "pressure" )
            {
                region_ids.insert( region_ids.end(),
                    condition.region_ids.begin(), condition.region_ids.end() );
            }
        }
        absl::c_sort( region_ids );
        region_ids.erase( std::unique( region_ids.begin(), region_ids.end() ),
            region_ids.end() );
        BoundaryCondition condition;
        condition.name = "boundaryTemperature_isothermal";
        condition.field_name = "temperature";
        condition.region_ids = std::move( region_ids );
        condition.scale = temperature;
        boundary_conditions_.push_back( std::move( condition ) );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_simulation_sections(
        pugi::xml_node& problem_node ) const
    {
        writer_->write_physics( problem_node );
        auto field_specifications =
            problem_node.append_child( "FieldSpecifications" );
        writer_->write_initial_conditions( field_specifications );
        writer_->write_source_terms( field_specifications );
        write_boundary_conditions( field_specifications );
        if( writer_->has_functions() || !boundary_tables_.empty() )
        {
            auto functions = problem_node.append_child( "Functions" );
            writer_->write_functions( functions );
            write_boundary_functions( functions );
        }
        writer_->write_events_and_outputs( problem_node );
        writer_->apply_solver_options( problem_node );
    }

    template class opengeode_geosciencesio_model_api GeosExporterImpl< BRep >;
    template class opengeode_geosciencesio_model_api
        GeosExporterImpl< StructuralModel >;
} // namespace geode::internal