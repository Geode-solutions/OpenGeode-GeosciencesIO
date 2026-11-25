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
#include <geode/geometry/basic_objects/tetrahedron.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/distance.hpp>
#include <geode/geometry/information.hpp>
#include <geode/geometry/nn_search.hpp>
#include <geode/geometry/point.hpp>
#include <geode/geometry/position.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/geosciences_io/model/common.hpp>

#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/polyhedral_solid.hpp>
#include <geode/mesh/core/regular_grid_solid.hpp>
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

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/solid_mechanics.hpp>
#include <geode/geosciences/simulation/helpers/simulation_report.hpp>
#include <geode/geosciences_io/model/internal/geos_boundary_tables.hpp>
#include <geode/geosciences_io/model/internal/geos_export.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_checker.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_rules.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_writer.hpp>
#include <geode/geosciences_io/model/internal/geos_xml.hpp>
#include <geode/geosciences_io/model/internal/simulation_boundary_values.hpp>
#include <geode/geosciences_io/model/internal/simulation_cell_properties.hpp>
#include <geode/geosciences_io/model/internal/simulation_export.hpp>
#include <geode/geosciences_io/model/internal/simulation_vtu.hpp>
#include <geode/geosciences_io/model/internal/simulation_xml.hpp>
#include <geode/model/representation/core/brep.hpp>

namespace
{
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

    // Mean of a double cell attribute over the first elements, in double
    // precision
    double double_attribute_mean( const geode::AttributeManager& manager,
        const geode::uuid& attribute_id,
        geode::index_t nb_elements )
    {
        const auto attribute =
            manager.find_read_only_attribute< double >( attribute_id );
        double sum{ 0 };
        for( const auto element : geode::Range{ nb_elements } )
        {
            sum += attribute->value( element );
        }
        return nb_elements == 0 ? 0. : sum / nb_elements;
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
    namespace
    {
        const SimulationNaming GEOS_NAMING{ "[GeosExporter]", "GEOS",
            &is_valid_geos_list_item, GEOS_LIST_SEPARATORS };
    } // namespace

    const absl::linked_hash_map< PHYSICAL_PROPERTY_NAME, std::string_view >
        PHYSICAL_PROPERTY_GEOS_NAMES{
            { PHYSICAL_PROPERTY_NAME::permeability, "rockPerm_permeability" },
            { PHYSICAL_PROPERTY_NAME::porosity,
                "rockPorosity_referencePorosity" },
        };

    constexpr std::array< std::pair< PHYSICAL_PROPERTY_NAME, std::string_view >,
        6 >
        MECHANICS_PROPERTIES{ {
            { PHYSICAL_PROPERTY_NAME::bulk_modulus, "bulkModulus" },
            { PHYSICAL_PROPERTY_NAME::shear_modulus, "shearModulus" },
            { PHYSICAL_PROPERTY_NAME::density, "density" },
            { PHYSICAL_PROPERTY_NAME::grain_density, "density" },
            { PHYSICAL_PROPERTY_NAME::biot_coefficient, "biotCoefficient" },
            { PHYSICAL_PROPERTY_NAME::thermal_expansion_coefficient,
                "thermalExpansionCoefficient" },
        } };

    std::optional< GeosCellProperty >& mechanics_context_property(
        GeosDeckContext& context, PHYSICAL_PROPERTY_NAME property )
    {
        switch( property )
        {
            case PHYSICAL_PROPERTY_NAME::bulk_modulus:
                return context.bulk_modulus;
            case PHYSICAL_PROPERTY_NAME::shear_modulus:
                return context.shear_modulus;
            case PHYSICAL_PROPERTY_NAME::biot_coefficient:
                return context.biot_coefficient;
            case PHYSICAL_PROPERTY_NAME::thermal_expansion_coefficient:
                return context.thermal_expansion;
            case PHYSICAL_PROPERTY_NAME::density:
            case PHYSICAL_PROPERTY_NAME::grain_density:
                return context.density;
            default:
                throw OpenGeodeGeosciencesIOModelException{ nullptr,
                    OpenGeodeException::TYPE::internal,
                    "[GeosExporter] Physical property ",
                    static_cast< int >( property ),
                    " is not a mechanics property of the GEOS deck" };
        }
    }

    constexpr auto INITIAL_STRESS_TENSOR = "initial_stress_tensor";

    constexpr auto GRAIN_BULK_MODULUS = "grain_bulk_modulus";

    // GEOS grain bulk modulus of incompressible grains (Biot coefficient 1)
    constexpr double INCOMPRESSIBLE_GRAIN_BULK_MODULUS{ 1e27 };

    template < typename Model >
    GeosExporterImpl< Model >::GeosExporterImpl(
        std::string_view directory, const Model& exported_model )
        : SimulationExporterImpl< Model >{ directory, exported_model,
              "[GeosExporter]" }
    {
    }

    template < typename Model >
    void GeosExporterImpl< Model >::convert_model()
    {
        bool all_regular_grids{ model().nb_blocks() != 0 };
        for( const auto& block : model().blocks() )
        {
            if( block.mesh().type_name() != RegularGrid3D::type_name_static() )
            {
                all_regular_grids = false;
                break;
            }
        }
        OpenGeodeGeosciencesIOModelException::check_exception(
            !all_regular_grids, nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot export a model whose Blocks are all "
            "RegularGrids: convert the Block meshes into HybridSolids "
            "first." );
        solid_.emplace( model(), region_definitions_, "[GeosExporter]" );
    }

    template < typename Model >
    GeosExporterImpl< Model >::~GeosExporterImpl() = default;

    template < typename Model >
    void GeosExporterImpl< Model >::write_files() const
    {
        pugi::xml_document doc_xml;
        auto pb_node = doc_xml.append_child( "Problem" );

        auto mesh_node = pb_node.append_child( "Mesh" );
        write_mesh_files( mesh_node );

        if( !well_perforations().empty() )
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
        OpenGeodeGeosciencesIOModelException::check_exception(
            doc_xml.save_file( filename_xml.c_str(), PUGIXML_TEXT( "    " ),
                pugi::format_indent_attributes ),
            nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot write file ", filename_xml );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::check_well_name(
        std::string_view name ) const
    {
        OpenGeodeGeosciencesIOModelException::check_exception(
            is_valid_geos_list_item( name ), nullptr,
            OpenGeodeException::TYPE::data,
            "[GeosExporter] Well perforations name \"", name,
            "\" must be non empty and contain no space, comma nor brace." );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::prepare_export()
    {
        clear_export();
        define_regions();
        if( problem() != nullptr )
        {
            check_problem();
        }
        convert_model();
        initialize_surface_cells( nb_solid_regions() );
        transfer_physical_properties();
        transfer_boundary_conditions();
        if( problem() != nullptr && problem()->solid_mechanics() )
        {
            transfer_mechanics_properties();
            transfer_mechanics_boundary_conditions();
        }
        compute_well_boxes();
        if( problem() != nullptr )
        {
            prepare_simulation();
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::clear_export()
    {
        surface_cells_.clear();
        surface_cells_region_.clear();
        surface_regions_.clear();
        imported_fields_.clear();
        boundary_conditions_.clear();
        boundary_grids_.clear();
        boundary_tables_.clear();
        porosity_attribute_id_.reset();
        permeability_attribute_id_.reset();
        mechanics_attribute_ids_.clear();
        initial_stress_attribute_id_.reset();
        grain_bulk_modulus_attribute_id_.reset();
        writer_.reset();
        solid_.reset();
    }

    template < typename Model >
    void GeosExporterImpl< Model >::define_regions()
    {
        region_definitions_ = create_region_definitions( model() );
        add_block_collection_regions( model(), region_definitions_ );
    }

    template < typename Model >
    index_t GeosExporterImpl< Model >::nb_solid_regions() const
    {
        return static_cast< index_t >(
            region_definitions_.solid_regions.size() );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::initialize_surface_cells(
        index_t first_surface_region_id )
    {
        const auto& polygons_mapping =
            solid_->mappings().surface_polygons_mapping;
        auto region_id = first_surface_region_id;
        for( const auto& surface : model().surfaces() )
        {
            index_t nb_skipped_polygons{ 0 };
            for( const auto polygon_id : Range{ surface.mesh().nb_polygons() } )
            {
                // Only polygons on a solid facet define GEOS faces
                if( !polygons_mapping.has_mapping_input(
                        { surface.id(), polygon_id } ) )
                {
                    nb_skipped_polygons++;
                    continue;
                }
                auto vertices =
                    polygon_unique_vertices( model(), surface, polygon_id );
                for( auto& vertex : vertices )
                {
                    vertex = solid_->mappings().unique_vertices_mapping.in2out(
                        vertex );
                }
                surface_cells_.emplace_back( std::move( vertices ) );
                surface_cells_region_.push_back( region_id );
            }
            if( nb_skipped_polygons != 0 )
            {
                Logger::warning( "[GeosExporter] ", nb_skipped_polygons,
                    " polygons of Surface ", surface.id().string(),
                    " are not on a Block polyhedron facet, they are not "
                    "exported." );
            }
            surface_regions_.emplace( surface.id(), region_id );
            region_id++;
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::compute_well_boxes()
    {
        well_boxes_.clear();
        if( well_perforations().empty() )
        {
            return;
        }
        const auto aabb = create_aabb_tree( solid_->mesh() );
        for( const auto well_id : Indices{ well_perforations() } )
        {
            const auto& well = well_perforations()[well_id];
            const auto& well_name = well_names()[well_id];
            BoundingBox3D perf_box;
            bool has_perforated_cell{ false };
            for( const auto point_id : Range( well->nb_vertices() ) )
            {
                const auto& point = well->point( point_id );
                const auto cell_id =
                    containing_cell( solid_->mesh(), aabb, point );
                if( cell_id == NO_ID )
                {
                    Logger::warning( "[GeosExporter] Perforation ", point_id,
                        " of well ", well_name, " (", point.string(),
                        ") is not inside the model, it is ignored." );
                    continue;
                }
                has_perforated_cell = true;
                for( const auto vertex_id :
                    solid_->mesh().polyhedron_vertices( cell_id ) )
                {
                    perf_box.add_point( solid_->mesh().point( vertex_id ) );
                }
            }
            if( !has_perforated_cell )
            {
                Logger::warning( "[GeosExporter] No perforation of well ",
                    well_name,
                    " is inside the model, no perforation Box is exported." );
                well_boxes_.emplace_back( std::nullopt );
                continue;
            }
            well_boxes_.emplace_back( std::move( perf_box ) );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::write_well_perforations_boxes(
        pugi::xml_node& root ) const
    {
        for( const auto well_id : Indices{ well_boxes_ } )
        {
            const auto& perf_box = well_boxes_[well_id];
            if( !perf_box )
            {
                continue;
            }
            auto box_node = root.append_child( "Box" );
            box_node.append_attribute( "name" ).set_value(
                well_names()[well_id].c_str() );
            static constexpr auto SAFETY_OFFSET = 100. * GLOBAL_EPSILON;
            const auto box_corner = []( const Point3D& corner, double offset ) {
                return absl::StrCat( "{",
                    geos_real( corner.value( 0 ) + offset ), ", ",
                    geos_real( corner.value( 1 ) + offset ), ", ",
                    geos_real( corner.value( 2 ) + offset ), "}" );
            };
            box_node.append_attribute( "xMin" ).set_value(
                box_corner( perf_box->min(), -SAFETY_OFFSET ).c_str() );
            box_node.append_attribute( "xMax" ).set_value(
                box_corner( perf_box->max(), SAFETY_OFFSET ).c_str() );
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
                specification_node.append_attribute( "direction" )
                    .set_value(
                        geos_array( condition.direction.value() ).c_str() );
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
                geos_real( condition.scale ).c_str() );
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
    void GeosExporterImpl< Model >::transfer_physical_properties()
    {
        for( const auto& [property, geos_name] : PHYSICAL_PROPERTY_GEOS_NAMES )
        {
            if( !model().has_physical_property( property ) )
            {
                continue;
            }
            const auto& property_info =
                model().physical_property_info( property );
            auto attribute_name = solid_->transfer_block_property(
                GEOS_NAMING, model(), property_info, geos_name );
            if( !attribute_name )
            {
                continue;
            }
            imported_fields_.emplace_back(
                std::move( attribute_name.value() ), to_string( geos_name ) );
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
        for( const auto& [property, geos_name] : MECHANICS_PROPERTIES )
        {
            if( !model().has_physical_property( property )
                || !is_physical_property_read( *problem(), property ) )
            {
                continue;
            }
            const auto& property_info =
                model().physical_property_info( property );
            if( !solid_->transfer_block_property(
                    GEOS_NAMING, model(), property_info, geos_name ) )
            {
                continue;
            }
            mechanics_attribute_ids_.emplace(
                property, property_info.attribute_id );
        }
        transfer_initial_stress();
        if( problem()->has_flow() )
        {
            compute_grain_bulk_modulus();
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_initial_stress()
    {
        const auto factor = stress_factor( *problem(), true );
        // Index of the first tensor component given by each property
        for( const auto& [property, name, first_component] :
            { std::tuple{ PHYSICAL_PROPERTY_NAME::initial_stress,
                  "initial_stress", local_index_t{ 0 } },
                std::tuple{ PHYSICAL_PROPERTY_NAME::initial_shear_stress,
                    "initial_shear_stress", local_index_t{ 3 } } } )
        {
            if( !model().has_physical_property( property ) )
            {
                continue;
            }
            transfer_stress_components( GEOS_NAMING, model(),
                solid_->mappings(), model().physical_property_info( property ),
                name, first_component, factor, initial_stress_tensor() );
        }
        // Created even without initial stress property: build_deck_context
        // imports it, then add_initial_pore_pressure shifts it
        if( problem()->has_flow() )
        {
            initial_stress_tensor();
        }
    }

    template < typename Model >
    VariableAttribute< StressTensor >&
        GeosExporterImpl< Model >::initial_stress_tensor()
    {
        auto& manager = solid_->mesh().polyhedron_attribute_manager();
        if( !initial_stress_attribute_id_ )
        {
            AttributeValues< StressTensor > values;
            values.default_value = { 0, 0, 0, 0, 0, 0 };
            initial_stress_attribute_id_ =
                manager.template create_attribute< VariableAttribute,
                    StressTensor >(
                    INITIAL_STRESS_TENSOR, values, AttributeProperties{} );
        }
        return *manager
                    .template find_attribute< VariableAttribute, StressTensor >(
                        initial_stress_attribute_id_.value() );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::compute_grain_bulk_modulus()
    {
        auto& manager = solid_->mesh().polyhedron_attribute_manager();
        const auto bulk_modulus =
            manager.template find_read_only_attribute< double >(
                mechanics_attribute_ids_.at(
                    PHYSICAL_PROPERTY_NAME::bulk_modulus ) );
        const auto biot_coefficient =
            manager.template find_read_only_attribute< double >(
                mechanics_attribute_ids_.at(
                    PHYSICAL_PROPERTY_NAME::biot_coefficient ) );
        const auto name = unique_attribute_name(
            GEOS_NAMING, manager, GRAIN_BULK_MODULUS, "grainBulkModulus" );
        AttributeValues< double > values;
        values.default_value = 0.;
        grain_bulk_modulus_attribute_id_ =
            manager.template create_attribute< VariableAttribute, double >(
                name, values, AttributeProperties{} );
        auto grain_bulk_modulus =
            manager.template find_attribute< VariableAttribute, double >(
                grain_bulk_modulus_attribute_id_.value() );
        for( const auto polyhedron : Range{ solid_->mesh().nb_polyhedra() } )
        {
            const auto biot = biot_coefficient->value( polyhedron );
            grain_bulk_modulus->set_value( polyhedron,
                biot >= 1.
                    ? INCOMPRESSIBLE_GRAIN_BULK_MODULUS
                    : bulk_modulus->value( polyhedron ) / ( 1. - biot ) );
        }
    }

    // GEOS poromechanics stores the stress of the solid, whose total stress
    // is shifted by the pore pressure: solid = total + b.p
    template < typename Model >
    void GeosExporterImpl< Model >::add_initial_pore_pressure(
        const RegionAssignment& initial_conditions )
    {
        const auto biot_coefficient =
            solid_->mesh()
                .polyhedron_attribute_manager()
                .template find_read_only_attribute< double >(
                    mechanics_attribute_ids_.at(
                        PHYSICAL_PROPERTY_NAME::biot_coefficient ) );
        std::vector< double > region_pressures;
        region_pressures.reserve( initial_conditions.region_items.size() );
        for( const auto& condition_id : initial_conditions.region_items )
        {
            region_pressures.push_back(
                dynamic_cast< const UniformInitialCondition& >(
                    problem()->initial_condition( condition_id.value() ) )
                    .pressure() );
        }
        auto& tensor = initial_stress_tensor();
        for( const auto polyhedron : Range{ solid_->mesh().nb_polyhedra() } )
        {
            const auto shift = biot_coefficient->value( polyhedron )
                               * region_pressures[solid_->region( polyhedron )];
            tensor.modify_value(
                polyhedron, [shift]( StressTensor& tensor_value ) {
                    for( const auto axis : LRange{ 3 } )
                    {
                        tensor_value[axis] += shift;
                    }
                } );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_boundary_conditions()
    {
        for( const auto& [property, field] : BOUNDARY_CONDITION_GEOS_FIELDS )
        {
            if( !model().has_physical_property( property ) )
            {
                continue;
            }
            BoundaryCondition prototype;
            prototype.field_name = to_string( field.field_name );
            transfer_boundary_condition(
                boundary_condition_attribute_id( model(), property ),
                field.name_prefix, prototype, 1. );
        }
        for( const auto& [property, fraction] : FLUID_FRACTION_GEOS_NAMES )
        {
            if( !model().has_physical_property( property ) )
            {
                continue;
            }
            const auto component =
                boundary_fraction_component( problem(), fraction );
            if( !component )
            {
                // Reported by GeosSimulationChecker
                continue;
            }
            BoundaryCondition prototype;
            prototype.field_name = "globalCompFraction";
            prototype.component = component;
            transfer_boundary_condition(
                boundary_condition_attribute_id( model(), property ),
                absl::StrCat( "boundaryComposition_", fraction.fluid ),
                prototype, 1. );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_mechanics_boundary_conditions()
    {
        if( model().has_physical_property(
                PHYSICAL_PROPERTY_NAME::boundary_displacement ) )
        {
            transfer_displacements( boundary_condition_attribute_id(
                model(), PHYSICAL_PROPERTY_NAME::boundary_displacement ) );
        }
        if( model().has_physical_property(
                PHYSICAL_PROPERTY_NAME::boundary_traction ) )
        {
            transfer_tractions( boundary_condition_attribute_id(
                model(), PHYSICAL_PROPERTY_NAME::boundary_traction ) );
        }
        if( model().has_physical_property(
                PHYSICAL_PROPERTY_NAME::boundary_normal_traction ) )
        {
            BoundaryCondition prototype;
            prototype.kind = BOUNDARY_KIND::normal_traction;
            // GEOS normal tractions pull the faces outwards
            transfer_boundary_condition(
                boundary_condition_attribute_id(
                    model(), PHYSICAL_PROPERTY_NAME::boundary_normal_traction ),
                "boundaryNormalTraction", prototype,
                stress_factor( *problem(), true ) );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_displacements(
        const uuid& attribute_id )
    {
        absl::linked_hash_map< std::pair< local_index_t, double >,
            std::vector< index_t > >
            value_regions;
        for( const auto& surface : model().surfaces() )
        {
            const auto values = surface_values< std::array< double, 3 > >(
                surface, attribute_id, GEOS_NAMING.context );
            if( !values )
            {
                continue;
            }
            // Uniform on the Surface: checked by
            // SimulationChecker::check_physical_problem_properties
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
            }
        }
        std::array< index_t, 3 > axis_condition_ids{ 0, 0, 0 };
        for( auto& [axis_value, region_ids] : value_regions )
        {
            const auto axis = axis_value.first;
            BoundaryCondition condition;
            condition.name = absl::StrCat( "boundaryDisplacement_",
                AXIS_NAMES[axis], "_", axis_condition_ids[axis]++ );
            condition.field_name = "totalDisplacement";
            condition.component = axis;
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
        auto groups = group_surface_values< std::array< double, 3 > >(
            model(), surface_regions_, attribute_id, GEOS_NAMING.context );
        for( const auto& variable : groups.variable_surfaces )
        {
            for( const auto axis : LRange{ 3 } )
            {
                const auto values = axis_values( variable.values, axis );
                if( absl::c_all_of( values, []( double value ) {
                        return value == 0;
                    } ) )
                {
                    continue;
                }
                BoundaryCondition condition;
                condition.name = absl::StrCat( "boundaryTraction_",
                    AXIS_NAMES[axis], "_surface", variable.region_id );
                condition.kind = BOUNDARY_KIND::vector_traction;
                condition.direction = std::array< double, 3 >{ 0, 0, 0 };
                condition.direction->at( axis ) = 1;
                condition.function_name = add_boundary_table(
                    *variable.surface, values, condition.name );
                condition.region_ids = { variable.region_id };
                // The table stores the actual values, hence a unit scale
                condition.scale = 1.;
                boundary_conditions_.push_back( std::move( condition ) );
            }
        }
        index_t condition_id{ 0 };
        for( auto& [traction, region_ids] : groups.uniform_regions )
        {
            const Vector3D vector{ traction };
            const auto norm = vector.length();
            if( norm == 0 )
            {
                continue;
            }
            BoundaryCondition condition;
            condition.name =
                absl::StrCat( "boundaryTraction_", condition_id++ );
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
        auto groups = group_surface_values< double >(
            model(), surface_regions_, attribute_id, GEOS_NAMING.context );
        for( auto& variable : groups.variable_surfaces )
        {
            for( auto& value : variable.values )
            {
                value *= factor;
            }
            auto condition = prototype;
            condition.name =
                absl::StrCat( name_prefix, "_surface", variable.region_id );
            condition.function_name = add_boundary_table(
                *variable.surface, variable.values, condition.name );
            condition.region_ids = { variable.region_id };
            // The table stores the actual values, hence a unit scale
            condition.scale = 1.;
            boundary_conditions_.push_back( std::move( condition ) );
        }
        index_t condition_id{ 0 };
        for( auto& [value, region_ids] : groups.uniform_regions )
        {
            auto condition = prototype;
            condition.name = absl::StrCat( name_prefix, "_", condition_id++ );
            condition.region_ids = std::move( region_ids );
            condition.scale = factor * value;
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
            centers.push_back( face_center( mesh, polygon_id ) );
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
        Logger::info( "[GeosExporter] Boundary condition table of Surface ",
            surface.id().string(), ": ", grid.coordinates[0].size(), " x ",
            grid.coordinates[1].size(), " x ", grid.coordinates[2].size(),
            " nodes" );
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
        const auto nb_polyhedra = solid_->mesh().nb_polyhedra();
        AttributeManager cell_attributes;
        cell_attributes.copy( solid_->mesh().polyhedron_attribute_manager() );
        cell_attributes.resize( nb_polyhedra + surface_cells_.size() );
        auto region_attribute =
            cell_attributes.find_attribute< VariableAttribute, index_t >(
                solid_->region_attribute_id() );
        for( const auto surface_cell : Indices{ surface_cells_region_ } )
        {
            region_attribute->set_value( nb_polyhedra + surface_cell,
                surface_cells_region_[surface_cell] );
        }
        duplicate_imported_fields( cell_attributes );
        if( !save_simulation_vtu(
                file_vtu, solid_->mesh(), surface_cells_, cell_attributes ) )
        {
            throw OpenGeodeGeosciencesIOModelException{ nullptr,
                OpenGeodeException::TYPE::data,
                "[GeosExporter] Cannot write the model solid of type ",
                solid_->mesh().type_name().get(),
                ": only TetrahedralSolid, HybridSolid and PolyhedralSolid are "
                "supported." };
        }
        return filename;
    }
    template < typename Model >
    void GeosExporterImpl< Model >::write_well_perforation_file() const
    {
        index_t well_id{ 0 };
        for( const auto& well : well_perforations() )
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
            OpenGeodeGeosciencesIOModelException::check_exception(
                !cell_attributes.attribute_ids_matching_name( field.vtu_name ),
                nullptr, OpenGeodeException::TYPE::data,
                "[GeosExporter] Cannot import a copy of ", field.source.name,
                " named ", field.vtu_name,
                ": another exported attribute has this name." );
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
        context.regions = region_definitions_;
        const auto& manager = solid_->mesh().polyhedron_attribute_manager();
        const auto nb_polyhedra = solid_->mesh().nb_polyhedra();
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
        const auto cell_property = [&manager, nb_polyhedra](
                                       const uuid& attribute_id ) {
            return GeosCellProperty{
                GeosCellArray{ attribute_id,
                    manager.find_generic_attribute( attribute_id )
                        ->name()
                        .value() },
                double_attribute_mean( manager, attribute_id, nb_polyhedra )
            };
        };
        for( const auto& [property, attribute_id] : mechanics_attribute_ids_ )
        {
            mechanics_context_property( context, property ) =
                cell_property( attribute_id );
        }
        if( grain_bulk_modulus_attribute_id_ )
        {
            context.grain_bulk_modulus =
                cell_property( grain_bulk_modulus_attribute_id_.value() );
        }
        if( initial_stress_attribute_id_ )
        {
            context.initial_stress =
                GeosCellArray{ initial_stress_attribute_id_.value(),
                    INITIAL_STRESS_TENSOR };
        }
        context.boundary_conditions = model_boundary_conditions( model() );
        // Only the wells perforating the model define a GEOS set
        for( const auto well_id : Indices{ well_boxes_ } )
        {
            if( well_boxes_[well_id] )
            {
                context.well_names.push_back( well_names()[well_id] );
            }
        }
        return context;
    }

    template < typename Model >
    void GeosExporterImpl< Model >::check_problem() const
    {
        const auto report = GeosSimulationChecker{
            *problem(), model(),
            GeosExportData{ region_definitions_,
                { well_names().begin(), well_names().end() } }
        }.check();
        OpenGeodeGeosciencesIOModelException::check_exception(
            report.is_simulable(), nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot write the GEOS simulation deck.\n",
            report.string() );
        for( const auto& issue : report.issues() )
        {
            Logger::warning( "[GeosExporter] ", issue.string() );
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::prepare_simulation()
    {
        context_ = build_deck_context();
        // Assignments are checked by GeosSimulationChecker
        // Solid mechanics data are cell properties: no material to assign
        RegionAssignment materials;
        RegionAssignment initial_conditions;
        if( problem()->has_flow() )
        {
            materials = assign_region_materials(
                *problem(), model(), context_.regions, "GEOS" );
            initial_conditions = assign_region_initial_conditions(
                *problem(), model(), context_.regions, "GEOS" );
        }
        if( problem()->has_flow() && problem()->solid_mechanics() )
        {
            add_initial_pore_pressure( initial_conditions );
        }
        if( const auto temperature = geos_default_boundary_temperature(
                *problem(), context_.boundary_conditions ) )
        {
            add_isothermal_boundary_temperature( temperature.value() );
        }
        writer_ = std::make_unique< GeosSimulationWriter >(
            *problem(), context_, materials, initial_conditions );
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