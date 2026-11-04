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

#include <geode/io/mesh/detail/vtu_solid_output_impl.hpp>

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/solid_mechanics.hpp>
#include <geode/geosciences/simulation/helpers/physical_problem_validation.hpp>
#include <geode/geosciences_io/model/internal/geos_boundary_tables.hpp>
#include <geode/geosciences_io/model/internal/geos_cell_properties.hpp>
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

    // Directory without trailing separator nor dot segments
    std::filesystem::path normalized_directory( std::string_view directory )
    {
        auto path = std::filesystem::path{ geode::to_string( directory ) }
                        .lexically_normal();
        if( !path.empty() && !path.has_filename() )
        {
            path = path.parent_path();
        }
        return path;
    }

    // The polyhedron is split in tetrahedra joining its barycenter to a fan
    // triangulation of its facets
    bool is_point_in_polyhedron( const geode::SolidMesh3D& solid,
        geode::index_t polyhedron_id,
        const geode::Point3D& point )
    {
        const auto barycenter = solid.polyhedron_barycenter( polyhedron_id );
        for( const auto facet :
            geode::LRange{ solid.nb_polyhedron_facets( polyhedron_id ) } )
        {
            const auto vertices =
                solid.polyhedron_facet_vertices( { polyhedron_id, facet } );
            const auto& first = solid.point( vertices[0] );
            for( const auto vertex : geode::LRange{ 1, vertices.size() - 1 } )
            {
                const geode::Tetrahedron tetrahedron{ first,
                    solid.point( vertices[vertex] ),
                    solid.point( vertices[vertex + 1] ), barycenter };
                const auto position =
                    geode::point_tetrahedron_position( point, tetrahedron );
                if( position != geode::POSITION::outside
                    && position != geode::POSITION::parallel )
                {
                    return true;
                }
            }
        }
        return false;
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

    // Mechanics properties of the other kind of problem are ignored
    bool is_read_by(
        const PhysicalProblem& problem, PHYSICAL_PROPERTY_NAME property )
    {
        switch( property )
        {
            case PHYSICAL_PROPERTY_NAME::density:
                return !problem.has_flow();
            case PHYSICAL_PROPERTY_NAME::grain_density:
            case PHYSICAL_PROPERTY_NAME::biot_coefficient:
                return problem.has_flow();
            case PHYSICAL_PROPERTY_NAME::thermal_expansion_coefficient:
                return problem.has_flow() && problem.thermal().has_value();
            default:
                return true;
        }
    }

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
        std::string_view files_directory, const Model& model )
        : model_( model ),
          files_directory_{ normalized_directory( files_directory ).string() },
          prefix_{ normalized_directory( files_directory ).stem().string() }
    {
        OpenGeodeGeosciencesIOModelException::check_exception( !prefix_.empty(),
            nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot export in directory \"", files_directory,
            "\": its name is used as file prefix and must not be empty." );
        std::filesystem::create_directories( files_directory_ );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::convert_model()
    {
        bool all_regular_grids{ model_.nb_blocks() != 0 };
        for( const auto& block : model_.blocks() )
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
        std::tie( model_solid_, model2solid_ ) =
            convert_brep_into_solid( model_ );
        compute_block_polyhedra_offsets();
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
    void GeosExporterImpl< Model >::compute_block_polyhedra_offsets()
    {
        const auto& polyhedra_mapping = model2solid_.solid_polyhedra_mapping;
        block_polyhedra_offsets_.clear();
        block_polyhedra_offsets_.reserve( model_.nb_blocks() + 1 );
        index_t offset{ 0 };
        for( const auto& block : model_.blocks() )
        {
            const auto nb_polyhedra = block.mesh().nb_polyhedra();
            OpenGeodeGeosciencesIOModelException::check_exception(
                nb_polyhedra == 0
                    || polyhedra_mapping.in2out( { block.id(), 0 } ) == offset,
                nullptr, OpenGeodeException::TYPE::data,
                "[GeosExporter] Polyhedra of Block ", block.id().string(),
                " are not contiguous in the model solid." );
            block_polyhedra_offsets_.push_back( offset );
            offset += nb_polyhedra;
        }
        OpenGeodeGeosciencesIOModelException::check_exception(
            offset == model_solid_->nb_polyhedra(), nullptr,
            OpenGeodeException::TYPE::data,
            "[GeosExporter] Wrong number of polyhedra in the model solid." );
        block_polyhedra_offsets_.push_back( offset );
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
        OpenGeodeGeosciencesIOModelException::check_exception(
            doc_xml.save_file( filename_xml.c_str(), PUGIXML_TEXT( "    " ),
                pugi::format_indent_attributes ),
            nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Cannot write file ", filename_xml );
    }

    template < typename Model >
    void GeosExporterImpl< Model >::add_well_perforations(
        const PointSet3D& perforations, std::string_view name )
    {
        OpenGeodeGeosciencesIOModelException::check_exception(
            is_valid_geos_list_item( name ), nullptr,
            OpenGeodeException::TYPE::data,
            "[GeosExporter] Well perforations name \"", name,
            "\" must be non empty and contain no space, comma nor brace." );
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
        clear_export();
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
        compute_well_boxes();
        if( problem_ != nullptr )
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
        boundary_fraction_phases_.clear();
        porosity_attribute_id_.reset();
        permeability_attribute_id_.reset();
        mechanics_attribute_ids_.clear();
        initial_stress_attribute_id_.reset();
        grain_bulk_modulus_attribute_id_.reset();
        writer_.reset();
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
        index_t block_index{ 0 };
        for( const auto& block : model_.blocks() )
        {
            const auto region = block_regions.find( block.id() );
            OpenGeodeGeosciencesIOModelException::check_exception(
                region != block_regions.end(), nullptr,
                OpenGeodeException::TYPE::internal, "[GeosExporter] Block ",
                block.id().string(), " belongs to no region." );
            for( const auto polyhedron_id :
                Range{ block_polyhedra_offsets_[block_index],
                    block_polyhedra_offsets_[block_index + 1] } )
            {
                region_attribute_->set_value( polyhedron_id, region->second );
            }
            block_index++;
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
        const auto& polygons_mapping = model2solid_.surface_polygons_mapping;
        auto region_id = first_surface_region_id;
        for( const auto& surface : model_.surfaces() )
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
                    polygon_unique_vertices( model_, surface, polygon_id );
                for( auto& vertex : vertices )
                {
                    vertex =
                        model2solid_.unique_vertices_mapping.in2out( vertex );
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
        if( well_perforations_.empty() )
        {
            return;
        }
        const auto aabb = create_aabb_tree( *model_solid_ );
        for( const auto well_id : Indices{ well_perforations_ } )
        {
            const auto& well = well_perforations_[well_id];
            const auto& well_name = well_names_[well_id];
            BoundingBox3D perf_box;
            bool has_perforated_cell{ false };
            for( const auto point_id : Range( well->nb_vertices() ) )
            {
                const auto& point = well->point( point_id );
                const auto cell_id = containing_cell( aabb, point );
                if( cell_id == NO_ID )
                {
                    Logger::warning( "[GeosExporter] Perforation ", point_id,
                        " of well ", well_name, " (", point.string(),
                        ") is not inside the model, it is ignored." );
                    continue;
                }
                has_perforated_cell = true;
                for( const auto vertex_id :
                    model_solid_->polyhedron_vertices( cell_id ) )
                {
                    perf_box.add_point( model_solid_->point( vertex_id ) );
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
                well_names_[well_id].c_str() );
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
    index_t GeosExporterImpl< Model >::containing_cell(
        const AABBTree3D& aabb, const Point3D& point ) const
    {
        double distance_to_nearest_cell_center{
            std::numeric_limits< double >::max()
        };
        index_t selected_cell_id{ NO_ID };
        for( const auto cell_id : aabb.containing_boxes( point ) )
        {
            if( !is_point_in_polyhedron( *model_solid_, cell_id, point ) )
            {
                continue;
            }
            const auto distance = point_point_distance< 3 >(
                point, model_solid_->polyhedron_barycenter( cell_id ) );
            if( distance < distance_to_nearest_cell_center )
            {
                distance_to_nearest_cell_center = distance;
                selected_cell_id = cell_id;
            }
        }
        return selected_cell_id;
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
            if( !model_.has_physical_property( property ) )
            {
                continue;
            }
            const auto& property_info =
                model_.physical_property_info( property );
            auto attribute_name =
                transfer_block_attribute( model_, *model_solid_,
                    block_polyhedra_offsets_, property_info, geos_name );
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
            if( !model_.has_physical_property( property )
                || !is_read_by( *problem_, property ) )
            {
                continue;
            }
            const auto& property_info =
                model_.physical_property_info( property );
            if( !transfer_block_attribute( model_, *model_solid_,
                    block_polyhedra_offsets_, property_info, geos_name ) )
            {
                continue;
            }
            mechanics_attribute_ids_.emplace(
                property, property_info.attribute_id );
        }
        transfer_initial_stress();
        if( problem_->has_flow() )
        {
            compute_grain_bulk_modulus();
        }
    }

    template < typename Model >
    void GeosExporterImpl< Model >::transfer_initial_stress()
    {
        const auto factor = geos_stress_factor( *problem_ );
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
            transfer_stress_components( model_, model2solid_,
                model_.physical_property_info( property ), name,
                first_component, factor, initial_stress_tensor() );
        }
        // Created even without initial stress property: build_deck_context
        // imports it, then add_initial_pore_pressure shifts it
        if( problem_->has_flow() )
        {
            initial_stress_tensor();
        }
    }

    template < typename Model >
    VariableAttribute< StressTensor >&
        GeosExporterImpl< Model >::initial_stress_tensor()
    {
        auto& manager = model_solid_->polyhedron_attribute_manager();
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
        auto& manager = model_solid_->polyhedron_attribute_manager();
        const auto bulk_modulus =
            manager.template find_read_only_attribute< double >(
                mechanics_attribute_ids_.at(
                    PHYSICAL_PROPERTY_NAME::bulk_modulus ) );
        const auto biot_coefficient =
            manager.template find_read_only_attribute< double >(
                mechanics_attribute_ids_.at(
                    PHYSICAL_PROPERTY_NAME::biot_coefficient ) );
        const auto name = unique_attribute_name(
            manager, GRAIN_BULK_MODULUS, "grainBulkModulus" );
        AttributeValues< double > values;
        values.default_value = 0.;
        grain_bulk_modulus_attribute_id_ =
            manager.template create_attribute< VariableAttribute, double >(
                name, values, AttributeProperties{} );
        auto grain_bulk_modulus =
            manager.template find_attribute< VariableAttribute, double >(
                grain_bulk_modulus_attribute_id_.value() );
        for( const auto polyhedron : Range{ model_solid_->nb_polyhedra() } )
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
        const GeosRegionAssignment& initial_conditions )
    {
        const auto biot_coefficient =
            model_solid_->polyhedron_attribute_manager()
                .template find_read_only_attribute< double >(
                    mechanics_attribute_ids_.at(
                        PHYSICAL_PROPERTY_NAME::biot_coefficient ) );
        std::vector< double > region_pressures;
        region_pressures.reserve( initial_conditions.region_items.size() );
        for( const auto& condition_id : initial_conditions.region_items )
        {
            region_pressures.push_back(
                dynamic_cast< const UniformInitialCondition& >(
                    problem_->initial_condition( condition_id.value() ) )
                    .pressure() );
        }
        auto& tensor = initial_stress_tensor();
        for( const auto polyhedron : Range{ model_solid_->nb_polyhedra() } )
        {
            const auto shift =
                biot_coefficient->value( polyhedron )
                * region_pressures[region_attribute_->value( polyhedron )];
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
                "boundaryNormalTraction", prototype,
                geos_stress_factor( *problem_ ) );
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
            // Uniform on the Surface: checked by
            // validate_physical_problem_properties
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
            model_, surface_regions_, attribute_id );
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
            model_, surface_regions_, attribute_id );
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
                "[GeosExporter] Cannot write the model solid of type ",
                model_solid_->type_name().get(),
                ": only TetrahedralSolid, HybridSolid and PolyhedralSolid are "
                "supported." };
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
        context.boundary_fraction_phases = boundary_fraction_phases_;
        context.has_boundary_pressure = model_.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_pressure );
        context.has_boundary_temperature = model_.has_physical_property(
            PHYSICAL_PROPERTY_NAME::boundary_temperature );
        // Only the wells perforating the model define a GEOS set
        for( const auto well_id : Indices{ well_boxes_ } )
        {
            if( well_boxes_[well_id] )
            {
                context.well_names.push_back( well_names_[well_id] );
            }
        }
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
        for( auto& issue :
            validate_physical_problem_properties( *problem_, model_ ) )
        {
            issues.push_back( std::move( issue ) );
        }
        for( auto& issue : validate_geos_problem( *problem_, well_names_ ) )
        {
            issues.push_back( std::move( issue ) );
        }
        check_no_issue( issues );
        for( const auto& warning :
            physical_problem_warnings( *problem_, model_ ) )
        {
            Logger::warning( "[GeosExporter] ", warning );
        }
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
        if( problem_->has_flow() && problem_->solid_mechanics() )
        {
            add_initial_pore_pressure( initial_conditions );
        }
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