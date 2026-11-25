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

#include <geode/geosciences_io/model/internal/ogs_export.hpp>

#include <limits>
#include <utility>

#include <absl/algorithm/container.h>
#include <absl/container/inlined_vector.h>
#include <absl/strings/str_cat.h>

#include <pugixml.hpp>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/aabb.hpp>
#include <geode/geometry/distance.hpp>
#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/builder/solid_mesh_builder.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/regular_grid_solid.hpp>
#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/surface_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/helpers/aabb_solid_helpers.hpp>
#include <geode/mesh/helpers/convert_solid_mesh.hpp>

#include <geode/model/helpers/component_mesh_polygons.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/builder/brep_builder.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/helpers/simulation_report.hpp>

#include <geode/geosciences_io/model/internal/ogs_simulation_checker.hpp>
#include <geode/geosciences_io/model/internal/ogs_simulation_rules.hpp>
#include <geode/geosciences_io/model/internal/simulation_boundary_values.hpp>
#include <geode/geosciences_io/model/internal/simulation_vtu.hpp>

namespace
{
    constexpr std::array< geode::PHYSICAL_PROPERTY_NAME, 2 >
        OGS_BOUNDARY_PROPERTIES{
            geode::PHYSICAL_PROPERTY_NAME::boundary_pressure,
            geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction
        };

    // Builds the mesh of Surface polygons lying on bulk facets, its points
    // being the bulk vertices of the polygons
    class OgsBoundaryMeshBuilder
    {
    public:
        OgsBoundaryMeshBuilder( const geode::BRep& model,
            const geode::SolidMesh3D& bulk,
            const geode::ModelToMeshMappings& mappings )
            : model_( model ), bulk_( bulk ), mappings_( mappings )
        {
        }

        // Adds the polygon if it lies on a bulk facet, with its value (NaN
        // for none) given to each of its vertices
        void add_polygon( const geode::Surface3D& surface,
            geode::index_t polygon,
            double value )
        {
            if( !mappings_.surface_polygons_mapping.has_mapping_input(
                    { surface.id(), polygon } ) )
            {
                nb_skipped_++;
                return;
            }
            auto vertices =
                geode::polygon_unique_vertices( model_, surface, polygon );
            geode::PolyhedronFacetVertices facet_vertices;
            std::vector< geode::index_t > cell;
            for( const auto unique_vertex : vertices )
            {
                const auto vertex =
                    mappings_.unique_vertices_mapping.in2out( unique_vertex );
                facet_vertices.push_back( vertex );
                cell.push_back( local_vertex( vertex, value ) );
            }
            const auto polyhedra =
                bulk_.polyhedra_from_facet_vertices( facet_vertices );
            mesh_.bulk_element_ids.push_back(
                polyhedra.empty() ? 0 : polyhedra.front().polyhedron_id );
            mesh_.cells.push_back( std::move( cell ) );
        }

        // The mesh, with a "values" point array if values were given
        [[nodiscard]] geode::internal::OgsMeshFile mesh( bool with_values )
        {
            if( with_values )
            {
                std::vector< double > means;
                means.reserve( sums_.size() );
                for( const auto vertex : geode::Indices{ sums_ } )
                {
                    means.push_back( sums_[vertex] / counts_[vertex] );
                }
                mesh_.point_arrays.emplace_back( "values", std::move( means ) );
            }
            return std::move( mesh_ );
        }

        [[nodiscard]] geode::index_t nb_skipped_polygons() const
        {
            return nb_skipped_;
        }

    private:
        geode::index_t local_vertex( geode::index_t vertex, double value )
        {
            const auto [iterator, inserted] = local_vertices_.try_emplace(
                vertex, static_cast< geode::index_t >( mesh_.points.size() ) );
            if( inserted )
            {
                mesh_.points.push_back( bulk_.point( vertex ) );
                mesh_.bulk_node_ids.push_back( vertex );
                sums_.push_back( 0. );
                counts_.push_back( 0 );
            }
            const auto local = iterator->second;
            sums_[local] += value;
            counts_[local]++;
            return local;
        }

    private:
        const geode::BRep& model_;
        const geode::SolidMesh3D& bulk_;
        const geode::ModelToMeshMappings& mappings_;
        geode::internal::OgsMeshFile mesh_;
        absl::flat_hash_map< geode::index_t, geode::index_t > local_vertices_;
        std::vector< double > sums_;
        std::vector< geode::index_t > counts_;
        geode::index_t nb_skipped_{ 0 };
    };

    // Copy of the model with HybridSolid Blocks if its Blocks are all
    // RegularGrids, nothing otherwise
    std::optional< geode::BRep > ogs_grid_model( const geode::BRep& model )
    {
        if( model.nb_blocks() == 0 )
        {
            return std::nullopt;
        }
        for( const auto& block : model.blocks() )
        {
            if( block.mesh().type_name()
                != geode::RegularGrid3D::type_name_static() )
            {
                return std::nullopt;
            }
        }
        auto copy = model.clone();
        geode::BRepBuilder builder{ copy };
        for( const auto& block : copy.blocks() )
        {
            // The original grid: a cloned RegularGrid cannot be converted
            const auto& grid = model.block( block.id() ).mesh();
            auto hybrid = geode::convert_solid_mesh_into_hybrid_solid( grid );
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                hybrid.has_value(), nullptr,
                geode::OpenGeodeException::TYPE::internal,
                "[OgsExporter] Cannot convert RegularGrid Block ",
                block.id().string(), " into a HybridSolid" );
            const auto nb_vertices = grid.nb_vertices();
            builder.update_block_mesh( block, std::move( hybrid.value() ) );
            for( const auto vertex : geode::Range{ nb_vertices } )
            {
                const geode::ComponentMeshVertex component_vertex{
                    block.component_id(), vertex
                };
                builder.set_unique_vertex(
                    component_vertex, model.unique_vertex( component_vertex ) );
            }
        }
        return copy;
    }

    // Vertices defining the orientation of a VTK cell, OGS needing
    // det(v[1] - v[0], v[2] - v[0], v[3] - v[0]) > 0, and the local vertex
    // pairs swapped to mirror a cell
    struct OgsCellOrientation
    {
        std::array< geode::local_index_t, 4 > vertices;
        absl::InlinedVector< std::array< geode::local_index_t, 2 >, 4 > swaps;
    };

    std::optional< OgsCellOrientation > ogs_cell_orientation(
        geode::local_index_t nb_vertices )
    {
        switch( nb_vertices )
        {
        case 4: // tetrahedron
            return OgsCellOrientation{ { 0, 1, 2, 3 }, { { 1, 2 } } };
        case 5: // pyramid, apex last
            return OgsCellOrientation{ { 0, 1, 3, 4 }, { { 1, 3 } } };
        case 6: // prism: OGS reads the VTK wedge faces swapped
            return OgsCellOrientation{ { 0, 2, 1, 3 },
                { { 0, 3 }, { 1, 4 }, { 2, 5 } } };
        case 8: // hexahedron
            return OgsCellOrientation{ { 0, 1, 3, 4 },
                { { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } } };
        default:
            return std::nullopt;
        }
    }

    // Mirrors the cells OGS would read with a negative Jacobian
    void orient_ogs_cells( geode::SolidMesh3D& mesh )
    {
        auto builder = geode::SolidMeshBuilder3D::create( mesh );
        for( const auto cell : geode::Range{ mesh.nb_polyhedra() } )
        {
            const auto orientation =
                ogs_cell_orientation( mesh.nb_polyhedron_vertices( cell ) );
            if( !orientation )
            {
                continue;
            }
            const auto point = [&mesh, cell]( geode::local_index_t local ) {
                return mesh.point( mesh.polyhedron_vertex( { cell, local } ) );
            };
            const auto& vertices = orientation->vertices;
            const auto& origin = point( vertices[0] );
            const geode::Vector3D first{ origin, point( vertices[1] ) };
            const geode::Vector3D second{ origin, point( vertices[2] ) };
            const geode::Vector3D third{ origin, point( vertices[3] ) };
            if( first.cross( second ).dot( third ) >= 0. )
            {
                continue;
            }
            const auto cell_vertices = mesh.polyhedron_vertices( cell );
            for( const auto& [local, other] : orientation->swaps )
            {
                builder->set_polyhedron_vertex(
                    { cell, local }, cell_vertices[other] );
                builder->set_polyhedron_vertex(
                    { cell, other }, cell_vertices[local] );
            }
        }
    }

    geode::index_t nearest_cell_vertex( const geode::SolidMesh3D& bulk,
        geode::index_t cell,
        const geode::Point3D& point )
    {
        auto nearest = bulk.polyhedron_vertex( { cell, 0 } );
        auto nearest_distance = std::numeric_limits< double >::max();
        for( const auto vertex : bulk.polyhedron_vertices( cell ) )
        {
            const auto distance =
                geode::point_point_distance( point, bulk.point( vertex ) );
            if( distance < nearest_distance )
            {
                nearest_distance = distance;
                nearest = vertex;
            }
        }
        return nearest;
    }

    template < typename T >
    std::optional< std::vector< T > > cell_values(
        const geode::SolidMesh3D& solid, std::string_view name )
    {
        const auto& manager = solid.polyhedron_attribute_manager();
        const auto ids = manager.attribute_ids_matching_name( name );
        if( !ids || ids->size() != 1
            || manager.attribute_type( ids->front() ) != typeid( T ).name() )
        {
            return std::nullopt;
        }
        const auto attribute =
            manager.find_read_only_attribute< T >( ids->front() );
        std::vector< T > values;
        values.reserve( solid.nb_polyhedra() );
        for( const auto cell : geode::Range{ solid.nb_polyhedra() } )
        {
            values.push_back( attribute->value( cell ) );
        }
        return values;
    }

    template < typename T >
    void add_cell_attribute( geode::AttributeManager& manager,
        std::string_view name,
        const std::vector< T >& values )
    {
        geode::AttributeValues< T > default_values;
        default_values.default_value = T{};
        const auto attribute_id =
            manager.create_attribute< geode::VariableAttribute, T >(
                name, default_values, geode::AttributeProperties{} );
        auto attribute = manager.find_attribute< geode::VariableAttribute, T >(
            attribute_id );
        for( const auto cell : geode::Indices{ values } )
        {
            attribute->set_value( cell, values[cell] );
        }
    }
} // namespace

namespace geode::internal
{
    template < typename Model >
    OgsExporterImpl< Model >::OgsExporterImpl( std::string_view files_directory,
        const Model& model,
        std::string_view context )
        : SimulationExporterImpl< Model >( files_directory, model, context ),
          context_{ to_string( context ) }
    {
    }

    template < typename Model >
    OgsExporterImpl< Model >::~OgsExporterImpl() = default;

    template < typename Model >
    void OgsExporterImpl< Model >::check_well_name(
        std::string_view name ) const
    {
        OpenGeodeGeosciencesIOModelException::check_exception(
            ogs_naming().is_valid_name( name ), nullptr,
            OpenGeodeException::TYPE::data, context_, " Well name \"", name,
            "\" must only contain letters, digits, '_' and '-': it is part "
            "of an OGS file name." );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::prepare_export()
    {
        clear_export();
        regions_ = create_region_definitions( this->model() );
        if( this->problem() != nullptr )
        {
            check_problem();
        }
        build_bulk();
        orient_bulk();
        compute_material_ids();
        transfer_cell_properties();
        build_boundary_meshes();
        build_well_meshes();
        if( this->problem() != nullptr )
        {
            compute_initial_values();
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::clear_export()
    {
        regions_ = SimulationRegions{};
        solid_.reset();
        grid_model_.reset();
        hybrid_.reset();
        material_ids_.clear();
        porosity_.reset();
        permeability_.reset();
        storage_.reset();
        cell_pressures_.clear();
        initial_pressure_ = OgsInitialValue{};
        initial_capillary_pressure_ = OgsInitialValue{};
        surfaces_.clear();
        surface_indices_.clear();
        boundary_meshes_.clear();
        well_meshes_.clear();
    }

    template < typename Model >
    void OgsExporterImpl< Model >::check_problem() const
    {
        const auto report = OgsSimulationChecker{
            *this->problem(), this->model(),
            OgsExportData{
                { this->well_names().begin(), this->well_names().end() } }
        }.check();
        OpenGeodeGeosciencesIOModelException::check_exception(
            report.is_simulable(), nullptr, OpenGeodeException::TYPE::data,
            context_, " Cannot write the OGS project.\n", report.string() );
        for( const auto& issue : report.issues() )
        {
            Logger::warning( context_, " ", issue.string() );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_bulk()
    {
        grid_model_ = ogs_grid_model( this->model() );
        solid_.emplace( solid_model(), regions_, context_ );
        const auto& mesh = solid_->mesh();
        if( dynamic_cast< const TetrahedralSolid3D* >( &mesh ) != nullptr
            || dynamic_cast< const HybridSolid3D* >( &mesh ) != nullptr )
        {
            return;
        }
        auto hybrid = convert_solid_mesh_into_hybrid_solid( mesh );
        OpenGeodeGeosciencesIOModelException::check_exception(
            hybrid.has_value(), nullptr, OpenGeodeException::TYPE::data,
            context_,
            " OGS only reads tetrahedra, hexahedra, prisms and pyramids: the "
            "model Blocks have other polyhedra." );
        hybrid_ = std::move( hybrid.value() );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::orient_bulk()
    {
        orient_ogs_cells( bulk_mesh() );
    }

    template < typename Model >
    const SolidMesh3D& OgsExporterImpl< Model >::bulk() const
    {
        if( hybrid_ )
        {
            return *hybrid_;
        }
        return solid_->mesh();
    }

    template < typename Model >
    SolidMesh3D& OgsExporterImpl< Model >::bulk_mesh()
    {
        if( hybrid_ )
        {
            return *hybrid_;
        }
        return solid_->mesh();
    }

    template < typename Model >
    const BRep& OgsExporterImpl< Model >::solid_model() const
    {
        if( grid_model_ )
        {
            return grid_model_.value();
        }
        return this->model();
    }

    template < typename Model >
    void OgsExporterImpl< Model >::compute_material_ids()
    {
        const auto nb_cells = bulk().nb_polyhedra();
        material_ids_.resize( nb_cells );
        if( this->problem() == nullptr )
        {
            for( const auto cell : Range{ nb_cells } )
            {
                material_ids_[cell] =
                    static_cast< std::int32_t >( solid_->region( cell ) );
            }
            return;
        }
        const auto assignment = assign_region_materials(
            *this->problem(), this->model(), regions_, "OGS" );
        for( const auto cell : Range{ nb_cells } )
        {
            const auto& material =
                assignment.region_items.at( solid_->region( cell ) );
            OpenGeodeGeosciencesIOModelException::check_exception(
                material.has_value(), nullptr, OpenGeodeException::TYPE::data,
                context_, " Cell ", cell, " has no material." );
            material_ids_[cell] =
                static_cast< std::int32_t >( material.value() );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::transfer_cell_properties()
    {
        const auto& model = solid_model();
        for( const auto property : { PHYSICAL_PROPERTY_NAME::porosity,
                 PHYSICAL_PROPERTY_NAME::permeability } )
        {
            if( !model.has_physical_property( property ) )
            {
                continue;
            }
            const auto& info = model.physical_property_info( property );
            const auto field = physical_property_string( property );
            const auto name = solid_->transfer_block_property(
                ogs_naming(), model, info, field );
            if( !name )
            {
                continue;
            }
            const auto& solid = solid_->mesh();
            if( property == PHYSICAL_PROPERTY_NAME::porosity )
            {
                porosity_ = cell_values< double >( solid, name.value() );
                continue;
            }
            permeability_ =
                cell_values< std::array< double, 3 > >( solid, name.value() );
            if( permeability_ )
            {
                continue;
            }
            if( const auto scalar =
                    cell_values< double >( solid, name.value() ) )
            {
                permeability_.emplace();
                for( const auto value : scalar.value() )
                {
                    permeability_->push_back( { value, value, value } );
                }
            }
        }
        if( permeability_ && is_two_phase() )
        {
            OpenGeodeGeosciencesIOModelException::check_exception(
                absl::c_all_of( permeability_.value(), is_isotropic ), nullptr,
                OpenGeodeException::TYPE::data, context_,
                " OGS two-phase flow needs an isotropic permeability: the "
                "model permeability differs along x, y and z." );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::compute_initial_values()
    {
        const auto& problem = *this->problem();
        if( porosity_ )
        {
            storage_.emplace();
            for( const auto cell : Indices{ porosity_.value() } )
            {
                const auto& material =
                    problem.materials()[static_cast< index_t >(
                        material_ids_[cell] )];
                storage_->push_back( ogs_storage(
                    porosity_.value()[cell], material.pore_compressibility ) );
            }
        }
        const auto assignment = assign_region_initial_conditions(
            problem, this->model(), regions_, "OGS" );
        std::vector< index_t > cell_conditions;
        cell_conditions.reserve( bulk().nb_polyhedra() );
        cell_pressures_.clear();
        for( const auto cell : Range{ bulk().nb_polyhedra() } )
        {
            const auto& condition_id =
                assignment.region_items.at( solid_->region( cell ) );
            OpenGeodeGeosciencesIOModelException::check_exception(
                condition_id.has_value(), nullptr,
                OpenGeodeException::TYPE::data, context_, " Cell ", cell,
                " has no initial condition." );
            const auto& condition =
                dynamic_cast< const UniformInitialCondition& >(
                    problem.initial_condition( condition_id.value() ) );
            cell_pressures_.push_back( condition.pressure() );
            cell_conditions.push_back( condition_id.value() );
        }
        initial_pressure_ =
            initial_value( "initial_pressure", cell_pressures_ );
        if( ogs_process( problem ) == OGS_PROCESS::two_phase_flow_pp )
        {
            compute_two_phase_values( cell_conditions );
        }
    }

    template < typename Model >
    double OgsExporterImpl< Model >::cell_capillary_pressure(
        index_t cell, double water_fraction ) const
    {
        const auto& problem = *this->problem();
        const auto& fluid =
            dynamic_cast< const DeadOilFluid& >( problem.fluid() );
        std::vector< double > composition( fluid.components().size(), 0. );
        composition[phase_component( fluid, PHASE::water ).value()] =
            water_fraction;
        composition[phase_component( fluid, PHASE::oil ).value()] =
            1. - water_fraction;
        const auto basis =
            dynamic_cast< const MultiphaseFlow& >( problem.flow() ).basis();
        const auto saturation = water_saturation(
            fluid, basis, composition, cell_pressures_[cell] );
        const auto& material =
            problem.materials()[static_cast< index_t >( material_ids_[cell] )];
        return capillary_pressure( *material.capillary_pressure, saturation );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::compute_two_phase_values(
        absl::Span< const index_t > cell_conditions )
    {
        const auto& problem = *this->problem();
        const auto& fluid =
            dynamic_cast< const DeadOilFluid& >( problem.fluid() );
        const auto water = phase_component( fluid, PHASE::water ).value();
        std::vector< double > capillary_pressures;
        capillary_pressures.reserve( cell_conditions.size() );
        for( const auto cell : Indices{ cell_conditions } )
        {
            const auto& condition =
                dynamic_cast< const UniformInitialCondition& >(
                    problem.initial_condition( cell_conditions[cell] ) );
            capillary_pressures.push_back( cell_capillary_pressure(
                cell, condition.composition()[water] ) );
        }
        initial_capillary_pressure_ =
            initial_value( "initial_capillary_pressure", capillary_pressures );
        convert_boundary_water_fractions();
    }

    template < typename Model >
    void OgsExporterImpl< Model >::convert_boundary_water_fractions()
    {
        for( auto& boundary : boundary_meshes_ )
        {
            if( boundary.property
                != PHYSICAL_PROPERTY_NAME::boundary_water_fraction )
            {
                continue;
            }
            auto& mesh = boundary.mesh;
            const auto nb_points = mesh.points.size();
            std::vector< double > fractions(
                nb_points, boundary.value.value_or( 0. ) );
            if( !boundary.value )
            {
                fractions = mesh.point_arrays.front().second;
            }
            std::vector< double > sums( nb_points, 0. );
            std::vector< index_t > counts( nb_points, 0 );
            for( const auto cell : Indices{ mesh.cells } )
            {
                const auto bulk_cell =
                    static_cast< index_t >( mesh.bulk_element_ids[cell] );
                for( const auto vertex : mesh.cells[cell] )
                {
                    sums[vertex] +=
                        cell_capillary_pressure( bulk_cell, fractions[vertex] );
                    counts[vertex]++;
                }
            }
            std::vector< double > pressures;
            pressures.reserve( nb_points );
            for( const auto vertex : Range{ nb_points } )
            {
                pressures.push_back(
                    sums[vertex] / static_cast< double >( counts[vertex] ) );
            }
            mesh.point_arrays.clear();
            if( absl::c_all_of( pressures, [&pressures]( double value ) {
                    return value == pressures.front();
                } ) )
            {
                boundary.value = pressures.front();
                continue;
            }
            boundary.value.reset();
            mesh.point_arrays.emplace_back( "values", std::move( pressures ) );
        }
    }

    template < typename Model >
    OgsInitialValue OgsExporterImpl< Model >::initial_value(
        std::string_view field, const std::vector< double >& cell_values )
    {
        if( absl::c_all_of( cell_values, [&cell_values]( double value ) {
                return value == cell_values.front();
            } ) )
        {
            return { cell_values.empty() ? 0. : cell_values.front(), "" };
        }
        auto& mesh = bulk_mesh();
        AttributeValues< double > values;
        values.default_value = 0.;
        const auto attribute_id =
            mesh.vertex_attribute_manager()
                .template create_attribute< VariableAttribute, double >(
                    field, values, AttributeProperties{} );
        auto attribute =
            mesh.vertex_attribute_manager()
                .template find_attribute< VariableAttribute, double >(
                    attribute_id );
        for( const auto vertex : Range{ mesh.nb_vertices() } )
        {
            const auto& polyhedra = mesh.polyhedra_around_vertex( vertex );
            double sum{ 0 };
            for( const auto& polyhedron : polyhedra )
            {
                sum += cell_values[polyhedron.polyhedron_id];
            }
            if( !polyhedra.empty() )
            {
                attribute->set_value(
                    vertex, sum / static_cast< double >( polyhedra.size() ) );
            }
        }
        return { std::nullopt, to_string( field ) };
    }

    template < typename Model >
    bool OgsExporterImpl< Model >::is_two_phase() const
    {
        return this->problem() != nullptr
               && ogs_process( *this->problem() )
                      == OGS_PROCESS::two_phase_flow_pp;
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_boundary_meshes()
    {
        for( const auto& surface : solid_model().surfaces() )
        {
            surface_indices_.emplace( surface.id(), surfaces_.size() );
            surfaces_.push_back( &surface );
        }
        const auto* problem = this->problem();
        for( const auto property : OGS_BOUNDARY_PROPERTIES )
        {
            if( property == PHYSICAL_PROPERTY_NAME::boundary_water_fraction
                && problem != nullptr
                && ogs_process( *problem ) != OGS_PROCESS::two_phase_flow_pp )
            {
                continue;
            }
            build_boundary_meshes( property );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_boundary_meshes(
        PHYSICAL_PROPERTY_NAME property )
    {
        const auto& model = solid_model();
        if( !model.has_physical_property( property ) )
        {
            return;
        }
        const auto& info = model.physical_property_info( property );
        if( info.component_type != Surface3D::component_type_static() )
        {
            return;
        }
        const auto groups = group_surface_values< double >(
            model, surface_indices_, info.attribute_id, context_ );
        const auto property_name = physical_property_string( property );
        const auto add_mesh = [this, property, property_name](
                                  OgsBoundaryMeshBuilder& builder,
                                  std::optional< double > value ) {
            if( builder.nb_skipped_polygons() != 0 )
            {
                Logger::warning( context_, " ", builder.nb_skipped_polygons(),
                    " polygons with ", property_name,
                    " are not on a Block polyhedron facet, they are not "
                    "exported." );
            }
            auto mesh = builder.mesh( !value.has_value() );
            if( mesh.cells.empty() )
            {
                return;
            }
            boundary_meshes_.push_back( { property,
                absl::StrCat( this->prefix(), "_", property_name, "_",
                    boundary_meshes_.size() ),
                value, std::move( mesh ) } );
        };
        for( const auto& [value, surface_ids] : groups.uniform_regions )
        {
            OgsBoundaryMeshBuilder builder{ model, bulk(), solid_->mappings() };
            for( const auto surface_id : surface_ids )
            {
                const auto& surface = *surfaces_[surface_id];
                for( const auto polygon :
                    Range{ surface.mesh().nb_polygons() } )
                {
                    builder.add_polygon( surface, polygon, value );
                }
            }
            add_mesh( builder, value );
        }
        for( const auto& variable : groups.variable_surfaces )
        {
            OgsBoundaryMeshBuilder builder{ model, bulk(), solid_->mappings() };
            for( const auto polygon : Indices{ variable.values } )
            {
                builder.add_polygon(
                    *variable.surface, polygon, variable.values[polygon] );
            }
            add_mesh( builder, std::nullopt );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_well_meshes()
    {
        if( this->well_perforations().empty() )
        {
            return;
        }
        const auto aabb = create_aabb_tree( bulk() );
        for( const auto well_id : Indices{ this->well_perforations() } )
        {
            const auto& well = *this->well_perforations()[well_id];
            const auto& name = this->well_names()[well_id];
            OgsWellMesh well_mesh{ name,
                absl::StrCat( this->prefix(), "_well_", name ), {} };
            for( const auto point_id : Range{ well.nb_vertices() } )
            {
                const auto& point = well.point( point_id );
                const auto cell = containing_cell( bulk(), aabb, point );
                if( cell == NO_ID )
                {
                    Logger::warning( context_, " Perforation ", point_id,
                        " of well ", name, " (", point.string(),
                        ") is not inside the model, it is ignored." );
                    continue;
                }
                const auto vertex = nearest_cell_vertex( bulk(), cell, point );
                auto& mesh = well_mesh.mesh;
                if( absl::c_linear_search( mesh.bulk_node_ids, vertex ) )
                {
                    continue;
                }
                mesh.cells.push_back(
                    { static_cast< index_t >( mesh.points.size() ) } );
                mesh.points.push_back( bulk().point( vertex ) );
                mesh.bulk_node_ids.push_back( vertex );
            }
            if( well_mesh.mesh.points.empty() )
            {
                Logger::warning( context_, " No perforation of well ", name,
                    " is inside the model, no well mesh is exported." );
                continue;
            }
            well_meshes_.push_back( std::move( well_mesh ) );
        }
    }

    template < typename Model >
    std::string OgsExporterImpl< Model >::bulk_file() const
    {
        return absl::StrCat( this->prefix(), ".vtu" );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::write_bulk() const
    {
        AttributeManager cells;
        cells.resize( bulk().nb_polyhedra() );
        add_cell_attribute( cells, "MaterialIDs", material_ids_ );
        if( porosity_ )
        {
            add_cell_attribute( cells,
                physical_property_string( PHYSICAL_PROPERTY_NAME::porosity ),
                porosity_.value() );
        }
        if( storage_ )
        {
            add_cell_attribute( cells, "storage", storage_.value() );
        }
        if( permeability_ )
        {
            const auto name = physical_property_string(
                PHYSICAL_PROPERTY_NAME::permeability );
            if( is_two_phase() )
            {
                // TWOPHASE_FLOW_PP reads a scalar permeability
                std::vector< double > scalars;
                scalars.reserve( permeability_->size() );
                for( const auto& value : permeability_.value() )
                {
                    scalars.push_back( value[0] );
                }
                add_cell_attribute( cells, name, scalars );
            }
            else
            {
                add_cell_attribute( cells, name, permeability_.value() );
            }
        }
        const auto written = save_simulation_vtu(
            absl::StrCat( this->files_directory(), "/", bulk_file() ), bulk(),
            {}, cells );
        OpenGeodeGeosciencesIOModelException::check_exception( written, nullptr,
            OpenGeodeException::TYPE::internal, context_,
            " Cannot write the bulk mesh of type ", bulk().type_name().get() );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::write_project() const
    {
        const auto& problem = *this->problem();
        OgsProjectContext context;
        context.process = ogs_process( problem ).value();
        context.prefix = this->prefix();
        context.has_porosity_field = porosity_.has_value();
        context.has_permeability_field = permeability_.has_value();
        context.boundaries = boundary_meshes_;
        context.wells = well_meshes_;
        context.initial_pressure = initial_pressure_;
        context.initial_capillary_pressure = initial_capillary_pressure_;
        context.reference_temperature = ogs_reference_temperature( problem );
        pugi::xml_document document;
        OgsProjectWriter{ problem, context }.write( document );
        const auto filename = absl::StrCat(
            this->files_directory(), "/", this->prefix(), ".prj" );
        OpenGeodeGeosciencesIOModelException::check_exception(
            document.save_file( filename.c_str(), "    " ), nullptr,
            OpenGeodeException::TYPE::data, context_, " Cannot write file ",
            filename );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::write_files() const
    {
        write_bulk();
        for( const auto& boundary : boundary_meshes_ )
        {
            save_ogs_mesh_file( absl::StrCat( this->files_directory(), "/",
                                    boundary.stem, ".vtu" ),
                boundary.mesh );
        }
        for( const auto& well : well_meshes_ )
        {
            save_ogs_mesh_file(
                absl::StrCat( this->files_directory(), "/", well.stem, ".vtu" ),
                well.mesh );
        }
        if( this->problem() != nullptr )
        {
            write_project();
        }
    }

    template class opengeode_geosciencesio_model_api OgsExporterImpl< BRep >;
    template class opengeode_geosciencesio_model_api
        OgsExporterImpl< StructuralModel >;
} // namespace geode::internal
