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

#include <cmath>
#include <limits>
#include <tuple>
#include <utility>

#include <absl/algorithm/container.h>
#include <absl/container/inlined_vector.h>
#include <absl/container/linked_hash_map.h>
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
#include <geode/geosciences/simulation/core/heat_source.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/helpers/simulation_report.hpp>

#include <geode/geosciences_io/model/internal/ogs_simulation_checker.hpp>
#include <geode/geosciences_io/model/internal/ogs_simulation_rules.hpp>
#include <geode/geosciences_io/model/internal/simulation_boundary_values.hpp>
#include <geode/geosciences_io/model/internal/simulation_cell_properties.hpp>
#include <geode/geosciences_io/model/internal/simulation_vtu.hpp>

namespace
{
    // Surface properties exported as boundary meshes: boundary_pressure and
    // boundary_water_fraction without problem, the ones of the primary
    // variables of the process otherwise
    absl::InlinedVector< geode::PHYSICAL_PROPERTY_NAME, 2 >
        ogs_boundary_properties( const geode::PhysicalProblem* problem )
    {
        const auto process = problem == nullptr
                                 ? std::nullopt
                                 : geode::internal::ogs_process( *problem );
        if( process == geode::internal::OGS_PROCESS::small_deformation )
        {
            return {};
        }
        absl::InlinedVector< geode::PHYSICAL_PROPERTY_NAME, 2 > properties{
            geode::PHYSICAL_PROPERTY_NAME::boundary_pressure
        };
        if( problem == nullptr
            || process == geode::internal::OGS_PROCESS::two_phase_flow_pp )
        {
            properties.push_back(
                geode::PHYSICAL_PROPERTY_NAME::boundary_water_fraction );
        }
        if( process && geode::internal::has_temperature( process.value() ) )
        {
            properties.push_back(
                geode::PHYSICAL_PROPERTY_NAME::boundary_temperature );
        }
        return properties;
    }

    // Initial temperature of a cell in a thermal problem
    double ogs_cell_temperature(
        const geode::UniformInitialCondition& condition,
        geode::index_t cell,
        std::string_view context )
    {
        const auto& temperature = condition.temperature();
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            temperature.has_value(), nullptr,
            geode::OpenGeodeException::TYPE::data, context, " Cell ", cell,
            " has no initial temperature: a thermal problem needs one." );
        return temperature.value();
    }

    // Where the values given to the polygons are written
    enum struct OGS_VALUES
    {
        none,
        // Mean of the polygon values at each vertex
        points,
        // Value of each polygon
        cells
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

        // Value given to the vertices of the next polygons
        void set_value( double value )
        {
            value_ = value;
        }

        // Orders the vertices of the next polygons with their right-hand
        // normal out of their bulk cell
        void set_outward()
        {
            outward_ = true;
        }

        // Adds the polygon if it lies on a bulk facet, the current value
        // being given to each of its vertices
        void add_polygon(
            const geode::Surface3D& surface, geode::index_t polygon )
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
                cell.push_back( local_vertex( vertex ) );
            }
            const auto polyhedra =
                bulk_.polyhedra_from_facet_vertices( facet_vertices );
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                !polyhedra.empty(), nullptr,
                geode::OpenGeodeException::TYPE::internal,
                "[OgsExporter] No bulk cell owns polygon ", polygon,
                " of Surface ", surface.id().string() );
            const auto owner = polyhedra.front().polyhedron_id;
            if( outward_ && !points_outward( facet_vertices, owner ) )
            {
                absl::c_reverse( cell );
            }
            mesh_.bulk_element_ids.push_back( owner );
            mesh_.cells.push_back( std::move( cell ) );
            cell_values_.push_back( value_ );
        }

        void add_surfaces( absl::Span< const geode::Surface3D* const > surfaces,
            absl::Span< const geode::index_t > surface_ids )
        {
            for( const auto surface_id : surface_ids )
            {
                add_surface( *surfaces[surface_id] );
            }
        }

        void add_surface( const geode::Surface3D& surface )
        {
            for( const auto polygon :
                geode::Range{ surface.mesh().nb_polygons() } )
            {
                add_polygon( surface, polygon );
            }
        }

        void add_surface_values(
            const geode::Surface3D& surface, absl::Span< const double > values )
        {
            for( const auto polygon : geode::Indices{ values } )
            {
                set_value( values[polygon] );
                add_polygon( surface, polygon );
            }
        }

        // The mesh, with a "values" point or cell array
        [[nodiscard]] geode::internal::OgsMeshFile mesh( OGS_VALUES values )
        {
            if( values == OGS_VALUES::cells )
            {
                mesh_.cell_arrays.emplace_back(
                    "values", std::move( cell_values_ ) );
            }
            if( values == OGS_VALUES::points )
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
        // The Newell normal of the facet points away from the polyhedron
        [[nodiscard]] bool points_outward(
            const geode::PolyhedronFacetVertices& vertices,
            geode::index_t polyhedron ) const
        {
            std::array< double, 3 > normal{ 0, 0, 0 };
            std::array< double, 3 > centroid{ 0, 0, 0 };
            for( const auto index : geode::Indices{ vertices } )
            {
                const auto& point = bulk_.point( vertices[index] );
                const auto& next =
                    bulk_.point( vertices[( index + 1 ) % vertices.size()] );
                for( const auto axis : geode::LRange{ 3 } )
                {
                    const auto first = ( axis + 1 ) % 3;
                    const auto second = ( axis + 2 ) % 3;
                    normal[axis] +=
                        ( point.value( first ) - next.value( first ) )
                        * ( point.value( second ) + next.value( second ) );
                    centroid[axis] +=
                        point.value( axis )
                        / static_cast< double >( vertices.size() );
                }
            }
            const geode::Vector3D outward{ bulk_.polyhedron_barycenter(
                                               polyhedron ),
                geode::Point3D{ centroid } };
            return geode::Vector3D{ normal }.dot( outward ) > 0.;
        }

        geode::index_t local_vertex( geode::index_t vertex )
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
            sums_[local] += value_;
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
        std::vector< double > cell_values_;
        geode::index_t nb_skipped_{ 0 };
        double value_{ 0 };
        bool outward_{ false };
    };

    // Surfaces (by index) of each (axis, value) of the uniform displacements;
    // NaN components are free
    absl::linked_hash_map< std::pair< geode::local_index_t, double >,
        std::vector< geode::index_t > >
        ogs_displacement_groups(
            absl::Span< const geode::Surface3D* const > surfaces,
            const geode::uuid& attribute_id,
            std::string_view context )
    {
        absl::linked_hash_map< std::pair< geode::local_index_t, double >,
            std::vector< geode::index_t > >
            groups;
        for( const auto surface_id : geode::Indices{ surfaces } )
        {
            const auto values =
                geode::internal::surface_values< std::array< double, 3 > >(
                    *surfaces[surface_id], attribute_id, context );
            if( !values )
            {
                continue;
            }
            // Uniform on the Surface: checked by
            // SimulationChecker::check_physical_problem_properties
            const auto& displacement = values->front();
            for( const auto axis : geode::LRange{ 3 } )
            {
                if( !std::isnan( displacement[axis] ) )
                {
                    groups[{ axis, displacement[axis] }].push_back(
                        surface_id );
                }
            }
        }
        return groups;
    }

    bool is_zero( absl::Span< const double > values )
    {
        return absl::c_all_of( values, []( double value ) {
            return value == 0.;
        } );
    }

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

    // Swaps of the base and top vertices of a prism or hexahedron
    absl::InlinedVector< std::array< geode::local_index_t, 2 >, 4 >
        ogs_base_top_swaps( geode::local_index_t nb_base_vertices )
    {
        absl::InlinedVector< std::array< geode::local_index_t, 2 >, 4 > swaps;
        for( const auto vertex : geode::LRange{ nb_base_vertices } )
        {
            swaps.push_back( { vertex, static_cast< geode::local_index_t >(
                                           vertex + nb_base_vertices ) } );
        }
        return swaps;
    }

    std::optional< OgsCellOrientation > ogs_cell_orientation(
        geode::local_index_t nb_vertices )
    {
        if( nb_vertices == geode::internal::NB_TETRAHEDRON_VERTICES )
        {
            return OgsCellOrientation{ { 0, 1, 2, 3 }, { { 1, 2 } } };
        }
        if( nb_vertices == geode::internal::NB_PYRAMID_VERTICES )
        {
            // Apex last
            return OgsCellOrientation{ { 0, 1, 3, 4 }, { { 1, 3 } } };
        }
        if( nb_vertices == geode::internal::NB_PRISM_VERTICES )
        {
            // OGS reads the VTK wedge faces swapped
            return OgsCellOrientation{ { 0, 2, 1, 3 },
                ogs_base_top_swaps( geode::internal::NB_PRISM_VERTICES / 2 ) };
        }
        if( nb_vertices == geode::internal::NB_HEXAHEDRON_VERTICES )
        {
            return OgsCellOrientation{ { 0, 1, 3, 4 },
                ogs_base_top_swaps(
                    geode::internal::NB_HEXAHEDRON_VERTICES / 2 ) };
        }
        return std::nullopt;
    }

    // Mirrors the cells OGS would read with a negative Jacobian
    void orient_ogs_cells( geode::SolidMesh3D& mesh )
    {
        auto builder = geode::SolidMeshBuilder3D::create( mesh );
        bool has_mirrored{ false };
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
            has_mirrored = true;
            const auto cell_vertices = mesh.polyhedron_vertices( cell );
            for( const auto& [local, other] : orientation->swaps )
            {
                builder->set_polyhedron_vertex(
                    { cell, local }, cell_vertices[other] );
                builder->set_polyhedron_vertex(
                    { cell, other }, cell_vertices[local] );
            }
        }
        // Mirroring moves the facets: their adjacencies are recomputed
        if( has_mirrored )
        {
            builder->compute_polyhedron_adjacencies();
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

    // Permeability of each cell, a scalar attribute giving the same value
    // along x, y and z
    std::optional< std::vector< std::array< double, 3 > > >
        ogs_permeability_values(
            const geode::SolidMesh3D& solid, std::string_view name )
    {
        if( auto values =
                cell_values< std::array< double, 3 > >( solid, name ) )
        {
            return values;
        }
        const auto scalars = cell_values< double >( solid, name );
        if( !scalars )
        {
            return std::nullopt;
        }
        std::vector< std::array< double, 3 > > values;
        values.reserve( scalars->size() );
        for( const auto value : scalars.value() )
        {
            values.push_back( { value, value, value } );
        }
        return values;
    }

    // The value shared by every item, none if they differ or if there is none
    std::optional< double > ogs_uniform_value(
        absl::Span< const double > values )
    {
        if( values.empty()
            || !absl::c_all_of( values, [&values]( double value ) {
                   return value == values.front();
               } ) )
        {
            return std::nullopt;
        }
        return values.front();
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

    // Permeability cell array: TWOPHASE_FLOW_PP reads a scalar
    void add_permeability_attribute( geode::AttributeManager& cells,
        const std::vector< std::array< double, 3 > >& permeability,
        bool scalar )
    {
        const auto name = geode::physical_property_string(
            geode::PHYSICAL_PROPERTY_NAME::permeability );
        if( !scalar )
        {
            add_cell_attribute( cells, name, permeability );
            return;
        }
        std::vector< double > scalars;
        scalars.reserve( permeability.size() );
        for( const auto& value : permeability )
        {
            scalars.push_back( value[0] );
        }
        add_cell_attribute( cells, name, scalars );
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
            build_heat_source_meshes();
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
        initial_temperature_ = OgsInitialValue{};
        mechanics_ = OgsMechanicsValues{};
        initial_stress_.reset();
        cell_arrays_.clear();
        surfaces_.clear();
        surface_indices_.clear();
        boundary_meshes_.clear();
        well_meshes_.clear();
        heat_source_meshes_.clear();
    }

    template < typename Model >
    void OgsExporterImpl< Model >::check_problem() const
    {
        const auto report = OgsSimulationChecker{
            *this->problem(), this->model(),
            OgsExportData{
                { this->well_names().begin(), this->well_names().end() },
                { this->heat_source_names().begin(),
                    this->heat_source_names().end() } }
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
        // Materials are ignored without flow: a single medium
        if( process() == OGS_PROCESS::small_deformation )
        {
            absl::c_fill( material_ids_, 0 );
            return;
        }
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
    std::optional< OGS_PROCESS > OgsExporterImpl< Model >::process() const
    {
        if( this->problem() == nullptr )
        {
            return std::nullopt;
        }
        return ogs_process( *this->problem() );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::transfer_cell_properties()
    {
        const auto ogs_process = process();
        if( ogs_process && has_displacement( ogs_process.value() ) )
        {
            transfer_mechanics_properties();
        }
        if( ogs_process == OGS_PROCESS::small_deformation )
        {
            return;
        }
        if( const auto name =
                transfer_property( PHYSICAL_PROPERTY_NAME::porosity ) )
        {
            porosity_ = cell_values< double >( solid_->mesh(), name.value() );
        }
        if( const auto name =
                transfer_property( PHYSICAL_PROPERTY_NAME::permeability ) )
        {
            permeability_ =
                ogs_permeability_values( solid_->mesh(), name.value() );
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
    std::optional< std::string > OgsExporterImpl< Model >::transfer_property(
        PHYSICAL_PROPERTY_NAME property )
    {
        const auto& model = solid_model();
        if( !model.has_physical_property( property ) )
        {
            return std::nullopt;
        }
        return solid_->transfer_block_property( ogs_naming(), model,
            model.physical_property_info( property ),
            physical_property_string( property ) );
    }

    template < typename Model >
    std::vector< double > OgsExporterImpl< Model >::property_cell_values(
        PHYSICAL_PROPERTY_NAME property )
    {
        if( const auto name = transfer_property( property ) )
        {
            if( auto values =
                    cell_values< double >( solid_->mesh(), name.value() ) )
            {
                return std::move( values.value() );
            }
        }
        return std::vector< double >( bulk().nb_polyhedra(), 0. );
    }

    template < typename Model >
    OgsCellValue OgsExporterImpl< Model >::cell_value(
        std::string_view field, std::vector< double > values )
    {
        if( const auto uniform = ogs_uniform_value( values ) )
        {
            return { uniform, "" };
        }
        cell_arrays_.emplace_back( to_string( field ), std::move( values ) );
        return { std::nullopt, to_string( field ) };
    }

    template < typename Model >
    void OgsExporterImpl< Model >::transfer_mechanics_properties()
    {
        const auto ogs_process = process().value();
        const auto bulk_moduli =
            property_cell_values( PHYSICAL_PROPERTY_NAME::bulk_modulus );
        const auto shear_moduli =
            property_cell_values( PHYSICAL_PROPERTY_NAME::shear_modulus );
        std::vector< double > youngs_moduli;
        std::vector< double > poissons_ratios;
        for( const auto cell : Indices{ bulk_moduli } )
        {
            const auto elasticity =
                ogs_elasticity( { bulk_moduli[cell], shear_moduli[cell] } );
            youngs_moduli.push_back( elasticity.youngs_modulus );
            poissons_ratios.push_back( elasticity.poissons_ratio );
        }
        mechanics_.youngs_modulus =
            cell_value( "youngs_modulus", std::move( youngs_moduli ) );
        mechanics_.poissons_ratio =
            cell_value( "poissons_ratio", std::move( poissons_ratios ) );
        const auto density_property =
            ogs_process == OGS_PROCESS::small_deformation
                ? PHYSICAL_PROPERTY_NAME::density
                : PHYSICAL_PROPERTY_NAME::grain_density;
        auto densities = property_cell_values( density_property );
        if( ogs_process == OGS_PROCESS::thermo_hydro_mechanics )
        {
            compute_solid_heat_capacity( densities );
            mechanics_.thermal_expansivity = cell_value(
                "thermal_expansion_coefficient",
                property_cell_values(
                    PHYSICAL_PROPERTY_NAME::thermal_expansion_coefficient ) );
        }
        mechanics_.solid_density =
            cell_value( physical_property_string( density_property ),
                std::move( densities ) );
        if( ogs_process != OGS_PROCESS::small_deformation )
        {
            mechanics_.biot_coefficient = cell_value( "biot_coefficient",
                property_cell_values(
                    PHYSICAL_PROPERTY_NAME::biot_coefficient ) );
        }
        transfer_initial_stress();
    }

    template < typename Model >
    void OgsExporterImpl< Model >::compute_solid_heat_capacity(
        absl::Span< const double > grain_densities )
    {
        const auto& materials = this->problem()->materials();
        std::vector< double > capacities;
        capacities.reserve( grain_densities.size() );
        for( const auto cell : Indices{ grain_densities } )
        {
            const auto& thermal =
                materials[static_cast< index_t >( material_ids_[cell] )]
                    .thermal;
            const auto volumetric =
                thermal ? thermal->volumetric_heat_capacity : 0.;
            // A zero density with heat capacity is refused by the checker
            capacities.push_back( grain_densities[cell] > 0.
                                      ? volumetric / grain_densities[cell]
                                      : 0. );
        }
        mechanics_.solid_heat_capacity =
            cell_value( "specific_heat_capacity", std::move( capacities ) );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::transfer_initial_stress()
    {
        const auto& model = solid_model();
        const std::array< std::tuple< PHYSICAL_PROPERTY_NAME, local_index_t >,
            2 >
            stress_properties{ { { PHYSICAL_PROPERTY_NAME::initial_stress, 0 },
                { PHYSICAL_PROPERTY_NAME::initial_shear_stress, 3 } } };
        const auto has_stress = absl::c_any_of(
            stress_properties, [&model]( const auto& stress_property ) {
                return model.has_physical_property(
                    std::get< 0 >( stress_property ) );
            } );
        if( !has_stress && process() == OGS_PROCESS::small_deformation )
        {
            return;
        }
        initial_stress_.emplace( bulk().nb_polyhedra(), OgsStressTensor{} );
        if( !has_stress )
        {
            return;
        }
        auto& manager = solid_->mesh().polyhedron_attribute_manager();
        AttributeValues< StressTensor > zero;
        zero.default_value = StressTensor{};
        const auto tensor_id =
            manager
                .template create_attribute< VariableAttribute, StressTensor >(
                    unique_attribute_name( ogs_naming(), manager,
                        "initial_stress_tensor", "initial_stress" ),
                    zero, AttributeProperties{} );
        auto tensor =
            manager.template find_attribute< VariableAttribute, StressTensor >(
                tensor_id );
        const auto factor = stress_factor( *this->problem(), true );
        for( const auto& [property, first_component] : stress_properties )
        {
            if( model.has_physical_property( property ) )
            {
                transfer_stress_components( ogs_naming(), model,
                    solid_->mappings(),
                    model.physical_property_info( property ),
                    physical_property_string( property ), first_component,
                    factor, *tensor );
            }
        }
        for( const auto cell : Range{ bulk().nb_polyhedra() } )
        {
            initial_stress_->at( cell ) =
                ogs_stress_tensor( tensor->value( cell ), 1. );
        }
        manager.delete_attribute( tensor_id );
    }

    template < typename Model >
    void OgsExporterImpl< Model >::compute_initial_values()
    {
        const auto& problem = *this->problem();
        if( process() == OGS_PROCESS::small_deformation )
        {
            return;
        }
        // Only LIQUID_FLOW and HT read a storage
        if( porosity_ && !has_displacement( process().value() ) )
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
        const auto is_thermal = has_temperature( process().value() );
        std::vector< index_t > cell_conditions;
        cell_conditions.reserve( bulk().nb_polyhedra() );
        std::vector< double > cell_temperatures;
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
            if( is_thermal )
            {
                cell_temperatures.push_back(
                    ogs_cell_temperature( condition, cell, context_ ) );
            }
        }
        initial_pressure_ =
            initial_value( "initial_pressure", cell_pressures_ );
        if( is_thermal )
        {
            initial_temperature_ =
                initial_value( "initial_temperature", cell_temperatures );
        }
        if( is_two_phase() )
        {
            compute_two_phase_values( cell_conditions );
        }
    }

    template < typename Model >
    double OgsExporterImpl< Model >::cell_capillary_pressure(
        const OgsCellWaterFraction& cell_fraction ) const
    {
        const auto cell = cell_fraction.cell;
        const auto water_fraction = cell_fraction.water_fraction;
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
                { cell, condition.composition()[water] } ) );
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
                == PHYSICAL_PROPERTY_NAME::boundary_water_fraction )
            {
                convert_water_fractions( boundary );
            }
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::convert_water_fractions(
        OgsBoundaryMesh& boundary ) const
    {
        auto& mesh = boundary.mesh;
        std::vector< double > fractions(
            mesh.points.size(), boundary.value.value_or( 0. ) );
        if( !boundary.value )
        {
            fractions = mesh.point_arrays.front().second;
        }
        auto pressures = vertex_capillary_pressures( mesh, fractions );
        mesh.point_arrays.clear();
        boundary.value = ogs_uniform_value( pressures );
        if( !boundary.value )
        {
            mesh.point_arrays.emplace_back( "values", std::move( pressures ) );
        }
    }

    template < typename Model >
    std::vector< double > OgsExporterImpl< Model >::vertex_capillary_pressures(
        const OgsMeshFile& mesh, absl::Span< const double > fractions ) const
    {
        const auto nb_points = mesh.points.size();
        std::vector< double > sums( nb_points, 0. );
        std::vector< index_t > counts( nb_points, 0 );
        for( const auto cell : Indices{ mesh.cells } )
        {
            const auto bulk_cell =
                static_cast< index_t >( mesh.bulk_element_ids[cell] );
            for( const auto vertex : mesh.cells[cell] )
            {
                sums[vertex] +=
                    cell_capillary_pressure( { bulk_cell, fractions[vertex] } );
                counts[vertex]++;
            }
        }
        for( const auto vertex : Indices{ sums } )
        {
            sums[vertex] /= static_cast< double >( counts[vertex] );
        }
        return sums;
    }

    template < typename Model >
    OgsInitialValue OgsExporterImpl< Model >::initial_value(
        std::string_view field, const std::vector< double >& cell_values )
    {
        if( const auto uniform = ogs_uniform_value( cell_values ) )
        {
            return { uniform, "" };
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
        for( const auto property : ogs_boundary_properties( this->problem() ) )
        {
            build_boundary_meshes( property );
        }
        const auto ogs_process = process();
        if( ogs_process && has_displacement( ogs_process.value() ) )
        {
            build_mechanics_boundary_meshes();
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_mechanics_boundary_meshes()
    {
        const auto& model = solid_model();
        const auto surface_attribute =
            [&model](
                PHYSICAL_PROPERTY_NAME property ) -> std::optional< uuid > {
            if( !model.has_physical_property( property ) )
            {
                return std::nullopt;
            }
            const auto& info = model.physical_property_info( property );
            if( info.component_type != Surface3D::component_type_static() )
            {
                return std::nullopt;
            }
            return info.attribute_id;
        };
        if( const auto attribute = surface_attribute(
                PHYSICAL_PROPERTY_NAME::boundary_displacement ) )
        {
            build_displacement_meshes( attribute.value() );
        }
        if( const auto attribute =
                surface_attribute( PHYSICAL_PROPERTY_NAME::boundary_traction ) )
        {
            build_traction_meshes( attribute.value() );
        }
        if( const auto attribute = surface_attribute(
                PHYSICAL_PROPERTY_NAME::boundary_normal_traction ) )
        {
            build_normal_traction_meshes( attribute.value() );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_displacement_meshes(
        const uuid& attribute_id )
    {
        const auto property = PHYSICAL_PROPERTY_NAME::boundary_displacement;
        for( const auto& [axis_value, surface_ids] :
            ogs_displacement_groups( surfaces_, attribute_id, context_ ) )
        {
            OgsBoundaryMeshBuilder builder{ solid_model(), bulk(),
                solid_->mappings() };
            builder.add_surfaces( surfaces_, surface_ids );
            add_boundary_mesh( property, builder.nb_skipped_polygons(),
                OgsBoundaryMesh{ property, "", axis_value.second,
                    builder.mesh( OGS_VALUES::none ), axis_value.first } );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_traction_meshes(
        const uuid& attribute_id )
    {
        const auto& model = solid_model();
        const auto property = PHYSICAL_PROPERTY_NAME::boundary_traction;
        const auto groups = group_surface_values< std::array< double, 3 > >(
            model, surface_indices_, attribute_id, context_ );
        // One mesh per axis of a non-zero uniform component
        for( const auto& [traction, surface_ids] : groups.uniform_regions )
        {
            for( const auto axis : LRange{ 3 } )
            {
                if( traction[axis] == 0. )
                {
                    continue;
                }
                OgsBoundaryMeshBuilder builder{ model, bulk(),
                    solid_->mappings() };
                builder.set_outward();
                builder.add_surfaces( surfaces_, surface_ids );
                add_boundary_mesh( property, builder.nb_skipped_polygons(),
                    OgsBoundaryMesh{ property, "", traction[axis],
                        builder.mesh( OGS_VALUES::none ), axis } );
            }
        }
        for( const auto& variable : groups.variable_surfaces )
        {
            build_variable_traction_meshes(
                *variable.surface, variable.values );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_variable_traction_meshes(
        const Surface3D& surface,
        absl::Span< const std::array< double, 3 > > tractions )
    {
        const auto property = PHYSICAL_PROPERTY_NAME::boundary_traction;
        for( const auto axis : LRange{ 3 } )
        {
            const auto values = axis_values( tractions, axis );
            if( is_zero( values ) )
            {
                continue;
            }
            OgsBoundaryMeshBuilder builder{ solid_model(), bulk(),
                solid_->mappings() };
            builder.set_outward();
            builder.add_surface_values( surface, values );
            add_boundary_mesh( property, builder.nb_skipped_polygons(),
                OgsBoundaryMesh{ property, "", std::nullopt,
                    builder.mesh( OGS_VALUES::cells ), axis } );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_normal_traction_meshes(
        const uuid& attribute_id )
    {
        const auto& model = solid_model();
        const auto property = PHYSICAL_PROPERTY_NAME::boundary_normal_traction;
        const auto factor = stress_factor( *this->problem(), true );
        auto groups = group_surface_values< double >(
            model, surface_indices_, attribute_id, context_ );
        for( const auto& [value, surface_ids] : groups.uniform_regions )
        {
            OgsBoundaryMeshBuilder builder{ model, bulk(), solid_->mappings() };
            builder.set_outward();
            for( const auto surface_id : surface_ids )
            {
                builder.add_surface( *surfaces_[surface_id] );
            }
            add_boundary_mesh( property, builder.nb_skipped_polygons(),
                OgsBoundaryMesh{ property, "", factor * value,
                    builder.mesh( OGS_VALUES::none ), std::nullopt } );
        }
        for( auto& variable : groups.variable_surfaces )
        {
            for( auto& value : variable.values )
            {
                value *= factor;
            }
            OgsBoundaryMeshBuilder builder{ model, bulk(), solid_->mappings() };
            builder.set_outward();
            builder.add_surface_values( *variable.surface, variable.values );
            add_boundary_mesh( property, builder.nb_skipped_polygons(),
                OgsBoundaryMesh{ property, "", std::nullopt,
                    builder.mesh( OGS_VALUES::cells ), std::nullopt } );
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
        for( const auto& [value, surface_ids] : groups.uniform_regions )
        {
            OgsBoundaryMeshBuilder builder{ model, bulk(), solid_->mappings() };
            builder.set_value( value );
            for( const auto surface_id : surface_ids )
            {
                builder.add_surface( *surfaces_[surface_id] );
            }
            add_boundary_mesh( property, builder.nb_skipped_polygons(),
                OgsBoundaryMesh{ property, "", value,
                    builder.mesh( OGS_VALUES::none ), std::nullopt } );
        }
        for( const auto& variable : groups.variable_surfaces )
        {
            OgsBoundaryMeshBuilder builder{ model, bulk(), solid_->mappings() };
            builder.add_surface_values( *variable.surface, variable.values );
            add_boundary_mesh( property, builder.nb_skipped_polygons(),
                OgsBoundaryMesh{ property, "", std::nullopt,
                    builder.mesh( OGS_VALUES::points ), std::nullopt } );
        }
    }

    template < typename Model >
    void OgsExporterImpl< Model >::add_boundary_mesh(
        PHYSICAL_PROPERTY_NAME property,
        index_t nb_skipped_polygons,
        OgsBoundaryMesh&& boundary )
    {
        const auto property_name = physical_property_string( property );
        if( nb_skipped_polygons != 0 )
        {
            Logger::warning( context_, " ", nb_skipped_polygons,
                " polygons with ", property_name,
                " are not on a Block polyhedron facet, they are not "
                "exported." );
        }
        if( boundary.mesh.cells.empty() )
        {
            return;
        }
        boundary.stem = absl::StrCat(
            this->prefix(), "_", property_name, "_", boundary_meshes_.size() );
        boundary_meshes_.push_back( std::move( boundary ) );
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
            auto well_mesh = build_well_mesh( well_id, aabb );
            if( well_mesh.mesh.points.empty() )
            {
                Logger::warning( context_, " No perforation of well ",
                    well_mesh.name,
                    " is inside the model, no well mesh is exported." );
                continue;
            }
            well_meshes_.push_back( std::move( well_mesh ) );
        }
    }

    template < typename Model >
    OgsWellMesh OgsExporterImpl< Model >::build_well_mesh(
        index_t well_id, const AABBTree3D& aabb ) const
    {
        const auto& name = this->well_names()[well_id];
        return { name, absl::StrCat( this->prefix(), "_well_", name ),
            point_vertices( *this->well_perforations()[well_id], aabb,
                { "Perforation", absl::StrCat( "well ", name ) } ) };
    }

    template < typename Model >
    OgsMeshFile OgsExporterImpl< Model >::point_vertices(
        const PointSet3D& points,
        const AABBTree3D& aabb,
        const OgsPointDescription& description ) const
    {
        OgsMeshFile mesh;
        for( const auto point_id : Range{ points.nb_vertices() } )
        {
            const auto& point = points.point( point_id );
            const auto cell = containing_cell( bulk(), aabb, point );
            if( cell == NO_ID )
            {
                Logger::warning( context_, " ", description.point, " ",
                    point_id, " of ", description.owner, " (", point.string(),
                    ") is not inside the model, it is ignored." );
                continue;
            }
            const auto vertex = nearest_cell_vertex( bulk(), cell, point );
            if( absl::c_linear_search( mesh.bulk_node_ids, vertex ) )
            {
                continue;
            }
            mesh.cells.push_back(
                { static_cast< index_t >( mesh.points.size() ) } );
            mesh.points.push_back( bulk().point( vertex ) );
            mesh.bulk_node_ids.push_back( vertex );
        }
        return mesh;
    }

    template < typename Model >
    void OgsExporterImpl< Model >::build_heat_source_meshes()
    {
        const auto heat_sources = this->problem()->heat_sources();
        for( const auto heat_source_id : Indices{ heat_sources } )
        {
            heat_source_meshes_.push_back(
                heat_sources[heat_source_id].block
                    ? build_block_heat_source_mesh( heat_source_id )
                    : build_point_heat_source_mesh( heat_source_id ) );
        }
    }

    // The points are given to the exporter: checked by the OGS checker
    template < typename Model >
    OgsHeatSourceMesh OgsExporterImpl< Model >::build_point_heat_source_mesh(
        index_t heat_source_id ) const
    {
        const auto& name = this->problem()->heat_sources()[heat_source_id].name;
        const auto names = this->heat_source_names();
        const auto points = static_cast< index_t >(
            absl::c_find( names, name ) - names.begin() );
        const auto aabb = create_aabb_tree( bulk() );
        auto mesh = point_vertices( *this->heat_source_points()[points], aabb,
            { "Point", absl::StrCat( "heat source ", name ) } );
        OpenGeodeGeosciencesIOModelException::check_exception(
            !mesh.points.empty(), nullptr, OpenGeodeException::TYPE::data,
            context_, " No point of heat source ", name,
            " is inside the model." );
        return { heat_source_id, absl::StrCat( this->prefix(), "_heat_", name ),
            std::move( mesh ), 0. };
    }

    // The Block is a Block of the model: checked by the regions checks
    template < typename Model >
    OgsHeatSourceMesh OgsExporterImpl< Model >::build_block_heat_source_mesh(
        index_t heat_source_id ) const
    {
        const auto& source = this->problem()->heat_sources()[heat_source_id];
        const auto offsets = solid_->block_polyhedra_offsets();
        index_t block_index{ 0 };
        for( const auto& block : solid_model().blocks() )
        {
            if( block.id() == source.block.value() )
            {
                break;
            }
            block_index++;
        }
        OgsHeatSourceMesh heat_mesh{ heat_source_id,
            absl::StrCat( this->prefix(), "_heat_", source.name ), {}, 0. };
        auto& mesh = heat_mesh.mesh;
        mesh.volumetric = true;
        absl::flat_hash_map< index_t, index_t > local_vertices;
        for( const auto cell :
            Range{ offsets[block_index], offsets[block_index + 1] } )
        {
            std::vector< index_t > local_cell;
            for( const auto vertex : bulk().polyhedron_vertices( cell ) )
            {
                const auto [local, inserted] = local_vertices.try_emplace(
                    vertex, static_cast< index_t >( mesh.points.size() ) );
                if( inserted )
                {
                    mesh.points.push_back( bulk().point( vertex ) );
                    mesh.bulk_node_ids.push_back( vertex );
                }
                local_cell.push_back( local->second );
            }
            mesh.cells.push_back( std::move( local_cell ) );
            mesh.bulk_element_ids.push_back( cell );
            heat_mesh.volume += ogs_cell_volume( bulk(), cell );
        }
        OpenGeodeGeosciencesIOModelException::check_exception(
            !mesh.cells.empty(), nullptr, OpenGeodeException::TYPE::data,
            context_, " Block ", source.block->string(), " of heat source ",
            source.name, " has no exported cell." );
        return heat_mesh;
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
            add_permeability_attribute(
                cells, permeability_.value(), is_two_phase() );
        }
        for( const auto& [name, values] : cell_arrays_ )
        {
            add_cell_attribute( cells, name, values );
        }
        if( initial_stress_ )
        {
            add_cell_attribute(
                cells, "initial_stress", initial_stress_.value() );
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
        context.heat_sources = heat_source_meshes_;
        context.initial_pressure = initial_pressure_;
        context.initial_capillary_pressure = initial_capillary_pressure_;
        context.initial_temperature = initial_temperature_;
        context.reference_temperature = ogs_reference_temperature( problem );
        context.youngs_modulus = mechanics_.youngs_modulus;
        context.poissons_ratio = mechanics_.poissons_ratio;
        context.solid_density = mechanics_.solid_density;
        context.biot_coefficient = mechanics_.biot_coefficient;
        context.thermal_expansivity = mechanics_.thermal_expansivity;
        context.solid_heat_capacity = mechanics_.solid_heat_capacity;
        if( initial_stress_ )
        {
            context.initial_stress_field = "initial_stress";
        }
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
        for( const auto& heat_source : heat_source_meshes_ )
        {
            save_ogs_mesh_file( absl::StrCat( this->files_directory(), "/",
                                    heat_source.stem, ".vtu" ),
                heat_source.mesh );
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
