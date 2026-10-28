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

#include <geode/geosciences_io/model/internal/geos_simulation_checker.hpp>

#include <cmath>
#include <functional>
#include <string>
#include <typeinfo>
#include <vector>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_set.h>
#include <absl/strings/str_cat.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/constant_attribute.hpp>
#include <geode/basic/range.hpp>

#include <geode/geometry/aabb.hpp>
#include <geode/geometry/point.hpp>

#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/polyhedral_solid.hpp>
#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/surface_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/helpers/aabb_solid_helpers.hpp>
#include <geode/mesh/helpers/detail/create_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/helpers/physical_problem_blocks.hpp>
#include <geode/geosciences/simulation/helpers/simulation_model_checks.hpp>

#include <geode/geosciences_io/model/internal/geos_boundary_conditions.hpp>
#include <geode/geosciences_io/model/internal/geos_regions.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_validation.hpp>
#include <geode/geosciences_io/model/internal/geos_simulation_writer.hpp>
#include <geode/geosciences_io/model/internal/geos_xml.hpp>

namespace
{
    using geode::PHYSICAL_PROPERTY_NAME;

    constexpr double FRACTION_SUM_TOLERANCE{ 1e-6 };

    constexpr std::array< PHYSICAL_PROPERTY_NAME, 5 > BOUNDARY_PROPERTIES{
        PHYSICAL_PROPERTY_NAME::boundary_pressure,
        PHYSICAL_PROPERTY_NAME::boundary_temperature,
        PHYSICAL_PROPERTY_NAME::boundary_oil_fraction,
        PHYSICAL_PROPERTY_NAME::boundary_gas_fraction,
        PHYSICAL_PROPERTY_NAME::boundary_water_fraction
    };

    std::vector< geode::uuid > surface_ids( const geode::BRep& model )
    {
        std::vector< geode::uuid > surfaces;
        for( const auto& surface : model.surfaces() )
        {
            surfaces.push_back( surface.id() );
        }
        return surfaces;
    }

    std::vector< double > polygon_values(
        const geode::Surface3D& surface, const geode::uuid& attribute_id )
    {
        const auto& mesh = surface.mesh();
        const auto attribute =
            mesh.polygon_attribute_manager().find_read_only_attribute< double >(
                attribute_id );
        std::vector< double > values;
        values.reserve( mesh.nb_polygons() );
        for( const auto polygon : geode::Range{ mesh.nb_polygons() } )
        {
            values.push_back( attribute->value( polygon ) );
        }
        return values;
    }

    bool are_equal( absl::Span< const double > values )
    {
        return absl::c_all_of( values, [&values]( double value ) {
            return value == values.front();
        } );
    }

    class GeosExportChecker
    {
    public:
        GeosExportChecker( const geode::BRep& model,
            geode::internal::GeosDeckContext regions,
            const geode::PhysicalProblem* problem,
            std::optional< absl::Span< const geode::internal::GeosWell > >
                wells )
            : model_( model ), problem_( problem ), wells_( wells )
        {
            result_.context = std::move( regions );
        }

        geode::internal::GeosExportCheck run()
        {
            if( problem_ != nullptr )
            {
                geode::check_physical_problem_consistency(
                    *problem_, model_, report() );
            }
            check_mesh();
            check_rock_properties();
            check_boundary_properties();
            if( wells_ )
            {
                check_wells();
            }
            if( problem_ != nullptr )
            {
                check_problem();
            }
            return std::move( result_ );
        }

    private:
        geode::SimulationReport& report()
        {
            return result_.report;
        }

        geode::internal::GeosDeckContext& context()
        {
            return result_.context;
        }

        [[nodiscard]] bool has_fluid() const
        {
            return problem_->has_flow() && problem_->flow().has_fluid();
        }

        [[nodiscard]] bool is_multiphase() const
        {
            return has_fluid()
                   && dynamic_cast< const geode::MultiphaseFlow* >(
                          &problem_->flow() )
                          != nullptr;
        }

        void check_mesh()
        {
            geode::check_simulation_blocks( model_, report() );
            geode::check_simulation_surfaces_on_blocks(
                model_, surface_ids( model_ ), report() );
            geode::check_simulation_unique_vertices_in_blocks(
                model_, report() );
            check_merged_mesh_type();
        }

        void check_merged_mesh_type()
        {
            if( model_.nb_blocks() == 0 )
            {
                return;
            }
            std::vector< std::reference_wrapper< const geode::SolidMesh3D > >
                meshes;
            for( const auto& block : model_.blocks() )
            {
                meshes.emplace_back( block.mesh() );
            }
            const auto merged =
                geode::detail::create_mesh< geode::SolidMesh3D >( meshes );
            if( dynamic_cast< const geode::TetrahedralSolid3D* >( merged.get() )
                    == nullptr
                && dynamic_cast< const geode::HybridSolid3D* >( merged.get() )
                       == nullptr
                && dynamic_cast< const geode::PolyhedralSolid3D* >(
                       merged.get() )
                       == nullptr )
            {
                report().add_error( "The Block meshes merge into a ",
                    merged->type_name().get(),
                    ": GEOS needs tetrahedral, hybrid or polyhedral Block "
                    "meshes" );
            }
        }

        void check_rock_properties()
        {
            geode::PhysicalPropertyRequirement porosity{
                PHYSICAL_PROPERTY_NAME::porosity
            };
            porosity.merged = true;
            porosity_ =
                geode::check_physical_property( model_, porosity, report() );
            geode::PhysicalPropertyRequirement permeability{
                PHYSICAL_PROPERTY_NAME::permeability
            };
            permeability.merged = true;
            permeability.value_types = {
                geode::PHYSICAL_PROPERTY_VALUE_TYPE::vector3
            };
            permeability_ = geode::check_physical_property(
                model_, permeability, report() );
            context().porosity_array = cell_array( *porosity_ );
            context().permeability_array = cell_array( *permeability_ );
            check_imported_names();
        }

        std::optional< geode::internal::GeosCellArray > cell_array(
            const geode::PhysicalPropertyCoverage& coverage )
        {
            if( !coverage.is_carried() )
            {
                return std::nullopt;
            }
            const auto& attribute_id =
                model_.physical_property_info( coverage.property ).attribute_id;
            for( const auto& block : model_.blocks() )
            {
                if( !coverage.carries( block.id() ) )
                {
                    continue;
                }
                const auto attribute =
                    block.mesh()
                        .polyhedron_attribute_manager()
                        .find_generic_attribute( attribute_id );
                if( attribute->properties().time )
                {
                    report().add_error( "The ",
                        geode::physical_property_string( coverage.property ),
                        " attribute is a time series step: GEOS imports "
                        "static properties only" );
                }
                return geode::internal::GeosCellArray{ attribute_id,
                    attribute->name().value_or( "" ) };
            }
            return std::nullopt;
        }

        void check_imported_names()
        {
            std::vector< std::string > names{
                geode::internal::GEOS_REGION_ATTRIBUTE_NAME
            };
            for( const auto* array :
                { &context().permeability_array, &context().porosity_array } )
            {
                if( !*array || ( *array )->name.empty() )
                {
                    continue;
                }
                const auto& name = ( *array )->name;
                if( !geode::internal::is_geos_name( name ) )
                {
                    report().add_error( "Block attribute name \"", name,
                        "\" must only contain letters, digits, '_' and '-' "
                        "to be imported by GEOS" );
                }
                if( absl::c_linear_search( names, name ) )
                {
                    report().add_error( "Block attribute name \"", name,
                        "\" is already used in the GEOS mesh: rename the "
                        "attribute" );
                }
                names.push_back( name );
            }
            imported_names_ = std::move( names );
        }

        void check_boundary_properties()
        {
            exported_boundaries_ =
                geode::internal::geos_boundary_properties( problem_ );
            for( const auto property : BOUNDARY_PROPERTIES )
            {
                // Properties GEOS does not write only get the ignored warning
                geode::SimulationReport ignored_issues;
                boundary_coverages_.push_back( geode::check_physical_property(
                    model_, geode::PhysicalPropertyRequirement{ property },
                    is_exported( property ) ? report() : ignored_issues ) );
                const auto& coverage = boundary_coverages_.back();
                if( property != PHYSICAL_PROPERTY_NAME::boundary_pressure
                    && model_.has_physical_property( property ) )
                {
                    geode::check_physical_property_presence(
                        model_, coverage, "it is not exported", report() );
                }
            }
            context().has_boundary_pressure =
                boundary_coverage( PHYSICAL_PROPERTY_NAME::boundary_pressure )
                    .is_carried();
            context().has_boundary_temperature = boundary_coverage(
                PHYSICAL_PROPERTY_NAME::boundary_temperature )
                                                     .is_carried();
            for( const auto& [property, fraction] :
                geode::internal::geos_fluid_fractions() )
            {
                if( is_exported( property )
                    && boundary_coverage( property ).is_carried() )
                {
                    context().boundary_fraction_phases.push_back(
                        fraction.phase );
                }
            }
            check_boundary_tables();
        }

        [[nodiscard]] bool is_exported( PHYSICAL_PROPERTY_NAME property ) const
        {
            return absl::c_linear_search( exported_boundaries_, property );
        }

        [[nodiscard]] const geode::PhysicalPropertyCoverage& boundary_coverage(
            PHYSICAL_PROPERTY_NAME property ) const
        {
            return boundary_coverages_[static_cast< std::size_t >(
                absl::c_find( BOUNDARY_PROPERTIES, property )
                - BOUNDARY_PROPERTIES.begin() )];
        }

        // Surfaces with non-uniform values are exported as GEOS tables
        void check_boundary_tables()
        {
            for( const auto& surface : model_.surfaces() )
            {
                if( !has_non_uniform_values( surface ) )
                {
                    continue;
                }
                const auto& mesh = surface.mesh();
                std::vector< geode::Point3D > centers;
                centers.reserve( mesh.nb_polygons() );
                for( const auto polygon : geode::Range{ mesh.nb_polygons() } )
                {
                    centers.push_back( mesh.polygon_barycenter( polygon ) );
                }
                if( !geode::internal::separating_table_coordinates( centers ) )
                {
                    report().add_error(
                        geode::simulation_component_description( surface ),
                        " has non-uniform boundary conditions but some of "
                        "its polygons share the same center: GEOS tables "
                        "cannot tell them apart" );
                }
            }
        }

        [[nodiscard]] bool has_non_uniform_values(
            const geode::Surface3D& surface ) const
        {
            for( const auto& coverage : boundary_coverages_ )
            {
                if( !is_exported( coverage.property )
                    || !coverage.covers( surface.id() ) )
                {
                    continue;
                }
                const auto& attribute_id =
                    model_.physical_property_info( coverage.property )
                        .attribute_id;
                if( !are_equal( polygon_values( surface, attribute_id ) ) )
                {
                    return true;
                }
            }
            return false;
        }

        void check_wells()
        {
            std::vector< geode::AABBTree< 3 > > trees;
            for( const auto& block : model_.blocks() )
            {
                trees.push_back( geode::create_aabb_tree( block.mesh() ) );
            }
            for( const auto& well : wells_.value() )
            {
                if( !geode::internal::is_geos_name( well.name ) )
                {
                    report().add_error( "Well name \"", well.name,
                        "\" must only contain letters, digits, '_' and '-'" );
                }
                const auto& points = well.perforations;
                geode::index_t nb_outside{ 0 };
                for( const auto vertex : geode::Range{ points.nb_vertices() } )
                {
                    if( absl::c_none_of(
                            trees, [&points, vertex]( const auto& tree ) {
                                return !tree.containing_boxes(
                                                points.point( vertex ) )
                                            .empty();
                            } ) )
                    {
                        nb_outside++;
                    }
                }
                if( nb_outside == points.nb_vertices() )
                {
                    report().add_error(
                        "Well ", well.name, " has no perforation in a Block" );
                }
                else if( nb_outside != 0 )
                {
                    report().add_warning( "Well ", well.name, ": ", nb_outside,
                        " of ", points.nb_vertices(),
                        " perforations are outside the Blocks and ignored" );
                }
            }
        }

        void check_problem()
        {
            std::optional< absl::Span< const std::string > > well_names;
            if( wells_ )
            {
                for( const auto& well : wells_.value() )
                {
                    well_names_.emplace_back( well.name );
                }
                well_names = absl::Span< const std::string >{ well_names_ };
            }
            geode::internal::validate_geos_problem(
                *problem_, well_names, report() );
            if( !problem_->materials().empty() )
            {
                check_material_rock_values();
            }
            check_region_assignments();
            if( !has_fluid() )
            {
                return;
            }
            check_problem_boundary_conditions();
            geode::internal::validate_geos_deck(
                *problem_, context(), report() );
            check_porosity_copies();
        }

        void check_material_rock_values()
        {
            geode::check_rock_property_availability(
                *problem_, model_, *porosity_, report() );
            geode::check_rock_property_availability(
                *problem_, model_, *permeability_, report() );
            check_unused_material_rock_values();
            for( const auto& material : problem_->materials() )
            {
                if( porosity_->is_carried() && material.porosity )
                {
                    report().add_warning( "Material ", material.name,
                        " porosity is only a GEOS default: the model porosity "
                        "is imported" );
                }
                if( permeability_->is_carried() && material.permeability )
                {
                    report().add_warning( "Material ", material.name,
                        " permeability is only a GEOS default: the model "
                        "permeability is imported" );
                }
            }
        }

        // GEOS writes every material, even those covering no Block
        void check_unused_material_rock_values()
        {
            absl::flat_hash_set< geode::index_t > used;
            for( const auto& [block, material] :
                geode::physical_problem_block_materials( *problem_, model_ ) )
            {
                used.insert( material );
            }
            const auto materials = problem_->materials();
            for( const auto material_id : geode::Indices{ materials } )
            {
                if( used.contains( material_id ) )
                {
                    continue;
                }
                const auto& material = materials[material_id];
                if( !context().porosity_array && !material.porosity )
                {
                    report().add_error( "Material ", material.name,
                        " needs a porosity: GEOS writes it although it covers "
                        "no Block, and the model has no porosity array" );
                }
                if( !context().permeability_array && !material.permeability )
                {
                    report().add_error( "Material ", material.name,
                        " needs a permeability: GEOS writes it although it "
                        "covers no Block, and the model has no permeability "
                        "array" );
                }
            }
        }

        void check_region_assignments()
        {
            result_.materials = geode::internal::assign_geos_materials(
                *problem_, model_, context() );
            result_.initial_conditions =
                geode::internal::assign_geos_initial_conditions(
                    *problem_, model_, context() );
            for( const auto* assignment : { &result_.materials.value(),
                     &result_.initial_conditions.value() } )
            {
                for( const auto& issue : assignment->issues )
                {
                    if( context().regions_are_inferred )
                    {
                        report().add_warning( issue );
                    }
                    else
                    {
                        report().add_error( issue );
                    }
                }
            }
        }

        void check_problem_boundary_conditions()
        {
            std::vector< PHYSICAL_PROPERTY_NAME > ignored;
            std::vector< geode::PhysicalPropertyCoverage > coinciding;
            for( const auto property : BOUNDARY_PROPERTIES )
            {
                if( !is_exported( property ) )
                {
                    ignored.push_back( property );
                }
                else if( boundary_coverage( property ).is_carried() )
                {
                    coinciding.push_back( boundary_coverage( property ) );
                }
            }
            geode::warn_ignored_physical_properties( model_, ignored,
                is_multiphase() ? "GEOS" : "GEOS single phase simulations",
                report() );
            geode::check_coinciding_physical_properties(
                model_, coinciding, report() );
            geode::check_physical_property_presence( model_,
                boundary_coverage( PHYSICAL_PROPERTY_NAME::boundary_pressure ),
                "nothing drives the flow in the GEOS deck", report() );
            check_boundary_ranges();
            if( is_multiphase() )
            {
                check_fraction_sums( coinciding );
            }
        }

        void check_boundary_ranges()
        {
            const auto& fluid = problem_->fluid();
            const auto& pressure =
                boundary_coverage( PHYSICAL_PROPERTY_NAME::boundary_pressure );
            const auto pressures =
                geode::internal::geos_fluid_pressure_range( fluid );
            if( pressures && pressure.value_bounds
                && ( pressure.value_bounds->at( 0 ) < pressures->at( 0 )
                     || pressure.value_bounds->at( 1 ) > pressures->at( 1 ) ) )
            {
                report().add_warning( "Boundary pressures are outside the "
                                      "fluid validity range" );
            }
            const auto& temperature = boundary_coverage(
                PHYSICAL_PROPERTY_NAME::boundary_temperature );
            const auto temperatures =
                geode::internal::geos_fluid_temperature_range( fluid );
            if( temperatures && temperature.value_bounds
                && ( temperature.value_bounds->at( 0 ) < temperatures->at( 0 )
                     || temperature.value_bounds->at( 1 )
                            > temperatures->at( 1 ) ) )
            {
                report().add_warning( "Boundary temperatures are outside the "
                                      "fluid validity range" );
            }
        }

        void check_fraction_sums(
            absl::Span< const geode::PhysicalPropertyCoverage > coverages )
        {
            std::vector< geode::uuid > fraction_ids;
            for( const auto& coverage : coverages )
            {
                if( geode::internal::geos_fluid_fractions().contains(
                        coverage.property ) )
                {
                    fraction_ids.push_back(
                        model_.physical_property_info( coverage.property )
                            .attribute_id );
                }
            }
            if( fraction_ids.empty() )
            {
                return;
            }
            for( const auto& surface : model_.surfaces() )
            {
                check_surface_fraction_sum( surface, fraction_ids );
            }
        }

        void check_surface_fraction_sum( const geode::Surface3D& surface,
            absl::Span< const geode::uuid > fraction_ids )
        {
            const auto& manager = surface.mesh().polygon_attribute_manager();
            std::vector< double > sums( surface.mesh().nb_polygons(), 0. );
            for( const auto& attribute_id : fraction_ids )
            {
                if( !manager.attribute_exists( attribute_id )
                    || manager.attribute_type( attribute_id )
                           != typeid( double ).name() )
                {
                    return;
                }
                const auto values = polygon_values( surface, attribute_id );
                for( const auto polygon : geode::Indices{ sums } )
                {
                    sums[polygon] += values[polygon];
                }
            }
            const auto wrong = absl::c_find_if( sums, []( double sum ) {
                return std::fabs( sum - 1 ) > FRACTION_SUM_TOLERANCE;
            } );
            if( wrong != sums.end() )
            {
                report().add_error(
                    geode::simulation_component_description( surface ),
                    " boundary fluid fractions do not sum to 1 (polygon ",
                    wrong - sums.begin(), ": ", *wrong, ")" );
            }
        }

        void check_porosity_copies()
        {
            if( !context().porosity_array )
            {
                return;
            }
            const auto nb_groups = static_cast< geode::index_t >(
                geode::internal::geos_porosity_groups( *problem_, context() )
                    .groups.size() );
            if( nb_groups < 2 )
            {
                return;
            }
            // The VTU cell arrays are the region attribute and the imported
            // rock properties
            const absl::flat_hash_set< std::string > vtu_names{
                imported_names_.begin(), imported_names_.end()
            };
            for( const auto group : geode::Range{ nb_groups } )
            {
                const auto copy = geode::internal::geos_porosity_vtu_name(
                    context().porosity_array->name, group, nb_groups );
                if( vtu_names.contains( copy ) )
                {
                    report().add_error( "Block attribute name \"", copy,
                        "\" is needed for a copy of the porosity of each GEOS "
                        "porosity model: rename the attribute" );
                }
            }
        }

    private:
        const geode::BRep& model_;
        const geode::PhysicalProblem* problem_;
        std::optional< absl::Span< const geode::internal::GeosWell > > wells_;
        geode::internal::GeosExportCheck result_;
        std::optional< geode::PhysicalPropertyCoverage > porosity_;
        std::optional< geode::PhysicalPropertyCoverage > permeability_;
        std::vector< geode::PhysicalPropertyCoverage > boundary_coverages_;
        /// Boundary properties the GEOS deck uses for this problem
        std::vector< PHYSICAL_PROPERTY_NAME > exported_boundaries_;
        std::vector< std::string > imported_names_;
        std::vector< std::string > well_names_;
    };
} // namespace

namespace geode::internal
{
    GeosExportCheck check_geos_export( const BRep& model,
        GeosDeckContext regions,
        const PhysicalProblem* problem,
        std::optional< absl::Span< const GeosWell > > wells )
    {
        return GeosExportChecker{ model, std::move( regions ), problem, wells }
            .run();
    }

    SimulationReport GeosSimulationChecker::check(
        const PhysicalProblem& problem, const BRep& model ) const
    {
        return check_geos_export( model, geos_inferred_region_context( model ),
            &problem, std::nullopt )
            .report;
    }
} // namespace geode::internal
