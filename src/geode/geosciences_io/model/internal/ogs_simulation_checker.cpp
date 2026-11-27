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

#include <geode/geosciences_io/model/internal/ogs_simulation_checker.hpp>

#include <string>
#include <string_view>
#include <utility>

#include <absl/algorithm/container.h>
#include <absl/container/btree_set.h>
#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/strip.h>
#include <absl/types/span.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/range.hpp>

#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/surface.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/simulation/core/flow_model.hpp>
#include <geode/geosciences/simulation/core/fluid.hpp>
#include <geode/geosciences/simulation/core/heat_source.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/saturation_functions.hpp>
#include <geode/geosciences/simulation/core/schedule.hpp>
#include <geode/geosciences/simulation/core/source_term.hpp>
#include <geode/geosciences/simulation/helpers/internal/physical_problem_blocks.hpp>
#include <geode/geosciences/simulation/helpers/internal/simulation_checker_helpers.hpp>

#include <geode/geosciences_io/model/internal/ogs_simulation_rules.hpp>
#include <geode/geosciences_io/model/internal/simulation_boundary_values.hpp>
#include <geode/geosciences_io/model/internal/simulation_checks.hpp>
#include <geode/geosciences_io/model/internal/simulation_xml.hpp>

namespace
{
    using geode::SIMULATION_CHECK_STAGE;

    constexpr std::string_view NOT_YET{ "not supported by the OGS exporter "
                                        "yet" };
    constexpr std::string_view OGS_NAME_RULE{
        "name must only contain letters, digits, '_' and '-'"
    };

    std::string ogs_surface_name( const geode::Surface3D& surface )
    {
        return surface.name().value_or( surface.id().string() );
    }

    std::string source_subject( const geode::SourceTerm& source )
    {
        return absl::StrCat( "Source term on well \"", source.well, "\"" );
    }

    std::string heat_source_subject( const geode::HeatSource& source )
    {
        return absl::StrCat( "Heat source \"", source.name, "\"" );
    }

    // Values of a double Block property, none if the Block has none
    std::optional< std::vector< double > > block_property_values(
        const geode::BRep& model,
        const geode::Block3D& block,
        geode::PHYSICAL_PROPERTY_NAME property )
    {
        if( !model.has_physical_property( property ) )
        {
            return std::nullopt;
        }
        const auto& info = model.physical_property_info( property );
        const auto& manager = block.mesh().polyhedron_attribute_manager();
        if( info.component_type != geode::Block3D::component_type_static()
            || !manager.attribute_exists( info.attribute_id )
            || manager.attribute_type( info.attribute_id )
                   != typeid( double ).name() )
        {
            return std::nullopt;
        }
        const auto attribute =
            manager.find_read_only_attribute< double >( info.attribute_id );
        std::vector< double > values;
        for( const auto polyhedron :
            geode::Range{ block.mesh().nb_polyhedra() } )
        {
            values.push_back( attribute->value( polyhedron ) );
        }
        return values;
    }

    // Every cell of the Block has a positive value of the property
    bool is_positive_everywhere( const geode::BRep& model,
        const geode::Block3D& block,
        geode::PHYSICAL_PROPERTY_NAME property )
    {
        const auto values = block_property_values( model, block, property );
        return values && absl::c_all_of( values.value(), []( double value ) {
            return value > 0.;
        } );
    }

    // Mean value of each Surface at each of its unique vertices
    using OgsVertexSurfaceValues = absl::flat_hash_map< geode::index_t,
        absl::flat_hash_map< geode::index_t, std::pair< double, double > > >;

    std::optional< std::vector< double > > ogs_surface_values(
        const geode::Surface3D& surface, const geode::uuid& attribute_id )
    {
        try
        {
            return geode::internal::surface_values< double >(
                surface, attribute_id, "[OgsExporter]" );
        }
        catch( const geode::OpenGeodeException& )
        {
            // Wrong attribute types are reported by the property checks
            return std::nullopt;
        }
    }

    struct OgsIndexedSurface
    {
        const geode::Surface3D& surface;
        geode::index_t index;
    };

    void add_ogs_surface_values( const geode::BRep& model,
        const OgsIndexedSurface& indexed,
        absl::Span< const double > values,
        OgsVertexSurfaceValues& vertex_values )
    {
        const auto& surface = indexed.surface;
        const auto& mesh = surface.mesh();
        for( const auto polygon : geode::Range{ mesh.nb_polygons() } )
        {
            const auto value =
                values.size() == 1 ? values.front() : values[polygon];
            for( const auto vertex : mesh.polygon_vertices( polygon ) )
            {
                auto& [sum, count] = vertex_values[model.unique_vertex(
                    { surface.component_id(), vertex } )][indexed.index];
                sum += value;
                count += 1.;
            }
        }
    }

    OgsVertexSurfaceValues ogs_vertex_surface_values(
        const geode::BRep& model, const geode::uuid& attribute_id )
    {
        OgsVertexSurfaceValues vertex_values;
        geode::index_t surface_index{ 0 };
        for( const auto& surface : model.surfaces() )
        {
            if( const auto values =
                    ogs_surface_values( surface, attribute_id ) )
            {
                add_ogs_surface_values( model, { surface, surface_index },
                    values.value(), vertex_values );
            }
            surface_index++;
        }
        return vertex_values;
    }

    using OgsSurfacePairs =
        absl::btree_set< std::pair< geode::index_t, geode::index_t > >;

    // Adds the pairs of Surfaces with different mean values at one vertex
    void add_ogs_conflicts( const absl::flat_hash_map< geode::index_t,
                                std::pair< double, double > >& surface_values,
        OgsSurfacePairs& conflicts )
    {
        for( const auto& [first, first_value] : surface_values )
        {
            for( const auto& [second, second_value] : surface_values )
            {
                if( first < second
                    && first_value.first / first_value.second
                           != second_value.first / second_value.second )
                {
                    conflicts.emplace( first, second );
                }
            }
        }
    }

    // Pairs of Surfaces (by index) with different mean values at a shared
    // vertex
    OgsSurfacePairs ogs_conflicting_surfaces(
        const OgsVertexSurfaceValues& vertex_values )
    {
        OgsSurfacePairs conflicts;
        for( const auto& [vertex, surface_values] : vertex_values )
        {
            geode_unused( vertex );
            add_ogs_conflicts( surface_values, conflicts );
        }
        return conflicts;
    }

    // Attribute of the boundary condition when registered on the Surfaces
    std::optional< geode::uuid > surface_boundary_attribute(
        const geode::BRep& model, geode::PHYSICAL_PROPERTY_NAME property )
    {
        if( !model.has_physical_property( property ) )
        {
            return std::nullopt;
        }
        const auto& info = model.physical_property_info( property );
        if( info.component_type != geode::Surface3D::component_type_static() )
        {
            return std::nullopt;
        }
        return info.attribute_id;
    }

    // The Surface has values of the attribute, as the boundary meshes read
    // them
    bool has_surface_values( const geode::Surface3D& surface,
        const std::optional< geode::uuid >& attribute_id )
    {
        const auto& mesh = surface.mesh();
        return attribute_id && mesh.nb_polygons() != 0
               && mesh.polygon_attribute_manager().attribute_exists(
                   attribute_id.value() );
    }

    // A Surface has a boundary pressure and no boundary temperature
    bool has_pressure_without_temperature( const geode::BRep& model )
    {
        const auto pressure = surface_boundary_attribute(
            model, geode::PHYSICAL_PROPERTY_NAME::boundary_pressure );
        const auto temperature = surface_boundary_attribute(
            model, geode::PHYSICAL_PROPERTY_NAME::boundary_temperature );
        for( const auto& surface : model.surfaces() )
        {
            if( has_surface_values( surface, pressure )
                && !has_surface_values( surface, temperature ) )
            {
                return true;
            }
        }
        return false;
    }

    class OgsProblemChecker
    {
    public:
        OgsProblemChecker( const geode::PhysicalProblem& problem,
            geode::SimulationReport& report )
            : problem_( problem ), report_( report )
        {
        }

        // Wells and heat source points are only checked with the export
        // data
        void check_problem( const geode::internal::OgsExportData* export_data )
        {
            std::optional< absl::Span< const std::string > > well_names;
            std::optional< absl::Span< const std::string > > heat_source_names;
            if( export_data != nullptr )
            {
                well_names = export_data->well_names;
                heat_source_names = export_data->heat_source_names;
            }
            check_solver_options();
            check_material_names();
            check_restart();
            if( !check_process() )
            {
                return;
            }
            check_source_terms( well_names );
            check_heat_sources( heat_source_names );
            if( geode::internal::ogs_process( problem_ )
                == geode::internal::OGS_PROCESS::two_phase_flow_pp )
            {
                check_two_phase();
            }
        }

        void check_boundary_conditions(
            const geode::internal::ModelBoundaryConditions& conditions )
        {
            if( geode::internal::ogs_process( problem_ )
                != geode::internal::OGS_PROCESS::two_phase_flow_pp )
            {
                return;
            }
            if( absl::c_linear_search(
                    conditions.fraction_phases, geode::PHASE::gas ) )
            {
                add_issue( "boundary_gas_fraction",
                    "OGS two-phase flow has no gas phase" );
            }
        }

    private:
        template < typename... Args >
        void add_issue( std::string subject, const Args&... message )
        {
            report_.add_error( SIMULATION_CHECK_STAGE::simulator,
                std::move( subject ), message... );
        }

        // Adds an error and returns false if OGS has no process for the
        // problem
        bool check_process()
        {
            if( geode::internal::ogs_process( problem_ ) )
            {
                return true;
            }
            if( problem_.solid_mechanics() && problem_.has_flow()
                && dynamic_cast< const geode::MultiphaseFlow* >(
                       &problem_.flow() )
                       != nullptr )
            {
                add_issue( "Solid mechanics", "OGS has no multiphase "
                                              "poromechanics process in the "
                                              "OGS exporter" );
                return false;
            }
            if( !problem_.has_flow() || !problem_.flow().has_fluid() )
            {
                add_issue( "Flow model", "the OGS exporter needs a flow model "
                                         "with a fluid" );
                return false;
            }
            if( problem_.thermal()
                && dynamic_cast< const geode::MultiphaseFlow* >(
                       &problem_.flow() )
                       != nullptr )
            {
                add_issue( "Thermal parameters",
                    "OGS thermal two-phase flow is not supported by the OGS "
                    "exporter" );
                return false;
            }
            report_unsupported_fluid();
            return false;
        }

        void report_unsupported_fluid()
        {
            const auto& fluid = problem_.fluid();
            if( dynamic_cast< const geode::Co2BrineFluid* >( &fluid )
                != nullptr )
            {
                add_issue( "Fluid", "Co2BrineFluid is not supported by OGS" );
                return;
            }
            if( dynamic_cast< const geode::DeadOilFluid* >( &fluid )
                == nullptr )
            {
                add_issue( "Flow model", problem_.flow().type_name(), " with ",
                    fluid.type_name(), " is not supported by OGS" );
                return;
            }
            if( absl::c_linear_search( fluid.phases(), geode::PHASE::gas ) )
            {
                add_issue( "Fluid",
                    "the gas phase of DeadOilFluid is not supported by OGS, "
                    "which has no black-oil model: OGS two-phase flow needs "
                    "water and oil only" );
                return;
            }
            add_issue( "Fluid",
                "OGS two-phase flow needs a MultiphaseFlow of water and oil" );
        }

        void check_material_names()
        {
            for( const auto& material : problem_.materials() )
            {
                if( !geode::internal::ogs_naming().is_valid_name(
                        material.name ) )
                {
                    add_issue(
                        geode::internal::material_subject( material.name ),
                        OGS_NAME_RULE );
                }
            }
        }

        void check_restart()
        {
            const auto& output = problem_.output();
            if( output.restart || output.restart_frequency )
            {
                report_.add_warning( SIMULATION_CHECK_STAGE::simulator,
                    "Output settings",
                    "restart files are not written by the OGS exporter" );
            }
        }

        void check_source_terms(
            std::optional< absl::Span< const std::string > > well_names )
        {
            const auto is_thermal = geode::internal::has_temperature(
                geode::internal::ogs_process( problem_ ).value() );
            for( const auto& source : problem_.source_terms() )
            {
                if( well_names )
                {
                    geode::internal::check_source_term_well(
                        source, well_names.value(), report_ );
                }
                if( !geode::internal::ogs_naming().is_valid_name(
                        source.well ) )
                {
                    add_issue(
                        source_subject( source ), "well ", OGS_NAME_RULE );
                }
                // The Nodal source term is on the pressure only
                if( is_thermal && source.rate > 0 )
                {
                    report_.add_warning( SIMULATION_CHECK_STAGE::simulator,
                        source_subject( source ),
                        "OGS injects the fluid at the local temperature" );
                }
            }
        }

        // Point heat sources are located by the points given to the exporter
        void check_heat_sources(
            std::optional< absl::Span< const std::string > > names )
        {
            for( const auto& source : problem_.heat_sources() )
            {
                if( !geode::internal::ogs_naming().is_valid_name(
                        source.name ) )
                {
                    add_issue( heat_source_subject( source ), OGS_NAME_RULE );
                }
                if( source.block || !names
                    || absl::c_linear_search( names.value(), source.name ) )
                {
                    continue;
                }
                geode::SimulationIssue issue{
                    geode::SIMULATION_ISSUE_SEVERITY::error,
                    SIMULATION_CHECK_STAGE::simulator,
                    heat_source_subject( source ),
                    "the heat source has no points given to the exporter"
                };
                issue.fix = "add the heat source points to the exporter";
                report_.add_issue( std::move( issue ) );
            }
        }

        void check_two_phase()
        {
            if( !problem_.source_terms().empty() )
            {
                add_issue( "Source terms",
                    "source terms in two-phase flow are ", NOT_YET );
            }
            for( const auto& material : problem_.materials() )
            {
                const auto subject =
                    geode::internal::material_subject( material.name );
                if( material.permeability
                    && !geode::internal::is_isotropic(
                        material.permeability.value() ) )
                {
                    add_issue( subject, "OGS two-phase flow needs an "
                                        "isotropic permeability" );
                }
                if( material.relative_permeability )
                {
                    check_law( subject, [&material] {
                        geode_unused(
                            geode::internal::relative_permeability_curves(
                                *material.relative_permeability ) );
                    } );
                }
                if( !material.capillary_pressure )
                {
                    geode::SimulationIssue issue{
                        geode::SIMULATION_ISSUE_SEVERITY::error,
                        SIMULATION_CHECK_STAGE::simulator, subject,
                        "OGS two-phase flow needs a capillary pressure: the "
                        "capillary pressure is an unknown of OGS"
                    };
                    issue.fix = "add a capillary pressure to the material";
                    report_.add_issue( std::move( issue ) );
                    continue;
                }
                check_law( subject, [&material] {
                    geode_unused( geode::internal::saturation_curve(
                        *material.capillary_pressure ) );
                } );
            }
        }

        template < typename Conversion >
        void check_law( const std::string& subject, Conversion conversion )
        {
            try
            {
                conversion();
            }
            catch( const geode::OpenGeodeException& exception )
            {
                add_issue( subject,
                    absl::StripPrefix( exception.what(), "[OgsExporter] " ) );
            }
        }

        void check_solver_options()
        {
            geode::internal::check_solver_option_keys( problem_,
                geode::internal::ogs_solver_options(),
                [this]( std::string_view message ) {
                    add_issue( "Solver options", message );
                } );
        }

    private:
        const geode::PhysicalProblem& problem_;
        geode::SimulationReport& report_;
    };
} // namespace

namespace geode::internal
{
    OgsSimulationChecker::OgsSimulationChecker(
        const PhysicalProblem& problem, const BRep& model )
        : SimulationChecker( problem, model )
    {
    }

    OgsSimulationChecker::OgsSimulationChecker( const PhysicalProblem& problem,
        const BRep& model,
        OgsExportData export_data )
        : SimulationChecker( problem, model ),
          export_data_( std::move( export_data ) )
    {
    }

    void OgsSimulationChecker::check_simulator()
    {
        report().set_simulator( "OGS" );
        if( export_data_ )
        {
            for( const auto& name : export_data_->well_names )
            {
                if( !ogs_naming().is_valid_name( name ) )
                {
                    report().add_error( SIMULATION_CHECK_STAGE::simulator,
                        absl::StrCat( "Well \"", name, "\"" ), "well ",
                        OGS_NAME_RULE );
                }
            }
        }
        OgsProblemChecker problem_checker{ problem(), report() };
        problem_checker.check_problem(
            export_data_ ? &export_data_.value() : nullptr );
        check_blocks();
        check_unique_vertices_in_blocks();
        check_cells();
        const auto process = ogs_process( problem() );
        if( !process )
        {
            return;
        }
        check_properties();
        problem_checker.check_boundary_conditions(
            model_boundary_conditions( model() ) );
        check_conflicting_boundaries( process.value() );
    }

    void OgsSimulationChecker::check_conflicting_boundaries(
        OGS_PROCESS process )
    {
        if( process != OGS_PROCESS::small_deformation )
        {
            check_conflicting_boundaries(
                PHYSICAL_PROPERTY_NAME::boundary_pressure );
        }
        if( process == OGS_PROCESS::two_phase_flow_pp )
        {
            check_conflicting_boundaries(
                PHYSICAL_PROPERTY_NAME::boundary_water_fraction );
        }
        if( has_temperature( process ) )
        {
            check_conflicting_boundaries(
                PHYSICAL_PROPERTY_NAME::boundary_temperature );
        }
    }

    void OgsSimulationChecker::check_conflicting_boundaries(
        PHYSICAL_PROPERTY_NAME property )
    {
        const auto& model = this->model();
        if( !model.has_physical_property( property ) )
        {
            return;
        }
        const auto& info = model.physical_property_info( property );
        if( info.component_type != Surface3D::component_type_static() )
        {
            return;
        }
        const auto conflicts = ogs_conflicting_surfaces(
            ogs_vertex_surface_values( model, info.attribute_id ) );
        std::vector< const Surface3D* > surfaces;
        for( const auto& surface : model.surfaces() )
        {
            surfaces.push_back( &surface );
        }
        for( const auto& [first, second] : conflicts )
        {
            report().add_warning( SIMULATION_CHECK_STAGE::simulator,
                to_string( physical_property_string( property ) ), "Surfaces ",
                ogs_surface_name( *surfaces[first] ), " and ",
                ogs_surface_name( *surfaces[second] ),
                " share vertices with different values: OGS keeps one of "
                "them at each shared vertex" );
        }
    }

    void OgsSimulationChecker::check_cells()
    {
        for( const auto& block : model().blocks() )
        {
            const auto& mesh = block.mesh();
            std::vector< index_t > unsupported;
            for( const auto polyhedron : Range{ mesh.nb_polyhedra() } )
            {
                if( !is_ogs_polyhedron( mesh, polyhedron ) )
                {
                    unsupported.push_back( polyhedron );
                }
            }
            if( unsupported.empty() )
            {
                continue;
            }
            SimulationIssue issue{ SIMULATION_ISSUE_SEVERITY::error,
                SIMULATION_CHECK_STAGE::mesh,
                absl::StrCat(
                    "Block ", block.name().value_or( block.id().string() ) ),
                absl::StrCat( unsupported.size(),
                    " polyhedra are not tetrahedra, hexahedra, prisms or "
                    "pyramids, the only cells OGS reads" ) };
            auto& location = issue.add_location( block );
            location.element_type = SIMULATION_ELEMENT_TYPE::polyhedron;
            location.elements = std::move( unsupported );
            location.nb_elements = mesh.nb_polyhedra();
            issue.fix = "remesh the Block with tetrahedra or hexahedra";
            report().add_issue( std::move( issue ) );
        }
    }

    // OGS reads a model porosity or permeability on every cell, merged from
    // the Block meshes. Without such a property, the project uses the value
    // of each material.
    void OgsSimulationChecker::check_properties()
    {
        const auto process = ogs_process( problem() ).value();
        if( has_displacement( process ) )
        {
            check_internal_tractions();
        }
        if( process == OGS_PROCESS::thermo_hydro_mechanics )
        {
            check_heat_capacity_density();
        }
        check_boundary_temperature();
        if( process == OGS_PROCESS::small_deformation )
        {
            return;
        }
        std::vector< PhysicalPropertyCoverage > coverages;
        for( const auto property : { PHYSICAL_PROPERTY_NAME::porosity,
                 PHYSICAL_PROPERTY_NAME::permeability } )
        {
            PhysicalPropertyRequirement requirement{ property };
            requirement.merged = true;
            coverages.push_back( check_physical_property( requirement ) );
        }
        check_rock_property_availability( coverages );
        const auto pressure =
            check_physical_property( PhysicalPropertyRequirement{
                PHYSICAL_PROPERTY_NAME::boundary_pressure } );
        if( problem().source_terms().empty() )
        {
            check_physical_property_presence(
                pressure, "nothing drives the flow in the OGS project" );
        }
        if( ogs_process( problem() ) == OGS_PROCESS::two_phase_flow_pp )
        {
            check_physical_property( PhysicalPropertyRequirement{
                PHYSICAL_PROPERTY_NAME::boundary_water_fraction } );
            warn_ignored_physical_properties(
                { PHYSICAL_PROPERTY_NAME::boundary_oil_fraction } );
        }
    }

    void OgsSimulationChecker::check_internal_tractions()
    {
        for( const auto property : { PHYSICAL_PROPERTY_NAME::boundary_traction,
                 PHYSICAL_PROPERTY_NAME::boundary_normal_traction } )
        {
            const auto attribute =
                surface_boundary_attribute( model(), property );
            for( const auto& surface : model().surfaces() )
            {
                if( !has_surface_values( surface, attribute )
                    || model().nb_incidences( surface.id() ) < 2 )
                {
                    continue;
                }
                SimulationIssue issue{ SIMULATION_ISSUE_SEVERITY::error,
                    SIMULATION_CHECK_STAGE::simulator,
                    to_string( physical_property_string( property ) ),
                    absl::StrCat( "Surface ", ogs_surface_name( surface ),
                        " is between two Blocks: OGS loads a traction on the "
                        "outer side of a model boundary" ) };
                issue.add_location( surface );
                issue.fix = "load the traction on model boundary Surfaces only";
                report().add_issue( std::move( issue ) );
            }
        }
    }

    // The Solid density of THERMO_HYDRO_MECHANICS is the grain density, and
    // its specific heat capacity C_v / grain density
    void OgsSimulationChecker::check_heat_capacity_density()
    {
        const auto block_materials =
            physical_problem_block_materials( problem(), model() );
        for( const auto& block : model().blocks() )
        {
            const auto material = block_materials.find( block.id() );
            if( material == block_materials.end() )
            {
                continue;
            }
            const auto& thermal =
                problem().materials()[material->second].thermal;
            if( !thermal || thermal->volumetric_heat_capacity <= 0.
                || is_positive_everywhere(
                    model(), block, PHYSICAL_PROPERTY_NAME::grain_density ) )
            {
                continue;
            }
            SimulationIssue issue{ SIMULATION_ISSUE_SEVERITY::error,
                SIMULATION_CHECK_STAGE::simulator,
                absl::StrCat(
                    "Block ", block.name().value_or( block.id().string() ) ),
                "OGS needs a positive grain_density to carry the rock heat "
                "capacity" };
            issue.add_location( block );
            issue.fix = "give the Block a positive grain_density";
            report().add_issue( std::move( issue ) );
        }
    }

    // Only HT reads boundary_temperature. On a pressure boundary without
    // one, OGS keeps a zero conductive flux and inflowing fluid takes the
    // temperature of the boundary vertex: warned once.
    void OgsSimulationChecker::check_boundary_temperature()
    {
        if( !has_temperature( ogs_process( problem() ).value() ) )
        {
            warn_ignored_physical_properties(
                { PHYSICAL_PROPERTY_NAME::boundary_temperature } );
            return;
        }
        check_physical_property( PhysicalPropertyRequirement{
            PHYSICAL_PROPERTY_NAME::boundary_temperature } );
        if( has_pressure_without_temperature( model() ) )
        {
            report().add_warning( SIMULATION_CHECK_STAGE::simulator,
                "Boundary conditions",
                "OGS lets the temperature evolve on pressure boundaries "
                "without boundary_temperature" );
        }
    }
} // namespace geode::internal
