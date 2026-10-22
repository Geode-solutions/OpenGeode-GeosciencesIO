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

#include <geode/geosciences_io/model/internal/brep_fem_output.hpp>

#include <limits>
#include <string_view>

#include <absl/container/flat_hash_map.h>

#include <geode/basic/attribute.hpp>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/io.hpp>
#include <geode/basic/types.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/solid_edges.hpp>
#include <geode/mesh/core/solid_facets.hpp>
#include <geode/mesh/core/surface_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>
#include <geode/mesh/io/tetrahedral_solid_output.hpp>

#include <geode/model/helpers/convert_model_meshes.hpp>
#include <geode/model/helpers/convert_to_mesh.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/line.hpp>

namespace
{

    class BRepFemOutputImpl
    {
        static constexpr auto APERTURE_ATTRIBUTE_NAME =
            "diagres_discontinuity_aperture";
        static constexpr auto CONDUIT_AREA_ATTRIBUTE_NAME =
            "diagres_conduit_area";
        static constexpr auto CONDUCTIVITY_ATTRIBUTE_NAME =
            "diagres_conductivity";
        static constexpr auto SURFACE_NAME_ATTRIBUTE = "surface_name";
        static constexpr auto LINE_NAME_ATTRIBUTE = "line_name";

    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };
        BRepFemOutputImpl( std::string_view filename, const geode::BRep& brep )
            : filename_{ filename }, brep_( brep )
        {
        }

        void write_file()
        {
            geode::Logger::info(
                "[BRepSolidFemOutput::write] Writing fem file." );
            auto [solid, model_to_mesh_mapping] =
                geode::convert_brep_into_solid( brep_ );
            const auto* tet_solid_ptr =
                dynamic_cast< const geode::TetrahedralSolid3D* >( solid.get() );
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                tet_solid_ptr != nullptr && solid->nb_polyhedra() != 0, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[BRepFemOutput::write] Only BReps whose Blocks are all "
                "non-empty TetrahedralSolid meshes can be saved in FEM "
                "format" );
            const auto& tet_solid = *tet_solid_ptr;
            copy_tetrahedra_attributes_in_solid(
                tet_solid, model_to_mesh_mapping );
            paint_polyhedron_block_id_attributes(
                tet_solid, model_to_mesh_mapping );
            paint_vertices_block_id_attributes(
                tet_solid, model_to_mesh_mapping );
            paint_surface_infos( tet_solid, model_to_mesh_mapping );
            paint_line_infos( tet_solid, model_to_mesh_mapping );
            geode::save_tetrahedral_solid( tet_solid, filename_ );
        }

    private:
        static double attribute_value( const geode::AttributeBase& attribute,
            const geode::ReadOnlyAttribute< double >* double_attribute,
            geode::index_t element )
        {
            if( double_attribute )
            {
                return double_attribute->value( element );
            }
            return attribute.generic_value( element );
        }

        void paint_facets_attribute( const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping,
            std::string_view attribute_name,
            double default_value ) const
        {
            geode::AttributeValues< double > solid_attribute_values;
            solid_attribute_values.default_value = default_value;
            solid_attribute_values.no_value = default_value;
            geode::AttributeProperties solid_attribute_properties;
            solid_attribute_properties.assignable = false;
            solid_attribute_properties.interpolable = false;
            solid_attribute_properties.transferable = true;
            auto solid_attribute_id =
                tet_solid.facets()
                    .facet_attribute_manager()
                    .create_attribute< geode::VariableAttribute, double >(
                        attribute_name, solid_attribute_values,
                        solid_attribute_properties );
            auto solid_attribute =
                tet_solid.facets()
                    .facet_attribute_manager()
                    .find_attribute< geode::VariableAttribute, double >(
                        solid_attribute_id );
            geode::AttributeValues< std::string >
                solid_surface_name_attribute_values;
            solid_surface_name_attribute_values.default_value = "No_name";
            solid_surface_name_attribute_values.no_value = "No_name";
            auto solid_surface_name_attribute_id =
                tet_solid.facets()
                    .facet_attribute_manager()
                    .create_attribute< geode::VariableAttribute, std::string >(
                        SURFACE_NAME_ATTRIBUTE,
                        solid_surface_name_attribute_values,
                        solid_attribute_properties );
            auto solid_surface_name_attribute =
                tet_solid.facets()
                    .facet_attribute_manager()
                    .find_attribute< geode::VariableAttribute, std::string >(
                        solid_surface_name_attribute_id );
            for( const auto& surface : brep_.surfaces() )
            {
                const auto& mesh = surface.mesh();
                const auto attribute_ids =
                    mesh.polygon_attribute_manager()
                        .attribute_ids_matching_name( attribute_name );
                if( !attribute_ids.has_value() )
                {
                    continue;
                }
                const auto surface_attribute =
                    mesh.polygon_attribute_manager().find_generic_attribute(
                        attribute_ids.value().front() );
                const auto* surface_double_attribute =
                    dynamic_cast< const geode::ReadOnlyAttribute< double >* >(
                        surface_attribute.get() );
                geode::index_t nb_skipped{ 0 };
                for( const auto polygon : geode::Range{ mesh.nb_polygons() } )
                {
                    const geode::MeshElement element{ surface.id(), polygon };
                    if( !model_to_mesh_mapping.surface_polygons_mapping
                            .has_mapping_input( element ) )
                    {
                        nb_skipped++;
                        continue;
                    }
                    const auto facet_in_solid =
                        model_to_mesh_mapping.surface_polygons_mapping.in2out(
                            element );
                    solid_attribute->set_value( facet_in_solid,
                        attribute_value( *surface_attribute,
                            surface_double_attribute, polygon ) );
                    if( const auto name = surface.name() )
                    {
                        solid_surface_name_attribute->set_value(
                            facet_in_solid, name.value() );
                    }
                }
                if( nb_skipped != 0 )
                {
                    geode::Logger::warning( "[BRepFemOutput::write] ",
                        nb_skipped, " polygons of Surface ",
                        surface.id().string(),
                        " are not facets of the tetrahedral mesh, their ",
                        attribute_name, " values are ignored" );
                }
            }
        }

        void paint_surface_infos( const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            for( const auto& surface : brep_.surfaces() )
            {
                if( surface.mesh()
                        .polygon_attribute_manager()
                        .attribute_ids_matching_name( APERTURE_ATTRIBUTE_NAME )
                        .has_value() )
                {
                    paint_facets_attribute( tet_solid, model_to_mesh_mapping,
                        APERTURE_ATTRIBUTE_NAME, -1.0 );
                    break;
                }
            }
            paint_facets_attribute( tet_solid, model_to_mesh_mapping,
                CONDUCTIVITY_ATTRIBUTE_NAME, -1.0 );
        }

        void paint_line_attribute( const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping,
            std::string_view attribute_name,
            double default_value ) const
        {
            geode::AttributeValues< double > solid_attribute_values;
            solid_attribute_values.default_value = default_value;
            solid_attribute_values.no_value = default_value;
            geode::AttributeProperties solid_attribute_properties;
            solid_attribute_properties.assignable = false;
            solid_attribute_properties.interpolable = false;
            solid_attribute_properties.transferable = true;
            const auto solid_attribute_id =
                tet_solid.edges()
                    .edge_attribute_manager()
                    .create_attribute< geode::VariableAttribute, double >(
                        attribute_name, solid_attribute_values,
                        solid_attribute_properties );
            auto solid_attribute =
                tet_solid.edges()
                    .edge_attribute_manager()
                    .find_attribute< geode::VariableAttribute, double >(
                        solid_attribute_id );
            geode::AttributeValues< std::string >
                solid_line_name_attribute_values;
            solid_line_name_attribute_values.default_value = "No_name";
            solid_line_name_attribute_values.no_value = "No_name";
            const auto solid_line_name_attribute_id =
                tet_solid.edges()
                    .edge_attribute_manager()
                    .create_attribute< geode::VariableAttribute, std::string >(
                        LINE_NAME_ATTRIBUTE, solid_line_name_attribute_values,
                        solid_attribute_properties );
            auto solid_line_name_attribute =
                tet_solid.edges()
                    .edge_attribute_manager()
                    .find_attribute< geode::VariableAttribute, std::string >(
                        solid_line_name_attribute_id );
            for( const auto& line : brep_.lines() )
            {
                const auto& mesh = line.mesh();
                const auto attribute_ids =
                    mesh.edge_attribute_manager().attribute_ids_matching_name(
                        attribute_name );
                if( !attribute_ids.has_value() )
                {
                    continue;
                }
                const auto line_attribute =
                    mesh.edge_attribute_manager().find_generic_attribute(
                        attribute_ids.value().front() );
                const auto* line_double_attribute =
                    dynamic_cast< const geode::ReadOnlyAttribute< double >* >(
                        line_attribute.get() );
                geode::index_t nb_skipped{ 0 };
                for( const auto edge : geode::Range{ mesh.nb_edges() } )
                {
                    const geode::MeshElement element{ line.id(), edge };
                    if( !model_to_mesh_mapping.line_edges_mapping
                            .has_mapping_input( element ) )
                    {
                        nb_skipped++;
                        continue;
                    }
                    const auto edge_in_solid =
                        model_to_mesh_mapping.line_edges_mapping.in2out(
                            element );
                    solid_attribute->set_value(
                        edge_in_solid, attribute_value( *line_attribute,
                                           line_double_attribute, edge ) );
                    if( const auto name = line.name() )
                    {
                        solid_line_name_attribute->set_value(
                            edge_in_solid, name.value() );
                    }
                }
                if( nb_skipped != 0 )
                {
                    geode::Logger::warning( "[BRepFemOutput::write] ",
                        nb_skipped, " edges of Line ", line.id().string(),
                        " are not edges of the tetrahedral mesh, their ",
                        attribute_name, " values are ignored" );
                }
            }
        }

        void paint_line_infos( const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            for( const auto& line : brep_.lines() )
            {
                if( line.mesh()
                        .edge_attribute_manager()
                        .attribute_ids_matching_name(
                            CONDUIT_AREA_ATTRIBUTE_NAME )
                        .has_value() )
                {
                    paint_line_attribute( tet_solid, model_to_mesh_mapping,
                        CONDUIT_AREA_ATTRIBUTE_NAME, -1.0 );
                    break;
                }
            }
            paint_line_attribute( tet_solid, model_to_mesh_mapping,
                CONDUCTIVITY_ATTRIBUTE_NAME, -1.0 );
        }

        void paint_polyhedron_block_id_attributes(
            const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            geode::AttributeValues< std::string > block_id_attribute_values;
            block_id_attribute_values.default_value = "No_name";
            block_id_attribute_values.no_value = "No_name";
            geode::AttributeProperties block_id_attribute_properties;
            block_id_attribute_properties.assignable = false;
            block_id_attribute_properties.interpolable = false;
            block_id_attribute_properties.transferable = true;
            auto attribute_p_id =
                tet_solid.polyhedron_attribute_manager()
                    .create_attribute< geode::VariableAttribute, std::string >(
                        "Block_ID_polyhedron", block_id_attribute_values,
                        block_id_attribute_properties );
            auto attribute_p =
                tet_solid.polyhedron_attribute_manager()
                    .find_attribute< geode::VariableAttribute, std::string >(
                        attribute_p_id );
            for( const auto& block : brep_.blocks() )
            {
                for( const auto polyhedron_id :
                    geode::Range{ block.mesh().nb_polyhedra() } )
                {
                    if( const auto name = block.name() )
                    {
                        attribute_p->set_value(
                            model_to_mesh_mapping.solid_polyhedra_mapping
                                .in2out( { block.id(), polyhedron_id } ),
                            name.value() );
                    }
                }
            }
        }

        void paint_vertices_block_id_attributes(
            const geode::TetrahedralSolid3D& tet_solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            auto attribute_v_id =
                tet_solid.vertex_attribute_manager()
                    .create_attribute< geode::VariableAttribute,
                        std::vector< std::string > >(
                        "Block_ID_vertex", {}, geode::AttributeProperties{} );
            auto attribute_v =
                tet_solid.vertex_attribute_manager()
                    .find_attribute< geode::VariableAttribute,
                        std::vector< std::string > >( attribute_v_id );
            for( const auto& block : brep_.blocks() )
            {
                for( const auto vertex_id :
                    geode::Range{ block.mesh().nb_vertices() } )
                {
                    const auto unique_vertex = brep_.unique_vertex(
                        { block.component_id(), vertex_id } );
                    const auto vertex_out =
                        model_to_mesh_mapping.unique_vertices_mapping.in2out(
                            unique_vertex );
                    auto vertex_blocks = attribute_v->value( vertex_out );
                    vertex_blocks.push_back(
                        block.name().value_or( block.id().string() ) );
                    attribute_v->set_value( vertex_out, vertex_blocks );
                }
            }
        }

        template < typename T >
        static void copy_block_attribute_values(
            const geode::AttributeBase& block_attribute,
            const geode::Block3D& block,
            const geode::ModelToMeshMappings& model_to_mesh_mapping,
            geode::VariableAttribute< double >& solid_attribute )
        {
            const auto& typed_attribute =
                dynamic_cast< const geode::ReadOnlyAttribute< T >& >(
                    block_attribute );
            for( const auto polyhedron :
                geode::Range{ block.mesh().nb_polyhedra() } )
            {
                solid_attribute.set_value(
                    model_to_mesh_mapping.solid_polyhedra_mapping.in2out(
                        { block.id(), polyhedron } ),
                    static_cast< double >(
                        typed_attribute.value( polyhedron ) ) );
            }
        }

        void copy_tetrahedra_attributes_in_solid(
            const geode::TetrahedralSolid3D& solid,
            const geode::ModelToMeshMappings& model_to_mesh_mapping ) const
        {
            absl::flat_hash_map< std::string,
                std::shared_ptr< geode::VariableAttribute< double > > >
                solid_attributes;
            for( const auto& block : brep_.blocks() )
            {
                const auto& manager =
                    block.mesh().polyhedron_attribute_manager();
                for( const auto& attribute_id : manager.attribute_ids() )
                {
                    const auto block_attribute =
                        manager.find_generic_attribute( attribute_id );
                    if( !block_attribute
                        || !block_attribute->properties().transferable )
                    {
                        continue;
                    }
                    const auto type = block_attribute->type();
                    const auto is_double = type == typeid( double ).name();
                    if( !is_double && type != typeid( float ).name() )
                    {
                        continue;
                    }
                    auto& solid_attribute =
                        solid_attributes[block_attribute->name().value()];
                    if( !solid_attribute )
                    {
                        geode::AttributeValues< double > solid_attribute_values;
                        solid_attribute_values.default_value =
                            std::numeric_limits< double >::quiet_NaN();
                        solid_attribute_values.no_value =
                            solid_attribute_values.default_value;
                        const auto solid_attribute_id =
                            solid.polyhedron_attribute_manager()
                                .create_attribute< geode::VariableAttribute,
                                    double >( block_attribute->name().value(),
                                    solid_attribute_values,
                                    block_attribute->properties() );
                        solid_attribute =
                            solid.polyhedron_attribute_manager()
                                .find_attribute< geode::VariableAttribute,
                                    double >( solid_attribute_id );
                    }
                    if( is_double )
                    {
                        copy_block_attribute_values< double >( *block_attribute,
                            block, model_to_mesh_mapping, *solid_attribute );
                    }
                    else
                    {
                        copy_block_attribute_values< float >( *block_attribute,
                            block, model_to_mesh_mapping, *solid_attribute );
                    }
                }
            }
        }

    private:
        std::string_view filename_;
        const geode::BRep& brep_;
    };
} // namespace

namespace geode::internal
{
    std::vector< std::string > BRepFemOutput::write( const BRep& brep ) const
    {
        BRepFemOutputImpl impl{ filename(), brep };
        impl.write_file();
        return { to_string( filename() ) };
    }

    bool BRepFemOutput::is_saveable( const BRep& brep ) const
    {
        if( brep.nb_blocks() == 0 )
        {
            return false;
        }
        for( const auto& block : brep.blocks() )
        {
            if( block.mesh().type_name()
                    != TetrahedralSolid3D::type_name_static()
                || block.mesh().nb_polyhedra() == 0 )
            {
                return false;
            }
        }
        for( const auto unique_vertex : Range{ brep.nb_unique_vertices() } )
        {
            if( !brep.has_component_mesh_vertices(
                    unique_vertex, Block3D::component_type_static() ) )
            {
                return false;
            }
        }
        return true;
    }
} // namespace geode::internal
