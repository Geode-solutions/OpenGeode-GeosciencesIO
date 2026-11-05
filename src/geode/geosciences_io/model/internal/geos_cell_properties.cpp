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

#include <geode/geosciences_io/model/internal/geos_cell_properties.hpp>

#include <typeinfo>

#include <absl/algorithm/container.h>
#include <absl/strings/str_cat.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/mapping.hpp>
#include <geode/basic/range.hpp>

#include <geode/mesh/core/solid_mesh.hpp>

#include <geode/model/helpers/convert_to_mesh.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/solid_mechanics.hpp>
#include <geode/geosciences_io/model/internal/geos_xml.hpp>

namespace geode::internal
{
    double geos_stress_factor( const PhysicalProblem& problem )
    {
        return problem.solid_mechanics()->sign_convention
                       == STRESS_SIGN_CONVENTION::compression_positive
                   ? -1.
                   : 1.;
    }

    void transfer_stress_components( const BRep& model,
        const ModelToMeshMappings& model2solid,
        const PhysicalProperties::Info& property_info,
        std::string_view name,
        local_index_t first_component,
        double factor,
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
                    [&value, first_component, factor](
                        StressTensor& tensor_value ) {
                        for( const auto axis : LRange{ 3 } )
                        {
                            tensor_value[first_component + axis] =
                                factor * value[axis];
                        }
                    } );
            }
        }
    }

    template < typename... Types >
    struct BlockAttributeTypes
    {
    };

    using SupportedBlockAttributeTypes = BlockAttributeTypes< double,
        float,
        int,
        unsigned int,
        long,
        unsigned long,
        long long,
        unsigned long long,
        std::array< double, 2 >,
        std::array< double, 3 >,
        std::array< double, 6 >,
        std::array< double, 9 >,
        std::array< float, 3 > >;

    // Copies the Block values in a VariableAttribute of the solid, the
    // polyhedra of each Block being contiguous in the solid
    template < typename T, typename Model >
    bool copy_typed_block_attribute( const Model& model,
        absl::Span< const index_t > block_polyhedra_offsets,
        std::string_view attribute_type,
        std::string_view attribute_name,
        const uuid& attribute_id,
        AttributeManager& solid_manager )
    {
        if( attribute_type != typeid( T ).name() )
        {
            return false;
        }
        AttributeProperties properties;
        properties.transferable = true;
        solid_manager.create_attribute< VariableAttribute, T >(
            attribute_name, attribute_id, AttributeValues< T >{}, properties );
        auto solid_attribute =
            solid_manager.find_attribute< VariableAttribute, T >(
                attribute_id );
        index_t block_index{ 0 };
        for( const auto& block : model.blocks() )
        {
            const auto& block_mesh = block.mesh();
            const auto block_attribute =
                block_mesh.polyhedron_attribute_manager()
                    .template find_read_only_attribute< T >( attribute_id );
            const auto offset = block_polyhedra_offsets[block_index++];
            for( const auto polyhedron_id : Range{ block_mesh.nb_polyhedra() } )
            {
                solid_attribute->set_value( offset + polyhedron_id,
                    block_attribute->value( polyhedron_id ) );
            }
        }
        return true;
    }

    template < typename Model, typename... Types >
    bool copy_block_attribute( const Model& model,
        absl::Span< const index_t > block_polyhedra_offsets,
        std::string_view attribute_type,
        std::string_view attribute_name,
        const uuid& attribute_id,
        AttributeManager& solid_manager,
        BlockAttributeTypes< Types... > /*unused*/ )
    {
        return (
            copy_typed_block_attribute< Types >( model, block_polyhedra_offsets,
                attribute_type, attribute_name, attribute_id, solid_manager )
            || ... );
    }

    // The VTU writer renames the attributes sharing a name, which would not
    // match the field names given to GEOS
    std::string unique_attribute_name( const AttributeManager& manager,
        std::string_view attribute_name,
        std::string_view geos_name )
    {
        auto valid_name = to_string( attribute_name );
        if( !is_valid_geos_list_item( valid_name ) )
        {
            absl::c_replace_if(
                valid_name,
                []( char character ) {
                    return GEOS_LIST_SEPARATORS.find( character )
                           != std::string_view::npos;
                },
                '_' );
            if( valid_name.empty() )
            {
                valid_name = to_string( geos_name );
            }
            Logger::warning( "[GeosExporter] Attribute name \"", attribute_name,
                "\" cannot be used in GEOS lists, it is "
                "exported as ",
                valid_name, "." );
        }
        const std::string_view name{ valid_name };
        if( !manager.attribute_ids_matching_name( name ) )
        {
            return valid_name;
        }
        auto unique_name = absl::StrCat( name, "_", geos_name );
        for( index_t suffix{ 0 };
            manager.attribute_ids_matching_name( unique_name ); suffix++ )
        {
            unique_name = absl::StrCat( name, "_", geos_name, "_", suffix );
        }
        Logger::warning( "[GeosExporter] Attribute name ", name,
            " is already exported, it is exported as ", unique_name, "." );
        return unique_name;
    }

    std::optional< std::string > transfer_block_attribute( const BRep& model,
        SolidMesh3D& solid,
        absl::Span< const index_t > block_polyhedra_offsets,
        const PhysicalProperties::Info& property_info,
        std::string_view geos_name )
    {
        const auto& attribute_id = property_info.attribute_id;
        OpenGeodeGeosciencesIOModelException::check_exception(
            property_info.component_type == Block3D::component_type_static(),
            nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Physical property attribute ",
            attribute_id.string(),
            " must be defined on Blocks to be exported." );
        std::optional< std::string > attribute_type;
        std::optional< std::string > attribute_name;
        for( const auto& block : model.blocks() )
        {
            const auto& block_manager =
                block.mesh().polyhedron_attribute_manager();
            if( !block_manager.attribute_exists( attribute_id ) )
            {
                Logger::warning( "[GeosExporter] Physical property attribute ",
                    attribute_id.string(), " is missing on Block ",
                    block.id().string(), ", it is not exported." );
                return std::nullopt;
            }
            const auto type = block_manager.attribute_type( attribute_id );
            if( !attribute_type )
            {
                attribute_type.emplace( type );
                attribute_name =
                    block_manager.find_generic_attribute( attribute_id )
                        ->name();
            }
            OpenGeodeGeosciencesIOModelException::check_exception(
                type == attribute_type.value(), nullptr,
                OpenGeodeException::TYPE::data,
                "[GeosExporter] Physical property attribute ",
                attribute_id.string(),
                " has different value types on Blocks." );
        }
        OpenGeodeGeosciencesIOModelException::check_exception(
            attribute_name.has_value(), nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Physical property attribute ",
            attribute_id.string(), " has no name." );
        auto& solid_manager = solid.polyhedron_attribute_manager();
        if( const auto attribute =
                solid_manager.find_generic_attribute( attribute_id ) )
        {
            // Attribute already exported for another physical property
            return attribute->name().value();
        }
        const auto vtu_name = unique_attribute_name(
            solid_manager, attribute_name.value(), geos_name );
        OpenGeodeGeosciencesIOModelException::check_exception(
            copy_block_attribute( model, block_polyhedra_offsets,
                attribute_type.value(), vtu_name, attribute_id, solid_manager,
                SupportedBlockAttributeTypes{} ),
            nullptr, OpenGeodeException::TYPE::data,
            "[GeosExporter] Physical property attribute ",
            attribute_name.value(), " has an unsupported value type." );
        return vtu_name;
    }
} // namespace geode::internal
