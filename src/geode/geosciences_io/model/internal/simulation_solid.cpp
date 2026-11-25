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

#include <geode/geosciences_io/model/internal/simulation_solid.hpp>

#include <tuple>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/mapping.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/mesh/core/solid_mesh.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences_io/model/internal/simulation_cell_properties.hpp>
#include <geode/geosciences_io/model/internal/simulation_regions.hpp>

namespace
{
    constexpr auto REGION_ID_ATTRIBUTE_NAME = "attribute";

    std::vector< geode::index_t > compute_block_polyhedra_offsets(
        const geode::BRep& model,
        const geode::ModelToMeshMappings& mappings,
        const geode::SolidMesh3D& solid,
        std::string_view context )
    {
        const auto& polyhedra_mapping = mappings.solid_polyhedra_mapping;
        std::vector< geode::index_t > offsets;
        offsets.reserve( model.nb_blocks() + 1 );
        geode::index_t offset{ 0 };
        for( const auto& block : model.blocks() )
        {
            const auto nb_polyhedra = block.mesh().nb_polyhedra();
            geode::OpenGeodeGeosciencesIOModelException::check_exception(
                nb_polyhedra == 0
                    || polyhedra_mapping.in2out( { block.id(), 0 } ) == offset,
                nullptr, geode::OpenGeodeException::TYPE::data, context,
                " Polyhedra of Block ", block.id().string(),
                " are not contiguous in the model solid." );
            offsets.push_back( offset );
            offset += nb_polyhedra;
        }
        geode::OpenGeodeGeosciencesIOModelException::check_exception(
            offset == solid.nb_polyhedra(), nullptr,
            geode::OpenGeodeException::TYPE::data, context,
            " Wrong number of polyhedra in the model solid." );
        offsets.push_back( offset );
        return offsets;
    }
} // namespace

namespace geode::internal
{
    SimulationSolid::SimulationSolid( const BRep& model,
        const SimulationRegions& regions,
        std::string_view context )
    {
        std::tie( mesh_, mappings_ ) = convert_brep_into_solid( model );
        block_polyhedra_offsets_ = compute_block_polyhedra_offsets(
            model, mappings_, *mesh_, context );
        AttributeValues< index_t > region_attribute_values;
        region_attribute_values.default_value = NO_ID;
        region_attribute_values.no_value = NO_ID;
        AttributeProperties region_attribute_properties;
        region_attribute_properties.assignable = false;
        region_attribute_properties.interpolable = false;
        region_attribute_properties.transferable = true;
        region_attribute_id_ =
            mesh_->polyhedron_attribute_manager()
                .create_attribute< VariableAttribute, index_t >(
                    REGION_ID_ATTRIBUTE_NAME, region_attribute_values,
                    region_attribute_properties );
        region_attribute_ = mesh_->polyhedron_attribute_manager()
                                .find_attribute< VariableAttribute, index_t >(
                                    region_attribute_id_ );
        index_t block_index{ 0 };
        for( const auto& block : model.blocks() )
        {
            const auto region = regions.block_regions.find( block.id() );
            OpenGeodeGeosciencesIOModelException::check_exception(
                region != regions.block_regions.end(), nullptr,
                OpenGeodeException::TYPE::internal, context, " Block ",
                block.id().string(), " belongs to no region." );
            for( const auto polyhedron_id :
                Range{ block_polyhedra_offsets_[block_index],
                    block_polyhedra_offsets_[block_index + 1] } )
            {
                region_attribute_->set_value( polyhedron_id, region->second );
            }
            block_index++;
        }
    }

    SimulationSolid::~SimulationSolid() = default;

    SolidMesh3D& SimulationSolid::mesh()
    {
        return *mesh_;
    }

    const SolidMesh3D& SimulationSolid::mesh() const
    {
        return *mesh_;
    }

    const ModelToMeshMappings& SimulationSolid::mappings() const
    {
        return mappings_;
    }

    absl::Span< const index_t > SimulationSolid::block_polyhedra_offsets() const
    {
        return block_polyhedra_offsets_;
    }

    const uuid& SimulationSolid::region_attribute_id() const
    {
        return region_attribute_id_;
    }

    index_t SimulationSolid::region( index_t polyhedron ) const
    {
        return region_attribute_->value( polyhedron );
    }

    std::optional< std::string > SimulationSolid::transfer_block_property(
        const SimulationNaming& naming,
        const BRep& model,
        const PhysicalProperties::Info& info,
        std::string_view simulator_field_name )
    {
        return transfer_block_attribute( naming, model, *mesh_,
            block_polyhedra_offsets_, info, simulator_field_name );
    }
} // namespace geode::internal
