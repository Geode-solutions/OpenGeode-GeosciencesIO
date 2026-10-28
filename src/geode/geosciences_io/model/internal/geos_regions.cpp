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

#include <geode/geosciences_io/model/internal/geos_regions.hpp>

#include <absl/container/linked_hash_map.h>
#include <absl/strings/str_cat.h>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.hpp>
#include <geode/geosciences/explicit/representation/core/structural_model.hpp>

namespace
{
    void add_region( geode::internal::GeosDeckContext& context,
        std::string_view component_name,
        const geode::uuid& component_id,
        std::vector< geode::uuid > blocks )
    {
        const auto region_id =
            static_cast< geode::index_t >( context.solid_regions.size() );
        for( const auto& block : blocks )
        {
            context.block_regions.emplace( block, region_id );
        }
        context.solid_regions.push_back(
            { region_id, absl::StrCat( "region_", region_id ),
                geode::to_string( component_name ), component_id,
                std::move( blocks ) } );
    }

    void add_remaining_blocks(
        const geode::BRep& model, geode::internal::GeosDeckContext& context )
    {
        for( const auto& block : model.blocks() )
        {
            if( context.block_regions.contains( block.id() ) )
            {
                continue;
            }
            add_region( context, block.name().value_or( block.id().string() ),
                block.id(), { block.id() } );
        }
    }
} // namespace

namespace geode::internal
{
    GeosDeckContext geos_region_context( const BRep& model )
    {
        GeosDeckContext context;
        add_remaining_blocks( model, context );
        return context;
    }

    GeosDeckContext geos_region_context( const StructuralModel& model )
    {
        GeosDeckContext context;
        for( const auto& unit : model.stratigraphic_units() )
        {
            std::vector< uuid > blocks;
            for( const auto& block : model.stratigraphic_unit_items( unit ) )
            {
                blocks.push_back( block.id() );
            }
            if( !blocks.empty() )
            {
                add_region( context, unit.name().value_or( unit.id().string() ),
                    unit.id(), std::move( blocks ) );
            }
        }
        add_remaining_blocks( model, context );
        return context;
    }

    GeosDeckContext geos_inferred_region_context( const BRep& model )
    {
        GeosDeckContext context;
        context.regions_are_inferred = true;
        absl::linked_hash_map< uuid, std::vector< uuid > > unit_blocks;
        for( const auto& block : model.blocks() )
        {
            for( const auto& collection : model.collections( block.id() ) )
            {
                if( collection.type
                    == StratigraphicUnit3D::component_type_static() )
                {
                    unit_blocks[collection.id].push_back( block.id() );
                    break;
                }
            }
        }
        for( auto& [unit, blocks] : unit_blocks )
        {
            add_region( context,
                absl::StrCat( "StratigraphicUnit ", unit.string() ), unit,
                std::move( blocks ) );
        }
        add_remaining_blocks( model, context );
        return context;
    }
} // namespace geode::internal
