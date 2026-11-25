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

#include <geode/geosciences_io/model/internal/simulation_regions.hpp>

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <absl/strings/str_cat.h>
#include <absl/types/span.h>

#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/uuid.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/geosciences/explicit/mixin/core/stratigraphic_units.hpp>
#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences/simulation/core/initial_condition.hpp>
#include <geode/geosciences/simulation/core/material.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/helpers/internal/physical_problem_blocks.hpp>
#include <geode/geosciences/simulation/helpers/internal/simulation_checker_helpers.hpp>
#include <geode/geosciences/simulation/helpers/simulation_report.hpp>

namespace
{
    // Assigns items (materials or initial conditions) to the regions:
    // items with regions first, then the default item (no region) to the
    // remaining regions
    class RegionAssigner
    {
    public:
        RegionAssigner( const geode::BRep& model,
            const geode::internal::SimulationRegions& regions,
            std::string_view simulator_name )
            : model_( model ),
              regions_( regions ),
              simulator_name_( simulator_name )
        {
            assignment_.region_items.resize( regions_.solid_regions.size() );
        }

        void assign( geode::index_t item,
            const std::string& subject,
            absl::Span< const geode::uuid > components )
        {
            for( const auto& component : components )
            {
                const auto regions = regions_.regions_by_uuid.find( component );
                if( regions == regions_.regions_by_uuid.end() )
                {
                    report_unknown_region( subject, component );
                    continue;
                }
                for( const auto region : regions->second )
                {
                    assignment_.region_items[region] = item;
                }
            }
        }

        geode::internal::RegionAssignment finish(
            std::optional< geode::index_t > default_item )
        {
            for( auto& item : assignment_.region_items )
            {
                if( !item )
                {
                    item = default_item;
                }
            }
            return std::move( assignment_ );
        }

    private:
        void add_issue( geode::SimulationIssue issue )
        {
            assignment_.issues.push_back( std::move( issue ) );
        }

        // Regions which are not Blocks nor collections of Blocks are
        // reported by SimulationChecker::check_physical_problem_regions
        void report_unknown_region(
            const std::string& subject, const geode::uuid& component )
        {
            if( const auto merged = regions_.merged_blocks.find( component );
                merged != regions_.merged_blocks.end() )
            {
                report_merged_block( subject, component, merged->second );
                return;
            }
            const auto blocks = geode::internal::physical_problem_region_blocks(
                model_, component );
            if( !blocks )
            {
                return;
            }
            geode::SimulationIssue issue{
                geode::SIMULATION_ISSUE_SEVERITY::error,
                geode::SIMULATION_CHECK_STAGE::regions, subject,
                absl::StrCat( "region ", component.string(),
                    " contains Blocks merged into other ", simulator_name_,
                    " regions" )
            };
            issue.fix = "assign the StratigraphicUnits, or Blocks of no "
                        "StratigraphicUnit";
            for( const auto& block : blocks.value() )
            {
                issue.add_location( model_.block( block ) );
            }
            add_issue( std::move( issue ) );
        }

        void report_merged_block( const std::string& subject,
            const geode::uuid& block,
            geode::index_t region )
        {
            const auto& component_name =
                regions_.solid_regions[region].component_name;
            geode::SimulationIssue issue{
                geode::SIMULATION_ISSUE_SEVERITY::error,
                geode::SIMULATION_CHECK_STAGE::regions, subject,
                absl::StrCat( "the Block is part of the ", simulator_name_,
                    " region of ", component_name,
                    ", which merges several Blocks" )
            };
            issue.fix = absl::StrCat(
                "assign ", component_name, " instead of its Blocks" );
            issue.add_location( model_.block( block ) );
            add_issue( std::move( issue ) );
        }

    private:
        const geode::BRep& model_;
        const geode::internal::SimulationRegions& regions_;
        std::string simulator_name_;
        geode::internal::RegionAssignment assignment_;
    };

    geode::internal::RegionAssignment assign_regions( RegionAssigner assigner,
        absl::Span< const absl::Span< const geode::uuid > > item_regions,
        absl::Span< const std::string > item_subjects )
    {
        std::optional< geode::index_t > default_item;
        for( const auto item : geode::Indices{ item_regions } )
        {
            if( item_regions[item].empty() )
            {
                default_item = item;
                continue;
            }
            assigner.assign( item, item_subjects[item], item_regions[item] );
        }
        return assigner.finish( default_item );
    }

    geode::index_t add_region( geode::internal::SimulationRegions& definitions,
        std::string_view component_name,
        const geode::uuid& component_id )
    {
        const auto region_id =
            static_cast< geode::index_t >( definitions.solid_regions.size() );
        definitions.solid_regions.push_back(
            { region_id, absl::StrCat( "region_", region_id ),
                geode::to_string( component_name ), component_id } );
        definitions.regions_by_uuid[component_id].push_back( region_id );
        return region_id;
    }
} // namespace

namespace geode::internal
{
    SimulationRegions brep_block_regions( const BRep& model )
    {
        SimulationRegions definitions;
        for( const auto& block : model.blocks() )
        {
            const auto region_id =
                static_cast< index_t >( definitions.solid_regions.size() );
            definitions.solid_regions.push_back( { region_id,
                absl::StrCat( "region_", region_id ),
                block.name().value_or( block.id().string() ), block.id() } );
            definitions.block_regions.emplace( block.id(), region_id );
            definitions.regions_by_uuid[block.id()].push_back( region_id );
        }
        return definitions;
    }

    SimulationRegions structural_model_unit_regions(
        const StructuralModel& model, std::string_view context )
    {
        SimulationRegions definitions;
        for( const auto& unit : model.stratigraphic_units() )
        {
            const auto region_id = add_region( definitions,
                unit.name().value_or( unit.id().string() ), unit.id() );
            for( const auto& block : model.stratigraphic_unit_items( unit ) )
            {
                if( !definitions.block_regions.emplace( block.id(), region_id )
                         .second )
                {
                    Logger::warning( context, " Block ", block.id().string(),
                        " belongs to several StratigraphicUnits, only the "
                        "first one is used as its region." );
                    continue;
                }
                definitions.merged_blocks.emplace( block.id(), region_id );
            }
        }
        for( const auto& block : model.blocks() )
        {
            if( definitions.block_regions.contains( block.id() ) )
            {
                continue;
            }
            const auto region_id = add_region( definitions,
                block.name().value_or( block.id().string() ), block.id() );
            definitions.block_regions.emplace( block.id(), region_id );
        }
        return definitions;
    }

    SimulationRegions structural_model_block_regions(
        const StructuralModel& model )
    {
        auto definitions = brep_block_regions( model );
        for( const auto& unit : model.stratigraphic_units() )
        {
            auto& unit_regions = definitions.regions_by_uuid[unit.id()];
            for( const auto& block : model.stratigraphic_unit_items( unit ) )
            {
                unit_regions.push_back(
                    definitions.block_regions.at( block.id() ) );
            }
        }
        return definitions;
    }

    void add_block_collection_regions(
        const BRep& model, SimulationRegions& regions )
    {
        for( const auto& collection : model.block_collections() )
        {
            std::vector< index_t > collection_regions;
            bool has_merged_block{ false };
            for( const auto& block :
                model.block_collection_items( collection ) )
            {
                if( regions.merged_blocks.contains( block.id() ) )
                {
                    has_merged_block = true;
                    break;
                }
                collection_regions.push_back(
                    regions.block_regions.at( block.id() ) );
            }
            if( !has_merged_block )
            {
                regions.regions_by_uuid.emplace(
                    collection.id(), std::move( collection_regions ) );
            }
        }
    }

    RegionAssignment assign_region_materials( const PhysicalProblem& problem,
        const BRep& model,
        const SimulationRegions& regions,
        std::string_view simulator_name )
    {
        std::vector< absl::Span< const uuid > > material_regions;
        std::vector< std::string > subjects;
        for( const auto& material : problem.materials() )
        {
            material_regions.emplace_back( material.regions );
            subjects.push_back( material_subject( material.name ) );
        }
        return assign_regions( RegionAssigner{ model, regions, simulator_name },
            material_regions, subjects );
    }

    RegionAssignment assign_region_initial_conditions(
        const PhysicalProblem& problem,
        const BRep& model,
        const SimulationRegions& regions,
        std::string_view simulator_name )
    {
        std::vector< absl::Span< const uuid > > condition_regions;
        std::vector< std::string > subjects;
        for( const auto condition_id :
            Range{ problem.nb_initial_conditions() } )
        {
            condition_regions.emplace_back(
                problem.initial_condition( condition_id ).regions() );
            subjects.push_back( initial_condition_subject( condition_id ) );
        }
        return assign_regions( RegionAssigner{ model, regions, simulator_name },
            condition_regions, subjects );
    }
} // namespace geode::internal
