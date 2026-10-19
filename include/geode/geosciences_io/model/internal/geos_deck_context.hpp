#pragma once

#include <array>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <absl/container/flat_hash_map.h>
#include <absl/types/span.h>

#include <geode/basic/uuid.hpp>

#include <geode/geosciences/simulation/core/fluid.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace geode::internal
{
    /// GEOS cell region made of the model cells whose region attribute is id
    struct GeosRegion
    {
        index_t id{ NO_ID };
        std::string name;
        std::string component_name;
        uuid component_id;
    };

    /// Model data needed to translate a PhysicalProblem into a GEOS deck
    struct GeosDeckContext
    {
        std::vector< GeosRegion > solid_regions;
        /// Block, StratigraphicUnit or collection uuid -> GEOS region ids
        absl::flat_hash_map< uuid, std::vector< index_t > > regions_by_uuid;
        /// VTU array name and mean value of the porosity PhysicalProperty
        std::optional< std::string > porosity_attribute;
        std::optional< double > mean_porosity;
        /// VTU array name and mean value of the permeability PhysicalProperty
        std::optional< std::string > permeability_attribute;
        std::optional< std::array< double, 3 > > mean_permeability;
        /// Phases of the boundary fluid fraction PhysicalProperties
        std::vector< PHASE > boundary_fraction_phases;
        bool has_boundary_pressure{ false };
        bool has_boundary_temperature{ false };
        /// Names of the well perforation sets given to the exporter
        std::vector< std::string > well_names;
    };

    /// Model array imported in GEOS, possibly duplicated under a new name
    struct GeosImportedField
    {
        std::string source_attribute;
        std::string vtu_name;
        std::string geos_name;
    };

    /// Item (material or initial condition) assigned to each solid region
    struct GeosRegionAssignment
    {
        std::vector< std::optional< index_t > > region_items;
        std::vector< std::string > issues;
    };

    /*!
     * Assigns items to regions: items with regions first, then the default
     * item (no region) to the remaining regions.
     */
    [[nodiscard]] GeosRegionAssignment assign_geos_regions(
        const GeosDeckContext& context,
        absl::Span< const absl::Span< const uuid > > item_regions,
        std::string_view item_kind );
} // namespace geode::internal
