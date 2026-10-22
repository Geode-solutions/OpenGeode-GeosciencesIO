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
#include <string_view>

#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_cat.h>

#include <geode/geosciences_io/model/helpers/brep_geos_export.hpp>

#include <geode/geosciences_io/model/internal/geos_export.hpp>
#include <geode/model/mixin/core/block.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/basic/pimpl_impl.hpp>

namespace geode
{
    class BRepGeosExporter::Impl : public internal::GeosExporterImpl< BRep >
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( Impl );

    public:
        Impl() = delete;
        Impl( const BRep& brep, std::string_view files_directory )
            : internal::GeosExporterImpl< BRep >( files_directory, brep )
        {
        }
        virtual ~Impl() = default;

    protected:
        [[nodiscard]] RegionDefinitions create_region_definitions(
            const BRep& model ) const final
        {
            RegionDefinitions definitions;
            for( const auto& block : model.blocks() )
            {
                const auto region_id =
                    static_cast< index_t >( definitions.regions.size() );
                definitions.regions.push_back(
                    { region_id, absl::StrCat( "region_", region_id ),
                        block.name().value_or( block.id().string() ),
                        block.id() } );
                definitions.block_regions.emplace( block.id(), region_id );
                definitions.regions_by_uuid[block.id()].push_back( region_id );
            }
            return definitions;
        }
    };

    BRepGeosExporter::BRepGeosExporter(
        const BRep& brep, std::string_view files_directory )
        : impl_( brep, files_directory )
    {
    }

    BRepGeosExporter::~BRepGeosExporter() = default;

    void BRepGeosExporter::add_well_perforations(
        const PointSet3D& well_perforation, std::string_view name )
    {
        impl_->add_well_perforations( well_perforation, name );
    }

    void BRepGeosExporter::set_physical_problem(
        const PhysicalProblem& problem )
    {
        impl_->set_physical_problem( problem );
    }

    void BRepGeosExporter::run()
    {
        impl_->prepare_export();
        impl_->write_files();
    }

} // namespace geode
