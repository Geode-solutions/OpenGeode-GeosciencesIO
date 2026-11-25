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

#include <geode/geosciences_io/model/helpers/brep_ogs_export.hpp>

#include <geode/geosciences_io/model/internal/ogs_export.hpp>
#include <geode/geosciences_io/model/internal/simulation_regions.hpp>
#include <geode/model/representation/core/brep.hpp>

#include <geode/basic/pimpl_impl.hpp>

namespace geode
{
    class BRepOgsExporter::Impl : public internal::OgsExporterImpl< BRep >
    {
        OPENGEODE_DISABLE_COPY_AND_MOVE( Impl );

    public:
        Impl() = delete;
        Impl( const BRep& brep, std::string_view files_directory )
            : internal::OgsExporterImpl< BRep >(
                files_directory, brep, "[OgsExporter]" )
        {
        }
        ~Impl() override = default;

    protected:
        [[nodiscard]] internal::SimulationRegions create_region_definitions(
            const BRep& model ) const final
        {
            auto regions = internal::brep_block_regions( model );
            internal::add_block_collection_regions( model, regions );
            return regions;
        }
    };

    BRepOgsExporter::BRepOgsExporter(
        const BRep& model, std::string_view files_directory )
        : impl_( model, files_directory )
    {
    }

    BRepOgsExporter::~BRepOgsExporter() = default;

    void BRepOgsExporter::add_well_perforations(
        const PointSet3D& well_perforations, std::string_view name )
    {
        impl_->add_well_perforations( well_perforations, name );
    }

    void BRepOgsExporter::set_physical_problem( const PhysicalProblem& problem )
    {
        impl_->set_physical_problem( problem );
    }

    void BRepOgsExporter::run()
    {
        impl_->prepare_export();
        impl_->write_files();
    }
} // namespace geode
