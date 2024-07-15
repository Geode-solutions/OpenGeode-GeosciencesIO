/*
 * Copyright (c) 2019 - 2024 Geode-solutions
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

#include <geode/geosciences_io/model/helpers/brep_geos_export.h>

#include <geode/geosciences_io/model/private/geos_export.h>
#include <geode/model/mixin/core/block.h>
#include <geode/model/representation/core/brep.h>

#include <geode/basic/pimpl_impl.h>

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
        absl::flat_hash_map< uuid, index_t > create_region_attribute_map(
            const BRep& model ) const final
        {
            auto region_id = 0;
            absl::flat_hash_map< uuid, index_t > region_map_id;
            for( const auto& block : model.blocks() )
            {
                region_map_id.emplace( block.id(), region_id++ );
            }
            return region_map_id;
        }
    };

    BRepGeosExporter::BRepGeosExporter(
        const BRep& brep, std::string_view files_directory )
        : impl_( brep, files_directory )
    {
    }

    BRepGeosExporter::~BRepGeosExporter() = default;

    void BRepGeosExporter::add_well_perforations(
        const PointSet3D& well_perforation )
    {
        impl_->add_well_perforations( well_perforation );
    }
    void BRepGeosExporter::add_cell_property_1d( std::string_view name )
    {
        impl_->add_cell_property1d( name );
    }
    void BRepGeosExporter::add_cell_property_2d( std::string_view name )
    {
        impl_->add_cell_property2d( name );
    }
    void BRepGeosExporter::add_cell_property_3d( std::string_view name )
    {
        impl_->add_cell_property3d( name );
    }

    void BRepGeosExporter::run()
    {
        impl_->prepare_export();
        impl_->write_files();
    }

} // namespace geode
