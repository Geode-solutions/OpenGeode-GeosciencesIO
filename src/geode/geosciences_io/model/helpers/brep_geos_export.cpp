/*
 * Copyright (c) 2019 - 2023 Geode-solutions
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

#include <geode/geosciences_io/model/helpers/brep_geos_export.h>

#include <geode/geosciences_io/model/helpers/geos_export.h>

#include <geode/basic/pimpl_impl.h>

namespace geode
{
    class BRepGeosExporter::Impl : public geode::GeosExporterImpl< geode::BRep >
    {
    public:
        Impl( const BRep& brep, absl::string_view files_directory )
            : geode::GeosExporterImpl< geode::BRep >( files_directory, brep )
        {
        }
        virtual ~Impl() = default;

    protected:
        absl::flat_hash_map< geode::uuid, geode::index_t >
            create_region_attribute_map( const geode::BRep& model ) const final
        {
            auto region_id{ 0 };
            absl::flat_hash_map< geode::uuid, geode::index_t > region_map_id;
            for( const auto& block : model.blocks() )
            {
                region_map_id.emplace( block.id(), region_id++ );
            }
            return region_map_id;
        }
    };

} // namespace geode

namespace geode
{
    BRepGeosExporter::BRepGeosExporter(
        const BRep& brep, absl::string_view files_directory )
        : impl_( brep, files_directory )
    {
    }

    BRepGeosExporter::~BRepGeosExporter() = default;

    void BRepGeosExporter::add_well_perforations(
        const PointSet3D& well_perforation )
    {
        impl_->add_well_perforations( well_perforation );
    }
    void BRepGeosExporter::add_cell_property( absl::string_view name )
    {
        impl_->add_cell_property( name );
    }

    void BRepGeosExporter::run()
    {
        impl_->prepare_export();
        impl_->write_files();
    }

} // namespace geode
