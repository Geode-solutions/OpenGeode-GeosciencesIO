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
#include <absl/container/flat_hash_map.h>
#include <absl/strings/string_view.h>

#include <geode/geosciences_io/model/helpers/structural_model_geos_export.h>

#include <geode/geosciences_io/model/helpers/geos_export.h>

// #include <geode/geosciences/explicit/mixin/core/stratigraphic_unit.h>
#include <geode/geosciences/explicit/mixin/core/stratigraphic_units.h>

#include <geode/geosciences/explicit/representation/core/structural_model.h>

#include <geode/basic/pimpl_impl.h>

namespace geode
{
    class StructuralModelGeosExporter::Impl
        : public geode::GeosExporterImpl< geode::StructuralModel >
    {
    public:
        Impl( const StructuralModel& model, absl::string_view files_directory )
            : geode::GeosExporterImpl< geode::StructuralModel >(
                files_directory, model )
        {
        }
        virtual ~Impl() = default;

    protected:
        absl::flat_hash_map< geode::uuid, geode::index_t >
            create_region_attribute_map(
                const geode::StructuralModel& model ) const final
        {
            auto region_id = 0;
            absl::flat_hash_map< geode::uuid, geode::index_t > region_map_id;
            for( const auto& strat_unit : model.stratigraphic_units() )
            {
                for( const auto& strat_unit_item :
                    model.stratigraphic_unit_items( strat_unit ) )
                {
                    region_map_id.emplace( strat_unit_item.id(), region_id );
                }
                region_id += 1;
            }
            return region_map_id;
        }
    };

} // namespace geode

namespace geode
{
    StructuralModelGeosExporter::StructuralModelGeosExporter(
        const StructuralModel& model, absl::string_view files_directory )
        : impl_( model, files_directory )
    {
    }

    StructuralModelGeosExporter::~StructuralModelGeosExporter() = default;

    void StructuralModelGeosExporter::add_well_perforations(
        const PointSet3D& well_perforation )
    {
        impl_->add_well_perforations( well_perforation );
    }
    void StructuralModelGeosExporter::add_cell_property_1d(
        absl::string_view name )
    {
        impl_->add_cell_property1d( name );
    }
    void StructuralModelGeosExporter::add_cell_property_2d(
        absl::string_view name )
    {
        impl_->add_cell_property2d( name );
    }
    void StructuralModelGeosExporter::add_cell_property_3d(
        absl::string_view name )
    {
        impl_->add_cell_property3d( name );
    }

    void StructuralModelGeosExporter::run()
    {
        impl_->prepare_export();
        impl_->write_files();
    }

} // namespace geode
