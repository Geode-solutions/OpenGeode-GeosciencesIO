/*
 * Copyright (c) 2019 - 2021 Geode-solutions
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

#include <geode/geosciences/private/ml_output_structural_model.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <queue>

#include <geode/basic/filename.h>

#include <geode/geometry/basic_objects.h>
#include <geode/geometry/bounding_box.h>
#include <geode/geometry/signed_mensuration.h>

#include <geode/mesh/core/polygonal_surface.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/corner.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/mixin/core/vertex_identifier.h>

#include <geode/geosciences/private/gocad_common.h>
#include <geode/geosciences/private/ml_output_impl.h>
#include <geode/geosciences/representation/core/structural_model.h>

namespace
{
    class MLOutputImplSM
        : public geode::detail::MLOutputImpl< geode::StructuralModel >
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };

        MLOutputImplSM(
            absl::string_view filename, const geode::StructuralModel& model )
            : geode::detail::MLOutputImpl< geode::StructuralModel >(
                filename, model ),
              model_( model )
        {
        }

    private:
        void write_geological_tfaces() override
        {
            for( const auto& fault : model_.faults() )
            {
                for( const auto& item : model_.fault_items( fault ) )
                {
                    if( components().find( item.id() ) != components().end() )
                    {
                        geode::Logger::warn( "[MLOutput] A Surface from ",
                            fault.name(),
                            " belongs to several collections. It has been "
                            "exported only once" );
                        continue;
                    }
                    this->file() << "TFACE " << component_id() << SPACE
                                 << fault_map_.at( fault.type() ) << SPACE
                                 << fault.name() << EOL;
                    write_key_triangle( item );
                    components().emplace( item.id(), component_id()++ );
                }
            }
            for( const auto& horizon : model_.horizons() )
            {
                for( const auto& item : model_.horizon_items( horizon ) )
                {
                    if( components().find( item.id() ) != components().end() )
                    {
                        geode::Logger::warn( "[MLOutput] A Surface from ",
                            horizon.name(),
                            " belongs to several collections. It has been "
                            "exported only once" );
                        continue;
                    }
                    this->file() << "TFACE " << component_id() << SPACE
                                 << horizon_map_.at( horizon.type() ) << SPACE
                                 << horizon.name() << EOL;
                    write_key_triangle( item );
                    components().emplace( item.id(), component_id()++ );
                }
            }
        }

        void write_geological_tsurfs() override
        {
            for( const auto& fault : model_.faults() )
            {
                this->file() << "TSURF " << fault.name() << EOL;
            }
            for( const auto& horizon : model_.horizons() )
            {
                this->file() << "TSURF " << horizon.name() << EOL;
            }
        }

        void write_geological_regions() override
        {
            for( const auto& stratigraphic_unit : model_.stratigraphic_units() )
            {
                this->file() << "LAYER " << stratigraphic_unit.name() << EOL
                             << SPACE << SPACE;
                geode::index_t counter{ 0 };
                for( const auto& item :
                    model_.stratigraphic_unit_items( stratigraphic_unit ) )
                {
                    this->file()
                        << components().at( item.id() ) << SPACE << SPACE;
                    counter++;
                    if( counter % 5 == 0 )
                    {
                        this->file() << EOL << SPACE << SPACE;
                    }
                }
                this->file() << 0 << EOL;
            }

            for( const auto& fault_block : model_.fault_blocks() )
            {
                this->file() << "FAULT_BLOCK " << fault_block.name() << EOL
                             << SPACE << SPACE;
                geode::index_t counter{ 0 };
                for( const auto& item :
                    model_.fault_block_items( fault_block ) )
                {
                    this->file()
                        << components().at( item.id() ) << SPACE << SPACE;
                    counter++;
                    if( counter % 5 == 0 )
                    {
                        this->file() << EOL << SPACE << SPACE;
                    }
                }
                this->file() << 0 << EOL;
            }
        }

        void write_geological_model_surfaces()
        {
            for( const auto& fault : model_.faults() )
            {
                this->file() << "GOCAD TSurf 1" << EOL;
                geode::detail::HeaderData header;
                header.name = fault.name().data();
                geode::detail::write_header( this->file(), header );
                geode::detail::write_CRS( this->file(), {} );
                this->file() << "GEOLOGICAL_FEATURE " << fault.name() << EOL;
                this->file() << "GEOLOGICAL_TYPE "
                             << fault_map_.at( fault.type() ) << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.fault_items( fault ) )
                {
                    this->file() << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;

                find_corners_and_line_starts(
                    model_.fault_items( fault ), line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                this->file() << "END" << EOL;
            }
            for( const auto& horizon : model_.horizons() )
            {
                this->file() << "GOCAD TSurf 1" << EOL;
                geode::detail::HeaderData header;
                header.name = horizon.name().data();
                geode::detail::write_header( this->file(), header );
                geode::detail::write_CRS( this->file(), {} );
                this->file() << "GEOLOGICAL_FEATURE " << horizon.name() << EOL;
                this->file() << "GEOLOGICAL_TYPE "
                             << horizon_map_.at( horizon.type() ) << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.horizon_items( horizon ) )
                {
                    this->file() << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;

                find_corners_and_line_starts(
                    model_.horizon_items( horizon ), line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                this->file() << "END" << EOL;
            }
        }

    private:
        const geode::StructuralModel& model_;
        const absl::flat_hash_map< geode::Fault3D::FAULT_TYPE, std::string >
            fault_map_ = { { geode::Fault3D::FAULT_TYPE::NO_TYPE, "fault" },
                { geode::Fault3D::FAULT_TYPE::NORMAL, "normal_fault" },
                { geode::Fault3D::FAULT_TYPE::REVERSE, "reverse_fault" },
                { geode::Fault3D::FAULT_TYPE::STRIKE_SLIP, "fault" },
                { geode::Fault3D::FAULT_TYPE::LISTRIC, "fault" },
                { geode::Fault3D::FAULT_TYPE::DECOLLEMENT, "fault" } };
        const absl::flat_hash_map< geode::Horizon3D::HORIZON_TYPE, std::string >
            horizon_map_ = { { geode::Horizon3D::HORIZON_TYPE::NO_TYPE,
                                 "none" },
                { geode::Horizon3D::HORIZON_TYPE::CONFORMAL, "top" },
                { geode::Horizon3D::HORIZON_TYPE::TOPOGRAPHY, "topographic" },
                { geode::Horizon3D::HORIZON_TYPE::INTRUSION, "intrusive" },
                { geode::Horizon3D::HORIZON_TYPE::NON_CONFORMAL,
                    "unconformity" } };
    };
} // namespace

namespace geode
{
    namespace detail
    {
        void MLOutputStructuralModel::write() const
        {
            const auto only_triangles =
                check_structural_model_polygons( structural_model() );
            if( !only_triangles )
            {
                geode::Logger::info(
                    "[MLOutput::write] Can not export into .ml a "
                    "StructuralModel with non triangular surface polygons." );
                return;
            }
            MLOutputImplSM impl{ filename(), structural_model() };
            impl.determine_surface_to_regions_signs();
            impl.write_this->file();
        }
    } // namespace detail
} // namespace geode
