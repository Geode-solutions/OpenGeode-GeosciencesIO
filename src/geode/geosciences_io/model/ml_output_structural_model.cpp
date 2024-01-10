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

#include <geode/geosciences_io/model/private/ml_output_structural_model.h>

#include <string>
#include <vector>

#include <geode/mesh/core/surface_mesh.h>

#include <geode/geosciences/explicit/representation/builder/structural_model_builder.h>
#include <geode/geosciences/explicit/representation/core/structural_model.h>
#include <geode/geosciences_io/model/private/gocad_common.h>
#include <geode/geosciences_io/model/private/ml_output_impl.h>

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
                    file() << "TFACE " << component_id() << SPACE
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
                    file() << "TFACE " << component_id() << SPACE
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
                file() << "TSURF " << fault.name() << EOL;
            }
            for( const auto& horizon : model_.horizons() )
            {
                file() << "TSURF " << horizon.name() << EOL;
            }
        }

        std::vector< geode::uuid > unclassified_tsurfs() const
        {
            std::vector< geode::uuid > result;
            for( const auto& surface : model_.surfaces() )
            {
                if( model_.nb_collections( surface.id() ) == 0 )
                {
                    result.push_back( surface.id() );
                }
            }
            return result;
        }

        void write_geological_regions() override
        {
            for( const auto& stratigraphic_unit : model_.stratigraphic_units() )
            {
                file() << "LAYER " << stratigraphic_unit.name() << EOL << SPACE
                       << SPACE;
                geode::index_t counter{ 0 };
                for( const auto& item :
                    model_.stratigraphic_unit_items( stratigraphic_unit ) )
                {
                    file() << components().at( item.id() ) << SPACE << SPACE;
                    counter++;
                    if( counter % 5 == 0 )
                    {
                        file() << EOL << SPACE << SPACE;
                    }
                }
                file() << 0 << EOL;
            }

            for( const auto& fault_block : model_.fault_blocks() )
            {
                file() << "FAULT_BLOCK " << fault_block.name() << EOL << SPACE
                       << SPACE;
                geode::index_t counter{ 0 };
                for( const auto& item :
                    model_.fault_block_items( fault_block ) )
                {
                    file() << components().at( item.id() ) << SPACE << SPACE;
                    counter++;
                    if( counter % 5 == 0 )
                    {
                        file() << EOL << SPACE << SPACE;
                    }
                }
                file() << 0 << EOL;
            }
        }

        void write_geological_model_surfaces() override
        {
            for( const auto& fault : model_.faults() )
            {
                file() << "GOCAD TSurf 1" << EOL;
                geode::detail::HeaderData header;
                header.name = geode::to_string( fault.name() );
                geode::detail::write_header( file(), header );
                geode::detail::write_CRS( file(), {} );
                file() << "GEOLOGICAL_FEATURE " << fault.name() << EOL;
                file() << "GEOLOGICAL_TYPE " << fault_map_.at( fault.type() )
                       << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.fault_items( fault ) )
                {
                    file() << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;

                find_corners_and_line_starts(
                    model_.fault_items( fault ), line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                file() << "END" << EOL;
            }
            for( const auto& horizon : model_.horizons() )
            {
                file() << "GOCAD TSurf 1" << EOL;
                geode::detail::HeaderData header;
                header.name = geode::to_string( horizon.name() );
                geode::detail::write_header( file(), header );
                geode::detail::write_CRS( file(), {} );
                file() << "GEOLOGICAL_FEATURE " << horizon.name() << EOL;
                file() << "GEOLOGICAL_TYPE "
                       << horizon_map_.at( horizon.type() ) << EOL;
                geode::index_t current_offset{ OFFSET_START };
                for( const auto& item : model_.horizon_items( horizon ) )
                {
                    file() << "TFACE" << EOL;
                    current_offset = write_surface( item, current_offset );
                }
                std::vector< std::array< geode::index_t, 2 > > line_starts;

                find_corners_and_line_starts(
                    model_.horizon_items( horizon ), line_starts );
                write_corners( line_starts );
                write_line_starts( current_offset, line_starts );
                file() << "END" << EOL;
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
        std::vector< std::string > MLOutputStructuralModel::write(
            const StructuralModel& structural_model ) const
        {
            OPENGEODE_EXCEPTION( check_brep_polygons( structural_model ),
                "[MLOutput::write] Can not export into .ml a "
                "StructuralModel with non triangular surface polygons." );
            if( structural_model.nb_model_boundaries() > 0 )
            {
                MLOutputImplSM impl{ filename(), structural_model };
                impl.write_file();
            }
            else
            {
                const auto new_structural_model =
                    clone_with_model_boundaries( structural_model );
                MLOutputImplSM impl{ filename(), new_structural_model };
                impl.write_file();
            }
            return { to_string( filename() ) };
        }

        bool MLOutputStructuralModel::is_saveable(
            const StructuralModel& structural_model ) const
        {
            return check_brep_polygons( structural_model );
        }
    } // namespace detail
} // namespace geode
