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

#include <geode/geosciences_io/model/internal/ml_output_structural_model.hpp>

#include <string>
#include <vector>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/geosciences/explicit/representation/builder/structural_model_builder.hpp>
#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences_io/model/internal/gocad_common.hpp>
#include <geode/geosciences_io/model/internal/ml_output_impl.hpp>

namespace
{
    class MLOutputImplSM
        : public geode::internal::MLOutputImpl< geode::StructuralModel >
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };

        MLOutputImplSM(
            std::string_view filename, const geode::StructuralModel& model )
            : geode::internal::MLOutputImpl< geode::StructuralModel >(
                  filename, model ),
              model_( model )
        {
        }

    private:
        void assign_geological_surfaces() override
        {
            for( const auto& fault : model_.faults() )
            {
                assign_surfaces( fault, model_.fault_items( fault ) );
            }
            for( const auto& horizon : model_.horizons() )
            {
                assign_surfaces( horizon, model_.horizon_items( horizon ) );
            }
        }

        void write_geological_tfaces() override
        {
            for( const auto& fault : model_.faults() )
            {
                write_collection_tfaces( owned_surfaces( fault.id() ),
                    fault_type( fault ), tsurf_name( fault ) );
            }
            for( const auto& horizon : model_.horizons() )
            {
                write_collection_tfaces( owned_surfaces( horizon.id() ),
                    horizon_type( horizon ), tsurf_name( horizon ) );
            }
        }

        void write_geological_tsurfs() override
        {
            for( const auto& fault : model_.faults() )
            {
                if( !owned_surfaces( fault.id() ).empty() )
                {
                    file() << "TSURF " << tsurf_name( fault ) << EOL;
                }
            }
            for( const auto& horizon : model_.horizons() )
            {
                if( !owned_surfaces( horizon.id() ).empty() )
                {
                    file() << "TSURF " << tsurf_name( horizon ) << EOL;
                }
            }
        }

        void write_geological_regions() override
        {
            for( const auto& stratigraphic_unit : model_.stratigraphic_units() )
            {
                file() << "LAYER " << component_name( stratigraphic_unit )
                       << EOL << SPACE << SPACE;
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
                file() << "FAULT_BLOCK " << component_name( fault_block ) << EOL
                       << SPACE << SPACE;
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
                write_collection_surface( owned_surfaces( fault.id() ),
                    fault_type( fault ), tsurf_name( fault ) );
            }
            for( const auto& horizon : model_.horizons() )
            {
                write_collection_surface( owned_surfaces( horizon.id() ),
                    horizon_type( horizon ), tsurf_name( horizon ) );
            }
        }

        std::string_view fault_type( const geode::Fault3D& fault ) const
        {
            const auto it = fault_map_.find( fault.type() );
            if( it == fault_map_.end() )
            {
                geode::Logger::warning( "[MLOutput] Fault ",
                    component_name( fault ),
                    " has an unknown type. It has been exported as fault" );
                return "fault";
            }
            return it->second;
        }

        std::string_view horizon_type( const geode::Horizon3D& horizon ) const
        {
            const auto it = horizon_map_.find( horizon.contact_type() );
            if( it == horizon_map_.end() )
            {
                geode::Logger::warning( "[MLOutput] Horizon ",
                    component_name( horizon ),
                    " has an unknown contact type. It has been exported as "
                    "top" );
                return "top";
            }
            return it->second;
        }

    private:
        const geode::StructuralModel& model_;
        const absl::flat_hash_map< geode::Fault3D::FAULT_TYPE, std::string >
            fault_map_ = { { geode::Fault3D::FAULT_TYPE::no_type, "fault" },
                { geode::Fault3D::FAULT_TYPE::normal, "normal_fault" },
                { geode::Fault3D::FAULT_TYPE::reverse, "reverse_fault" },
                { geode::Fault3D::FAULT_TYPE::strike_slip, "fault" },
                { geode::Fault3D::FAULT_TYPE::listric, "fault" },
                { geode::Fault3D::FAULT_TYPE::decollement, "fault" } };
        // Erosion and baselap contacts have no dedicated .ml type: they are
        // exported as unconformities
        const absl::flat_hash_map< geode::Horizon3D::CONTACT_TYPE, std::string >
            horizon_map_ = { { geode::Horizon3D::CONTACT_TYPE::conformal,
                                 "top" },
                { geode::Horizon3D::CONTACT_TYPE::topography, "topographic" },
                { geode::Horizon3D::CONTACT_TYPE::intrusion, "intrusive" },
                { geode::Horizon3D::CONTACT_TYPE::discontinuity,
                    "unconformity" },
                { geode::Horizon3D::CONTACT_TYPE::erosion, "unconformity" },
                { geode::Horizon3D::CONTACT_TYPE::baselap, "unconformity" } };
    };
} // namespace

namespace geode::internal
{
    std::vector< std::string > MLOutputStructuralModel::write(
        const StructuralModel& structural_model ) const
    {
        OpenGeodeGeosciencesIOModelException::check_exception(
            check_brep_polygons( structural_model ), nullptr,
            OpenGeodeException::TYPE::data,
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
} // namespace geode::internal
