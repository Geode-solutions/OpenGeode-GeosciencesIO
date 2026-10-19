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

#include <functional>
#include <string>
#include <vector>

#include <absl/container/flat_hash_set.h>

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
            assign_surfaces_to_tsurfs();
        }

    private:
        using TSurfSurfaces =
            std::vector< std::reference_wrapper< const geode::Surface3D > >;

        /*!
         * In a .ml file, each TFACE belongs to a single TSURF. A Surface
         * shared by several collections (e.g. a topographic Horizon that is
         * also a ModelBoundary) is assigned to the first one, in the export
         * order: ModelBoundaries, Faults, then Horizons. A collection
         * without any Surface of its own is not exported.
         */
        void assign_surfaces_to_tsurfs()
        {
            absl::flat_hash_set< geode::uuid > assigned;
            for( const auto& boundary : model_.model_boundaries() )
            {
                for( const auto& item :
                    model_.model_boundary_items( boundary ) )
                {
                    assigned.insert( item.id() );
                }
            }
            for( const auto& fault : model_.faults() )
            {
                assign_surfaces( fault, model_.fault_items( fault ), assigned );
            }
            for( const auto& horizon : model_.horizons() )
            {
                assign_surfaces(
                    horizon, model_.horizon_items( horizon ), assigned );
            }
        }

        template < typename Collection, typename ItemRange >
        void assign_surfaces( const Collection& collection,
            const ItemRange& items,
            absl::flat_hash_set< geode::uuid >& assigned )
        {
            auto& surfaces = tsurf_surfaces_[collection.id()];
            for( const auto& item : items )
            {
                if( !assigned.insert( item.id() ).second )
                {
                    geode::Logger::warning( "[MLOutput] A Surface from ",
                        component_name( collection ),
                        " belongs to several collections. It has been "
                        "exported only once" );
                    continue;
                }
                surfaces.emplace_back( item );
            }
            if( surfaces.empty() )
            {
                geode::Logger::warning( "[MLOutput] ",
                    component_name( collection ),
                    " has no Surface of its own. It has not been exported" );
            }
        }

        const TSurfSurfaces& tsurf_surfaces( const geode::uuid& id ) const
        {
            return tsurf_surfaces_.at( id );
        }

        void write_tfaces( const TSurfSurfaces& surfaces,
            std::string_view type,
            std::string_view name )
        {
            for( const auto& surface : surfaces )
            {
                file() << "TFACE " << component_id() << SPACE << type << SPACE
                       << name << EOL;
                write_key_triangle( surface.get() );
                components().emplace( surface.get().id(), component_id()++ );
            }
        }

        void write_geological_tfaces() override
        {
            for( const auto& fault : model_.faults() )
            {
                write_tfaces( tsurf_surfaces( fault.id() ),
                    fault_map_.at( fault.type() ), component_name( fault ) );
            }
            for( const auto& horizon : model_.horizons() )
            {
                write_tfaces( tsurf_surfaces( horizon.id() ),
                    horizon_type( horizon ), component_name( horizon ) );
            }
        }

        void write_geological_tsurfs() override
        {
            for( const auto& fault : model_.faults() )
            {
                if( !tsurf_surfaces( fault.id() ).empty() )
                {
                    file() << "TSURF " << component_name( fault ) << EOL;
                }
            }
            for( const auto& horizon : model_.horizons() )
            {
                if( !tsurf_surfaces( horizon.id() ).empty() )
                {
                    file() << "TSURF " << component_name( horizon ) << EOL;
                }
            }
        }

        std::vector< geode::uuid > unclassified_tsurfs() const override
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

        void write_geological_surface( const TSurfSurfaces& surfaces,
            std::string_view type,
            std::string_view name )
        {
            if( surfaces.empty() )
            {
                return;
            }
            file() << "GOCAD TSurf 1" << EOL;
            geode::internal::HeaderData header;
            header.name = geode::to_string( name );
            geode::internal::write_header( file(), header );
            geode::internal::write_CRS( file(), {} );
            file() << "GEOLOGICAL_FEATURE " << name << EOL;
            file() << "GEOLOGICAL_TYPE " << type << EOL;
            geode::index_t current_offset{ OFFSET_START };
            for( const auto& surface : surfaces )
            {
                file() << "TFACE" << EOL;
                current_offset = write_surface( surface.get(), current_offset );
            }
            std::vector< std::array< geode::index_t, 2 > > line_starts;
            find_corners_and_line_starts( surfaces, line_starts );
            write_corners( line_starts );
            write_line_starts( current_offset, line_starts );
            file() << "END" << EOL;
        }

        void write_geological_model_surfaces() override
        {
            for( const auto& fault : model_.faults() )
            {
                write_geological_surface( tsurf_surfaces( fault.id() ),
                    fault_map_.at( fault.type() ), component_name( fault ) );
            }
            for( const auto& horizon : model_.horizons() )
            {
                write_geological_surface( tsurf_surfaces( horizon.id() ),
                    horizon_type( horizon ), component_name( horizon ) );
            }
        }

        std::string_view horizon_type( const geode::Horizon3D& horizon ) const
        {
            const auto it = horizon_map_.find( horizon.contact_type() );
            if( it == horizon_map_.end() )
            {
                return "top";
            }
            return it->second;
        }

    private:
        const geode::StructuralModel& model_;
        absl::flat_hash_map< geode::uuid, TSurfSurfaces > tsurf_surfaces_;
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
