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
#include <geode/basic/pimpl_impl.h>

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <absl/container/flat_hash_map.h>
#include <absl/strings/string_view.h>

#include <ghc/filesystem.hpp>

#include <geode/basic/attribute.h>
#include <geode/basic/attribute_manager.h>
#include <geode/basic/filename.h>

#include <geode/geosciences_io/model/common.h>

#include <geode/mesh/core/hybrid_solid.h>
#include <geode/mesh/core/polyhedral_solid.h>
#include <geode/mesh/core/solid_mesh.h>
#include <geode/mesh/core/surface_mesh.h>
#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/mesh/io/hybrid_solid_output.h>
#include <geode/mesh/io/polyhedral_solid_output.h>
#include <geode/mesh/io/tetrahedral_solid_output.h>

#include <geode/geosciences_io/model/helpers/brep_geos_export.h>
#include <geode/model/helpers/convert_to_mesh.h>
#include <geode/model/mixin/core/block.h>
#include <geode/model/representation/core/brep.h>
namespace
{
    /* absl::flat_hash_map< geode::uuid, geode::index_t >
         create_region_attribute_map( const geode::StructuralModel& model )
     {
         auto region_id{ 0 };
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
     }*/
    absl::flat_hash_map< geode::uuid, geode::index_t >
        create_region_attribute_map( const geode::BRep& model )
    {
        auto region_id{ 0 };
        absl::flat_hash_map< geode::uuid, geode::index_t > region_map_id;
        for( const auto& block : model.blocks() )
        {
            region_map_id.emplace( block.id(), region_id );
            region_id += 1;
        }
        return region_map_id;
    }
} // namespace
namespace geode
{
    class BRepGeosExporter::Impl
    {
    public:
        Impl( const BRep& brep, absl::string_view files_directory )
            : model_( brep ),
              files_directory_{ ghc::filesystem::path{
                  geode::to_string( files_directory ) }
                                    .string() },
              filename_prefix_{ geode::to_string( files_directory ) }
        {
            if( ghc::filesystem::path{ geode::to_string( files_directory ) }
                    .is_relative() )
            {
                ghc::filesystem::create_directory(
                    ghc::filesystem::current_path() / files_directory_ );
            }
            else
            {
                ghc::filesystem::create_directory( files_directory_ );
            }
        }
        absl::string_view files_directory() const
        {
            return files_directory_;
        }
        absl::string_view prefix() const
        {
            return filename_prefix_;
        }
        void export_meshes() const
        {
            std::string REGION_ID_ATTRIBUTE_NAME{ "attribute" };
            std::unique_ptr< geode::SurfaceMesh3D > surface;
            std::unique_ptr< geode::SolidMesh3D > solid;
            std::tie( surface, solid ) =
                geode::convert_brep_into_surface_and_solid( model_ );

            const auto brep_component_uuid_attribute =
                solid->polyhedron_attribute_manager()
                    .find_attribute<
                        geode::uuid_from_conversion_attribute_type >(
                        geode::uuid_from_conversion_attribute_name );

            const std::shared_ptr< geode::VariableAttribute< geode::index_t > >
                region_attribute{
                    solid->polyhedron_attribute_manager()
                        .template find_or_create_attribute<
                            geode::VariableAttribute, geode::index_t >(
                            REGION_ID_ATTRIBUTE_NAME, geode::NO_ID )
                };
            auto region_map_id = create_region_attribute_map( model_ );

            for( const auto polyhedron_id :
                geode::Range( solid->nb_polyhedra() ) )
            {
                region_attribute->set_value( polyhedron_id,
                    region_map_id
                        .find( brep_component_uuid_attribute->value(
                            polyhedron_id ) )
                        ->second );
            }

            solid->vertex_attribute_manager().delete_attribute(
                geode::unique_vertex_from_conversion_attribute_name );
            solid->polyhedron_attribute_manager().delete_attribute(
                geode::MESH_ELEMENT_ATTRIBUTE_NAME );
            solid->polyhedron_attribute_manager().delete_attribute(
                geode::uuid_from_conversion_attribute_name );

            const auto file =
                absl::StrCat( files_directory(), "/", prefix(), ".vtu" );
            if( const auto* tetra =
                    dynamic_cast< const geode::TetrahedralSolid3D* >(
                        solid.get() ) )
            {
                geode::save_tetrahedral_solid( *tetra, file );
            }
            else if( const auto* hybrid =
                         dynamic_cast< const geode::HybridSolid3D* >(
                             solid.get() ) )
            {
                geode::save_hybrid_solid( *hybrid, file );
            }
            else if( const auto* poly =
                         dynamic_cast< const geode::PolyhedralSolid3D* >(
                             solid.get() ) )
            {
                geode::save_polyhedral_solid( *poly, file );
            }
            else
            {
                throw geode::OpenGeodeException(
                    "[Blocks::save_geos] Cannot find the explicit "
                    "SolidMesh type" );
            }
        }

    private:
        const BRep& model_;
        DEBUG_CONST std::string files_directory_;
        DEBUG_CONST std::string filename_prefix_;
    };

} // namespace geode

namespace geode
{
    BRepGeosExporter::BRepGeosExporter(
        const BRep& brep, absl::string_view files_directory )
        : impl_( brep, files_directory )
    {
    }

    void BRepGeosExporter::export_meshes() const
    {
        impl_->export_meshes();
    }

} // namespace geode
