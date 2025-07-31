/*
 * Copyright (c) 2019 - 2025 Geode-solutions. All rights reserved.
 */

#include <geode/geosciences_io/model/internal/strati_lso_input.hpp>

#include <fstream>
#include <optional>
#include <string>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/file.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/solid_mesh.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>

#include <geode/model/mixin/core/block.hpp>

#include <geode/geosciences/explicit/representation/io/structural_model_input.hpp>
#include <geode/geosciences/implicit/geometry/stratigraphic_point.hpp>
#include <geode/geosciences/implicit/representation/builder/stratigraphic_model_builder.hpp>
#include <geode/geosciences/implicit/representation/core/stratigraphic_model.hpp>

namespace
{
    class StratigraphicLSOInputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr auto GEOL_ATTRIBUTE_NAME =
            "SnS!ingridStructuralWorkflow/skua_model_Geology";
        static constexpr auto STRATI_ATTRIBUTE_NAME =
            "SnS!ingridStructuralWorkflow/skua_model_stratigraphy";

        StratigraphicLSOInputImpl( std::string_view filename )
            : filename_( filename )
        {
        }

        geode::StratigraphicModel read_file()
        {
            geode::StratigraphicModel model{
                geode::StructuralModelInputFactory::create( "lso", filename_ )
                    ->read()
            };
            geode::StratigraphicModelBuilder builder{ model };

            for( const auto& block : model.blocks() )
            {
                auto& manager = block.mesh().vertex_attribute_manager();
                OPENGEODE_DATA_EXCEPTION(
                    manager.attribute_exists( STRATI_ATTRIBUTE_NAME )
                        && manager.attribute_exists( GEOL_ATTRIBUTE_NAME ),
                    "[ImplicitLSOInput] Could not find the properties "
                    "associated to StratigraphicModeling in the file, "
                    "named '",
                    GEOL_ATTRIBUTE_NAME, "' and '", STRATI_ATTRIBUTE_NAME,
                    "'." );
                OPENGEODE_DATA_EXCEPTION(
                    ( block.mesh().type_name()
                        == geode::TetrahedralSolid3D::type_name_static() ),
                    "[ImplicitLSOInput] Blocks must be meshed as "
                    "TetrahedralSolids, which is not the case for block with "
                    "uuid '",
                    block.id().string(), "'." );
                const auto strati_attribute =
                    manager.find_attribute< double >( STRATI_ATTRIBUTE_NAME );
                const auto geol_attribute =
                    manager.find_attribute< std::array< double, 2 > >(
                        GEOL_ATTRIBUTE_NAME );
                for( const auto vertex_id :
                    geode::Range{ block.mesh().nb_vertices() } )
                {
                    auto geol_value = geol_attribute->value( vertex_id );
                    builder.set_stratigraphic_coordinates( block, vertex_id,
                        geode::StratigraphicPoint3D{
                            { geol_value[0], geol_value[1],
                                strati_attribute->value( vertex_id ) } } );
                }
            }
            return model;
        }

    private:
        std::string filename_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        StratigraphicModel StratigraphicLSOInput::read()
        {
            StratigraphicLSOInputImpl impl{ filename() };
            return impl.read_file();
        }

        Percentage StratigraphicLSOInput::is_loadable() const
        {
            const auto structural_percent =
                is_structural_model_loadable( filename() );
            if( structural_percent.value() != 1 )
            {
                return structural_percent;
            }
            std::ifstream file{ to_string( filename() ) };
            const auto line = goto_keyword_if_it_exists( file, "PROPERTIES" );
            if( !line )
            {
                return Percentage{ 0 };
            }
            if( line->find( StratigraphicLSOInputImpl::GEOL_ATTRIBUTE_NAME )
                == std::string::npos )
            {
                return Percentage{ 0 };
            }
            if( line->find( StratigraphicLSOInputImpl::STRATI_ATTRIBUTE_NAME )
                == std::string::npos )
            {
                return Percentage{ 0 };
            }
            return Percentage{ 1 };
        }
    } // namespace internal
} // namespace geode
