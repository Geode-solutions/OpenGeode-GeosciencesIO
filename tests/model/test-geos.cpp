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

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/uuid.hpp>
#include <geode/basic/variable_attribute.hpp>
#include <geode/tests_config.hpp>

#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>

#include <geode/geosciences_io/model/helpers/brep_geos_export.hpp>

#include <geode/geometry/point.hpp>
#include <geode/io/mesh/common.hpp>
#include <geode/io/model/common.hpp>
#include <geode/mesh/builder/point_set_builder.hpp>
#include <geode/mesh/core/point_set.hpp>

#include <geode/mesh/core/geode/geode_point_set.hpp>
#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/io/hybrid_solid_input.hpp>

#include <geode/model/mixin/core/block.hpp>
#include <geode/model/mixin/core/physical_properties.hpp>
#include <geode/model/representation/builder/brep_builder.hpp>
#include <geode/model/representation/core/brep.hpp>
#include <geode/model/representation/io/brep_input.hpp>
#include <geode/model/representation/io/brep_output.hpp>

void test_picasso()
{
    // Load structural model
    auto model =
        geode::load_brep( absl::StrCat( geode::DATA_PATH, "picasso.og_brep" ) );
    geode::BRepGeosExporter exporter( model, "picasso" );
    exporter.run();
}
void toy_model()
{
    auto model = geode::load_brep( absl::StrCat(
        geode::DATA_PATH, "adaptive_brep_perm_and_poro.og_brep" ) );
    geode::BRepGeosExporter exporter( model, "toy_model" );
    exporter.add_cell_property_1d( "permeability" );
    exporter.add_cell_property_1d( "porosity" );
    auto point_set = geode::PointSet3D::create(
        geode::OpenGeodePointSet3D::impl_name_static() );
    auto builder = geode::PointSetBuilder3D::create( *point_set );
    builder->create_point( geode::Point3D{ { 20., 20., 10. } } );
    exporter.add_well_perforations( *point_set );
    exporter.run();
}
namespace
{

    std::vector< double > read_geos_file( std::string_view filename )
    {
        std::ifstream file{ absl::StrCat( geode::DATA_PATH, filename ) };
        geode::OpenGeodeGeosciencesIOModelException::test(
            file.good(), "Cannot open ", filename );
        return { std::istream_iterator< double >{ file },
            std::istream_iterator< double >{} };
    }

    geode::index_t closest_index(
        const std::vector< double >& coordinates, double value )
    {
        const auto closest = std::min_element( coordinates.begin(),
            coordinates.end(), [value]( double lhs, double rhs ) {
                return std::fabs( lhs - value ) < std::fabs( rhs - value );
            } );
        return static_cast< geode::index_t >(
            std::distance( coordinates.begin(), closest ) );
    }

    struct Spe10Data
    {
        Spe10Data()
            : xlin{ read_geos_file( "xlin.geos" ) },
              ylin{ read_geos_file( "ylin.geos" ) },
              zlin{ read_geos_file( "zlin.geos" ) },
              permx{ read_geos_file( "permx.geos" ) },
              permy{ read_geos_file( "permy.geos" ) },
              permz{ read_geos_file( "permz.geos" ) },
              poro{ read_geos_file( "poro.geos" ) }
        {
            const auto nb_cells = xlin.size() * ylin.size() * zlin.size();
            geode::OpenGeodeGeosciencesIOModelException::test(
                permx.size() == nb_cells && permy.size() == nb_cells
                    && permz.size() == nb_cells && poro.size() == nb_cells,
                "SPE10 data sizes are inconsistent" );
        }

        geode::index_t cell_index( const geode::Point3D& center ) const
        {
            const auto i = closest_index( xlin, center.value( 0 ) );
            const auto j = closest_index( ylin, center.value( 1 ) );
            const auto k = closest_index( zlin, center.value( 2 ) );
            return i + j * xlin.size() + k * xlin.size() * ylin.size();
        }

        std::array< double, 3 > permeability( geode::index_t cell ) const
        {
            return { permx[cell] * MILLIDARCY_TO_M2,
                permy[cell] * MILLIDARCY_TO_M2,
                permz[cell] * MILLIDARCY_TO_M2 };
        }

        std::vector< double > xlin;
        std::vector< double > ylin;
        std::vector< double > zlin;
        std::vector< double > permx;
        std::vector< double > permy;
        std::vector< double > permz;
        std::vector< double > poro;
    };

    void create_block_unique_vertices( geode::BRep& model )
    {
        if( model.nb_unique_vertices() != 0 )
        {
            return;
        }
        geode::Logger::info( "Creating missing unique vertices" );
        geode::BRepBuilder builder{ model };
        for( const auto& block : model.blocks() )
        {
            const auto nb_vertices = block.mesh().nb_vertices();
            const auto first = builder.create_unique_vertices( nb_vertices );
            for( const auto vertex : geode::Range{ nb_vertices } )
            {
                builder.set_unique_vertex(
                    { block.component_id(), vertex }, first + vertex );
            }
        }
    }

    void add_spe10_physical_properties(
        geode::BRep& model, const Spe10Data& spe10 )
    {
        const geode::uuid permeability_id;
        const geode::uuid porosity_id;
        geode::AttributeProperties properties;
        properties.assignable = false;
        properties.interpolable = false;
        properties.transferable = true;
        geode::index_t nb_cells{ 0 };
        for( const auto& block : model.blocks() )
        {
            const auto& mesh = block.mesh();
            auto& manager = mesh.polyhedron_attribute_manager();
            manager.create_attribute< geode::VariableAttribute,
                std::array< double, 3 > >( "permeability", permeability_id,
                { { 0., 0., 0. }, { 0., 0., 0. } }, properties );
            manager.create_attribute< geode::VariableAttribute, double >(
                "porosity", porosity_id, { 0., 0. }, properties );
            auto permeability =
                manager.find_attribute< geode::VariableAttribute,
                    std::array< double, 3 > >( permeability_id );
            auto porosity =
                manager.find_attribute< geode::VariableAttribute, double >(
                    porosity_id );
            for( const auto polyhedron : geode::Range{ mesh.nb_polyhedra() } )
            {
                const auto cell = spe10.cell_index(
                    mesh.polyhedron_barycenter( polyhedron ) );
                permeability->set_value(
                    polyhedron, spe10.permeability( cell ) );
                porosity->set_value( polyhedron, spe10.poro[cell] );
            }
            nb_cells += mesh.nb_polyhedra();
        }
        geode::OpenGeodeGeosciencesIOModelException::test(
            nb_cells == spe10.poro.size(),
            "Wrong number of cells in grid: ", nb_cells, " instead of ",
            spe10.poro.size() );
        geode::BRepBuilder builder{ model };
        builder.set_physical_property(
            geode::PHYSICAL_PROPERTY_NAME::permeability,
            geode::Block3D::component_type_static(), permeability_id );
        builder.set_physical_property( geode::PHYSICAL_PROPERTY_NAME::porosity,
            geode::Block3D::component_type_static(), porosity_id );
    }

} // namespace

void test_grid_geos()
{
    auto model = geode::load_brep(
        absl::StrCat( geode::DATA_PATH, "grid_geos.og_brep" ) );
    create_block_unique_vertices( model );
    const Spe10Data spe10;
    add_spe10_physical_properties( model, spe10 );
    geode::save_brep( model, "grid_geos_with_physical_properties.og_brep" );
    geode::BRepGeosExporter exporter( model, "grid_geos" );
    exporter.run();
}

int main()
{
    try
    {
        geode::OpenGeodeGeosciencesIOModelLibrary::initialize();
        geode::OpenGeodeIOMeshLibrary::initialize();
        geode::OpenGeodeIOModelLibrary::initialize();
        // test_picasso();
        // toy_model();
        test_grid_geos();
        geode::Logger::info( "TEST SUCCESS" );

        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}
