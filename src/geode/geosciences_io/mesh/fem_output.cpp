/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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
#include <fstream>
#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/io.hpp>
#include <geode/basic/types.hpp>

#include <geode/basic/filename.hpp>

#include <geode/geosciences_io/mesh/internal/fem_output.hpp>

#include <geode/geometry/point.hpp>

namespace
{
    std::string format_ranges( absl::Span< const geode::index_t > elements )
    {
        if( elements.empty() )
        {
            return "";
        }
        std::string result;
        geode::index_t start = elements[0];
        geode::index_t prev = elements[0];
        for( const auto element : geode::Range{ 1, elements.size() } )
        {
            if( elements[element] == prev + 1 )
            {
                prev = elements[element];
                continue;
            }
            if( start == prev )
            {
                absl::StrAppendFormat( &result, "%d ", start );
            }
            else
            {
                absl::StrAppendFormat( &result, "%d-%d ", start, prev );
            }
            start = prev = elements[element];
        }
        if( start == prev )
        {
            absl::StrAppendFormat( &result, "%d", start );
        }
        else
        {
            absl::StrAppendFormat( &result, "%d-%d", start, prev );
        }

        return result;
    }

    class SolidFemOutputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };
        SolidFemOutputImpl(
            std::string_view filename, const geode::TetrahedralSolid3D& solid )
            : file_{ geode::to_string( filename ) }, solid_( solid )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void write_file()
        {
            geode::Logger::info( "[SolidFemOutput::write] Writing fem file." );
            const auto directories =
                geode::filepath_without_filename( "filename" );
            if( !directories.empty() )
            {
                std::filesystem::create_directories( directories );
            }
            write_problem();
            write_class();
            write_dimension();
            write_scale();
            write_var_node();
            write_node_coordinates();
            write_ref_nodal_dist();
            write_ref_element_dist();
            // write_mass_transport_properties();
            // write_flow_properties();
            // write_heat_transport_properties();
            write_nodal_sets();
            write_element_sets();
            write_gravity();
            write_end();
        }

        void write_problem()
        {
            file_ << "PROBLEM:" << EOL;
        }

        void write_class()
        {
            file_ << "CLASS (v.7.202.18152)" << EOL;
            std::string four_spaces = "    ";
            file_ << "   0    0    0    3    0    0    8    8    0    0" << EOL;
        }

        void write_dimension()
        {
            file_ << "DIMENS" << EOL;
            file_ << SPACE << solid_.nb_vertices() << SPACE
                  << solid_.nb_polyhedra() << " 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0"
                  << EOL;
        }

        void write_scale()
        {
            file_ << "SCALE" << EOL << EOL;
        }

        void write_var_node()
        {
            file_ << "VARNODE" << EOL;

            file_ << SPACE << solid_.nb_polyhedra() << " 4 4" << EOL;
            for( const auto polyhedron : geode::Range{ solid_.nb_polyhedra() } )
            {
                std::string poly_line = " 6 ";
                for( const auto vertex :
                    solid_.polyhedron_vertices( polyhedron ) )
                {
                    poly_line += std::to_string( vertex + 1 ) + SPACE;
                }
                file_ << poly_line << EOL;
            }
        }

        void write_node_coordinates()
        {
            file_ << "XYZCOOR" << EOL;
            for( const auto vertex : geode::Range{ solid_.nb_vertices() } )
            {
                file_ << SPACE << solid_.point( vertex ).value( 0 ) << ", "
                      << solid_.point( vertex ).value( 1 ) << ", "
                      << solid_.point( vertex ).value( 2 ) << EOL;
            }
        }

        void write_property( geode::AttributeBase& attribute )
        {
            file_ << attribute.name() << EOL;

            absl::flat_hash_map< double, std::vector< geode::index_t > >
                attribute_distribution =
                    create_attribute_distribution( attribute );
            std::vector< double > values;
            values.reserve( attribute_distribution.size() );
            for( auto val : attribute_distribution )
            {
                values.push_back( val.first );
            }
            absl::c_sort( values );
            for( const auto value : values )
            {
                std::string line = "";
                const auto ranges =
                    format_ranges( attribute_distribution[value] );
                line += ranges + SPACE;
                file_ << "  " << value << "  " << line << EOL;
            }
        }

        absl::flat_hash_map< double, std::vector< geode::index_t > >
            create_attribute_distribution( geode::AttributeBase& attribute )
        {
            absl::flat_hash_map< double, std::vector< geode::index_t > >
                att_dist;
            for( const auto elem : geode::Range{ solid_.nb_polyhedra() } )
            {
                const auto value = attribute.generic_value( elem );
                if( att_dist.contains( value ) )
                {
                    att_dist[value].push_back( elem );
                }
                else
                {
                    att_dist[value] = { elem };
                }
            }
            return att_dist;
        }

        void write_ref_nodal_dist()
        {
            geode::index_t nb_vtx_att{ 3 };
            std::vector< std::shared_ptr< geode::AttributeBase > > attributes;
            attributes.reserve(
                solid_.vertex_attribute_manager().attribute_names().size() );
            for( const auto name :
                solid_.vertex_attribute_manager().attribute_names() )
            {
                const auto attribute =
                    solid_.vertex_attribute_manager().find_generic_attribute(
                        name );
                if( !attribute || !attribute->is_genericable()
                    || name == "points" )
                {
                    continue;
                }
                attributes.push_back( attribute );
                nb_vtx_att++;
            }
            file_ << "REF_DIS_I" << EOL;
            file_ << SPACE << nb_vtx_att << "," << SPACE << solid_.nb_vertices()
                  << ","
                  << " 0" << EOL;
            write_xyz_dist();
            for( const auto& attribute : attributes )
            {
                write_property( *attribute );
            }
        }

        absl::flat_hash_map< double, std::vector< geode::index_t > >
            create_porosity_region_map(
                const geode::ReadOnlyAttribute< double >& porosity_attribute )
        {
            absl::flat_hash_map< double, std::vector< geode::index_t > >
                porosity_region_map;
            for( const auto polyhedron : geode::Range{ solid_.nb_polyhedra() } )
            {
                const auto value = porosity_attribute.value( polyhedron );
                if( !porosity_region_map.contains( value ) )
                {
                    porosity_region_map[value] = { polyhedron };
                    continue;
                }
                porosity_region_map[value].push_back( polyhedron );
            }
            return porosity_region_map;
        }

        void write_porosity()
        {
            if( !solid_.polyhedron_attribute_manager().attribute_exists(
                    "Porosity" ) )
            {
                return;
            }
            const auto& porosity_attribute =
                solid_.polyhedron_attribute_manager().find_attribute< double >(
                    "Porosity" );
            if( !porosity_attribute )
            {
                return;
            }
            file_ << "MAT_I_TRAN" << EOL;
            file_ << "201 0.30 \"Porosity\"" << EOL;
            const auto porosity_map =
                create_porosity_region_map( *porosity_attribute );
            for( const auto& [porosity_value, tetrahedra] : porosity_map )
            {
                std::string line = "";
                const auto ranges = format_ranges( tetrahedra );
                line += ranges + SPACE;
                file_ << "  " << porosity_value << "  " << line << EOL;
            }
        }

        void write_properties() {}

        void write_mass_transport_properties()
        {
            write_porosity();
            // TODO: which other properties to add?
        }

        void write_flow_properties()
        {
            // TODO: which properties to add?
        }

        void write_heat_transport_properties()
        { // TODO: which properties to add?
        }

        absl::flat_hash_map< double, std::vector< geode::index_t > >
            create_coord_dist( int dim )
        {
            absl::flat_hash_map< double, std::vector< geode::index_t > >
                coord_dist;
            for( const auto vertex : geode::Range{ solid_.nb_vertices() } )
            {
                if( coord_dist.contains( solid_.point( vertex ).value( dim ) ) )
                {
                    coord_dist[solid_.point( vertex ).value( dim )].push_back(
                        vertex );
                }
                else
                {
                    coord_dist[solid_.point( vertex ).value( dim )] = {
                        vertex
                    };
                }
            }
            return coord_dist;
        }

        void write_one_coord_component_dist(
            const std::string& dim_name, const geode::index_t dim )
        {
            file_ << dim_name << EOL;
            absl::flat_hash_map< double, std::vector< geode::index_t > > dist =
                create_coord_dist( dim );
            std::vector< double > values;
            values.reserve( dist.size() );
            for( auto val = dist.begin(); val != dist.end(); ++val )
            {
                values.push_back( val->first );
            }
            std::sort( values.begin(), values.end() );
            for( const auto value : values )
            {
                std::string line = "";
                for( const auto vtx : dist[value] )
                {
                    line += std::to_string( vtx + 1 ) + SPACE;
                }
                file_ << "     " << value << "  " << line << EOL;
            }
        }

        void write_xyz_dist()
        {
            write_one_coord_component_dist( "X", 0 );
            write_one_coord_component_dist( "Y", 1 );
            write_one_coord_component_dist( "Z", 2 );
        }

        void write_ref_element_dist()
        {
            file_ << "REF_DISE_I" << EOL;
            int nb_elem_att =
                solid_.polyhedron_attribute_manager().attribute_names().size()
                - 3;
            file_ << SPACE << nb_elem_att << "," << SPACE
                  << solid_.nb_polyhedra() << ","
                  << " 0" << EOL;
            for( const auto name :
                solid_.polyhedron_attribute_manager().attribute_names() )
            {
                const auto attribute = solid_.polyhedron_attribute_manager()
                                           .find_generic_attribute( name );
                if( !attribute || !attribute->is_genericable()
                    || name == "tetrahedron_vertices"
                    || name == "tetrahedron_adjacents"
                    || name == "geode_active" )
                {
                    continue;
                }
                write_property( *attribute );
            }
        }

        void write_nodal_sets()
        {
            if( solid_.vertex_attribute_manager().attribute_exists(
                    "Block_ID_vertex" ) )
            {
                auto attribute_v =
                    solid_.vertex_attribute_manager()
                        .find_attribute< std::vector< std::string_view > >(
                            "Block_ID_vertex" );
                file_ << "NODALSETS" << EOL;
                const auto vertex_regions =
                    create_region_map( *attribute_v, false );
                for( const auto& vertex_region : vertex_regions )
                {
                    std::string line = "";
                    const auto ranges = format_ranges( vertex_region.second );
                    line += ranges + SPACE;
                    file_ << "  " << vertex_region.first << "  " << line << EOL;
                }
            }
        }

        void write_element_sets()
        {
            if( solid_.polyhedron_attribute_manager().attribute_exists(
                    "Block_ID_polyhedron" ) )
            {
                file_ << "ELEMENTALSETS" << EOL;
                auto attribute_p = solid_.polyhedron_attribute_manager()
                                       .find_attribute< std::string_view >(
                                           "Block_ID_polyhedron" );
                const auto elem_regions =
                    create_region_map( *attribute_p, true );
                for( const auto& elem_region : elem_regions )
                {
                    std::string line = "";
                    const auto ranges = format_ranges( elem_region.second );
                    line += ranges + SPACE;
                    file_ << "  " << elem_region.first << "  " << line << EOL;
                }
            }
        }

        absl::flat_hash_map< std::string_view, std::vector< geode::index_t > >
            create_region_map(
                geode::ReadOnlyAttribute< std::string_view >& attribute,
                bool create_element_region )
        {
            absl::flat_hash_map< std::string_view,
                std::vector< geode::index_t > >
                region_map;
            geode::index_t nb_obj =
                ( create_element_region ? solid_.nb_polyhedra()
                                        : solid_.nb_vertices() );
            for( const auto obj : geode::Range{ nb_obj } )
            {
                const auto value = attribute.value( obj );
                if( region_map.contains( value ) )
                {
                    region_map[value].push_back( obj + 1 );
                }
                else
                {
                    region_map[value] = { obj + 1 };
                }
            }
            return region_map;
        }

        absl::flat_hash_map< std::string_view, std::vector< geode::index_t > >
            create_region_map(
                geode::ReadOnlyAttribute< std::vector< std::string_view > >&
                    attribute,
                bool create_element_region )
        {
            absl::flat_hash_map< std::string_view,
                std::vector< geode::index_t > >
                region_map;
            int nb_obj = ( create_element_region ? solid_.nb_polyhedra()
                                                 : solid_.nb_vertices() );
            for( const auto obj : geode::Range{ nb_obj } )
            {
                const auto& values = attribute.value( obj );
                if( values.empty() )
                {
                    continue;
                }
                for( const auto& value : values )
                {
                    if( region_map.contains( value ) )
                    {
                        region_map[value].push_back( obj + 1 );
                    }
                    else
                    {
                        region_map[value] = { obj + 1 };
                    }
                }
            }
            return region_map;
        }

        void write_gravity()
        {
            file_ << "GRAVITY" << EOL;
            file_ << " 0 0 -1" << EOL;
        }

        void write_end()
        {
            file_ << "END" << EOL;
        }

    private:
        std::ofstream file_;
        const geode::TetrahedralSolid3D& solid_;
        std::vector< std::shared_ptr< geode::AttributeBase > > generic_att_{};
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::vector< std::string > SolidFemOutput::write(
            const TetrahedralSolid3D& solid ) const
        {
            SolidFemOutputImpl impl{ filename(), solid };
            impl.write_file();
            return { to_string( filename() ) };
        }
    } // namespace internal
} // namespace geode
