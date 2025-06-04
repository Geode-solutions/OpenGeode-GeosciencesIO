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
    static constexpr char EOL{ '\n' };
    static constexpr char SPACE{ ' ' };

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

    std::string add_spaces( geode::index_t n )
    {
        std::string result;
        for( geode::index_t i = 0; i < n; ++i )
        {
            result += SPACE;
        }
        return result;
    }

    class SolidFemOutputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
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
            write_discrete_features();
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
                    att_dist[value].push_back( elem + 1 );
                }
                else
                {
                    att_dist[value] = { elem + 1 };
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

        struct DiscreteFeature1D
        {
            // std::string name;
            std::vector< geode::index_t > nodes;
            geode::index_t feature_id{ 1 };
        };

        struct DiscreteFeature2D
        {
            // std::string name;
            std::vector< geode::index_t > nodes;
            geode::index_t feature_id{ 1 };
        };

        template < typename Feature >
        struct FeatureGroup
        {
            std::string name;
            std::vector< Feature > features;
            std::vector< geode::index_t > feature_ids;

            FeatureGroup( std::string group_name ) : name( group_name ) {}

            geode::index_t nb_features() const
            {
                return features.size();
            }
        };

        struct DiscreteFeatures
        {
            std::vector< FeatureGroup< DiscreteFeature1D > > features1D_groups;
            std::vector< FeatureGroup< DiscreteFeature2D > > features2D_groups;

            geode::index_t nb_features() const
            {
                geode::index_t nb_features = 0;
                for( const auto& feature1D_group : features1D_groups )
                {
                    nb_features += feature1D_group.nb_features();
                }
                for( const auto& feature2D_group : features2D_groups )
                {
                    nb_features += feature2D_group.nb_features();
                }
                return nb_features;
            }

            geode::index_t nb_groups() const
            {
                return features1D_groups.size() + features2D_groups.size();
            }

            void build_feature_ids()
            {
                geode::index_t feature_id = 1;
                for( auto& feature_group : features2D_groups )
                {
                    for( auto& feature : feature_group.features )
                    {
                        feature.feature_id = feature_id;
                        feature_group.feature_ids.push_back( feature_id );
                        feature_id++;
                    }
                }
                for( auto& feature_group : features1D_groups )
                {
                    for( auto& feature : feature_group.features )
                    {
                        feature.feature_id = feature_id;
                        feature_group.feature_ids.push_back( feature_id );
                        feature_id++;
                    }
                }
            }
        };

        DiscreteFeatures build_discrete_features()
        {
            DiscreteFeature1D feature1D;
            feature1D.nodes = { 1, 2 };
            DiscreteFeature1D feature1D_b;
            feature1D_b.nodes = { 3, 4 };
            DiscreteFeature1D feature1D_c;
            feature1D_c.nodes = { 5, 6 };
            DiscreteFeature2D feature2D;
            feature2D.nodes = { 10, 11, 12 };
            DiscreteFeature2D feature2D_b;
            feature2D_b.nodes = { 13, 14, 15 };
            FeatureGroup< DiscreteFeature1D > feature1D_group_1{ "well_1" };
            feature1D_group_1.features = { feature1D, feature1D_b };
            FeatureGroup< DiscreteFeature1D > feature1D_group_2{ "well_2" };
            feature1D_group_2.features = { feature1D_c };
            FeatureGroup< DiscreteFeature2D > feature2D_group_1{ "surface_1" };
            feature2D_group_1.features = { feature2D, feature2D_b };
            DiscreteFeatures features;
            features.features1D_groups = { feature1D_group_1,
                feature1D_group_2 };
            features.features2D_groups = { feature2D_group_1 };

            return features;
        }

        void write_discrete_feature_header()
        {
            file_ << "DFE" << EOL;
            file_ << add_spaces( 1 )
                  << "<?xml version=\"1.0\" encoding=\"utf-8\" "
                     "standalone=\"no\" ?>"
                  << EOL;
            file_ << add_spaces( 1 ) << "<fractures>" << EOL;
        }

        void write_2d_feature_signature(
            const std::vector< DiscreteFeature2D >& features )
        {
            for( const auto& feature : features )
            {
                std::string line =
                    "  " + std::to_string( feature.nodes.size() );
                for( const auto node : feature.nodes )
                {
                    line += " " + std::to_string( node );
                }
                file_ << line << EOL;
            }
        }

        void write_feature_signature(
            geode::index_t nb_features, std::string signature )
        {
            geode::index_t iterator = 0;
            while( iterator < nb_features )
            {
                file_ << add_spaces( 2 ) << signature << EOL;
                iterator++;
            }
        }

        void write_discrete_feature_signatures(
            const DiscreteFeatures& features )
        {
            const auto nb_features = features.nb_features();
            file_ << add_spaces( 2 )
                  << "<fep count=\"" + std::to_string( nb_features ) + "\">"
                  << EOL;
            file_ << add_spaces( 3 ) << "<![CDATA[" << EOL;
            for( const auto& feature2D_group : features.features2D_groups )
            {
                write_feature_signature(
                    feature2D_group.nb_features(), "c2d3,darcy" );
            }
            for( const auto& feature1D_group : features.features1D_groups )
            {
                write_feature_signature(
                    feature1D_group.nb_features(), "c1d2,darcy" );
            }
            file_ << add_spaces( 3 ) << "]]>" << EOL;
            file_ << add_spaces( 2 ) << "</fep>" << EOL;
        }

        template < typename Feature >
        void write_feature_nodal_incidence_matrix(
            const std::vector< Feature >& features )
        {
            for( const auto& feature : features )
            {
                std::string line =
                    "  " + std::to_string( feature.nodes.size() );
                for( const auto node : feature.nodes )
                {
                    line += " " + std::to_string( node );
                }
                file_ << line << EOL;
            }
        }

        void write_discrete_feature_nodal_incidence_matrix(
            const DiscreteFeatures& features )
        {
            const auto nb_features = features.nb_features();
            file_ << add_spaces( 2 )
                  << "<nop count=\"" + std::to_string( nb_features ) + "\">"
                  << EOL;
            file_ << add_spaces( 3 ) << "<![CDATA[" << EOL;
            for( const auto& feature2D_group : features.features2D_groups )
            {
                write_feature_nodal_incidence_matrix< DiscreteFeature2D >(
                    feature2D_group.features );
            }
            for( const auto& feature1D_group : features.features1D_groups )
            {
                write_feature_nodal_incidence_matrix< DiscreteFeature1D >(
                    feature1D_group.features );
            }
            file_ << add_spaces( 3 ) << "]]>" << EOL;
            file_ << add_spaces( 2 ) << "</nop>" << EOL;
        }

        template < typename Feature >
        void write_discrete_feature_group(
            std::vector< FeatureGroup< Feature > >& feature_groups )
        {
            for( const auto& feature_group : feature_groups )
            {
                file_ << add_spaces( 3 ) << "<group name=\""
                      << feature_group.name
                      << "\" law=\"darcy\" mode=\"unstructured\">" << EOL;
                file_ << add_spaces( 4 ) << "<elements count=\""
                      << std::to_string( feature_group.nb_features() ) << "\">"
                      << EOL;
                file_ << add_spaces( 5 ) << "<![CDATA[";
                file_ << format_ranges( feature_group.feature_ids );
                file_ << "]]>" << EOL;
                file_ << add_spaces( 3 ) << "</group>" << EOL;
            }
        }

        void write_discrete_feature_groups( DiscreteFeatures& features )
        {
            features.build_feature_ids();
            file_ << add_spaces( 2 ) << "<groups count=\""
                  << features.nb_groups() << "\">" << EOL;
            write_discrete_feature_group< DiscreteFeature2D >(
                features.features2D_groups );
            write_discrete_feature_group< DiscreteFeature1D >(
                features.features1D_groups );
            file_ << add_spaces( 2 ) << "</groups>" << EOL;
        }

        void write_discrete_feature_properties_header()
        {
            file_ << add_spaces( 2 ) << "<properties>" << EOL;
            file_ << add_spaces( 3 ) << "<flow>" << EOL;
            file_ << add_spaces( 4 ) << "<materials>" << EOL;
        }

        void write_discrete_feature_properties_tail()
        {
            file_ << add_spaces( 4 ) << "</materials>" << EOL;
            file_ << add_spaces( 3 ) << "</flow>" << EOL;
            file_ << add_spaces( 2 ) << "</properties>" << EOL;
        }

        void write_discrete_feature_properties(
            const DiscreteFeatures& features )
        {
            geode_unused( features );
            write_discrete_feature_properties_header();
            // todo
            write_discrete_feature_properties_tail();
        }

        void write_discrete_feature_tail()
        {
            file_ << " </fractures>" << EOL;
        }

        void write_discrete_features()
        {
            auto features = build_discrete_features();
            write_discrete_feature_header();
            write_discrete_feature_signatures( features );
            write_discrete_feature_nodal_incidence_matrix( features );
            write_discrete_feature_groups( features );
            write_discrete_feature_properties( features );
            write_discrete_feature_tail();
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
