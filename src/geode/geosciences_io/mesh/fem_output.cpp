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
#include <geode/basic/filename.hpp>
#include <geode/basic/io.hpp>
#include <geode/basic/types.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/solid_edges.hpp>
#include <geode/mesh/core/solid_facets.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

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
        return std::string( n, ' ' );
    }

    class SolidFemOutputImpl
    {
        static constexpr const char CDATA_TAG_START[]{ "<![CDATA[" };
        static constexpr const char CDATA_TAG_END[]{ "]]>" };
        static constexpr const char APERTURE_ATTRIBUTE_NAME[]{
            "diagres_discontinuity_aperture"
        };
        static constexpr const char CONDUIT_AREA_ATTRIBUTE_NAME[]{
            "diagres_conduit_area"
        };
        static constexpr const char CONDUCTIVITY_ATTRIBUTE_NAME[]{
            "diagres_conductivity"
        };
        static constexpr const char SURFACE_NAME_ATTRIBUTE[]{ "surface_name" };
        static constexpr const char LINE_NAME_ATTRIBUTE[]{ "line_name" };

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
                    const auto vertex_id = vertex + 1;
                    absl::StrAppend( &poly_line, vertex_id, add_spaces( 1 ) );
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
                absl::StrAppend( &line, ranges, add_spaces( 1 ) );
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
                absl::StrAppend( &line, ranges, add_spaces( 1 ) );
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
                    const auto vertex_id = vtx + 1;
                    absl::StrAppend( &line, vertex_id, add_spaces( 1 ) );
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
                    absl::StrAppend( &line, ranges, add_spaces( 1 ) );
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
                    absl::StrAppend( &line, ranges, add_spaces( 1 ) );
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

        struct Property
        {
            std::string name;
            absl::flat_hash_map< double, std::vector< geode::index_t > >
                values_to_features;
        };

        struct Properties
        {
            std::vector< Property > properties;
        };

        struct DiscreteFeature1D
        {
            std::vector< geode::index_t > nodes;
            geode::index_t feature_id{ 1 };
        };

        struct DiscreteFeature2D
        {
            std::vector< geode::index_t > nodes;
            geode::index_t feature_id{ 1 };
        };

        template < typename Feature >
        struct FeatureGroup
        {
            std::string name;
            std::vector< Feature > features;
            std::vector< geode::index_t > feature_ids;
            std::vector< Property > properties;

            FeatureGroup( std::string group_name ) : name( group_name ) {}

            geode::index_t nb_features() const
            {
                return features.size();
            }

            std::optional< geode::index_t > get_property(
                const std::string& other_property_name )
            {
                geode::index_t property_index = 0;
                for( auto& property : properties )
                {
                    if( property.name == other_property_name )
                    {
                        return property_index;
                    }
                    property_index++;
                }
                return std::nullopt;
            }
        };

        struct DiscreteFeatures
        {
            std::vector< FeatureGroup< DiscreteFeature1D > > features1D_groups;
            std::vector< FeatureGroup< DiscreteFeature2D > > features2D_groups;

            Properties properties;

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
        };

        class DiscreteFeatureBuilder
        {
            const std::string FEATURE1D_GROUP_NAME = "Edge_Feature_LDS";
            const std::string FEATURE2D_GROUP_NAME = "Surface_Feature_LDS";

        public:
            DiscreteFeatureBuilder( DiscreteFeatures& features,
                const geode::TetrahedralSolid3D& solid )
                : features_( features ), solid_( solid ) {};

            void build_discrete_features()
            {
                build_2d_feature_groups();
                build_1d_feature_groups();
            }

        private:
            void build_2d_feature_groups()
            {
                auto& feature2D_group =
                    features_.features2D_groups.emplace_back(
                        FEATURE2D_GROUP_NAME );
                if( solid_.facets().facet_attribute_manager().attribute_exists(
                        APERTURE_ATTRIBUTE_NAME ) )
                {
                    const auto aperture_attribute =
                        solid_.facets()
                            .facet_attribute_manager()
                            .find_attribute< double >(
                                APERTURE_ATTRIBUTE_NAME );
                    const auto conductivity_attribute =
                        solid_.facets()
                            .facet_attribute_manager()
                            .find_or_create_attribute< geode::VariableAttribute,
                                double >( CONDUCTIVITY_ATTRIBUTE_NAME, -1.0 );
                    for( const auto facet :
                        geode::Range{ solid_.facets().nb_facets() } )
                    {
                        const auto aperture_value =
                            aperture_attribute->value( facet );
                        if( aperture_value < 0.0 )
                        {
                            continue;
                        }
                        std::vector< geode::index_t > facet_vertices;
                        for( const auto& vertex :
                            solid_.facets().facet_vertices( facet ) )
                        {
                            facet_vertices.push_back( vertex + 1 );
                        }
                        auto& feature2D =
                            feature2D_group.features.emplace_back();
                        feature2D.nodes = facet_vertices;
                        feature2D.feature_id = feature_id_;
                        feature2D_group.feature_ids.push_back( feature_id_ );
                        auto aperture_property = feature2D_group.get_property(
                            APERTURE_ATTRIBUTE_NAME );
                        if( !aperture_property )
                        {
                            Property new_property;
                            new_property.name = APERTURE_ATTRIBUTE_NAME;
                            absl::flat_hash_map< double,
                                std::vector< geode::index_t > >
                                property_map;
                            property_map[aperture_value] = { feature_id_ };
                            new_property.values_to_features = property_map;
                            feature2D_group.properties.push_back(
                                new_property );
                        }
                        else
                        {
                            feature2D_group
                                .properties[aperture_property.value()]
                                .values_to_features[aperture_value]
                                .push_back( feature_id_ );
                        }
                        const auto conductivity_value =
                            conductivity_attribute->value( facet );
                        if( conductivity_value > 0.0 )
                        {
                            auto conductivity_property =
                                feature2D_group.get_property(
                                    CONDUCTIVITY_ATTRIBUTE_NAME );
                            if( !conductivity_property )
                            {
                                Property new_property;
                                new_property.name = CONDUCTIVITY_ATTRIBUTE_NAME;
                                absl::flat_hash_map< double,
                                    std::vector< geode::index_t > >
                                    property_map;
                                property_map[conductivity_value] = {
                                    feature_id_
                                };
                                new_property.values_to_features = property_map;
                                feature2D_group.properties.push_back(
                                    new_property );
                            }
                            else
                            {
                                feature2D_group
                                    .properties[conductivity_property.value()]
                                    .values_to_features[conductivity_value]
                                    .push_back( feature_id_ );
                            }
                        }
                        feature_id_++;
                    }
                }
            }

            void build_1d_feature_groups()
            {
                auto& feature1D_group =
                    features_.features1D_groups.emplace_back(
                        FEATURE1D_GROUP_NAME );
                if( solid_.edges().edge_attribute_manager().attribute_exists(
                        CONDUIT_AREA_ATTRIBUTE_NAME ) )
                {
                    const auto conduit_area_attribute =
                        solid_.edges()
                            .edge_attribute_manager()
                            .find_attribute< double >(
                                CONDUIT_AREA_ATTRIBUTE_NAME );
                    const auto conductivity_attribute =
                        solid_.edges()
                            .edge_attribute_manager()
                            .find_or_create_attribute< geode::VariableAttribute,
                                double >( CONDUCTIVITY_ATTRIBUTE_NAME, -1.0 );
                    for( const auto edge :
                        geode::Range{ solid_.edges().nb_edges() } )
                    {
                        const auto conduit_area_value =
                            conduit_area_attribute->value( edge );
                        if( conduit_area_value < 0.0 )
                        {
                            continue;
                        }
                        std::vector< geode::index_t > edge_vertices;
                        for( const auto& vertex :
                            solid_.edges().edge_vertices( edge ) )
                        {
                            edge_vertices.push_back( vertex + 1 );
                        }
                        auto& feature1D =
                            feature1D_group.features.emplace_back();
                        feature1D.nodes = edge_vertices;
                        feature1D.feature_id = feature_id_;
                        feature1D_group.feature_ids.push_back( feature_id_ );
                        auto conduit_area_property =
                            feature1D_group.get_property(
                                CONDUIT_AREA_ATTRIBUTE_NAME );
                        if( !conduit_area_property )
                        {
                            Property new_property;
                            new_property.name = CONDUIT_AREA_ATTRIBUTE_NAME;
                            absl::flat_hash_map< double,
                                std::vector< geode::index_t > >
                                property_map;
                            property_map[conduit_area_value] = { feature_id_ };
                            new_property.values_to_features = property_map;
                            feature1D_group.properties.push_back(
                                new_property );
                        }
                        else
                        {
                            feature1D_group
                                .properties[conduit_area_property.value()]
                                .values_to_features[conduit_area_value]
                                .push_back( feature_id_ );
                        }
                        const auto conductivity_value =
                            conductivity_attribute->value( edge );
                        if( conductivity_value > 0.0 )
                        {
                            auto conductivity_property =
                                feature1D_group.get_property(
                                    CONDUCTIVITY_ATTRIBUTE_NAME );
                            if( !conductivity_property )
                            {
                                Property new_property;
                                new_property.name = CONDUCTIVITY_ATTRIBUTE_NAME;
                                absl::flat_hash_map< double,
                                    std::vector< geode::index_t > >
                                    property_map;
                                property_map[conductivity_value] = {
                                    feature_id_
                                };
                                new_property.values_to_features = property_map;
                                feature1D_group.properties.push_back(
                                    new_property );
                            }
                            else
                            {
                                feature1D_group
                                    .properties[conductivity_property.value()]
                                    .values_to_features[conductivity_value]
                                    .push_back( feature_id_ );
                            }
                        }
                        feature_id_++;
                    }
                }
            }

        private:
            DiscreteFeatures& features_;
            const geode::TetrahedralSolid3D& solid_;
            geode::index_t feature_id_ = 1;
        };

        std::string xml_start_tag( std::string tag )
        {
            return absl::StrCat( "<", tag, ">" );
        }

        std::string xml_start_tag( std::string tag, std::string value )
        {
            return absl::StrCat( "<", tag, "=", value, ">" );
        }

        std::string xml_end_tag( std::string tag )
        {
            return absl::StrCat( "</", tag, ">" );
        }

        void write_discrete_feature_header()
        {
            file_ << "DFE" << EOL;
            file_ << add_spaces( 2 )
                  << "<?xml version=\"1.0\" encoding=\"utf-8\" "
                     "standalone=\"no\" ?>"
                  << EOL;
            file_ << add_spaces( 2 ) << xml_start_tag( "fractures" ) << EOL;
        }

        void write_feature_signature(
            geode::index_t nb_features, std::string signature )
        {
            geode::index_t iterator = 0;
            while( iterator < nb_features )
            {
                file_ << add_spaces( 1 ) << signature;
                if( iterator < nb_features - 1 )
                {
                    file_ << EOL;
                }
                iterator++;
            }
        }

        void write_discrete_feature_signatures(
            const DiscreteFeatures& features )
        {
            const auto nb_features = features.nb_features();
            const auto fep_count_value_str =
                absl::StrCat( "\"", nb_features, "\"" );
            file_ << add_spaces( 4 )
                  << xml_start_tag( "fep count", fep_count_value_str ) << EOL;
            file_ << add_spaces( 6 ) << CDATA_TAG_START << EOL;
            for( const auto& feature2D_group : features.features2D_groups )
            {
                write_feature_signature(
                    feature2D_group.nb_features(), "c2d3,darcy" );
            }
            file_ << EOL;
            for( const auto& feature1D_group : features.features1D_groups )
            {
                write_feature_signature(
                    feature1D_group.nb_features(), "c1d2,darcy" );
            }
            file_ << CDATA_TAG_END << EOL;
            file_ << add_spaces( 4 ) << xml_end_tag( "fep" ) << EOL;
        }

        template < typename Feature >
        void write_feature_nodal_incidence_matrix(
            const std::vector< Feature >& features )
        {
            geode::index_t iterator = 0;
            for( const auto& feature : features )
            {
                std::string line =
                    absl::StrCat( add_spaces( 5 ), feature.nodes.size() );
                for( const auto node : feature.nodes )
                {
                    absl::StrAppend( &line, ",", add_spaces( 1 ), node );
                }
                file_ << line;
                if( iterator < features.size() - 1 )
                {
                    file_ << EOL;
                }
                iterator++;
            }
        }

        void write_discrete_feature_nodal_incidence_matrix(
            const DiscreteFeatures& features )
        {
            const auto nb_features = features.nb_features();
            const auto nop_count_value_str =
                absl::StrCat( "\"", nb_features, "\"" );
            file_ << add_spaces( 4 )
                  << xml_start_tag( "nop count", nop_count_value_str ) << EOL;
            file_ << add_spaces( 6 ) << CDATA_TAG_START << EOL;
            for( const auto& feature2D_group : features.features2D_groups )
            {
                write_feature_nodal_incidence_matrix< DiscreteFeature2D >(
                    feature2D_group.features );
            }
            file_ << EOL;
            for( const auto& feature1D_group : features.features1D_groups )
            {
                write_feature_nodal_incidence_matrix< DiscreteFeature1D >(
                    feature1D_group.features );
            }
            file_ << CDATA_TAG_END << EOL;
            file_ << add_spaces( 4 ) << xml_end_tag( "nop" ) << EOL;
        }

        template < typename Feature >
        void write_discrete_feature_group(
            std::vector< FeatureGroup< Feature > >& feature_groups )
        {
            for( const auto& feature_group : feature_groups )
            {
                geode::Logger::info( "herreeeeeee" );
                geode::Logger::info( feature_group.feature_ids.size() );
                const auto group_name_value_str =
                    absl::StrCat( "\"", feature_group.name, "\"",
                        " law=\"darcy\" mode=\"unstructured\"" );
                file_ << add_spaces( 6 )
                      << xml_start_tag( "group name", group_name_value_str )
                      << EOL;
                file_ << add_spaces( 8 ) << "<elements count=\""
                      << std::to_string( feature_group.nb_features() ) << "\">"
                      << EOL;
                file_ << add_spaces( 10 ) << CDATA_TAG_START;
                file_ << format_ranges( feature_group.feature_ids );
                file_ << CDATA_TAG_END << EOL;
                file_ << add_spaces( 8 ) << xml_end_tag( "elements" ) << EOL;
                file_ << add_spaces( 6 ) << xml_end_tag( "group" ) << EOL;
            }
        }

        void write_discrete_feature_groups( DiscreteFeatures& features )
        {
            const auto groups_count_value_str =
                absl::StrCat( "\"", features.nb_groups(), "\"" );
            file_ << add_spaces( 4 )
                  << xml_start_tag( "groups count", groups_count_value_str )
                  << EOL;
            write_discrete_feature_group< DiscreteFeature2D >(
                features.features2D_groups );
            write_discrete_feature_group< DiscreteFeature1D >(
                features.features1D_groups );
            file_ << add_spaces( 4 ) << xml_end_tag( "groups" ) << EOL;
        }

        void write_discrete_feature_properties_header()
        {
            file_ << add_spaces( 4 ) << xml_start_tag( "properties" ) << EOL;
            file_ << add_spaces( 6 ) << xml_start_tag( "flow" ) << EOL;
            file_ << add_spaces( 8 ) << xml_start_tag( "materials" ) << EOL;
        }

        void write_discrete_feature_properties_tail()
        {
            file_ << add_spaces( 8 ) << xml_end_tag( "materials" ) << EOL;
            file_ << add_spaces( 6 ) << xml_end_tag( "flow" ) << EOL;
            file_ << add_spaces( 4 ) << xml_end_tag( "properties" ) << EOL;
        }

        void write_material_property( std::string material_id,
            absl::flat_hash_map< double, std::vector< geode::index_t > >&
                property_values_to_features )
        {
            file_ << add_spaces( 10 )
                  << xml_start_tag( "material id", material_id ) << EOL;
            file_ << add_spaces( 12 ) << CDATA_TAG_START << EOL;
            geode::index_t iterator = 0;
            for( const auto& [value, feature_ids] :
                property_values_to_features )
            {
                file_ << add_spaces( 14 ) << value << add_spaces( 1 )
                      << format_ranges( feature_ids );
                if( iterator < property_values_to_features.size() - 1 )
                {
                    file_ << EOL;
                }
                iterator++;
            }
            file_ << add_spaces( 1 ) << CDATA_TAG_END << EOL;
            file_ << add_spaces( 10 ) << xml_end_tag( "material" ) << EOL;
        }

        void combine_maps(
            absl::flat_hash_map< double, std::vector< geode::index_t > >& map1,
            const absl::flat_hash_map< double, std::vector< geode::index_t > >&
                map2 )
        {
            for( const auto& [value, feature_ids] : map2 )
            {
                if( map1.find( value ) == map1.end() )
                {
                    map1[value] = feature_ids;
                    continue;
                }
                map1[value].insert(
                    map1[value].end(), feature_ids.begin(), feature_ids.end() );
            }
        }

        void write_discrete_feature_properties( DiscreteFeatures& features )
        {
            write_discrete_feature_properties_header();
            const auto aperture_property_index =
                features.features2D_groups[0].get_property(
                    APERTURE_ATTRIBUTE_NAME );
            const auto area_property_index =
                features.features1D_groups[0].get_property(
                    CONDUIT_AREA_ATTRIBUTE_NAME );
            if( aperture_property_index )
            {
                const auto& property2D =
                    features.features2D_groups[0]
                        .properties[aperture_property_index.value()];
                auto values_to_features = property2D.values_to_features;
                if( area_property_index )
                {
                    const auto& property1D =
                        features.features1D_groups[0]
                            .properties[aperture_property_index.value()];
                    combine_maps(
                        values_to_features, property1D.values_to_features );
                }
                write_material_property( "\"AREA\"", values_to_features );
            }
            const auto conductivity_property_index =
                features.features2D_groups[0].get_property(
                    CONDUCTIVITY_ATTRIBUTE_NAME );
            const auto conductivity_1D_property_index =
                features.features1D_groups[0].get_property(
                    CONDUCTIVITY_ATTRIBUTE_NAME );
            if( conductivity_property_index )
            {
                const auto& property =
                    features.features2D_groups[0]
                        .properties[conductivity_property_index.value()];
                auto values_to_features = property.values_to_features;
                if( conductivity_1D_property_index )
                {
                    const auto& property1D =
                        features.features1D_groups[0]
                            .properties[conductivity_1D_property_index.value()];
                    combine_maps(
                        values_to_features, property1D.values_to_features );
                }
                write_material_property( "\"COND\"", values_to_features );
            }
            write_discrete_feature_properties_tail();
        }

        void write_discrete_feature_tail()
        {
            file_ << add_spaces( 1 ) << xml_end_tag( "fractures" ) << EOL;
        }

        void write_discrete_features()
        {
            DiscreteFeatures features{};
            DiscreteFeatureBuilder builder{ features, solid_ };
            builder.build_discrete_features();
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
