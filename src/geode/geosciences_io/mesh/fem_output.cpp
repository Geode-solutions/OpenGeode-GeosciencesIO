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
#include <cmath>
#include <fstream>
#include <limits>

#include <absl/algorithm/container.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/io.hpp>
#include <geode/basic/types.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/solid_edges.hpp>
#include <geode/mesh/core/solid_facets.hpp>
#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/geosciences_io/mesh/internal/fem_output.hpp>

#include <geode/geometry/point.hpp>

#include <geode/basic/detail/disable_debug_logger.hpp>

namespace
{
    static constexpr char EOL{ '\n' };
    static constexpr char SPACE{ ' ' };
    static constexpr geode::index_t CHUNK_SIZE{ 1500 };

    std::vector< std::string > format_range_chunks(
        absl::Span< const geode::index_t > elements,
        geode::index_t max_chunk_size )
    {
        std::vector< std::string > results;
        if( elements.empty() )
        {
            return results;
        }
        std::string current;
        geode::index_t start = elements[0];
        geode::index_t prev = elements[0];
        auto append_value = [&]( geode::index_t s, geode::index_t p ) {
            std::string tmp;
            if( s == p )
            {
                absl::StrAppendFormat( &tmp, "%d", s );
            }
            else
            {
                absl::StrAppendFormat( &tmp, "%d-%d", s, p );
            }
            if( current.size() + tmp.size() + 1 > max_chunk_size )
            {
                results.push_back( current );
                current.clear();
            }

            if( !current.empty() )
            {
                current += " ";
            }
            current += tmp;
        };
        for( const auto i : geode::Range{ 1, elements.size() } )
        {
            if( elements[i] == prev + 1 )
            {
                prev = elements[i];
                continue;
            }
            append_value( start, prev );
            start = prev = elements[i];
        }
        append_value( start, prev );
        if( !current.empty() )
        {
            results.push_back( std::move( current ) );
        }
        return results;
    }

    template < typename T >
    std::shared_ptr< geode::ReadOnlyAttribute< T > > find_typed_attribute(
        const geode::AttributeManager& manager, std::string_view name )
    {
        const auto attribute_ids = manager.attribute_ids_matching_name( name );
        if( !attribute_ids )
        {
            return nullptr;
        }
        for( const auto& attribute_id : attribute_ids.value() )
        {
            if( manager.attribute_type( attribute_id ) == typeid( T ).name() )
            {
                return manager.find_read_only_attribute< T >( attribute_id );
            }
        }
        geode::Logger::warning( "[SolidFemOutput::write] Attribute ", name,
            " is ignored: its type is not the expected one" );
        return nullptr;
    }

    double attribute_value( const geode::ReadOnlyAttribute< double >& attribute,
        geode::index_t element )
    {
        return attribute.value( element );
    }

    double attribute_value(
        const geode::AttributeBase& attribute, geode::index_t element )
    {
        return attribute.generic_value( element );
    }

    // Sorts the (value, id) pairs and calls func( value, ids ) for each value,
    // ids being sorted
    template < typename Value, typename Func >
    void for_each_value_group(
        std::vector< std::pair< Value, geode::index_t > >& values_ids,
        const Func& func )
    {
        absl::c_sort( values_ids );
        std::vector< geode::index_t > ids;
        geode::index_t group_start{ 0 };
        for( const auto i : geode::Indices{ values_ids } )
        {
            if( values_ids[i].first != values_ids[group_start].first )
            {
                func( values_ids[group_start].first, ids );
                ids.clear();
                group_start = i;
            }
            ids.push_back( values_ids[i].second );
        }
        if( !ids.empty() )
        {
            func( values_ids[group_start].first, ids );
        }
    }

    std::string add_spaces( geode::index_t n )
    {
        return std::string( n, SPACE );
    }

    class SolidFemOutputImpl
    {
        static constexpr auto CDATA_TAG_START = "<![CDATA[";
        static constexpr auto CDATA_TAG_END = "]]>";
        static constexpr auto APERTURE_ATTRIBUTE_NAME =
            "diagres_discontinuity_aperture";
        static constexpr auto CONDUIT_AREA_ATTRIBUTE_NAME =
            "diagres_conduit_area";
        static constexpr auto CONDUCTIVITY_ATTRIBUTE_NAME =
            "diagres_conductivity";
        static constexpr auto SURFACE_NAME_ATTRIBUTE = "surface_name";
        static constexpr auto LINE_NAME_ATTRIBUTE = "line_name";

    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        SolidFemOutputImpl(
            std::string_view filename, const geode::TetrahedralSolid3D& solid )
            : file_{ geode::to_string( filename ) }, solid_( solid )
        {
            file_.precision( std::numeric_limits< double >::max_digits10 );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                file_.good(), nullptr, geode::OpenGeodeException::TYPE::data,
                "Error while opening file: ", filename );
        }

        void write_file()
        {
            geode::Logger::info( "[SolidFemOutput::write] Writing fem file." );
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

        template < typename Value >
        void write_chunks( absl::Span< const geode::index_t > elements,
            const Value& attribute_value,
            geode::index_t max_chunk_size )
        {
            for( const auto& range_chunk :
                format_range_chunks( elements, max_chunk_size ) )
            {
                std::string line = "";
                absl::StrAppend( &line, range_chunk, add_spaces( 1 ) );
                file_ << "  " << attribute_value << "  " << line << EOL;
            }
        }

        void write_property(
            geode::AttributeBase& attribute, geode::index_t nb_elements )
        {
            file_ << attribute.name().value() << EOL;
            auto values_ids =
                attribute.type() == typeid( double ).name()
                    ? attribute_values_ids(
                          dynamic_cast<
                              const geode::ReadOnlyAttribute< double >& >(
                              attribute ),
                          nb_elements )
                    : attribute_values_ids( attribute, nb_elements );
            for_each_value_group( values_ids,
                [this]( double value, absl::Span< const geode::index_t > ids ) {
                    write_chunks( ids, value, CHUNK_SIZE );
                } );
        }

        template < typename Attribute >
        std::vector< std::pair< double, geode::index_t > > attribute_values_ids(
            const Attribute& attribute, geode::index_t nb_elements )
        {
            std::vector< std::pair< double, geode::index_t > > values_ids;
            values_ids.reserve( nb_elements );
            for( const auto elem : geode::Range{ nb_elements } )
            {
                const auto value = attribute_value( attribute, elem );
                if( !std::isnan( value ) )
                {
                    values_ids.emplace_back( value, elem + 1 );
                }
            }
            if( const auto nb_nan_values = nb_elements - values_ids.size();
                nb_nan_values != 0 )
            {
                geode::Logger::warning( "[SolidFemOutput::write] ",
                    nb_nan_values, " NaN values of attribute ",
                    attribute.name().value(), " are not exported" );
            }
            return values_ids;
        }

        static std::vector< std::shared_ptr< geode::AttributeBase > >
            exported_attributes( const geode::AttributeManager& manager )
        {
            std::vector< std::shared_ptr< geode::AttributeBase > > attributes;
            for( const auto& id : manager.attribute_ids() )
            {
                auto attribute = manager.find_generic_attribute( id );
                if( !attribute || !attribute->is_genericable()
                    || !attribute->properties().transferable )
                {
                    continue;
                }
                if( attribute->nb_items() != 1 )
                {
                    geode::Logger::warning( "[SolidFemOutput::write] "
                                            "Attribute ",
                        attribute->name().value(),
                        " is not exported: only scalar attributes are "
                        "supported" );
                    continue;
                }
                attributes.push_back( std::move( attribute ) );
            }
            return attributes;
        }

        void write_ref_nodal_dist()
        {
            static constexpr geode::index_t NB_COORDINATES{ 3 };
            const auto attributes =
                exported_attributes( solid_.vertex_attribute_manager() );
            file_ << "REF_DIS_I" << EOL;
            file_ << SPACE << NB_COORDINATES + attributes.size() << "," << SPACE
                  << solid_.nb_vertices() << ","
                  << " 0" << EOL;
            write_xyz_dist();
            for( const auto& attribute : attributes )
            {
                write_property( *attribute, solid_.nb_vertices() );
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
            const auto porosity_attribute_ids =
                solid_.polyhedron_attribute_manager()
                    .attribute_ids_matching_name( "Porosity" );
            if( !porosity_attribute_ids.has_value() )
            {
                return;
            }
            const auto& porosity_attribute =
                solid_.polyhedron_attribute_manager()
                    .find_read_only_attribute< double >(
                        porosity_attribute_ids.value().front() );
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
                write_chunks( tetrahedra, porosity_value, CHUNK_SIZE );
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

        void write_one_coord_component_dist(
            const std::string& dim_name, const geode::local_index_t dim )
        {
            file_ << dim_name << EOL;
            std::vector< std::pair< double, geode::index_t > > values_ids;
            values_ids.reserve( solid_.nb_vertices() );
            for( const auto vertex : geode::Range{ solid_.nb_vertices() } )
            {
                values_ids.emplace_back(
                    solid_.point( vertex ).value( dim ), vertex + 1 );
            }
            std::string line;
            for_each_value_group(
                values_ids, [this, &line]( double value,
                                absl::Span< const geode::index_t > ids ) {
                    line.clear();
                    for( const auto vertex_id : ids )
                    {
                        absl::StrAppend( &line, vertex_id, add_spaces( 1 ) );
                    }
                    file_ << "     " << value << "  " << line << EOL;
                } );
        }

        void write_xyz_dist()
        {
            write_one_coord_component_dist( "X", 0 );
            write_one_coord_component_dist( "Y", 1 );
            write_one_coord_component_dist( "Z", 2 );
        }

        void write_ref_element_dist()
        {
            const auto attributes =
                exported_attributes( solid_.polyhedron_attribute_manager() );
            file_ << "REF_DISE_I" << EOL;
            file_ << SPACE << attributes.size() << "," << SPACE
                  << solid_.nb_polyhedra() << ","
                  << " 0" << EOL;
            for( const auto& attribute : attributes )
            {
                write_property( *attribute, solid_.nb_polyhedra() );
            }
        }

        void write_nodal_sets()
        {
            const auto attribute_v =
                find_typed_attribute< std::vector< std::string > >(
                    solid_.vertex_attribute_manager(), "Block_ID_vertex" );
            if( !attribute_v )
            {
                return;
            }
            file_ << "NODALSETS" << EOL;
            std::vector< std::pair< std::string_view, geode::index_t > >
                names_ids;
            names_ids.reserve( solid_.nb_vertices() );
            for( const auto vertex : geode::Range{ solid_.nb_vertices() } )
            {
                for( const auto& name : attribute_v->value( vertex ) )
                {
                    names_ids.emplace_back( name, vertex + 1 );
                }
            }
            write_sets( names_ids );
        }

        void write_element_sets()
        {
            const auto attribute_p = find_typed_attribute< std::string >(
                solid_.polyhedron_attribute_manager(), "Block_ID_polyhedron" );
            if( !attribute_p )
            {
                return;
            }
            file_ << "ELEMENTALSETS" << EOL;
            std::vector< std::pair< std::string_view, geode::index_t > >
                names_ids;
            names_ids.reserve( solid_.nb_polyhedra() );
            for( const auto polyhedron : geode::Range{ solid_.nb_polyhedra() } )
            {
                names_ids.emplace_back(
                    attribute_p->value( polyhedron ), polyhedron + 1 );
            }
            write_sets( names_ids );
        }

        void write_sets(
            std::vector< std::pair< std::string_view, geode::index_t > >&
                names_ids )
        {
            for_each_value_group(
                names_ids, [this]( std::string_view name,
                               absl::Span< const geode::index_t > ids ) {
                    write_chunks( ids, name, CHUNK_SIZE );
                } );
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
                if( solid_.are_facets_enabled() )
                {
                    const auto& manager =
                        solid_.facets().facet_attribute_manager();
                    if( const auto aperture_attribute =
                            find_typed_attribute< double >(
                                manager, APERTURE_ATTRIBUTE_NAME ) )
                    {
                        build_2d_feature_groups( *aperture_attribute,
                            find_typed_attribute< double >(
                                manager, CONDUCTIVITY_ATTRIBUTE_NAME )
                                .get() );
                    }
                }
                if( solid_.are_edges_enabled() )
                {
                    const auto& manager =
                        solid_.edges().edge_attribute_manager();
                    if( const auto conduit_area_attribute =
                            find_typed_attribute< double >(
                                manager, CONDUIT_AREA_ATTRIBUTE_NAME ) )
                    {
                        build_1d_feature_groups( *conduit_area_attribute,
                            find_typed_attribute< double >(
                                manager, CONDUCTIVITY_ATTRIBUTE_NAME )
                                .get() );
                    }
                }
            }

        private:
            void build_2d_feature_groups(
                const geode::ReadOnlyAttribute< double >& aperture_attribute,
                const geode::ReadOnlyAttribute< double >*
                    conductivity_attribute )
            {
                auto& feature2D_group =
                    features_.features2D_groups.emplace_back(
                        FEATURE2D_GROUP_NAME );
                bool some_features_built{ false };
                for( const auto facet :
                    geode::Range{ solid_.facets().nb_facets() } )
                {
                    const auto aperture_value =
                        aperture_attribute.value( facet );
                    if( !( aperture_value >= 0.0 ) )
                    {
                        continue;
                    }
                    std::vector< geode::index_t > facet_vertices;
                    for( const auto& vertex :
                        solid_.facets().facet_vertices( facet ) )
                    {
                        facet_vertices.push_back( vertex + 1 );
                    }
                    auto& feature2D = feature2D_group.features.emplace_back();
                    feature2D.nodes = facet_vertices;
                    feature2D.feature_id = feature_id_;
                    feature2D_group.feature_ids.push_back( feature_id_ );
                    auto aperture_property =
                        feature2D_group.get_property( APERTURE_ATTRIBUTE_NAME );
                    if( !aperture_property )
                    {
                        Property new_property;
                        new_property.name = APERTURE_ATTRIBUTE_NAME;
                        absl::flat_hash_map< double,
                            std::vector< geode::index_t > >
                            property_map;
                        property_map[aperture_value] = { feature_id_ };
                        new_property.values_to_features = property_map;
                        feature2D_group.properties.push_back( new_property );
                    }
                    else
                    {
                        feature2D_group.properties[aperture_property.value()]
                            .values_to_features[aperture_value]
                            .push_back( feature_id_ );
                    }
                    const auto conductivity_value =
                        conductivity_attribute
                            ? conductivity_attribute->value( facet )
                            : -1.0;
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
                            property_map[conductivity_value] = { feature_id_ };
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
                    some_features_built = true;
                    feature_id_++;
                }
                if( !some_features_built )
                {
                    features_.features2D_groups.pop_back();
                    geode::Logger::warning(
                        "[SolidFemOutput::write] No 2D discrete feature "
                        "written: no facet has a non-negative ",
                        APERTURE_ATTRIBUTE_NAME, " value" );
                }
            }

            void build_1d_feature_groups(
                const geode::ReadOnlyAttribute< double >&
                    conduit_area_attribute,
                const geode::ReadOnlyAttribute< double >*
                    conductivity_attribute )
            {
                auto& feature1D_group =
                    features_.features1D_groups.emplace_back(
                        FEATURE1D_GROUP_NAME );
                bool some_features_built{ false };
                for( const auto edge :
                    geode::Range{ solid_.edges().nb_edges() } )
                {
                    const auto conduit_area_value =
                        conduit_area_attribute.value( edge );
                    if( !( conduit_area_value >= 0.0 ) )
                    {
                        continue;
                    }
                    std::vector< geode::index_t > edge_vertices;
                    for( const auto& vertex :
                        solid_.edges().edge_vertices( edge ) )
                    {
                        edge_vertices.push_back( vertex + 1 );
                    }
                    auto& feature1D = feature1D_group.features.emplace_back();
                    feature1D.nodes = edge_vertices;
                    feature1D.feature_id = feature_id_;
                    feature1D_group.feature_ids.push_back( feature_id_ );
                    auto conduit_area_property = feature1D_group.get_property(
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
                        feature1D_group.properties.push_back( new_property );
                    }
                    else
                    {
                        feature1D_group
                            .properties[conduit_area_property.value()]
                            .values_to_features[conduit_area_value]
                            .push_back( feature_id_ );
                    }
                    const auto conductivity_value =
                        conductivity_attribute
                            ? conductivity_attribute->value( edge )
                            : -1.0;
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
                            property_map[conductivity_value] = { feature_id_ };
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
                    some_features_built = true;
                    feature_id_++;
                }
                if( !some_features_built )
                {
                    features_.features1D_groups.pop_back();
                    geode::Logger::warning(
                        "[SolidFemOutput::write] No 1D discrete feature "
                        "written: no edge has a non-negative ",
                        CONDUIT_AREA_ATTRIBUTE_NAME, " value" );
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
                const auto chunks = format_range_chunks(
                    feature_group.feature_ids, CHUNK_SIZE );
                for( const auto chunk_id : geode::Range{ chunks.size() } )
                {
                    const auto& chunk = chunks[chunk_id];
                    file_ << chunk;
                    if( chunk_id < chunks.size() - 1 )
                    {
                        file_ << EOL;
                    }
                }
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
            std::vector< double > values;
            values.reserve( property_values_to_features.size() );
            for( const auto& value_features : property_values_to_features )
            {
                values.push_back( value_features.first );
            }
            absl::c_sort( values );
            geode::index_t iterator = 0;
            for( const auto value : values )
            {
                const auto range_chunks = format_range_chunks(
                    property_values_to_features.at( value ), CHUNK_SIZE );
                for( const auto chunk : geode::Indices{ range_chunks } )
                {
                    file_ << add_spaces( 14 ) << value << add_spaces( 1 )
                          << range_chunks[chunk];
                    if( chunk + 1 < range_chunks.size()
                        || iterator + 1 < property_values_to_features.size() )
                    {
                        file_ << EOL;
                    }
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
            if( features.features2D_groups.empty()
                && features.features1D_groups.empty() )
            {
                return;
            }
            write_discrete_feature_properties_header();
            std::optional< geode::index_t > aperture_property_index =
                std::nullopt;
            if( !features.features2D_groups.empty() )
            {
                aperture_property_index =
                    features.features2D_groups.front().get_property(
                        APERTURE_ATTRIBUTE_NAME );
            }
            std::optional< geode::index_t > area_property_index = std::nullopt;
            if( !features.features1D_groups.empty() )
            {
                area_property_index =
                    features.features1D_groups.front().get_property(
                        CONDUIT_AREA_ATTRIBUTE_NAME );
            }
            if( aperture_property_index )
            {
                const auto& property2D =
                    features.features2D_groups.front()
                        .properties[aperture_property_index.value()];
                auto values_to_features = property2D.values_to_features;
                if( area_property_index )
                {
                    const auto& property1D =
                        features.features1D_groups.front()
                            .properties[area_property_index.value()];
                    combine_maps(
                        values_to_features, property1D.values_to_features );
                }
                write_material_property( "\"AREA\"", values_to_features );
            }
            else
            {
                if( area_property_index )
                {
                    const auto& property1D =
                        features.features1D_groups.front()
                            .properties[area_property_index.value()];
                    auto values_to_features = property1D.values_to_features;
                    write_material_property( "\"AREA\"", values_to_features );
                }
            }
            std::optional< geode::index_t > conductivity_property_index =
                std::nullopt;
            if( !features.features2D_groups.empty() )
            {
                conductivity_property_index =
                    features.features2D_groups.front().get_property(
                        CONDUCTIVITY_ATTRIBUTE_NAME );
            }
            std::optional< geode::index_t > conductivity_1D_property_index =
                std::nullopt;
            if( !features.features1D_groups.empty() )
            {
                conductivity_1D_property_index =
                    features.features1D_groups.front().get_property(
                        CONDUCTIVITY_ATTRIBUTE_NAME );
            }
            if( conductivity_property_index )
            {
                const auto& property =
                    features.features2D_groups.front()
                        .properties[conductivity_property_index.value()];
                auto values_to_features = property.values_to_features;
                if( conductivity_1D_property_index )
                {
                    const auto& property1D =
                        features.features1D_groups.front()
                            .properties[conductivity_1D_property_index.value()];
                    combine_maps(
                        values_to_features, property1D.values_to_features );
                }
                write_material_property( "\"COND\"", values_to_features );
            }
            else
            {
                if( conductivity_1D_property_index )
                {
                    const auto& property1D =
                        features.features1D_groups.front()
                            .properties[conductivity_1D_property_index.value()];
                    auto values_to_features = property1D.values_to_features;
                    write_material_property( "\"COND\"", values_to_features );
                }
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
            if( features.nb_groups() == 0 )
            {
                return;
            }
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

namespace geode::internal
{
    std::vector< std::string > SolidFemOutput::write(
        const TetrahedralSolid3D& solid ) const
    {
        SolidFemOutputImpl impl{ filename(), solid };
        impl.write_file();
        return { to_string( filename() ) };
    }
} // namespace geode::internal

#include <geode/basic/detail/enable_debug_logger.hpp>