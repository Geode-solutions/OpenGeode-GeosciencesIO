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

#include <geode/geosciences_io/mesh/internal/json_grid_input.hpp>

#include <cmath>
#include <fstream>

#include <nlohmann/json.hpp>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_set.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/identifier_builder.hpp>
#include <geode/basic/percentage.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/core/light_regular_grid.hpp>

namespace
{
    constexpr auto GRID_KEY = "grid_df";
    constexpr auto NAME_KEY = "name";
    constexpr std::array< std::string_view, 3 > COORDINATE_KEYS{ "x", "y",
        "z" };
    constexpr std::array< std::string_view, 2 > BOOL_ATTRIBUTE_KEYS{
        "grid_topo_name", "grid_site_name"
    };

    nlohmann::json load_json( std::string_view filename )
    {
        std::ifstream file{ geode::to_string( filename ), std::ios::binary };
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            file.good(), nullptr, geode::OpenGeodeException::TYPE::data,
            "[JSONGridInput] Error while opening file: ", filename );
        nlohmann::json json;
        file >> json;
        return json;
    }

    double coordinate( const nlohmann::json& record, geode::local_index_t axis )
    {
        return record.at( COORDINATE_KEYS[axis] ).get< double >();
    }

    bool same_coordinate( const nlohmann::json& record0,
        const nlohmann::json& record1,
        geode::local_index_t axis )
    {
        return std::fabs(
                   coordinate( record0, axis ) - coordinate( record1, axis ) )
               <= geode::GLOBAL_EPSILON;
    }

    template < typename T >
    using NamedAttributes = std::vector< std::pair< std::string,
        std::shared_ptr< geode::VariableAttribute< T > > > >;

    struct GridAttributes
    {
        NamedAttributes< bool > bools;
        NamedAttributes< double > doubles;
    };

    template < typename T >
    std::shared_ptr< geode::VariableAttribute< T > > create_attribute(
        geode::AttributeManager& manager,
        std::string_view name,
        geode::AttributeValues< T > values )
    {
        geode::AttributeProperties properties;
        properties.assignable = false;
        properties.interpolable = false;
        properties.transferable = true;
        const auto attribute_id =
            manager.create_attribute< geode::VariableAttribute, T >(
                name, std::move( values ), std::move( properties ) );
        return manager.find_attribute< geode::VariableAttribute, T >(
            attribute_id );
    }

    void set_bool_value( geode::VariableAttribute< bool >& attribute,
        geode::index_t vertex,
        const nlohmann::json& record,
        const std::string& name )
    {
        const auto value_it = record.find( name );
        if( value_it == record.end() )
        {
            return;
        }
        if( value_it->is_boolean() )
        {
            attribute.set_value( vertex, value_it->get< bool >() );
        }
        else if( value_it->is_number() )
        {
            attribute.set_value( vertex, value_it->get< double >() != 0 );
        }
    }

    void set_double_value( geode::VariableAttribute< double >& attribute,
        geode::index_t vertex,
        const nlohmann::json& record,
        const std::string& name )
    {
        const auto value_it = record.find( name );
        if( value_it != record.end() && value_it->is_number() )
        {
            attribute.set_value( vertex, value_it->get< double >() );
        }
    }

    class JSONGridInputImpl
    {
    public:
        explicit JSONGridInputImpl( std::string_view filename )
            : json_( load_json( filename ) )
        {
            const auto grid_it = json_.find( GRID_KEY );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                grid_it != json_.end() && grid_it->is_array()
                    && !grid_it->empty(),
                nullptr, geode::OpenGeodeException::TYPE::data,
                "[JSONGridInput] Missing or empty ", GRID_KEY, " array" );
            records_ = std::move( *grid_it );
        }

        geode::LightRegularGrid3D read_file()
        {
            compute_grid_dimensions();
            geode::LightRegularGrid3D grid{ origin(), cells_number(),
                cells_length() };
            // TODO: coordinate_system (EPSG code) is not stored since
            // LightRegularGrid does not support coordinate reference systems
            set_grid_name( grid );
            read_attributes( grid );
            return grid;
        }

    private:
        [[nodiscard]] const nlohmann::json& record( geode::index_t index ) const
        {
            return records_[index];
        }

        /*
         * Records are ordered with z varying fastest, then x, then y
         */
        void compute_grid_dimensions()
        {
            const auto nb_records =
                static_cast< geode::index_t >( records_.size() );
            const auto& first = record( 0 );
            auto& nb_x = nb_vertices_[0];
            auto& nb_y = nb_vertices_[1];
            auto& nb_z = nb_vertices_[2];
            nb_z = 1;
            while( nb_z < nb_records
                   && same_coordinate( first, record( nb_z ), 0 )
                   && same_coordinate( first, record( nb_z ), 1 ) )
            {
                nb_z++;
            }
            nb_x = 1;
            while( nb_x * nb_z < nb_records
                   && same_coordinate( first, record( nb_x * nb_z ), 1 ) )
            {
                nb_x++;
            }
            nb_y = nb_records / ( nb_x * nb_z );
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                nb_x * nb_y * nb_z == nb_records, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[JSONGridInput] Number of records (", nb_records,
                ") does not match a regular grid of ", nb_x, "x", nb_z,
                " vertices per y layer" );
        }

        [[nodiscard]] geode::index_t nb_records() const
        {
            return nb_vertices_[0] * nb_vertices_[1] * nb_vertices_[2];
        }

        /*
         * Number of records between two consecutive vertices along the axis
         */
        [[nodiscard]] geode::index_t stride( geode::local_index_t axis ) const
        {
            const std::array< geode::index_t, 3 > strides{ nb_vertices_[2],
                nb_vertices_[0] * nb_vertices_[2], 1 };
            return strides[axis];
        }

        [[nodiscard]] geode::LightRegularGrid3D::VertexIndices vertex_indices(
            geode::index_t record_id ) const
        {
            return { ( record_id / stride( 0 ) ) % nb_vertices_[0],
                record_id / stride( 1 ), record_id % nb_vertices_[2] };
        }

        [[nodiscard]] geode::Point3D record_point(
            geode::index_t record_id ) const
        {
            const auto& current = record( record_id );
            return geode::Point3D{ { coordinate( current, 0 ),
                coordinate( current, 1 ), coordinate( current, 2 ) } };
        }

        [[nodiscard]] geode::Point3D origin() const
        {
            return record_point( 0 );
        }

        [[nodiscard]] std::array< geode::index_t, 3 > cells_number() const
        {
            return { nb_vertices_[0] - 1, nb_vertices_[1] - 1,
                nb_vertices_[2] - 1 };
        }

        [[nodiscard]] double cell_length( geode::local_index_t axis ) const
        {
            if( nb_vertices_[axis] < 2 )
            {
                return 1.;
            }
            return coordinate( record( stride( axis ) ), axis )
                   - coordinate( record( 0 ), axis );
        }

        [[nodiscard]] std::array< double, 3 > cells_length() const
        {
            return { cell_length( 0 ), cell_length( 1 ), cell_length( 2 ) };
        }

        [[nodiscard]] absl::flat_hash_set< std::string >
            bool_attribute_names() const
        {
            absl::flat_hash_set< std::string > names;
            for( const auto& key : BOOL_ATTRIBUTE_KEYS )
            {
                const auto name_it = json_.find( key );
                if( name_it != json_.end() && name_it->is_string() )
                {
                    names.emplace( name_it->get< std::string >() );
                }
            }
            return names;
        }

        [[nodiscard]] GridAttributes create_attributes(
            geode::LightRegularGrid3D& grid ) const
        {
            const auto bool_names = bool_attribute_names();
            auto& manager = grid.grid_vertex_attribute_manager();
            GridAttributes attributes;
            for( const auto& item : record( 0 ).items() )
            {
                const auto& name = item.key();
                if( absl::c_find( COORDINATE_KEYS, name )
                    != COORDINATE_KEYS.end() )
                {
                    continue;
                }
                if( bool_names.contains( name ) )
                {
                    geode::AttributeValues< bool > values;
                    values.default_value = false;
                    attributes.bools.emplace_back(
                        name, create_attribute( manager, name, values ) );
                    continue;
                }
                geode::AttributeValues< double > values;
                values.default_value = std::nan( "" );
                values.no_value = std::nan( "" );
                attributes.doubles.emplace_back(
                    name, create_attribute( manager, name, values ) );
            }
            return attributes;
        }

        void set_grid_name( geode::LightRegularGrid3D& grid ) const
        {
            const auto name_it = json_.find( NAME_KEY );
            if( name_it != json_.end() && name_it->is_string() )
            {
                geode::IdentifierBuilder{ grid }.set_name(
                    name_it->get< std::string >() );
            }
        }

        void read_attributes( geode::LightRegularGrid3D& grid ) const
        {
            const auto attributes = create_attributes( grid );
            for( const auto record_id : geode::Range{ nb_records() } )
            {
                const auto vertex =
                    grid.vertex_index( vertex_indices( record_id ) );
                geode::OpenGeodeGeosciencesIOMeshException::check_assertion(
                    grid.point( vertex ).inexact_equal(
                        record_point( record_id ) ),
                    "[JSONGridInput] Record ", record_id,
                    " does not match the expected grid vertex" );
                const auto& current = record( record_id );
                for( const auto& [name, attribute] : attributes.bools )
                {
                    set_bool_value( *attribute, vertex, current, name );
                }
                for( const auto& [name, attribute] : attributes.doubles )
                {
                    set_double_value( *attribute, vertex, current, name );
                }
            }
        }

    private:
        nlohmann::json json_;
        nlohmann::json records_;
        std::array< geode::index_t, 3 > nb_vertices_{ 0, 0, 0 };
    };
} // namespace

namespace geode::internal
{
    LightRegularGrid3D JSONGridInput::read()
    {
        JSONGridInputImpl reader{ filename() };
        return reader.read_file();
    }

    Percentage JSONGridInput::is_loadable() const
    {
        try
        {
            const auto json = load_json( filename() );
            const auto grid_it = json.find( GRID_KEY );
            if( grid_it == json.end() || !grid_it->is_array()
                || grid_it->empty() || !grid_it->front().is_object() )
            {
                return Percentage{ 0 };
            }
            for( const auto& key : COORDINATE_KEYS )
            {
                if( !grid_it->front().contains( key ) )
                {
                    return Percentage{ 0 };
                }
            }
            return Percentage{ 1 };
        }
        catch( ... )
        {
            return Percentage{ 0 };
        }
    }
} // namespace geode::internal
