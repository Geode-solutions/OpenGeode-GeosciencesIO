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
#include <limits>
#include <optional>

#include <nlohmann/json.hpp>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/identifier_builder.hpp>
#include <geode/basic/percentage.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/core/light_regular_grid.hpp>

/*
 * Expected file content:
 * {
 *   "grid_df": [ { "x": 2.5, "y": 2.5, "z": 10.0, "<attribute>": 1.0, ... },
 *                ... ],
 *   "name": "<grid name>",           (optional)
 *   "grid_topo_name": "<attribute>", (optional, attribute read as bool)
 *   "grid_site_name": "<attribute>", (optional, attribute read as bool)
 *   ...                              (other keys are ignored)
 * }
 * Each record of grid_df is a grid vertex, records can be in any order.
 */
namespace
{
    constexpr auto GRID_KEY = "grid_df";
    constexpr auto NAME_KEY = "name";
    constexpr std::array< std::string_view, 3 > COORDINATE_KEYS{ "x", "y",
        "z" };
    constexpr std::array< std::string_view, 2 > BOOL_ATTRIBUTE_KEYS{
        "grid_topo_name", "grid_site_name"
    };
    constexpr auto NO_DATA = std::numeric_limits< double >::quiet_NaN();

    /*
     * Values of one key across all records, NaN where the record has no
     * numeric value for this key
     */
    using Column = std::vector< double >;

    /*
     * SAX handler storing grid_df records column by column, so that the file
     * is never held in memory as a JSON document
     */
    class JSONGridSax
    {
    public:
        using number_integer_t = nlohmann::json::number_integer_t;
        using number_unsigned_t = nlohmann::json::number_unsigned_t;
        using number_float_t = nlohmann::json::number_float_t;
        using string_t = nlohmann::json::string_t;
        using binary_t = nlohmann::json::binary_t;

        /*
         * Parsing stops once max_records records have been read
         */
        explicit JSONGridSax( geode::index_t max_records )
            : max_records_( max_records )
        {
        }

        bool null()
        {
            return true;
        }

        bool boolean( bool value )
        {
            return number( value ? 1. : 0. );
        }

        bool number_integer( number_integer_t value )
        {
            return number( static_cast< double >( value ) );
        }

        bool number_unsigned( number_unsigned_t value )
        {
            return number( static_cast< double >( value ) );
        }

        bool number_float( number_float_t value, const string_t& /*unused*/ )
        {
            return number( value );
        }

        bool string( string_t& value )
        {
            if( depth_ == 1 )
            {
                top_level_strings_[current_key_] = std::move( value );
            }
            return true;
        }

        bool binary( binary_t& /*unused*/ )
        {
            return true;
        }

        bool start_object( std::size_t /*unused*/ )
        {
            depth_++;
            if( depth_ == 1 )
            {
                has_root_object_ = true;
            }
            else if( in_record() )
            {
                current_key_.clear();
            }
            return true;
        }

        bool end_object()
        {
            const auto record_ended = in_record();
            depth_--;
            if( record_ended )
            {
                nb_records_++;
                return nb_records_ < max_records_;
            }
            return true;
        }

        bool start_array( std::size_t /*unused*/ )
        {
            depth_++;
            if( depth_ == 2 && current_key_ == GRID_KEY )
            {
                in_grid_ = true;
            }
            return true;
        }

        bool end_array()
        {
            if( depth_ == 2 )
            {
                in_grid_ = false;
            }
            depth_--;
            return true;
        }

        bool key( string_t& value )
        {
            if( depth_ == 1 || in_record() )
            {
                current_key_ = std::move( value );
            }
            return true;
        }

        bool parse_error( std::size_t position,
            const std::string& /*unused*/,
            const nlohmann::detail::exception& exception )
        {
            error_ =
                absl::StrCat( "at byte ", position, ": ", exception.what() );
            return false;
        }

        [[nodiscard]] bool has_root_object() const
        {
            return has_root_object_;
        }

        [[nodiscard]] geode::index_t nb_records() const
        {
            return nb_records_;
        }

        [[nodiscard]] const std::string& error() const
        {
            return error_;
        }

        [[nodiscard]] const std::array< Column, 3 >& coordinates() const
        {
            return coordinates_;
        }

        [[nodiscard]] const std::vector< std::pair< std::string, Column > >&
            attributes() const
        {
            return attributes_;
        }

        [[nodiscard]] std::optional< std::string > top_level_string(
            std::string_view key ) const
        {
            const auto it = top_level_strings_.find( key );
            if( it == top_level_strings_.end() )
            {
                return std::nullopt;
            }
            return it->second;
        }

        /*
         * Pads every column with NaN up to the number of records
         */
        void finalize()
        {
            for( auto& coordinate : coordinates_ )
            {
                coordinate.resize( nb_records_, NO_DATA );
            }
            for( auto& [name, values] : attributes_ )
            {
                values.resize( nb_records_, NO_DATA );
            }
        }

    private:
        [[nodiscard]] bool in_record() const
        {
            return in_grid_ && depth_ == 3;
        }

        bool number( double value )
        {
            if( in_record() )
            {
                auto& values = column( current_key_ );
                values.resize( nb_records_ + 1, NO_DATA );
                values.back() = value;
            }
            return true;
        }

        Column& column( const std::string& key )
        {
            for( const auto axis : geode::LRange{ 3 } )
            {
                if( key == COORDINATE_KEYS[axis] )
                {
                    return coordinates_[axis];
                }
            }
            const auto [it, inserted] = attribute_ids_.try_emplace(
                key, static_cast< geode::index_t >( attributes_.size() ) );
            if( inserted )
            {
                attributes_.emplace_back( key, Column{} );
            }
            return attributes_[it->second].second;
        }

    private:
        geode::index_t max_records_;
        geode::index_t depth_{ 0 };
        bool has_root_object_{ false };
        bool in_grid_{ false };
        std::string current_key_;
        geode::index_t nb_records_{ 0 };
        std::array< Column, 3 > coordinates_;
        std::vector< std::pair< std::string, Column > > attributes_;
        absl::flat_hash_map< std::string, geode::index_t > attribute_ids_;
        absl::flat_hash_map< std::string, std::string > top_level_strings_;
        std::string error_;
    };

    JSONGridSax parse_json( std::string_view filename,
        geode::index_t max_records =
            std::numeric_limits< geode::index_t >::max() )
    {
        std::ifstream file{ geode::to_string( filename ), std::ios::binary };
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            file.good(), nullptr, geode::OpenGeodeException::TYPE::data,
            "[JSONGridInput] Error while opening file: ", filename );
        JSONGridSax sax{ max_records };
        const auto success = nlohmann::json::sax_parse( file, &sax );
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            success || sax.nb_records() == max_records, nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[JSONGridInput] Error while parsing file ", filename, " ",
            sax.error() );
        sax.finalize();
        return sax;
    }

    /*
     * Distinct values of a grid axis, regularly spaced
     */
    struct AxisDiscretization
    {
        double origin{ 0 };
        double cell_length{ 1 };
        geode::index_t nb_vertices{ 1 };
    };

    AxisDiscretization compute_axis_discretization(
        const Column& coordinates, geode::local_index_t axis )
    {
        auto values = coordinates;
        absl::c_sort( values );
        values.erase( std::unique( values.begin(), values.end(),
                          []( double lhs, double rhs ) {
                              return rhs - lhs <= geode::GLOBAL_EPSILON;
                          } ),
            values.end() );
        AxisDiscretization discretization;
        discretization.origin = values.front();
        discretization.nb_vertices =
            static_cast< geode::index_t >( values.size() );
        if( discretization.nb_vertices < 2 )
        {
            return discretization;
        }
        discretization.cell_length = ( values.back() - values.front() )
                                     / ( discretization.nb_vertices - 1 );
        for( const auto vertex : geode::Indices{ values } )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                std::fabs( values[vertex] - discretization.origin
                           - vertex * discretization.cell_length )
                    <= geode::GLOBAL_EPSILON,
                nullptr, geode::OpenGeodeException::TYPE::data,
                "[JSONGridInput] Coordinates along ", COORDINATE_KEYS[axis],
                " are not regularly spaced" );
        }
        return discretization;
    }

    class JSONGridInputImpl
    {
    public:
        explicit JSONGridInputImpl( std::string_view filename )
            : sax_( parse_json( filename ) )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                sax_.nb_records() > 0, nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[JSONGridInput] Missing or empty ", GRID_KEY, " array" );
            for( const auto axis : geode::LRange{ 3 } )
            {
                geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                    absl::c_none_of( sax_.coordinates()[axis],
                        []( double value ) {
                            return std::isnan( value );
                        } ),
                    nullptr, geode::OpenGeodeException::TYPE::data,
                    "[JSONGridInput] Missing numeric ", COORDINATE_KEYS[axis],
                    " coordinate in ", GRID_KEY, " records" );
                axes_[axis] = compute_axis_discretization(
                    sax_.coordinates()[axis], axis );
            }
        }

        geode::LightRegularGrid3D read_file()
        {
            geode::LightRegularGrid3D grid{ origin(), cells_number(),
                cells_length() };
            // TODO: coordinate_system (EPSG code) is not stored since
            // LightRegularGrid does not support coordinate reference systems
            set_grid_name( grid );
            read_attributes( grid, record_vertices( grid ) );
            return grid;
        }

    private:
        [[nodiscard]] geode::Point3D origin() const
        {
            return geode::Point3D{ { axes_[0].origin, axes_[1].origin,
                axes_[2].origin } };
        }

        [[nodiscard]] std::array< geode::index_t, 3 > cells_number() const
        {
            return { axes_[0].nb_vertices - 1, axes_[1].nb_vertices - 1,
                axes_[2].nb_vertices - 1 };
        }

        [[nodiscard]] std::array< double, 3 > cells_length() const
        {
            return { axes_[0].cell_length, axes_[1].cell_length,
                axes_[2].cell_length };
        }

        [[nodiscard]] geode::index_t vertex_index_along_axis(
            geode::index_t record_id, geode::local_index_t axis ) const
        {
            const auto& discretization = axes_[axis];
            return static_cast< geode::index_t >( std::lround(
                ( sax_.coordinates()[axis][record_id] - discretization.origin )
                / discretization.cell_length ) );
        }

        /*
         * Grid vertex of each record, checking that records and grid vertices
         * are in one-to-one correspondence
         */
        [[nodiscard]] std::vector< geode::index_t > record_vertices(
            const geode::LightRegularGrid3D& grid ) const
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                sax_.nb_records() == grid.nb_grid_vertices(), nullptr,
                geode::OpenGeodeException::TYPE::data,
                "[JSONGridInput] Number of records (", sax_.nb_records(),
                ") does not match the number of grid vertices (",
                grid.nb_grid_vertices(), ")" );
            std::vector< geode::index_t > vertices( sax_.nb_records() );
            std::vector< bool > visited( grid.nb_grid_vertices(), false );
            for( const auto record_id : geode::Range{ sax_.nb_records() } )
            {
                const auto vertex = grid.vertex_index(
                    { vertex_index_along_axis( record_id, 0 ),
                        vertex_index_along_axis( record_id, 1 ),
                        vertex_index_along_axis( record_id, 2 ) } );
                geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                    !visited[vertex], nullptr,
                    geode::OpenGeodeException::TYPE::data,
                    "[JSONGridInput] Record ", record_id,
                    " duplicates grid vertex ", vertex );
                visited[vertex] = true;
                vertices[record_id] = vertex;
            }
            return vertices;
        }

        [[nodiscard]] absl::flat_hash_set< std::string >
            bool_attribute_names() const
        {
            absl::flat_hash_set< std::string > names;
            for( const auto& key : BOOL_ATTRIBUTE_KEYS )
            {
                if( auto name = sax_.top_level_string( key ) )
                {
                    names.emplace( std::move( name.value() ) );
                }
            }
            return names;
        }

        void set_grid_name( geode::LightRegularGrid3D& grid ) const
        {
            if( auto name = sax_.top_level_string( NAME_KEY ) )
            {
                geode::IdentifierBuilder{ grid }.set_name(
                    std::move( name.value() ) );
            }
        }

        void read_attributes( geode::LightRegularGrid3D& grid,
            absl::Span< const geode::index_t > vertices ) const
        {
            const auto bool_names = bool_attribute_names();
            auto& manager = grid.grid_vertex_attribute_manager();
            for( const auto& [name, values] : sax_.attributes() )
            {
                if( bool_names.contains( name ) )
                {
                    geode::AttributeValues< bool > default_values;
                    default_values.default_value = false;
                    const auto attribute = create_attribute< bool >(
                        manager, name, std::move( default_values ) );
                    for( const auto record_id : geode::Indices{ values } )
                    {
                        const auto value = values[record_id];
                        if( !std::isnan( value ) )
                        {
                            attribute->set_value(
                                vertices[record_id], value != 0 );
                        }
                    }
                    continue;
                }
                geode::AttributeValues< double > default_values;
                default_values.default_value = NO_DATA;
                default_values.no_value = NO_DATA;
                const auto attribute = create_attribute< double >(
                    manager, name, std::move( default_values ) );
                for( const auto record_id : geode::Indices{ values } )
                {
                    attribute->set_value(
                        vertices[record_id], values[record_id] );
                }
            }
        }

        template < typename T >
        static std::shared_ptr< geode::VariableAttribute< T > >
            create_attribute( geode::AttributeManager& manager,
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

    private:
        JSONGridSax sax_;
        std::array< AxisDiscretization, 3 > axes_;
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
            const auto sax = parse_json( filename(), 1 );
            if( !sax.has_root_object() || sax.nb_records() == 0 )
            {
                return Percentage{ 0 };
            }
            for( const auto& coordinate : sax.coordinates() )
            {
                if( std::isnan( coordinate.front() ) )
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
