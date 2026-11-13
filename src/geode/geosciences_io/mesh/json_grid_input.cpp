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
#include <cstdint>
#include <fstream>
#include <limits>
#include <optional>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_set.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/identifier_builder.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/percentage.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>
#include <geode/geometry/vector.hpp>

#include <geode/mesh/core/light_regular_grid.hpp>

#include <geode/geosciences/explicit/geometry/geographic_coordinate_system_helper.hpp>

#include <geode/geosciences_io/mesh/internal/json_common.hpp>

/*
 * Expected file content:
 * {
 *   "grid_df": [ { "x": 2.5, "y": 2.5, "z": 10.0, "<attribute>": 1.0, ... },
 *                ... ],
 *   "name": "<grid name>",           (optional)
 *   "coordinate_system": "<code>",   (optional, e.g. "3949" or "EPSG:3949")
 *   "grid_topo_name": "<attribute>", (optional, attribute read as bool)
 *   "grid_site_name": "<attribute>", (optional, attribute read as bool)
 *   ...                              (other keys are ignored)
 * }
 * Each record of grid_df is a grid vertex, records can be in any order.
 * String values of the records are ignored.
 */
namespace
{
    constexpr auto GRID_KEY = "grid_df";
    constexpr auto GRID_LOG_PREFIX = "[JSONGridInput]";
    constexpr auto NAME_KEY = "name";
    constexpr std::array< std::string_view, 2 > BOOL_ATTRIBUTE_KEYS{
        "grid_topo_name", "grid_site_name"
    };
    constexpr auto& COORDINATE_KEYS = geode::internal::JSON_COORDINATE_KEYS;

    using Column = geode::internal::JSONColumn;

    geode::internal::JSONRecords parse_json( std::string_view filename,
        geode::index_t max_records =
            std::numeric_limits< geode::index_t >::max() )
    {
        return geode::internal::parse_json_records(
            filename, GRID_KEY, GRID_LOG_PREFIX, max_records );
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
                nullptr, geode::OpenGeodeException::TYPE::data, GRID_LOG_PREFIX,
                " Coordinates along ", COORDINATE_KEYS[axis],
                " are not regularly spaced" );
        }
        return discretization;
    }

    class JSONGridInputImpl
    {
    public:
        explicit JSONGridInputImpl( std::string_view filename )
            : records_( parse_json( filename ) )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                records_.nb_records > 0, nullptr,
                geode::OpenGeodeException::TYPE::data, GRID_LOG_PREFIX,
                " Missing or empty ", GRID_KEY, " array" );
            geode::internal::warn_ignored_json_strings(
                records_, GRID_KEY, GRID_LOG_PREFIX, "" );
            for( const auto axis : geode::LRange{ 3 } )
            {
                geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                    absl::c_none_of( records_.coordinates[axis],
                        []( double value ) {
                            return std::isnan( value );
                        } ),
                    nullptr, geode::OpenGeodeException::TYPE::data,
                    GRID_LOG_PREFIX, " Missing numeric ", COORDINATE_KEYS[axis],
                    " coordinate in ", GRID_KEY, " records" );
                axes_[axis] = compute_axis_discretization(
                    records_.coordinates[axis], axis );
            }
        }

        geode::LightRegularGrid3D read_file()
        {
            check_number_of_vertices();
            geode::LightRegularGrid3D grid{ origin(), cells_number(),
                cells_length() };
            set_grid_name( grid );
            set_grid_coordinate_system( grid );
            read_attributes( grid, record_vertices( grid ) );
            return grid;
        }

    private:
        void set_grid_coordinate_system( geode::LightRegularGrid3D& grid ) const
        {
            auto info = geode::internal::read_json_coordinate_system(
                records_, GRID_LOG_PREFIX );
            if( !info )
            {
                return;
            }
            const auto crs_name = info->authority_code();
            geode::assign_light_regular_grid_geographic_coordinate_system_info(
                grid, crs_name, std::move( info.value() ) );
        }

        void check_number_of_vertices() const
        {
            std::uint64_t nb_vertices{ 1 };
            for( const auto& axis : axes_ )
            {
                // Each factor is at most the number of records, stopping
                // as soon as the product exceeds it prevents any overflow
                nb_vertices *= axis.nb_vertices;
                if( nb_vertices > records_.nb_records )
                {
                    break;
                }
            }
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                nb_vertices == records_.nb_records, nullptr,
                geode::OpenGeodeException::TYPE::data, GRID_LOG_PREFIX,
                " Number of records (", records_.nb_records,
                ") does not match the number of grid vertices (",
                axes_[0].nb_vertices, " x ", axes_[1].nb_vertices, " x ",
                axes_[2].nb_vertices, ")" );
        }

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
            return static_cast< geode::index_t >(
                std::lround( ( records_.coordinates[axis][record_id]
                                 - discretization.origin )
                             / discretization.cell_length ) );
        }

        /*
         * Grid vertex of each record, checking that records and grid vertices
         * are in one-to-one correspondence
         */
        [[nodiscard]] std::vector< geode::index_t > record_vertices(
            const geode::LightRegularGrid3D& grid ) const
        {
            std::vector< geode::index_t > vertices( records_.nb_records );
            std::vector< bool > visited( grid.nb_grid_vertices(), false );
            for( const auto record_id : geode::Range{ records_.nb_records } )
            {
                const auto vertex = grid.vertex_index(
                    { vertex_index_along_axis( record_id, 0 ),
                        vertex_index_along_axis( record_id, 1 ),
                        vertex_index_along_axis( record_id, 2 ) } );
                geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                    !visited[vertex], nullptr,
                    geode::OpenGeodeException::TYPE::data, GRID_LOG_PREFIX,
                    " Record ", record_id, " duplicates grid vertex ", vertex );
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
                if( auto name = records_.top_level_string( key ) )
                {
                    names.emplace( std::move( name.value() ) );
                }
            }
            return names;
        }

        void set_grid_name( geode::LightRegularGrid3D& grid ) const
        {
            if( auto name = records_.top_level_string( NAME_KEY ) )
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
            for( const auto& [name, values] : records_.attributes )
            {
                if( bool_names.contains( name ) )
                {
                    geode::AttributeValues< bool > default_values;
                    default_values.default_value = false;
                    const auto attribute =
                        geode::internal::create_json_attribute< bool >(
                            manager, name, std::move( default_values ) );
                    for( const auto& [record_id, value] : values )
                    {
                        attribute->set_value( vertices[record_id], value != 0 );
                    }
                    continue;
                }
                const auto attribute =
                    geode::internal::create_json_double_attribute(
                        manager, name );
                for( const auto& [record_id, value] : values )
                {
                    attribute->set_value( vertices[record_id], value );
                }
            }
        }

    private:
        geode::internal::JSONRecords records_;
        std::array< AxisDiscretization, 3 > axes_;
    };
} // namespace

namespace geode::internal
{
    std::vector< LightRegularGrid3D > JSONGridInput::read()
    {
        JSONGridInputImpl reader{ filename() };
        std::vector< LightRegularGrid3D > objects;
        objects.emplace_back( reader.read_file() );
        return objects;
    }

    Percentage JSONGridInput::is_loadable() const
    {
        try
        {
            const auto sax = parse_json( filename(), 1 );
            if( !sax.has_root_object || sax.nb_records == 0 )
            {
                return Percentage{ 0 };
            }
            for( const auto& coordinate : sax.coordinates )
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
