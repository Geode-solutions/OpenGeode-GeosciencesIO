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

#include <geode/geosciences_io/mesh/internal/json_edged_curve_input.hpp>

#include <cmath>
#include <limits>
#include <optional>

#include <absl/algorithm/container.h>
#include <absl/container/flat_hash_map.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/percentage.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geometry/point.hpp>

#include <geode/mesh/builder/edged_curve_builder.hpp>
#include <geode/mesh/core/edged_curve.hpp>

#include <geode/geosciences/explicit/geometry/geographic_coordinate_system_helper.hpp>

#include <geode/geosciences_io/mesh/internal/json_common.hpp>

/*
 * Expected file content:
 * {
 *   "drillholes_df": [ { "drillhole": "<name>", "x": 59.1, "y": 52.5,
 *                        "z": 18.0, "<attribute>": 1.0, ... },
 *                      ... ],
 *   "drillhole_col_name": "<key>",   (optional, key of the records holding
 *                                     the drillhole name, "drillhole" by
 *                                     default)
 *   "coordinate_system": "<code>",   (optional, e.g. "3949" or "EPSG:3949")
 *   ...                              (other keys are ignored)
 * }
 * Each record of drillholes_df is a vertex of the drillhole named by the
 * value of its drillhole key. Each drillhole is an EdgedCurve3D named after
 * the drillhole, its vertices are the records of the drillhole in file order
 * and its edges link consecutive vertices. The numeric values of the records
 * are vertex attributes, the string values are ignored. Records without
 * drillhole or without numeric x, y and z are skipped.
 */
namespace
{
    constexpr auto DRILLHOLES_KEY = "drillholes_df";
    constexpr auto DRILLHOLE_COLUMN_KEY = "drillhole_col_name";
    constexpr auto DEFAULT_DRILLHOLE_COLUMN = "drillhole";
    constexpr auto CURVE_LOG_PREFIX = "[JSONEdgedCurveInput]";

    geode::internal::JSONRecords parse_drillholes( std::string_view filename,
        geode::index_t max_records =
            std::numeric_limits< geode::index_t >::max() )
    {
        return geode::internal::parse_json_records(
            filename, DRILLHOLES_KEY, CURVE_LOG_PREFIX, max_records );
    }

    /*
     * Records of one drillhole, in file order
     */
    struct Drillhole
    {
        std::string name;
        std::vector< geode::index_t > records;
    };

    class JSONEdgedCurveInputImpl
    {
    public:
        explicit JSONEdgedCurveInputImpl( std::string_view filename )
            : records_( parse_drillholes( filename ) ),
              drillhole_column_(
                  records_.top_level_string( DRILLHOLE_COLUMN_KEY )
                      .value_or( DEFAULT_DRILLHOLE_COLUMN ) )
        {
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                records_.nb_records > 0, nullptr,
                geode::OpenGeodeException::TYPE::data, CURVE_LOG_PREFIX,
                " Missing or empty ", DRILLHOLES_KEY, " array" );
            geode::internal::warn_ignored_json_strings(
                records_, DRILLHOLES_KEY, CURVE_LOG_PREFIX, drillhole_column_ );
            group_records();
            geode::OpenGeodeGeosciencesIOMeshException::check_exception(
                !drillholes_.empty(), nullptr,
                geode::OpenGeodeException::TYPE::data, CURVE_LOG_PREFIX,
                " No valid record with a drillhole and numeric x, y and z in ",
                DRILLHOLES_KEY );
        }

        std::vector< std::unique_ptr< geode::EdgedCurve3D > > read_curves(
            const geode::MeshImpl& impl ) const
        {
            const auto info = geode::internal::read_json_coordinate_system(
                records_, CURVE_LOG_PREFIX );
            std::vector< std::unique_ptr< geode::EdgedCurve3D > > curves;
            // Vertex of each record in the curve of its drillhole
            std::vector< geode::index_t > record_vertices(
                records_.nb_records, geode::NO_ID );
            for( const auto& drillhole : drillholes_ )
            {
                curves.push_back( create_curve( impl, drillhole, info ) );
                for( const auto vertex : geode::Indices{ drillhole.records } )
                {
                    record_vertices[drillhole.records[vertex]] = vertex;
                }
            }
            read_attributes( curves, record_vertices );
            return curves;
        }

    private:
        [[nodiscard]] bool has_numeric_coordinates(
            geode::index_t record_id ) const
        {
            return absl::c_none_of( records_.coordinates,
                [record_id]( const geode::internal::JSONColumn& coordinate ) {
                    return std::isnan( coordinate[record_id] );
                } );
        }

        [[nodiscard]] geode::Point3D point( geode::index_t record_id ) const
        {
            return geode::Point3D{ { records_.coordinates[0][record_id],
                records_.coordinates[1][record_id],
                records_.coordinates[2][record_id] } };
        }

        void group_records()
        {
            const auto column = absl::c_find_if(
                records_.strings, [this]( const auto& string_column ) {
                    return string_column.first == drillhole_column_;
                } );
            absl::flat_hash_map< std::string, geode::index_t > drillhole_ids;
            geode::index_t nb_without_drillhole{ records_.nb_records };
            geode::index_t nb_without_coordinates{ 0 };
            if( column != records_.strings.end() )
            {
                nb_without_drillhole -=
                    static_cast< geode::index_t >( column->second.size() );
                for( const auto& [record_id, name] : column->second )
                {
                    if( !has_numeric_coordinates( record_id ) )
                    {
                        nb_without_coordinates++;
                        continue;
                    }
                    const auto [it, inserted] = drillhole_ids.try_emplace( name,
                        static_cast< geode::index_t >( drillholes_.size() ) );
                    if( inserted )
                    {
                        drillholes_.push_back( { name, {} } );
                    }
                    drillholes_[it->second].records.push_back( record_id );
                }
            }
            if( nb_without_drillhole > 0 )
            {
                geode::Logger::warning( CURVE_LOG_PREFIX, " ",
                    nb_without_drillhole, " records of ", DRILLHOLES_KEY,
                    " without string value for key \"", drillhole_column_,
                    "\" are skipped" );
            }
            if( nb_without_coordinates > 0 )
            {
                geode::Logger::warning( CURVE_LOG_PREFIX, " ",
                    nb_without_coordinates, " records of ", DRILLHOLES_KEY,
                    " without numeric x, y and z are skipped" );
            }
        }

        [[nodiscard]] std::unique_ptr< geode::EdgedCurve3D > create_curve(
            const geode::MeshImpl& impl,
            const Drillhole& drillhole,
            const std::optional< geode::GeographicCoordinateSystemInfo >& info )
            const
        {
            auto curve = geode::EdgedCurve3D::create( impl );
            auto builder = geode::EdgedCurveBuilder3D::create( *curve );
            builder->set_name( drillhole.name );
            const auto nb_vertices =
                static_cast< geode::index_t >( drillhole.records.size() );
            builder->create_vertices( nb_vertices );
            for( const auto vertex : geode::Range{ nb_vertices } )
            {
                builder->set_point(
                    vertex, point( drillhole.records[vertex] ) );
                if( vertex > 0 )
                {
                    builder->create_edge( vertex - 1, vertex );
                }
            }
            if( info )
            {
                geode::assign_edged_curve_geographic_coordinate_system_info(
                    *curve, *builder, info->authority_code(), info.value() );
            }
            return curve;
        }

        void read_attributes(
            const std::vector< std::unique_ptr< geode::EdgedCurve3D > >& curves,
            absl::Span< const geode::index_t > record_vertices ) const
        {
            // Drillhole of each record
            std::vector< geode::index_t > record_curves(
                records_.nb_records, geode::NO_ID );
            for( const auto curve : geode::Indices{ drillholes_ } )
            {
                for( const auto record_id : drillholes_[curve].records )
                {
                    record_curves[record_id] = curve;
                }
            }
            for( const auto& [name, values] : records_.attributes )
            {
                if( name == drillhole_column_ )
                {
                    continue;
                }
                std::vector<
                    std::shared_ptr< geode::VariableAttribute< double > > >
                    attributes( curves.size() );
                for( const auto& [record_id, value] : values )
                {
                    const auto curve = record_curves[record_id];
                    if( curve == geode::NO_ID )
                    {
                        continue;
                    }
                    if( !attributes[curve] )
                    {
                        attributes[curve] =
                            geode::internal::create_json_double_attribute(
                                curves[curve]->vertex_attribute_manager(),
                                name );
                    }
                    attributes[curve]->set_value(
                        record_vertices[record_id], value );
                }
            }
        }

    private:
        geode::internal::JSONRecords records_;
        std::string drillhole_column_;
        std::vector< Drillhole > drillholes_;
    };
} // namespace

namespace geode::internal
{
    std::vector< std::unique_ptr< EdgedCurve3D > > JSONEdgedCurveInput::read(
        const MeshImpl& impl )
    {
        const JSONEdgedCurveInputImpl reader{ filename() };
        return reader.read_curves( impl );
    }

    Percentage JSONEdgedCurveInput::is_loadable() const
    {
        try
        {
            const auto records = parse_drillholes( filename(), 1 );
            if( !records.has_root_object || records.nb_records == 0 )
            {
                return Percentage{ 0 };
            }
            for( const auto& coordinate : records.coordinates )
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
