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

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <absl/container/flat_hash_map.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/types.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/geosciences/explicit/geometry/geographic_coordinate_system.hpp>

namespace geode::internal
{
    /*!
     * Value of the numeric attributes of the records that do not have one
     */
    inline constexpr auto JSON_NO_DATA =
        std::numeric_limits< double >::quiet_NaN();

    inline constexpr std::array< std::string_view, 3 > JSON_COORDINATE_KEYS{
        "x", "y", "z"
    };

    /*!
     * Values of one key across all records, NaN where the record has no
     * numeric value for this key
     */
    using JSONColumn = std::vector< double >;

    /*!
     * (record, value) pairs of one key, only for the records having a
     * numeric value for this key
     */
    using JSONSparseColumn = std::vector< std::pair< index_t, double > >;

    /*!
     * (record, value) pairs of one key, only for the records having a
     * string value for this key
     */
    using JSONStringColumn = std::vector< std::pair< index_t, std::string > >;

    /*!
     * Records of one array of a JSON file, stored column by column so that
     * the file is never held in memory as a JSON document.
     * Records are the objects of the array, the "x", "y" and "z" numeric
     * values are the coordinates, the other numeric values are attributes and
     * the string values are kept apart.
     */
    struct JSONRecords
    {
        [[nodiscard]] std::optional< std::string > top_level_string(
            std::string_view key ) const;

        index_t nb_records{ 0 };
        bool has_root_object{ false };
        /*!
         * One value per record, padded with NaN
         */
        std::array< JSONColumn, 3 > coordinates;
        /*!
         * Numeric values of the keys that are not coordinates
         */
        std::vector< std::pair< std::string, JSONSparseColumn > > attributes;
        /*!
         * String values of the record keys
         */
        std::vector< std::pair< std::string, JSONStringColumn > > strings;
        /*!
         * String values of the keys of the root object
         */
        absl::flat_hash_map< std::string, std::string > top_level_strings;
    };

    /*!
     * Parse a JSON file made of a root object containing an array of records.
     * The file is parsed as a stream, only the records of the array named
     * array_key are stored.
     * @param[in] array_key Key of the root object holding the records.
     * @param[in] log_prefix Prefix of the messages, such as "[JSONGridInput]".
     * @param[in] max_records Parsing stops once max_records records have been
     * read.
     * @exception OpenGeodeException if the file cannot be opened or parsed,
     * unless the error occurs after the array of records or after the
     * max_records first records.
     */
    [[nodiscard]] JSONRecords parse_json_records( std::string_view filename,
        std::string_view array_key,
        std::string_view log_prefix,
        index_t max_records );

    /*!
     * Emit one warning per key having string values in the records, except
     * for kept_key.
     */
    void warn_ignored_json_strings( const JSONRecords& records,
        std::string_view array_key,
        std::string_view log_prefix,
        std::string_view kept_key );

    /*!
     * Read the coordinate system of the root object ("coordinate_system" key,
     * authority code such as "3949" or "EPSG:3949").
     * @return the coordinate system, or nothing if there is none or if it is
     * not valid (a warning is emitted in this case).
     */
    [[nodiscard]] std::optional< GeographicCoordinateSystemInfo >
        read_json_coordinate_system(
            const JSONRecords& records, std::string_view log_prefix );

    template < typename T >
    [[nodiscard]] std::shared_ptr< VariableAttribute< T > >
        create_json_attribute( AttributeManager& manager,
            std::string_view name,
            AttributeValues< T > values )
    {
        AttributeProperties properties;
        properties.assignable = false;
        properties.interpolable = false;
        properties.transferable = true;
        const auto attribute_id =
            manager.create_attribute< VariableAttribute, T >(
                name, std::move( values ), std::move( properties ) );
        return manager.find_attribute< VariableAttribute, T >( attribute_id );
    }

    /*!
     * Create a double attribute with NaN as default and no value
     */
    [[nodiscard]] inline std::shared_ptr< VariableAttribute< double > >
        create_json_double_attribute(
            AttributeManager& manager, std::string_view name )
    {
        AttributeValues< double > values;
        values.default_value = JSON_NO_DATA;
        values.no_value = JSON_NO_DATA;
        return create_json_attribute< double >(
            manager, name, std::move( values ) );
    }
} // namespace geode::internal
