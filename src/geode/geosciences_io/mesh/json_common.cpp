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

#include <geode/geosciences_io/mesh/internal/json_common.hpp>

#include <fstream>

#include <nlohmann/json.hpp>

#include <absl/container/flat_hash_set.h>

#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>

#include <geode/geosciences/explicit/geometry/geographic_coordinate_system.hpp>

#include <geode/geosciences_io/mesh/common.hpp>

namespace
{
    constexpr auto COORDINATE_SYSTEM_KEY = "coordinate_system";

    /*
     * SAX handler storing the records of the array named array_key
     */
    class JSONRecordsSax
    {
    public:
        using number_integer_t = nlohmann::json::number_integer_t;
        using number_unsigned_t = nlohmann::json::number_unsigned_t;
        using number_float_t = nlohmann::json::number_float_t;
        using string_t = nlohmann::json::string_t;
        using binary_t = nlohmann::json::binary_t;

        JSONRecordsSax( geode::internal::JSONRecords& records,
            std::string_view array_key,
            geode::index_t max_records )
            : records_( records ),
              array_key_( array_key ),
              max_records_( max_records )
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
                records_.top_level_strings[current_key_] = std::move( value );
            }
            else if( in_record() )
            {
                const auto [it, inserted] = string_ids_.try_emplace(
                    current_key_,
                    static_cast< geode::index_t >( records_.strings.size() ) );
                if( inserted )
                {
                    records_.strings.emplace_back(
                        current_key_, geode::internal::JSONStringColumn{} );
                }
                records_.strings[it->second].second.emplace_back(
                    records_.nb_records, std::move( value ) );
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
                records_.has_root_object = true;
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
                records_.nb_records++;
                return records_.nb_records < max_records_;
            }
            return true;
        }

        bool start_array( std::size_t /*unused*/ )
        {
            depth_++;
            if( depth_ == 2 && current_key_ == array_key_ )
            {
                in_array_ = true;
            }
            return true;
        }

        bool end_array()
        {
            if( depth_ == 2 )
            {
                array_complete_ = array_complete_ || in_array_;
                in_array_ = false;
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

        [[nodiscard]] bool array_complete() const
        {
            return array_complete_;
        }

        [[nodiscard]] const std::string& error() const
        {
            return error_;
        }

    private:
        [[nodiscard]] bool in_record() const
        {
            return in_array_ && depth_ == 3;
        }

        bool number( double value )
        {
            if( !in_record() )
            {
                return true;
            }
            for( const auto axis : geode::LRange{ 3 } )
            {
                if( current_key_
                    == geode::internal::JSON_COORDINATE_KEYS[axis] )
                {
                    auto& values = records_.coordinates[axis];
                    values.resize( records_.nb_records + 1,
                        geode::internal::JSON_NO_DATA );
                    values.back() = value;
                    return true;
                }
            }
            const auto [it, inserted] = attribute_ids_.try_emplace(
                current_key_,
                static_cast< geode::index_t >( records_.attributes.size() ) );
            if( inserted )
            {
                records_.attributes.emplace_back(
                    current_key_, geode::internal::JSONSparseColumn{} );
            }
            records_.attributes[it->second].second.emplace_back(
                records_.nb_records, value );
            return true;
        }

    private:
        geode::internal::JSONRecords& records_;
        std::string_view array_key_;
        geode::index_t max_records_;
        geode::index_t depth_{ 0 };
        bool in_array_{ false };
        bool array_complete_{ false };
        std::string current_key_;
        absl::flat_hash_map< std::string, geode::index_t > attribute_ids_;
        absl::flat_hash_map< std::string, geode::index_t > string_ids_;
        std::string error_;
    };
} // namespace

namespace geode::internal
{
    std::optional< std::string > JSONRecords::top_level_string(
        std::string_view key ) const
    {
        const auto it = top_level_strings.find( key );
        if( it == top_level_strings.end() )
        {
            return std::nullopt;
        }
        return it->second;
    }

    JSONRecords parse_json_records( std::string_view filename,
        std::string_view array_key,
        std::string_view log_prefix,
        index_t max_records )
    {
        std::ifstream file{ to_string( filename ), std::ios::binary };
        OpenGeodeGeosciencesIOMeshException::check_exception( file.good(),
            nullptr, OpenGeodeException::TYPE::data, log_prefix,
            " Error while opening file: ", filename );
        JSONRecords records;
        JSONRecordsSax sax{ records, array_key, max_records };
        const auto success = nlohmann::json::sax_parse( file, &sax );
        if( !success && records.nb_records != max_records
            && sax.array_complete() )
        {
            Logger::warning( log_prefix, " Error while parsing file ", filename,
                " after the ", array_key, " array, ", sax.error(),
                ": the rest of the file is ignored" );
        }
        else
        {
            OpenGeodeGeosciencesIOMeshException::check_exception(
                success || records.nb_records == max_records, nullptr,
                OpenGeodeException::TYPE::data, log_prefix,
                " Error while parsing file ", filename, " ", sax.error() );
        }
        for( auto& coordinate : records.coordinates )
        {
            coordinate.resize( records.nb_records, JSON_NO_DATA );
        }
        return records;
    }

    void warn_ignored_json_strings( const JSONRecords& records,
        std::string_view array_key,
        std::string_view log_prefix,
        std::string_view kept_key )
    {
        for( const auto& string_column : records.strings )
        {
            if( string_column.first != kept_key )
            {
                Logger::warning( log_prefix, " String values of key \"",
                    string_column.first, "\" in ", array_key,
                    " records are not supported and are ignored" );
            }
        }
    }

    std::optional< GeographicCoordinateSystemInfo > read_json_coordinate_system(
        const JSONRecords& records, std::string_view log_prefix )
    {
        const auto code = records.top_level_string( COORDINATE_SYSTEM_KEY );
        if( !code )
        {
            return std::nullopt;
        }
        try
        {
            return GeographicCoordinateSystemInfo::from_authority_code(
                code.value() );
        }
        catch( const OpenGeodeException& exception )
        {
            Logger::warning( log_prefix, " Coordinate system ", code.value(),
                " is ignored: ", exception.what() );
            return std::nullopt;
        }
    }
} // namespace geode::internal
