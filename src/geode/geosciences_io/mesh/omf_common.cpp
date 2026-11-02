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

#include <geode/geosciences_io/mesh/internal/omf_common.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <limits>
#include <typeinfo>

#include <zlib.h>

#include <absl/algorithm/container.h>
#include <absl/strings/str_cat.h>

#include <geode/basic/attribute_manager.hpp>
#include <geode/basic/file.hpp>
#include <geode/basic/logger.hpp>
#include <geode/basic/range.hpp>
#include <geode/basic/variable_attribute.hpp>

#include <geode/image/core/rgb_color.hpp>

namespace
{
    constexpr std::array< unsigned char, 4 > MAGIC_NUMBER{ 0x84, 0x83, 0x82,
        0x81 };
    constexpr std::string_view VERSION_PREFIX{ "OMF-v" };
    constexpr std::streamoff VERSION_OFFSET{ 4 };
    constexpr std::streamoff JSON_START_OFFSET{ 52 };
    constexpr std::streamoff HEADER_SIZE{ 60 };
    constexpr std::size_t INFLATE_CHUNK_SIZE{ 1 << 16 };
    constexpr std::string_view CLASS_KEY{ "__class__" };
    constexpr std::string_view ARRAY_KEY{ "array" };

    struct OMFHeader
    {
        uint64_t json_start{ 0 };
        uint64_t file_size{ 0 };
    };

    std::optional< OMFHeader > read_omf_header( std::ifstream& file )
    {
        std::array< char, HEADER_SIZE > header;
        file.seekg( 0, std::ios::end );
        const auto file_size = static_cast< uint64_t >( file.tellg() );
        file.seekg( 0, std::ios::beg );
        if( !file.read( header.data(), HEADER_SIZE ) )
        {
            file.clear();
            return std::nullopt;
        }
        if( std::memcmp(
                header.data(), MAGIC_NUMBER.data(), MAGIC_NUMBER.size() )
            != 0 )
        {
            return std::nullopt;
        }
        if( std::string_view{
                header.data() + VERSION_OFFSET, VERSION_PREFIX.size() }
            != VERSION_PREFIX )
        {
            return std::nullopt;
        }
        uint64_t json_start{ 0 };
        for( const auto byte : geode::Range{ 8 } )
        {
            json_start |= static_cast< uint64_t >( static_cast< unsigned char >(
                              header[JSON_START_OFFSET + byte] ) )
                          << ( 8 * byte );
        }
        return OMFHeader{ json_start, file_size };
    }

    std::string_view omf_class_name( const nlohmann::json& object )
    {
        const auto it = object.find( CLASS_KEY );
        if( it == object.end() || !it->is_string() )
        {
            return {};
        }
        return it->get_ref< const std::string& >();
    }

    /*
     * Value of a required key, with a clear error if it is missing
     */
    const nlohmann::json& omf_member(
        const nlohmann::json& object, std::string_view key )
    {
        const auto it = object.find( key );
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            it != object.end(), nullptr, geode::OpenGeodeException::TYPE::data,
            "[OMFInput] Missing key \"", key, "\" in ",
            omf_class_name( object ), " object" );
        return *it;
    }

    geode::Vector3D omf_vector( const nlohmann::json& object,
        std::string_view key,
        const geode::Vector3D& default_value )
    {
        const auto it = object.find( key );
        if( it == object.end() || !it->is_array() || it->size() != 3 )
        {
            return default_value;
        }
        geode::Vector3D result;
        for( const auto axis : geode::LRange{ 3 } )
        {
            const auto& value = ( *it )[axis];
            if( !value.is_number() )
            {
                return default_value;
            }
            result.set_value( axis, value.get< double >() );
        }
        return result;
    }

    /*
     * Number of values usable for a data: the minimum of the available and
     * expected numbers, with a warning when they differ
     */
    geode::index_t omf_usable_size( std::string_view data_name,
        std::size_t nb_available,
        geode::index_t nb_expected )
    {
        if( nb_available != nb_expected )
        {
            geode::Logger::warning( "[OMFInput] Data \"", data_name, "\" has ",
                nb_available, " values instead of ", nb_expected,
                ": missing values are left undefined and extra values are "
                "ignored" );
        }
        return static_cast< geode::index_t >(
            std::min< std::size_t >( nb_available, nb_expected ) );
    }

    geode::AttributeProperties omf_attribute_properties( bool interpolable )
    {
        geode::AttributeProperties properties;
        properties.assignable = false;
        properties.interpolable = interpolable;
        properties.transferable = true;
        return properties;
    }

    /*
     * Find the attribute with this name and type, or create it. If an
     * attribute with this name but another type exists (e.g. two merged
     * elements with data of the same name but different kinds), the type
     * suffix is appended to the name.
     */
    template < typename T >
    std::shared_ptr< geode::VariableAttribute< T > >
        find_or_create_omf_attribute( geode::AttributeManager& manager,
            std::string_view name,
            std::string_view type_suffix,
            geode::AttributeValues< T > values,
            bool interpolable )
    {
        auto candidate = geode::to_string( name );
        while( true )
        {
            const auto ids = manager.attribute_ids_matching_name( candidate );
            if( !ids )
            {
                const auto id =
                    manager.create_attribute< geode::VariableAttribute, T >(
                        candidate, std::move( values ),
                        omf_attribute_properties( interpolable ) );
                return manager.find_attribute< geode::VariableAttribute, T >(
                    id );
            }
            for( const auto& id : ids.value() )
            {
                if( manager.attribute_type( id ) == typeid( T ).name() )
                {
                    return manager
                        .find_attribute< geode::VariableAttribute, T >( id );
                }
            }
            absl::StrAppend( &candidate, "_", type_suffix );
        }
    }

    geode::AttributeValues< double > omf_no_data_double()
    {
        geode::AttributeValues< double > values;
        values.default_value = std::numeric_limits< double >::quiet_NaN();
        values.no_value = values.default_value;
        return values;
    }

    template < geode::index_t dimension >
    geode::AttributeValues< geode::Vector< dimension > > omf_no_data_vector()
    {
        geode::AttributeValues< geode::Vector< dimension > > values;
        for( const auto axis : geode::LRange{ dimension } )
        {
            values.default_value.set_value(
                axis, std::numeric_limits< double >::quiet_NaN() );
        }
        values.no_value = values.default_value;
        return values;
    }

    /*
     * Black is the OpenGeode no value color (see
     * convert_raster_image_into_grid), it is also used for invalid colors
     */
    geode::AttributeValues< geode::RGBColor > omf_no_data_color()
    {
        geode::AttributeValues< geode::RGBColor > values;
        values.no_value = values.default_value;
        return values;
    }

    geode::AttributeValues< std::string > omf_no_data_string()
    {
        geode::AttributeValues< std::string > values;
        values.no_value = values.default_value;
        return values;
    }

    std::vector< std::string > read_omf_strings( const nlohmann::json& array )
    {
        const auto it = array.find( ARRAY_KEY );
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            it != array.end() && it->is_array(), nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[OMFInput] Expected an inline list of strings" );
        std::vector< std::string > strings;
        strings.reserve( it->size() );
        for( const auto& value : *it )
        {
            strings.emplace_back(
                value.is_string() ? value.get< std::string >() : "" );
        }
        return strings;
    }

    /*
     * Integer value of a real, -1 if it is not finite or out of the int64_t
     * range
     */
    int64_t omf_real_to_integer( double real )
    {
        // 2^63, exactly representable as a double
        constexpr auto LIMIT =
            -static_cast< double >( std::numeric_limits< int64_t >::min() );
        if( !std::isfinite( real ) || real < -LIMIT || real >= LIMIT )
        {
            return -1;
        }
        return static_cast< int64_t >( real );
    }

    /*
     * Integer value of a JSON number, -1 if it is not a number or out of the
     * int64_t range
     */
    int64_t omf_json_integer( const nlohmann::json& value )
    {
        if( value.is_number_unsigned() )
        {
            const auto integer = value.get< uint64_t >();
            return integer > static_cast< uint64_t >(
                       std::numeric_limits< int64_t >::max() )
                       ? -1
                       : static_cast< int64_t >( integer );
        }
        if( value.is_number_integer() )
        {
            return value.get< int64_t >();
        }
        if( value.is_number_float() )
        {
            return omf_real_to_integer( value.get< double >() );
        }
        return -1;
    }

    uint64_t omf_json_offset(
        const nlohmann::json& binary, std::string_view key )
    {
        const auto it = binary.find( key );
        if( it == binary.end() )
        {
            return 0;
        }
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            it->is_number_unsigned(), nullptr,
            geode::OpenGeodeException::TYPE::data, "[OMFInput] Array block ",
            key, " ", it->dump(), " is not a non negative integer" );
        return it->get< uint64_t >();
    }

    geode::RGBColor omf_color( int64_t red, int64_t green, int64_t blue )
    {
        const auto clamp = []( int64_t value ) {
            return static_cast< geode::local_index_t >(
                std::clamp< int64_t >( value, 0, 255 ) );
        };
        return { clamp( red ), clamp( green ), clamp( blue ) };
    }

    std::vector< geode::RGBColor > read_omf_inline_colors(
        const nlohmann::json& array )
    {
        const auto it = array.find( ARRAY_KEY );
        geode::OpenGeodeGeosciencesIOMeshException::check_exception(
            it != array.end() && it->is_array(), nullptr,
            geode::OpenGeodeException::TYPE::data,
            "[OMFInput] Expected an inline list of colors" );
        std::vector< geode::RGBColor > colors;
        colors.reserve( it->size() );
        for( const auto& value : *it )
        {
            if( value.is_array() && value.size() == 3
                && absl::c_all_of(
                    value, []( const nlohmann::json& component ) {
                        return component.is_number();
                    } ) )
            {
                colors.push_back( omf_color( omf_json_integer( value[0] ),
                    omf_json_integer( value[1] ),
                    omf_json_integer( value[2] ) ) );
            }
            else
            {
                colors.emplace_back();
            }
        }
        return colors;
    }

    class OMFDataReader
    {
    public:
        OMFDataReader( const geode::internal::OMFFile& file,
            const nlohmann::json& data,
            geode::index_t nb_values,
            geode::AttributeManager& manager,
            absl::FunctionRef< geode::index_t( geode::index_t ) >
                attribute_element )
            : file_( file ),
              data_( data ),
              name_( geode::internal::OMFFile::name( data ) ),
              nb_values_( nb_values ),
              manager_( manager ),
              attribute_element_( attribute_element )
        {
        }

        void read()
        {
            const auto data_class = omf_class_name( data_ );
            if( data_class == "ScalarData" )
            {
                read_scalar_data();
            }
            else if( data_class == "Vector2Data" )
            {
                read_vector_data< 2 >();
            }
            else if( data_class == "Vector3Data" )
            {
                read_vector_data< 3 >();
            }
            else if( data_class == "ColorData" )
            {
                read_color_data();
            }
            else if( data_class == "StringData"
                     || data_class == "DateTimeData" )
            {
                set_values( name_, "string", omf_no_data_string(), false,
                    read_omf_strings( array() ) );
            }
            else if( data_class == "MappedData" )
            {
                read_mapped_data();
            }
            else
            {
                geode::Logger::warning( "[OMFInput] Data \"", name_,
                    "\" of unsupported class ", data_class, " is skipped" );
            }
        }

    private:
        const nlohmann::json& array() const
        {
            return file_.object( omf_member( data_, ARRAY_KEY ) );
        }

        template < typename T >
        void set_values( std::string_view name,
            std::string_view type_suffix,
            geode::AttributeValues< T > default_values,
            bool interpolable,
            const std::vector< T >& values )
        {
            const auto attribute = find_or_create_omf_attribute< T >( manager_,
                name, type_suffix, std::move( default_values ), interpolable );
            const auto nb = omf_usable_size( name, values.size(), nb_values_ );
            for( const auto value : geode::Range{ nb } )
            {
                const auto element = attribute_element_( value );
                if( element != geode::NO_ID )
                {
                    attribute->set_value( element, values[value] );
                }
            }
        }

        void read_scalar_data()
        {
            set_values( name_, "scalar", omf_no_data_double(), true,
                file_.read_scalars( array(), nb_values_ ) );
        }

        template < geode::index_t dimension >
        void read_vector_data()
        {
            const auto components = file_.read_scalars(
                array(), uint64_t{ nb_values_ } * dimension );
            if( components.size() % dimension != 0 )
            {
                geode::Logger::warning( "[OMFInput] Data \"", name_,
                    "\" has a number of components not multiple of ", dimension,
                    ": the last incomplete vector is ignored" );
            }
            std::vector< geode::Vector< dimension > > vectors(
                components.size() / dimension );
            for( const auto vector : geode::Indices{ vectors } )
            {
                for( const auto axis : geode::LRange{ dimension } )
                {
                    vectors[vector].set_value(
                        axis, components[dimension * vector + axis] );
                }
            }
            set_values( name_, absl::StrCat( "vector", dimension ),
                omf_no_data_vector< dimension >(), true, vectors );
        }

        void read_color_data()
        {
            const auto& colors = array();
            if( omf_class_name( colors ) == "ColorArray" )
            {
                set_values( name_, "color", omf_no_data_color(), false,
                    read_omf_inline_colors( colors ) );
                return;
            }
            const auto components =
                file_.read_integers( colors, uint64_t{ nb_values_ } * 3 );
            std::vector< geode::RGBColor > values( components.size() / 3 );
            for( const auto color : geode::Indices{ values } )
            {
                values[color] = omf_color( components[3 * color],
                    components[3 * color + 1], components[3 * color + 2] );
            }
            set_values( name_, "color", omf_no_data_color(), false, values );
        }

        /*
         * The indices are stored in "<name>_index" (-1 for no data) and the
         * value of each legend is stored in an attribute named after the
         * legend
         */
        void read_mapped_data()
        {
            const auto indices = file_.read_integers( array(), nb_values_ );
            const auto nb =
                omf_usable_size( name_, indices.size(), nb_values_ );
            geode::AttributeValues< int > index_values;
            index_values.default_value = -1;
            index_values.no_value = -1;
            const auto index_attribute = find_or_create_omf_attribute< int >(
                manager_, absl::StrCat( name_, "_index" ), "index",
                std::move( index_values ), false );
            for( const auto value : geode::Range{ nb } )
            {
                const auto index = indices[value];
                const auto element = attribute_element_( value );
                if( element != geode::NO_ID && index >= 0
                    && index <= std::numeric_limits< int >::max() )
                {
                    index_attribute->set_value(
                        element, static_cast< int >( index ) );
                }
            }
            for( const auto& legend_uid :
                data_.value( "legends", nlohmann::json::array() ) )
            {
                read_legend( file_.object( legend_uid ), indices, nb );
            }
        }

        void read_legend( const nlohmann::json& legend,
            const std::vector< int64_t >& indices,
            geode::index_t nb )
        {
            auto legend_name = geode::internal::OMFFile::name( legend );
            if( legend_name.empty() )
            {
                legend_name = name_;
            }
            const auto& legend_values =
                file_.object( omf_member( legend, "values" ) );
            const auto legend_class = omf_class_name( legend_values );
            if( legend_class == "StringArray"
                || legend_class == "DateTimeArray" )
            {
                set_legend_values( legend_name, "string", omf_no_data_string(),
                    read_omf_strings( legend_values ), indices, nb );
            }
            else if( legend_class == "ColorArray" )
            {
                set_legend_values( legend_name, "color", omf_no_data_color(),
                    read_omf_inline_colors( legend_values ), indices, nb );
            }
            else if( legend_class == "ScalarArray" )
            {
                set_legend_values( legend_name, "scalar", omf_no_data_double(),
                    file_.read_scalars( legend_values ), indices, nb );
            }
            else
            {
                geode::Logger::warning( "[OMFInput] Legend \"", legend_name,
                    "\" of unsupported class ", legend_class, " is skipped" );
            }
        }

        template < typename T >
        void set_legend_values( std::string_view name,
            std::string_view type_suffix,
            geode::AttributeValues< T > default_values,
            const std::vector< T >& legend,
            const std::vector< int64_t >& indices,
            geode::index_t nb )
        {
            const auto attribute = find_or_create_omf_attribute< T >( manager_,
                name, type_suffix, std::move( default_values ), false );
            for( const auto value : geode::Range{ nb } )
            {
                const auto index = indices[value];
                const auto element = attribute_element_( value );
                if( element != geode::NO_ID && index >= 0
                    && static_cast< uint64_t >( index ) < legend.size() )
                {
                    attribute->set_value( element, legend[index] );
                }
            }
        }

    private:
        const geode::internal::OMFFile& file_;
        const nlohmann::json& data_;
        std::string_view name_;
        geode::index_t nb_values_;
        geode::AttributeManager& manager_;
        absl::FunctionRef< geode::index_t( geode::index_t ) >
            attribute_element_;
    };
} // namespace

namespace geode::internal
{
    OMFFile::OMFFile( std::string_view filename )
        : filename_{ to_string( filename ) },
          file_{ filename_, std::ios::binary }
    {
        OpenGeodeGeosciencesIOMeshException::check_exception( file_.good(),
            nullptr, OpenGeodeException::TYPE::data,
            "[OMFInput] Error while opening file: ", filename );
        const auto header = read_omf_header( file_ );
        OpenGeodeGeosciencesIOMeshException::check_exception(
            header.has_value(), nullptr, OpenGeodeException::TYPE::data,
            "[OMFInput] ", filename,
            " is not an OMF v1 file (wrong magic number or version)" );
        json_start_ = header->json_start;
        OpenGeodeGeosciencesIOMeshException::check_exception(
            json_start_ >= static_cast< uint64_t >( HEADER_SIZE )
                && json_start_ < header->file_size,
            nullptr, OpenGeodeException::TYPE::data, "[OMFInput] ", filename,
            " has an invalid JSON start (", json_start_, " for a file of ",
            header->file_size, " bytes)" );
        std::string json( header->file_size - json_start_, '\0' );
        file_.seekg( static_cast< std::streamoff >( json_start_ ) );
        file_.read(
            json.data(), static_cast< std::streamsize >( json.size() ) );
        index_ = nlohmann::json::parse( json, nullptr, false );
        OpenGeodeGeosciencesIOMeshException::check_exception(
            index_.is_object(), nullptr, OpenGeodeException::TYPE::data,
            "[OMFInput] ", filename, " has an invalid JSON index" );
        for( const auto& [uid, object] : index_.items() )
        {
            if( object.is_object() && omf_class_name( object ) == "Project" )
            {
                project_ = &object;
                break;
            }
        }
        OpenGeodeGeosciencesIOMeshException::check_exception(
            project_ != nullptr, nullptr, OpenGeodeException::TYPE::data,
            "[OMFInput] ", filename, " has no Project" );
        project_origin_ = omf_vector( *project_, "origin", Vector3D{} );
    }

    Percentage OMFFile::is_loadable( std::string_view filename,
        std::string_view element_class,
        std::string_view geometry_class )
    {
        return is_loadable(
            filename, [element_class, geometry_class]( const OMFFile& omf ) {
                return omf.elements( element_class, geometry_class ).size();
            } );
    }

    Percentage OMFFile::is_loadable( std::string_view filename,
        absl::FunctionRef< std::size_t( const OMFFile& ) >
            nb_loadable_elements )
    {
        try
        {
            const auto omf = loadability_file( filename );
            if( !omf )
            {
                return Percentage{ 0 };
            }
            const auto nb_elements = omf->nb_elements();
            if( nb_elements == 0 )
            {
                return Percentage{ 0 };
            }
            return Percentage{ static_cast< double >(
                                   nb_loadable_elements( *omf ) )
                               / static_cast< double >( nb_elements ) };
        }
        catch( ... )
        {
            return Percentage{ 0 };
        }
    }

    std::shared_ptr< const OMFFile > OMFFile::loadability_file(
        std::string_view filename )
    {
        struct CachedFile
        {
            std::string path;
            std::uintmax_t size{ 0 };
            std::filesystem::file_time_type time;
            std::shared_ptr< const OMFFile > file;
        };
        thread_local CachedFile cache;
        auto path = to_string( filename );
        std::error_code error;
        const auto size = std::filesystem::file_size( path, error );
        const auto time = error
                              ? std::filesystem::file_time_type{}
                              : std::filesystem::last_write_time( path, error );
        if( !error && cache.file && cache.path == path && cache.size == size
            && cache.time == time )
        {
            return cache.file;
        }
        cache.file.reset();
        {
            std::ifstream file{ path, std::ios::binary };
            if( !file.good() || !read_omf_header( file ) )
            {
                return nullptr;
            }
        }
        auto omf = std::make_shared< OMFFile >( filename );
        omf->warn_missing_elements_ = false;
        omf->file_.close();
        if( !error )
        {
            cache = { std::move( path ), size, time, omf };
        }
        return omf;
    }

    std::size_t OMFFile::nb_elements() const
    {
        std::size_t count{ 0 };
        for( const auto& uid :
            project_->value( "elements", nlohmann::json::array() ) )
        {
            const auto it = uid.is_string()
                                ? index_.find( uid.get< std::string >() )
                                : index_.end();
            if( it != index_.end() && it->is_object() )
            {
                count++;
            }
        }
        return count;
    }

    std::string_view OMFFile::project_name() const
    {
        return name( *project_ );
    }

    std::vector< const nlohmann::json* > OMFFile::elements(
        std::string_view element_class, std::string_view geometry_class ) const
    {
        std::vector< const nlohmann::json* > result;
        for( const auto& uid :
            project_->value( "elements", nlohmann::json::array() ) )
        {
            const auto it = uid.is_string()
                                ? index_.find( uid.get< std::string >() )
                                : index_.end();
            if( it == index_.end() || !it->is_object() )
            {
                if( warn_missing_elements_ )
                {
                    Logger::warning( "[OMFInput] Project element ", uid.dump(),
                        " is missing from the file and is skipped" );
                }
                continue;
            }
            if( omf_class_name( *it ) != element_class )
            {
                continue;
            }
            if( !geometry_class.empty() )
            {
                const auto geometry_it = it->find( "geometry" );
                if( geometry_it == it->end() || !geometry_it->is_string()
                    || !index_.contains( geometry_it->get< std::string >() )
                    || omf_class_name(
                           index_.at( geometry_it->get< std::string >() ) )
                           != geometry_class )
                {
                    continue;
                }
            }
            result.push_back( &*it );
        }
        return result;
    }

    const nlohmann::json& OMFFile::object( const nlohmann::json& uid ) const
    {
        const auto it = uid.is_string()
                            ? index_.find( uid.get< std::string >() )
                            : index_.end();
        OpenGeodeGeosciencesIOMeshException::check_exception(
            it != index_.end() && it->is_object(), nullptr,
            OpenGeodeException::TYPE::data, "[OMFInput] Object ", uid.dump(),
            " referenced in ", filename_, " does not exist" );
        return *it;
    }

    const nlohmann::json& OMFFile::geometry(
        const nlohmann::json& element ) const
    {
        return object( omf_member( element, "geometry" ) );
    }

    std::string_view OMFFile::name( const nlohmann::json& object )
    {
        const auto it = object.find( "name" );
        if( it == object.end() || !it->is_string() )
        {
            return {};
        }
        return it->get_ref< const std::string& >();
    }

    std::string_view OMFFile::mesh_name(
        absl::Span< const nlohmann::json* const > elements ) const
    {
        if( elements.size() == 1 )
        {
            return name( *elements.front() );
        }
        return project_name();
    }

    std::vector< double > OMFFile::read_tensor(
        const nlohmann::json& geometry, std::string_view key )
    {
        const auto it = geometry.find( key );
        OpenGeodeGeosciencesIOMeshException::check_exception(
            it != geometry.end() && it->is_array() && !it->empty(), nullptr,
            OpenGeodeException::TYPE::data, "[OMFInput] Missing grid ", key );
        std::vector< double > tensor;
        tensor.reserve( it->size() );
        for( const auto& value : *it )
        {
            OpenGeodeGeosciencesIOMeshException::check_exception(
                value.is_number() && value.get< double >() > 0
                    && std::isfinite( value.get< double >() ),
                nullptr, OpenGeodeException::TYPE::data, "[OMFInput] Grid ",
                key, " contains a non positive cell size" );
            tensor.push_back( value.get< double >() );
        }
        return tensor;
    }

    Vector3D OMFFile::axis( const nlohmann::json& geometry,
        std::string_view key,
        local_index_t default_axis )
    {
        Vector3D default_value;
        default_value.set_value( default_axis, 1 );
        const auto axis = omf_vector( geometry, key, default_value );
        OpenGeodeGeosciencesIOMeshException::check_exception(
            axis.length() > GLOBAL_EPSILON, nullptr,
            OpenGeodeException::TYPE::data, "[OMFInput] Grid ", key,
            " is a null vector" );
        return axis.normalize();
    }

    Vector3D OMFFile::origin( const nlohmann::json& geometry ) const
    {
        return project_origin_ + omf_vector( geometry, "origin", Vector3D{} );
    }

    std::vector< Point3D > OMFFile::read_points(
        const nlohmann::json& geometry, std::string_view key ) const
    {
        const auto coordinates =
            read_scalars( object( omf_member( geometry, key ) ) );
        OpenGeodeGeosciencesIOMeshException::check_exception(
            coordinates.size() % 3 == 0, nullptr,
            OpenGeodeException::TYPE::data, "[OMFInput] Number of ", key,
            " coordinates (", coordinates.size(), ") is not a multiple of 3" );
        OpenGeodeGeosciencesIOMeshException::check_exception(
            coordinates.size() / 3 < NO_ID, nullptr,
            OpenGeodeException::TYPE::data, "[OMFInput] Too many ", key );
        const auto translation = origin( geometry );
        std::vector< Point3D > points( coordinates.size() / 3 );
        for( const auto point : Indices{ points } )
        {
            points[point] =
                Point3D{ { coordinates[3 * point], coordinates[3 * point + 1],
                    coordinates[3 * point + 2] } }
                + translation;
        }
        return points;
    }

    template < index_t nb_vertices_per_element >
    std::vector< std::array< index_t, nb_vertices_per_element > >
        OMFFile::read_elements( const nlohmann::json& geometry,
            std::string_view key,
            index_t nb_vertices ) const
    {
        const auto indices =
            read_integers( object( omf_member( geometry, key ) ) );
        OpenGeodeGeosciencesIOMeshException::check_exception(
            indices.size() % nb_vertices_per_element == 0, nullptr,
            OpenGeodeException::TYPE::data, "[OMFInput] Number of ", key,
            " indices (", indices.size(), ") is not a multiple of ",
            nb_vertices_per_element );
        std::vector< std::array< index_t, nb_vertices_per_element > > elements;
        elements.reserve( indices.size() / nb_vertices_per_element );
        index_t nb_invalid{ 0 };
        index_t nb_degenerate{ 0 };
        for( std::size_t first = 0; first < indices.size();
            first += nb_vertices_per_element )
        {
            auto& element = elements.emplace_back();
            for( const auto vertex : LRange{ nb_vertices_per_element } )
            {
                const auto index = indices[first + vertex];
                if( index < 0 || index >= nb_vertices )
                {
                    element.fill( NO_ID );
                    nb_invalid++;
                    break;
                }
                element[vertex] = static_cast< index_t >( index );
                if( absl::c_count(
                        absl::MakeConstSpan( element.data(), vertex ),
                        element[vertex] )
                    != 0 )
                {
                    element.fill( NO_ID );
                    nb_degenerate++;
                    break;
                }
            }
        }
        if( nb_invalid > 0 )
        {
            Logger::warning( "[OMFInput] ", nb_invalid, " ", key,
                " referencing vertices outside [0, ", nb_vertices,
                ") are skipped" );
        }
        if( nb_degenerate > 0 )
        {
            Logger::warning( "[OMFInput] ", nb_degenerate, " ", key,
                " with repeated vertices are skipped" );
        }
        return elements;
    }

    template std::vector< std::array< index_t, 2 > >
        OMFFile::read_elements< 2 >(
            const nlohmann::json&, std::string_view, index_t ) const;
    template std::vector< std::array< index_t, 3 > >
        OMFFile::read_elements< 3 >(
            const nlohmann::json&, std::string_view, index_t ) const;

    std::vector< double > OMFFile::read_optional_scalars(
        const nlohmann::json& geometry,
        std::string_view key,
        uint64_t max_values ) const
    {
        const auto it = geometry.find( key );
        if( it == geometry.end() || it->is_null() )
        {
            return {};
        }
        return read_scalars( object( *it ), max_values );
    }

    void OMFFile::read_data( const nlohmann::json& element,
        std::string_view location,
        index_t nb_values,
        AttributeManager& manager,
        absl::FunctionRef< index_t( index_t ) > attribute_element ) const
    {
        for( const auto& data_uid :
            element.value( "data", nlohmann::json::array() ) )
        {
            const nlohmann::json* data{ nullptr };
            try
            {
                data = &object( data_uid );
                if( data->value( "location", "" ) != location )
                {
                    continue;
                }
                OMFDataReader{ *this, *data, nb_values, manager,
                    attribute_element }
                    .read();
            }
            catch( const std::exception& exception )
            {
                Logger::warning( "[OMFInput] Data ",
                    data ? absl::StrCat( "\"", name( *data ), "\"" )
                         : data_uid.dump(),
                    " of element \"", name( element ),
                    "\" is skipped: ", exception.what() );
            }
        }
    }

    std::vector< double > OMFFile::read_scalars(
        const nlohmann::json& array, uint64_t max_values ) const
    {
        const auto& content = omf_member( array, ARRAY_KEY );
        if( content.is_array() )
        {
            std::vector< double > values;
            values.reserve( content.size() );
            for( const auto& value : content )
            {
                values.push_back(
                    value.is_number()
                        ? value.get< double >()
                        : std::numeric_limits< double >::quiet_NaN() );
            }
            return values;
        }
        const auto dtype = content.value( "dtype", "" );
        if( dtype == "<f8" )
        {
            return read_binary< double >( content, max_values );
        }
        OpenGeodeGeosciencesIOMeshException::check_exception( dtype == "<i8",
            nullptr, OpenGeodeException::TYPE::data,
            "[OMFInput] Unsupported array dtype ", dtype );
        const auto integers = read_binary< int64_t >( content, max_values );
        std::vector< double > values( integers.size() );
        for( const auto value : Indices{ values } )
        {
            values[value] = static_cast< double >( integers[value] );
        }
        return values;
    }

    std::vector< int64_t > OMFFile::read_integers(
        const nlohmann::json& array, uint64_t max_values ) const
    {
        const auto& content = omf_member( array, ARRAY_KEY );
        if( content.is_array() )
        {
            std::vector< int64_t > values;
            values.reserve( content.size() );
            for( const auto& value : content )
            {
                if( !value.is_array() )
                {
                    values.push_back( omf_json_integer( value ) );
                    continue;
                }
                for( const auto& component : value )
                {
                    values.push_back( omf_json_integer( component ) );
                }
            }
            return values;
        }
        const auto dtype = content.value( "dtype", "" );
        if( dtype == "<i8" )
        {
            return read_binary< int64_t >( content, max_values );
        }
        OpenGeodeGeosciencesIOMeshException::check_exception( dtype == "<f8",
            nullptr, OpenGeodeException::TYPE::data,
            "[OMFInput] Unsupported array dtype ", dtype );
        const auto reals = read_binary< double >( content, max_values );
        std::vector< int64_t > values( reals.size() );
        for( const auto value : Indices{ values } )
        {
            values[value] = omf_real_to_integer( reals[value] );
        }
        return values;
    }

    template < typename T >
    std::vector< T > OMFFile::read_binary(
        const nlohmann::json& binary, uint64_t max_values ) const
    {
        static_assert( sizeof( T ) == 8 );
        const auto start = omf_json_offset( binary, "start" );
        const auto length = omf_json_offset( binary, "length" );
        OpenGeodeGeosciencesIOMeshException::check_exception(
            start >= static_cast< uint64_t >( HEADER_SIZE )
                && length <= json_start_ && start <= json_start_ - length
                && length <= std::numeric_limits< uInt >::max(),
            nullptr, OpenGeodeException::TYPE::data, "[OMFInput] Array block [",
            start, ", ", start + length, ") is outside the data section [",
            HEADER_SIZE, ", ", json_start_, ")" );
        std::vector< char > compressed( length );
        file_.clear();
        file_.seekg( static_cast< std::streamoff >( start ) );
        file_.read(
            compressed.data(), static_cast< std::streamsize >( length ) );
        OpenGeodeGeosciencesIOMeshException::check_exception( file_.good(),
            nullptr, OpenGeodeException::TYPE::data,
            "[OMFInput] Failed to read array block at ", start );

        z_stream stream{};
        OpenGeodeGeosciencesIOMeshException::check_exception(
            inflateInit( &stream ) == Z_OK, nullptr,
            OpenGeodeException::TYPE::internal,
            "[OMFInput] Failed to initialize zlib" );
        stream.next_in = reinterpret_cast< Bytef* >( compressed.data() );
        stream.avail_in = static_cast< uInt >( length );
        const auto max_size = max_values * sizeof( T );
        std::vector< T > result;
        std::size_t size{ 0 };
        int status{ Z_OK };
        bool truncated{ false };
        while( status == Z_OK )
        {
            result.resize( result.size() + INFLATE_CHUNK_SIZE / sizeof( T ) );
            const auto capacity = result.size() * sizeof( T );
            stream.next_out =
                reinterpret_cast< Bytef* >( result.data() ) + size;
            stream.avail_out = static_cast< uInt >( capacity - size );
            status = inflate( &stream, Z_NO_FLUSH );
            size = capacity - stream.avail_out;
            if( size > max_size )
            {
                size = max_size;
                truncated = true;
                break;
            }
            if( status == Z_BUF_ERROR && stream.avail_in > 0 )
            {
                status = Z_OK;
            }
        }
        inflateEnd( &stream );
        if( truncated )
        {
            Logger::warning( "[OMFInput] Array block at ", start,
                " has more than the ", max_values,
                " expected values: extra values are ignored" );
        }
        else
        {
            OpenGeodeGeosciencesIOMeshException::check_exception(
                status == Z_STREAM_END, nullptr, OpenGeodeException::TYPE::data,
                "[OMFInput] Corrupted compressed array block at ", start );
            OpenGeodeGeosciencesIOMeshException::check_exception(
                size % sizeof( T ) == 0, nullptr,
                OpenGeodeException::TYPE::data, "[OMFInput] Array block at ",
                start, " has ", size, " bytes, which is not a multiple of ",
                sizeof( T ) );
        }
        result.resize( size / sizeof( T ) );
        return result;
    }

    void set_omf_element_attribute( AttributeManager& manager,
        index_t element_id,
        index_t begin,
        index_t end )
    {
        AttributeValues< index_t > values;
        values.default_value = NO_ID;
        values.no_value = NO_ID;
        const auto attribute = find_or_create_omf_attribute< index_t >( manager,
            OMFFile::OMF_ELEMENT_ATTRIBUTE_NAME, "id", std::move( values ),
            false );
        for( const auto element : Range{ begin, end } )
        {
            attribute->set_value( element, element_id );
        }
    }
} // namespace geode::internal
