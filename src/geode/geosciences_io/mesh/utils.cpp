/*
 * Copyright (c) 2019 - 2023 Geode-solutions
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

#include <geode/geosciences_io/mesh/private/utils.h>

#include <fstream>

#include <absl/strings/match.h>

#include <geode/basic/logger.h>

namespace geode
{
    namespace detail
    {
        bool string_starts_with(
            absl::string_view string, absl::string_view check )
        {
            return absl::StartsWith( string, check );
        }

        bool line_starts_with( std::ifstream& file, absl::string_view check )
        {
            std::string line;
            std::getline( file, line );
            return string_starts_with( line, check );
        }

        void check_keyword( std::ifstream& file, absl::string_view keyword )
        {
            OPENGEODE_EXCEPTION( line_starts_with( file, keyword ),
                absl::StrCat( "Line should starts with \"", keyword, "\"" ) );
        }

        std::string goto_keyword( std::ifstream& file, absl::string_view word )
        {
            std::string line;
            while( std::getline( file, line ) )
            {
                if( string_starts_with( line, word ) )
                {
                    return line;
                }
            }
            throw geode::OpenGeodeException{
                "[goto_keyword] Cannot find the requested keyword: ", word
            };
            return "";
        }

        std::string goto_keywords(
            std::ifstream& file, absl::Span< const absl::string_view > words )
        {
            std::string line;
            while( std::getline( file, line ) )
            {
                for( const auto word : words )
                {
                    if( string_starts_with( line, word ) )
                    {
                        return line;
                    }
                }
            }
            throw geode::OpenGeodeException{
                "[goto_keywords] Cannot find one of the requested keywords"
            };
            return "";
        }

        absl::optional< std::string > goto_keyword_if_it_exists(
            std::ifstream& file, absl::string_view word )
        {
            std::string line;
            while( std::getline( file, line ) )
            {
                if( string_starts_with( line, word ) )
                {
                    return line;
                }
            }
            geode::Logger::debug(
                "[goto_keyword_if_it_exists] Couldn't find word ", word,
                " in the file, returning to file begin." );
            file.clear();
            file.seekg( std::ios::beg );
            return absl::nullopt;
        }
    } // namespace detail
} // namespace geode
