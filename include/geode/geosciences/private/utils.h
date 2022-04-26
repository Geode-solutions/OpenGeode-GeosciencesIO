/*
 * Copyright (c) 2019 Geode-solutions
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
#include <deque>
#include <string>

#include <geode/basic/bitsery_archive.h>
#include <geode/basic/uuid.h>

#include <geode/geosciences/detail/common.h>

namespace geode
{
    class BRep;
} // namespace geode

namespace geode
{
    namespace detail
    {
        bool string_starts_with(
            absl::string_view string, absl::string_view check );

        void check_keyword( std::ifstream& file, absl::string_view keyword );

        bool line_starts_with( std::ifstream& file, absl::string_view check );

        std::string goto_keyword( std::ifstream& file, absl::string_view word );

        std::string goto_keywords(
            std::ifstream& file, absl::Span< const absl::string_view > words );

        absl::optional< std::string > goto_keyword_if_it_exists(
            std::ifstream& file, absl::string_view word );

        double read_index_t( absl::string_view token );

        double read_double( absl::string_view token );
    } // namespace detail
} // namespace geode
