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

#include <geode/geosciences_io/model/private/geos_output.h>
#include <geode/geosciences_io/model/private/geos_output_brep.h>

namespace
{
    class GeosBRepOutputImpl
        : public geode::detail::GeosOutputImpl< geode::BRep >
    {
    public:
        GeosBRepOutputImpl(
            absl::string_view files_directory, const geode::BRep& model )
            : geode::detail::GeosOutputImpl< geode::BRep >(
                files_directory, model )
        {
        }
    };
} // namespace

namespace geode
{
    namespace detail
    {
        std::vector< std::string > GeosOutputBRep::write(
            const BRep& brep ) const
        {
            GeosBRepOutputImpl impl{ filename(), brep };
            impl.write_file();
            return { to_string( filename() ) };
        }

        bool GeosOutputBRep::is_saveable( const BRep& brep ) const
        {
            return false;
            // return check_brep_polygons( brep );
        }
    } // namespace detail
} // namespace geode
