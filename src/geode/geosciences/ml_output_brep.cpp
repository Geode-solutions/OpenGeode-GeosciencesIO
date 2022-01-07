/*
 * Copyright (c) 2019 - 2022 Geode-solutions
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

#include <geode/geosciences/private/ml_output_brep.h>

#include <geode/mesh/core/surface_mesh.h>

#include <geode/geosciences/private/gocad_common.h>
#include <geode/geosciences/private/ml_output_impl.h>

namespace
{
    class MLOutputImplBRep : public geode::detail::MLOutputImpl< geode::BRep >
    {
    public:
        MLOutputImplBRep( absl::string_view filename, const geode::BRep& model )
            : geode::detail::MLOutputImpl< geode::BRep >( filename, model )
        {
        }

    private:
        void write_geological_tfaces() override {}

        void write_geological_tsurfs() override {}

        void write_geological_regions() override {}

        void write_geological_model_surfaces() override {}
    };
} // namespace

namespace geode
{
    namespace detail
    {
        void MLOutputBRep::write() const
        {
            const auto only_triangles = check_brep_polygons( brep() );
            if( !only_triangles )
            {
                geode::Logger::info(
                    "[MLOutput::write] Can not export into .ml a "
                    "BRep with non triangular surface polygons." );
                return;
            }
            MLOutputImplBRep impl{ filename(), brep() };
            impl.write_file();
        }
    } // namespace detail
} // namespace geode
