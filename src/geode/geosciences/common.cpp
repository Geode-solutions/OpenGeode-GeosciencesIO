/*
 * Copyright (c) 2019 - 2021 Geode-solutions
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

#include <geode/geosciences/detail/common.h>

#include <geode/mesh/io/triangulated_surface_input.h>

#include <geode/geosciences/private/lso_input.h>
#include <geode/geosciences/private/ml_input.h>
#include <geode/geosciences/private/ml_output.h>
#include <geode/geosciences/private/ts_input.h>
#include <geode/geosciences/representation/io/structural_model_input.h>
#include <geode/geosciences/representation/io/structural_model_output.h>

namespace
{
    void register_structural_model_input()
    {
        geode::StructuralModelInputFactory::register_creator<
            geode::detail::MLInput >(
            geode::detail::MLInput::extension().data() );
        geode::StructuralModelInputFactory::register_creator<
            geode::detail::LSOInput >(
            geode::detail::LSOInput::extension().data() );
    }

    void register_structural_model_output()
    {
        geode::StructuralModelOutputFactory::register_creator<
            geode::detail::MLOutput >(
            geode::detail::MLOutput::extension().data() );
    }

    void register_triangulated_surface_input()
    {
        geode::TriangulatedSurfaceInputFactory3D::register_creator<
            geode::detail::TSInput >(
            geode::detail::TSInput::extension().data() );
    }

    OPENGEODE_LIBRARY_INITIALIZE( OpenGeode_GeosciencesIO_geosciences )
    {
        register_structural_model_input();
        register_structural_model_output();
        register_triangulated_surface_input();
    }
} // namespace

namespace geode
{
    namespace detail
    {
        void initialize_geosciences_io() {}
    } // namespace detail
} // namespace geode
