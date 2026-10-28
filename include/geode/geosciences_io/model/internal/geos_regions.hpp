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

#include <geode/geosciences_io/model/common.hpp>
#include <geode/geosciences_io/model/internal/geos_deck_context.hpp>

namespace geode
{
    class BRep;
    class StructuralModel;
} // namespace geode

namespace geode::internal
{
    /// Deck context holding one GEOS region per Block
    [[nodiscard]] GeosDeckContext geos_region_context( const BRep& model );

    /*!
     * Deck context holding one GEOS region per StratigraphicUnit (its Blocks
     * are merged), then one per remaining Block
     */
    [[nodiscard]] GeosDeckContext geos_region_context(
        const StructuralModel& model );

    /*!
     * Deck context of a model whose exporter is unknown: StratigraphicUnits
     * are found from the Block collections and merged as by the
     * StructuralModel exporter
     */
    [[nodiscard]] GeosDeckContext geos_inferred_region_context(
        const BRep& model );
} // namespace geode::internal
