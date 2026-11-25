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

#include <geode/geosciences_io/model/internal/simulation_checks.hpp>

#include <absl/algorithm/container.h>
#include <absl/strings/str_cat.h>

#include <geode/geosciences/simulation/core/source_term.hpp>

namespace geode::internal
{
    void check_source_term_well( const SourceTerm& source,
        absl::Span< const std::string > well_names,
        SimulationReport& report )
    {
        if( absl::c_find( well_names, source.well ) != well_names.end() )
        {
            return;
        }
        SimulationIssue issue{ SIMULATION_ISSUE_SEVERITY::error,
            SIMULATION_CHECK_STAGE::simulator,
            absl::StrCat( "Source term on well \"", source.well, "\"" ),
            "the well has no perforations given to the exporter" };
        issue.fix = "add the well perforations to the exporter";
        report.add_issue( std::move( issue ) );
    }
} // namespace geode::internal
