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

#include <string>
#include <string_view>

#include <absl/container/flat_hash_map.h>
#include <absl/functional/function_ref.h>

namespace pugi
{
    class xml_node;
} // namespace pugi

namespace geode
{
    class PhysicalProblem;
} // namespace geode

namespace geode::internal
{
    /// Solver option role -> XPath of the nodes it modifies
    using SolverOptionRoles = absl::flat_hash_map< std::string, std::string >;

    /// How the solver options of a simulator are read
    struct SimulatorSolverOptions
    {
        /// Key of the options in the problem solver options
        std::string_view key;
        /// Display name, used in issue messages
        std::string_view name;
        /// Prefix of the warnings
        std::string_view context;
        const SolverOptionRoles& roles;
    };

    /// Real value written with the fewest significant digits (at least 15)
    /// that read back to the same value
    [[nodiscard]] std::string xml_real( double value );

    void set_attribute( pugi::xml_node& node,
        std::string_view attribute,
        std::string_view value );

    void set_real_attribute(
        pugi::xml_node& node, std::string_view attribute, double value );

    void add_comment( pugi::xml_node& node, std::string_view comment );

    /*!
     * Applies the "<role>/<name>" solver options of a simulator: each option
     * is given to \p apply for every node selected by the XPath of its role.
     * Warns, with \p context as prefix, when no node matches a role.
     */
    void apply_solver_options( pugi::xml_node& root,
        const PhysicalProblem& problem,
        const SimulatorSolverOptions& simulator,
        absl::FunctionRef< void(
            pugi::xml_node&, std::string_view name, std::string_view value ) >
            apply );

    /// Reports through \p add_issue each solver option that is not
    /// "<role>/<attribute>" with a known role
    void check_solver_option_keys( const PhysicalProblem& problem,
        const SimulatorSolverOptions& simulator,
        absl::FunctionRef< void( std::string_view ) > add_issue );
} // namespace geode::internal
