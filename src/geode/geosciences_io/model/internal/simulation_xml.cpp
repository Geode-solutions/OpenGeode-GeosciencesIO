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

#include <geode/geosciences_io/model/internal/simulation_xml.hpp>

#include <cstdlib>
#include <vector>

#include <pugixml.hpp>

#include <absl/strings/str_cat.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_split.h>

#include <geode/basic/logger.hpp>

#include <geode/geosciences/simulation/core/physical_problem.hpp>
#include <geode/geosciences/simulation/core/solver_options.hpp>

namespace geode::internal
{
    std::string xml_real( double value )
    {
        static constexpr int MIN_DIGITS{ 15 };
        static constexpr int MAX_DIGITS{ 17 };
        for( auto digits = MIN_DIGITS; digits < MAX_DIGITS; digits++ )
        {
            auto text = absl::StrFormat( "%.*g", digits, value );
            if( std::strtod( text.c_str(), nullptr ) == value )
            {
                return text;
            }
        }
        return absl::StrFormat( "%.*g", MAX_DIGITS, value );
    }

    void set_attribute( pugi::xml_node& node,
        std::string_view attribute,
        std::string_view value )
    {
        const std::string attribute_name{ attribute };
        auto xml_attribute = node.attribute( attribute_name.c_str() );
        if( !xml_attribute )
        {
            xml_attribute = node.append_attribute( attribute_name.c_str() );
        }
        xml_attribute.set_value( std::string{ value }.c_str() );
    }

    void set_real_attribute(
        pugi::xml_node& node, std::string_view attribute, double value )
    {
        set_attribute( node, attribute, xml_real( value ) );
    }

    void add_comment( pugi::xml_node& node, std::string_view comment )
    {
        node.append_child( pugi::node_comment )
            .set_value( absl::StrCat( " ", comment, " " ).c_str() );
    }

    void apply_solver_options( pugi::xml_node& root,
        const PhysicalProblem& problem,
        const SimulatorSolverOptions& simulator,
        absl::FunctionRef< void(
            pugi::xml_node&, std::string_view name, std::string_view value ) >
            apply )
    {
        for( const auto& [key, value] :
            problem.solver_options().options( simulator.key ) )
        {
            const std::vector< std::string > parts =
                absl::StrSplit( key, absl::MaxSplits( '/', 1 ) );
            const auto& xpath = simulator.roles.at( parts[0] );
            const auto nodes = root.select_nodes( xpath.c_str() );
            if( nodes.empty() )
            {
                Logger::warning( simulator.context, " Solver option ", key,
                    " ignored: no ", parts[0], " node in the deck" );
            }
            for( const auto& selected : nodes )
            {
                auto node = selected.node();
                apply( node, parts[1], value );
            }
        }
    }

    void check_solver_option_keys( const PhysicalProblem& problem,
        const SimulatorSolverOptions& simulator,
        absl::FunctionRef< void( std::string_view ) > add_issue )
    {
        for( const auto& [key, value] :
            problem.solver_options().options( simulator.key ) )
        {
            const std::vector< std::string_view > parts =
                absl::StrSplit( key, absl::MaxSplits( '/', 1 ) );
            if( parts.size() != 2 || parts[1].empty()
                || !simulator.roles.contains( parts[0] ) )
            {
                add_issue( absl::StrCat( simulator.name, " solver option \"",
                    key,
                    R"(" must be "<role>/<attribute>" with a known role)" ) );
            }
        }
    }
} // namespace geode::internal
