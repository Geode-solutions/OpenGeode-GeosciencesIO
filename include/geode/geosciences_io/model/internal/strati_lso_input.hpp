/*
 * Copyright (c) 2019 - 2025 Geode-solutions. All rights reserved.
 */

#pragma once

#include <geode/geosciences/implicit/representation/io/stratigraphic_model_input.hpp>

#include <geode/geosciences_io/model/common.hpp>

namespace geode
{
    namespace internal
    {
        class StratigraphicLSOInput final : public StratigraphicModelInput
        {
        public:
            explicit StratigraphicLSOInput( std::string_view filename )
                : StratigraphicModelInput( filename )
            {
            }

            [[nodiscard]] static std::string_view extension()
            {
                static constexpr auto EXT = "lso";
                return EXT;
            }

            [[nodiscard]] StratigraphicModel read() final;

            AdditionalFiles additional_files() const final
            {
                return {};
            }

            index_t object_priority() const final
            {
                return 1;
            }

            Percentage is_loadable() const final;
        };
    } // namespace internal
} // namespace geode