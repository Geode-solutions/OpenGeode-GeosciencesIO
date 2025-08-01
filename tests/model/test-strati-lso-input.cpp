/*
 * Copyright (c) 2019 - 2025 Geode-solutions. All rights reserved.
 */

#include <geode/tests_config.hpp>

#include <geode/basic/assert.hpp>
#include <geode/basic/logger.hpp>

#include <geode/geosciences/implicit/representation/core/stratigraphic_model.hpp>
#include <geode/geosciences/implicit/representation/io/stratigraphic_model_input.hpp>

#include <geode/implicit/io/common.hpp>

void test_lso_input()
{
    geode::Logger::info( "Starting test" );

    auto implicit_model = geode::load_stratigraphic_model(
        absl::StrCat( geode::DATA_PATH, "vri.lso" ) );
}

int main()
{
    try
    {
        geode::ImplicitIOLibrary::initialize();
        test_lso_input();

        geode::Logger::info( "TEST SUCCESS" );
        return 0;
    }
    catch( ... )
    {
        return geode::geode_lippincott();
    }
}