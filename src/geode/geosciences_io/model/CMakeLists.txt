# Copyright (c) 2019 - 2024 Geode-solutions
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

add_geode_library(
    NAME model
    FOLDER "geode/geosciences_io/model"
    SOURCES
        "common.cpp"
        "helpers/brep_geos_export.cpp"
        "helpers/structural_model_geos_export.cpp"
        "internal/brep_fem_output.cpp"
        "internal/geos_export.cpp"
        "internal/gocad_common.cpp"
        "internal/lso_input.cpp"
        "internal/lso_output.cpp"
        "internal/ml_input.cpp"
        "internal/ml_output_brep.cpp"
        "internal/ml_output_structural_model.cpp"
        "internal/shp_input.cpp"
        "internal/horizons_stack_skua_input.cpp"
    PUBLIC_HEADERS
        "common.hpp"
    INTERNAL_HEADERS
        "helpers/brep_geos_export.hpp"
        "helpers/structural_model_geos_export.hpp"
        "internal/brep_fem_output.hpp"
        "internal/geos_export.hpp"
        "internal/gocad_common.hpp"
        "internal/lso_input.hpp"
        "internal/lso_output.hpp"
        "internal/ml_input.hpp"
        "internal/ml_output_brep.hpp"
        "internal/ml_output_impl.hpp"
        "internal/ml_output_structural_model.hpp"
        "internal/shp_input.hpp"
        "internal/horizons_stack_skua_input.hpp"
    PUBLIC_DEPENDENCIES
        OpenGeode::basic
    PRIVATE_DEPENDENCIES
        OpenGeode::geometry
        OpenGeode::mesh
        OpenGeode::model
        OpenGeode-Geosciences::explicit
        OpenGeode-Geosciences::implicit
        GDAL::GDAL
        pugixml::pugixml
        ${PROJECT_NAME}::mesh
)
