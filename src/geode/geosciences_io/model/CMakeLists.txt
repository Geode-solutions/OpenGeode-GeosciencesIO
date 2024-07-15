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
        "geos_export.cpp"
        "gocad_common.cpp"
        "lso_input.cpp"
        "lso_output.cpp"
        "ml_input.cpp"
        "ml_output_brep.cpp"
        "ml_output_structural_model.cpp"
        "shp_input.cpp"
        "horizons_stack_skua_input.cpp"
    PUBLIC_HEADERS
        "common.h"
    INTERNAL_HEADERS
        "helpers/brep_geos_export.h"
        "helpers/structural_model_geos_export.h"
        "internal/geos_export.h"
        "internal/gocad_common.h"
        "internal/lso_input.h"
        "internal/lso_output.h"
        "internal/ml_input.h"
        "internal/ml_output_brep.h"
        "internal/ml_output_impl.h"
        "internal/ml_output_structural_model.h"
        "internal/shp_input.h"
        "internal/horizons_stack_skua_input.h"
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
