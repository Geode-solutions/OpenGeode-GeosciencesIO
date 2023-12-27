# Copyright (c) 2019 - 2023 Geode-solutions
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
        "geos_output_brep.cpp"
        "geos_output_structural_model.cpp"
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
    PRIVATE_HEADERS
        "private/geos_output.h"
        "private/geos_output_brep.h"
        "private/geos_output_structural_model.h"
        "private/gocad_common.h"
        "private/lso_input.h"
        "private/lso_output.h"
        "private/ml_input.h"
        "private/ml_output_brep.h"
        "private/ml_output_impl.h"
        "private/ml_output_structural_model.h"
        "private/shp_input.h"
        "private/horizons_stack_skua_input.h"
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
