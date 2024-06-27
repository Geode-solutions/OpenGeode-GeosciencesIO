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
    NAME mesh
    FOLDER "geode/geosciences_io/mesh"
    SOURCES
        "common.cpp"
        "gocad_common.cpp"
        "grdecl_input.cpp"
        "pl_input.cpp"
        "pl_output.cpp"
        "ts_input.cpp"
        "ts_output.cpp"
        "vo_input.cpp"
        "wl_input.cpp"
        "well_dat_input.cpp"
        "well_dev_input.cpp"
        "well_txt_input.cpp"
    PUBLIC_HEADERS
        "common.h"
    PRIVATE_HEADERS
        "private/gocad_common.h"
        "private/grdecl_input.h"
        "private/pl_input.h"
        "private/pl_output.h"
        "private/ts_input.h"
        "private/ts_output.h"
        "private/vo_input.h"
        "private/well_input.h"
        "private/well_dat_input.h"
        "private/well_dev_input.h"
        "private/well_txt_input.h"
        "private/wl_input.h"
    PUBLIC_DEPENDENCIES
        OpenGeode::basic
    PRIVATE_DEPENDENCIES
        OpenGeode::geometry
        OpenGeode::mesh
)
