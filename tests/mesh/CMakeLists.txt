# Copyright (c) 2019 - 2025 Geode-solutions
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

find_package(GDAL REQUIRED CONFIG)
add_geode_test(
    SOURCE "test-solid-fem-output.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-dem.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        OpenGeode-IO::mesh
        GDAL::GDAL
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-geotiff.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        OpenGeode-IO::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-grdecl.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh

        ${PROJECT_NAME}::mesh
    ESSENTIAL
)
add_geode_test(
    SOURCE "test-pl.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-ts.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-vo.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-vs.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-well-dat.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-well-dev.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-well-txt.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)
add_geode_test(
    SOURCE "test-wl.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::mesh
        ${PROJECT_NAME}::mesh
)

