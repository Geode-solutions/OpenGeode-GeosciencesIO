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

add_geode_test(
    SOURCE "test-ml.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::model
        OpenGeode-Geosciences::explicit
        ${PROJECT_NAME}::model
    ESSENTIAL
)
add_geode_test(
    SOURCE "test-lso.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::model
        OpenGeode-Geosciences::explicit
        ${PROJECT_NAME}::model
)
add_geode_test(
    SOURCE "test-shp.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::model
        ${PROJECT_NAME}::model
)
add_geode_test(
    SOURCE "test-skua-horizons-stack.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::model
        OpenGeode-Geosciences::explicit
        OpenGeode-Geosciences::implicit
        ${PROJECT_NAME}::model
)
