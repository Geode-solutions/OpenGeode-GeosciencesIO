# Copyright (c) 2019 - 2021 Geode-solutions
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from setuptools import setup
from os import path

with open(path.join('${CMAKE_SOURCE_DIR}', 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
with open('${CMAKE_CURRENT_LIST_DIR}/requirements.txt') as f:
    install_requires = f.read().strip().split('\n')

setup(
    name='OpenGeode-GeosciencesIO',
    version='${CPACK_PACKAGE_VERSION}',
    description='Input/Output formats for OpenGeode-Geosciences',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Geode-solutions',
    author_email='contact@geode-solutions.com',
    url='https://github.com/Geode-solutions/OpenGeode-GeosciencesIO',
    packages=['opengeode_geosciencesio'],
    package_data={
        '': ['*.so', '*.dll', '*.pyd', '*.dylib']
    },
    install_requires=install_requires,
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License'
    ],
    zip_safe=False
)