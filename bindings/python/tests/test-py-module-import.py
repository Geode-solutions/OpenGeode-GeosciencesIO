# -*- coding: utf-8 -*-
# Copyright (c) 2019 - 2026 Geode-solutions
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

"""Each raw extension module imports alone, in a fresh interpreter: it imports
the Python packages of the types and libraries it uses by itself."""

import subprocess
import sys

import geode_test_env


def check(condition, message):
    if not condition:
        raise ValueError("[Test] " + message)


def import_alone(module, statements):
    code = f"import {module}\n" + statements
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    check(
        result.returncode == 0,
        f"import {module} alone failed:\n{result.stdout}\n{result.stderr}",
    )


if __name__ == "__main__":
    import_alone(
        "opengeode_geosciencesio_py_mesh",
        "import sys\n"
        "for name in ['opengeode', 'opengeode_io', 'opengeode_geosciences']:\n"
        "    assert name in sys.modules, name\n"
        "opengeode_geosciencesio_py_mesh.OpenGeodeGeosciencesIOMeshLibrary"
        ".initialize()\n",
    )
    import_alone(
        "opengeode_geosciencesio_py_model",
        "import sys\n"
        "for name in ['opengeode', 'opengeode_io', 'opengeode_geosciences']:\n"
        "    assert name in sys.modules, name\n"
        "m = opengeode_geosciencesio_py_model\n"
        "m.OpenGeodeGeosciencesIOModelLibrary.initialize()\n"
        "assert m.BRepGeosExporter.__doc__\n",
    )
