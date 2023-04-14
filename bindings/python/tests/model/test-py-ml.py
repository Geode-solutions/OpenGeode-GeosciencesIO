# -*- coding: utf-8 -*-
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

import os
import sys
import platform
if sys.version_info >= (3, 8, 0) and platform.system() == "Windows":
    for path in [x.strip() for x in os.environ['PATH'].split(';') if x]:
        os.add_dll_directory(path)

import opengeode
import opengeode_geosciences as geosciences
import opengeode_geosciencesio_py_model as geosciences_io

if __name__ == '__main__':
    geosciences_io.GeosciencesIOModelLibrary.initialize()
    test_dir = os.path.dirname(__file__)
    data_dir = os.path.abspath(os.path.join(
        test_dir, "../../../../tests/data"))

    model = geosciences.load_structural_model(
        os.path.join(data_dir, "modelA4.ml"))

    if model.nb_corners() != 52:
        raise ValueError(
            "[Test] Number of Corners in the loaded StructuralModel is not correct")
    if model.nb_lines() != 98:
        raise ValueError(
            "[Test] Number of Lines in the loaded StructuralModel is not correct")
    if model.nb_surfaces() != 55:
        raise ValueError(
            "[Test] Number of Surfaces in the loaded StructuralModel is not correct")
    if model.nb_blocks() != 8:
        raise ValueError(
            "[Test] Number of Blocks in the loaded StructuralModel is not correct")
    if model.nb_faults() != 2:
        raise ValueError(
            "[Test] Number of Faults in the loaded StructuralModel is not correct")
    if model.nb_horizons() != 3:
        raise ValueError(
            "[Test] Number of Horizons in the loaded StructuralModel is not correct")
    if model.nb_model_boundaries() != 6:
        raise ValueError(
            "[Test] Number of ModelBoundary in the loaded StructuralModel is not correct")

    geosciences.save_structural_model(model, "modelA4.og_strm")
    geosciences.save_structural_model(model, "modelA4_saved.ml")
    reload = geosciences.load_structural_model("modelA4_saved.ml")

    if reload.nb_corners() != 52:
        raise ValueError(
            "[Test] Number of Corners in the reloaded StructuralModel is not correct")
    if reload.nb_lines() != 98:
        raise ValueError(
            "[Test] Number of Lines in the reloaded StructuralModel is not correct")
    if reload.nb_surfaces() != 55:
        raise ValueError(
            "[Test] Number of Surfaces in the reloaded StructuralModel is not correct")
    if reload.nb_blocks() != 8:
        raise ValueError(
            "[Test] Number of Blocks in the reloaded StructuralModel is not correct")
    if reload.nb_faults() != 2:
        raise ValueError(
            "[Test] Number of Faults in the reloaded StructuralModel is not correct")
    if reload.nb_horizons() != 3:
        raise ValueError(
            "[Test] Number of Horizons in the reloaded StructuralModel is not correct")
    if reload.nb_model_boundaries() != 6:
        raise ValueError(
            "[Test] Number of ModelBoundary in the reloaded StructuralModel is not correct")
