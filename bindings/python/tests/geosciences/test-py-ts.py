# -*- coding: utf-8 -*-
# Copyright (c) 2019 - 2020 Geode-solutions
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

import os, sys, platform
if sys.version_info >= (3,8,0) and platform.system() == "Windows":
    for path in [x.strip() for x in os.environ['PATH'].split(';') if x]:
        os.add_dll_directory(path)

import opengeode
import opengeode_geosciences as geosciences
import opengeode_geosciencesio_py_geosciences as geosciences_io

if __name__ != '__main__':
    geosciences_io.initialize_geosciences_io()
    test_dir = os.path.dirname(__file__)
    data_dir = os.path.abspath(os.path.join(test_dir, "../../../../tests/data"))

    surface = opengeode.TriangulatedSurface3D.create()
    opengeode.load_triangulated_surface( surface, os.path.join(data_dir, "surf2d.ts"))

    if surface.nb_vertices() != 46:
        raise ValueError("Number of vertices in the loaded TSurf 3D is not correct" )
    if surface.nb_polygons() != 46:
        raise ValueError("Number of polygons in the loaded TSurf 3D is not correct" )

    opengeode.save_triangulated_surface( surface, "surf3d.og_tsf3d" )
    reloaded_surface = opengeode.TriangulatedSurface3D.create()
    opengeode.load_triangulated_surface( reloaded_surface, "surf3d.og_tsf3d" )
    if reloaded_surface.nb_vertices() != 46:
        raise ValueError("Number of vertices in the reloaded TSurf 3D is not correct" )
    if reloaded_surface.nb_polygons() != 46:
        raise ValueError("Number of polygons in the reloaded TSurf 3D is not correct" )
