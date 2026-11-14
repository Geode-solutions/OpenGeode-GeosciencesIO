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

"""The mesh formats of OpenGeode-GeosciencesIO (Gocad TS, VS, VO, PL, wells,
GeoTIFF, DEM, GRDECL, OMF, FEM) are read and written through the generic
functions of the opengeode package once opengeode_geosciencesio is imported,
and the errors are the OpenGeode ones."""

import os
import shutil
import tempfile

import geode_test_env

import opengeode
import opengeode_geosciencesio

DATA = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../../../../tests/data")
)


def check(condition, message):
    if not condition:
        raise ValueError("[Test] " + message)


def expect_data_error(function, *args):
    try:
        function(*args)
    except opengeode.OpenGeodeDataError as error:
        check(isinstance(error, RuntimeError), "Data errors are RuntimeError")
        return error
    raise ValueError("[Test] OpenGeodeDataError expected for " + repr(args))


def check_extensions(cls, read, written):
    inputs = opengeode.input_extensions(cls)
    outputs = opengeode.output_extensions(cls)
    for extension in read:
        check(extension in inputs, f"{cls.__name__} should read {extension}")
    for extension in written:
        check(extension in outputs, f"{cls.__name__} should write {extension}")


def test_extensions():
    check_extensions(opengeode.PointSet3D, ["vs", "omf"], ["vs"])
    check_extensions(
        opengeode.EdgedCurve3D,
        ["wl", "dat", "txt", "dev", "pl", "omf", "json"],
        ["pl"],
    )
    check_extensions(opengeode.TriangulatedSurface3D, ["ts", "omf"], ["ts"])
    check_extensions(opengeode.PolygonalSurface3D, ["dem", "tif", "tiff", "omf"], [])
    check_extensions(opengeode.TetrahedralSolid3D, [], ["fem"])
    check_extensions(opengeode.HybridSolid3D, ["grdecl"], [])
    check_extensions(opengeode.RegularGrid3D, ["vo"], [])
    check_extensions(opengeode.LightRegularGrid2D, ["tif", "tiff"], [])
    check_extensions(opengeode.LightRegularGrid3D, ["json", "omf"], [])
    # The formats of OpenGeode and OpenGeode-IO stay available
    check_extensions(opengeode.TriangulatedSurface3D, ["og_tsf3d", "stl"], ["obj"])
    # Formats are not shared between the types
    for cls in [opengeode.PointSet3D, opengeode.TetrahedralSolid3D]:
        check("ts" not in opengeode.input_extensions(cls), f"{cls.__name__} reads no ts")
    check("vs" not in opengeode.input_extensions(opengeode.EdgedCurve3D), "No vs curve")


def test_is_loadable():
    surface = opengeode.TriangulatedSurface3D
    check(
        opengeode.is_loadable(surface, os.path.join(DATA, "surf2d.ts")).value() == 1,
        "surf2d.ts should be loadable",
    )
    check(
        opengeode.is_loadable(surface, os.path.join(DATA, "missing.ts")).value() == 0,
        "A file that does not exist is not loadable",
    )
    check(
        opengeode.is_loadable(surface, os.path.join(DATA, "points.vs")).value() == 0,
        "A vs is not loadable as a surface",
    )
    check(
        opengeode.is_loadable(
            opengeode.PointSet3D, os.path.join(DATA, "surf2d.ts")
        ).value()
        == 0,
        "A ts is not loadable as a point set",
    )


def load(cls, filename):
    return opengeode.load(cls, os.path.join(DATA, filename))


def test_loading():
    cases = [
        (opengeode.TriangulatedSurface3D, "surf2d.ts", "nb_polygons", 46, 46),
        (opengeode.TriangulatedSurface3D, "omf_surfaces.omf", "nb_polygons", 2, 4),
        (opengeode.PointSet3D, "points.vs", "nb_vertices", 6, 6),
        (opengeode.PointSet3D, "omf_points.omf", "nb_vertices", 4, 4),
        (opengeode.EdgedCurve3D, "test.wl", "nb_edges", 54, 55),
        (opengeode.EdgedCurve3D, "test_well.dat", "nb_edges", 10, 11),
        (opengeode.EdgedCurve3D, "test_well.txt", "nb_edges", 10, 11),
        (opengeode.EdgedCurve3D, "test_well.dev", "nb_edges", 103, 104),
        (opengeode.PolygonalSurface3D, "bathy_IrishSea_DEM.dem", "nb_polygons", 453226, 456427),
        (opengeode.PolygonalSurface3D, "cea.tiff", "nb_polygons", 263682, 264710),
        (opengeode.PolygonalSurface3D, "omf_surfaces.omf", "nb_polygons", 2, 6),
        (opengeode.HybridSolid3D, "EclipseGridTest.grdecl", "nb_polyhedra", 24, 60),
        (opengeode.LightRegularGrid2D, "cea.tiff", "nb_cells", 264710, 264710),
        (opengeode.LightRegularGrid3D, "omf_volume.omf", "nb_cells", 24, 24),
    ]
    for cls, filename, name, nb_elements, nb_vertices in cases:
        path = os.path.join(DATA, filename)
        check(
            opengeode.is_loadable(cls, path).value() > 0,
            f"{filename} should be loadable as {cls.__name__}",
        )
        mesh = opengeode.load(cls, path)
        check(isinstance(mesh, cls), f"{filename} should be a {cls.__name__}")
        check(
            getattr(mesh, name)() == nb_elements,
            f"Wrong {name} for {filename}: {getattr(mesh, name)()}",
        )
        if hasattr(mesh, "nb_vertices"):
            check(mesh.nb_vertices() == nb_vertices, f"Wrong nb_vertices for {filename}")
    grid = load(opengeode.RegularGrid3D, "test.vo")
    check(grid.nb_cells() == 1000, "Wrong cells of test.vo")
    grid = load(opengeode.LightRegularGrid3D, "data.json")
    check(grid.nb_cells() == 4693, "Wrong cells of data.json")
    check(
        [grid.nb_cells_in_direction(i) for i in range(3)] == [19, 19, 13],
        "Wrong cells per direction of data.json",
    )
    # The .tif extension is registered as well as .tiff
    with tempfile.TemporaryDirectory() as directory:
        tif = os.path.join(directory, "cea.tif")
        shutil.copy(os.path.join(DATA, "cea.tiff"), tif)
        surface = opengeode.load(opengeode.PolygonalSurface3D, tif)
        check(
            (surface.nb_polygons(), surface.nb_vertices()) == (263682, 264710),
            "Wrong surface of cea.tif",
        )
        grid = opengeode.load(opengeode.LightRegularGrid2D, tif)
        check(
            [grid.nb_cells_in_direction(i) for i in range(2)] == [514, 515],
            "Wrong grid of cea.tif",
        )
    # A GeoTIFF without geotransform cannot be placed
    path = os.path.join(DATA, "bathy_IrishSea_GeoTIFF.tif")
    for cls in [opengeode.PolygonalSurface3D, opengeode.LightRegularGrid2D]:
        error = expect_data_error(opengeode.load, cls, path)
        check("geotransform" in error.message, f"geotransform error for {cls.__name__}")


def test_load_all():
    cases = [
        (opengeode.TriangulatedSurface3D, "2tsurfs.ts", [4, 4]),
        (opengeode.PointSet3D, "2vsets.vs", [6, 3]),
        (opengeode.EdgedCurve3D, "2plines.pl", [4, 3]),
        (opengeode.EdgedCurve3D, "2wells.wl", [55, 10]),
        (opengeode.EdgedCurve3D, "omf_lines.omf", [3, 2]),
        (opengeode.EdgedCurve3D, "data.json", [6, 5, 14, 8, 10, 6, 6, 7]),
    ]
    for cls, filename, sizes in cases:
        path = os.path.join(DATA, filename)
        meshes = opengeode.load_all(cls, path)
        check(
            [mesh.nb_vertices() for mesh in meshes] == sizes,
            f"Wrong objects of {filename}",
        )
        # load needs exactly one object
        expect_data_error(opengeode.load, cls, path)
    grids = opengeode.load_all(opengeode.RegularGrid3D, os.path.join(DATA, "2voxets.vo"))
    check([grid.nb_cells() for grid in grids] == [1000, 500], "Wrong voxets")


def test_saving(directory):
    surface = load(opengeode.TriangulatedSurface3D, "surf2d.ts")
    # The extension is case insensitive
    filename = os.path.join(directory, "UPPER.TS")
    check(opengeode.is_saveable(surface, filename), "ts should be saveable")
    check(opengeode.save(surface, filename) == [filename], "Wrong ts files")
    reloaded = opengeode.load(opengeode.TriangulatedSurface3D, filename)
    check(reloaded.nb_vertices() == 46, "ts reload: vertices")
    check(reloaded.nb_polygons() == 46, "ts reload: polygons")

    points = load(opengeode.PointSet3D, "points.vs")
    filename = os.path.join(directory, "points.vs")
    check(opengeode.save(points, filename) == [filename], "Wrong vs files")
    reloaded = opengeode.load(opengeode.PointSet3D, filename)
    check(reloaded.nb_vertices() == 6, "vs reload: vertices")
    for index in range(6):
        check(reloaded.point(index) == points.point(index), f"vs reload: point {index}")

    curve = load(opengeode.EdgedCurve3D, "test.wl")
    filename = os.path.join(directory, "curve.pl")
    check(opengeode.is_saveable(curve, filename), "pl should be saveable")
    check(opengeode.save(curve, filename) == [filename], "Wrong pl files")
    reloaded = opengeode.load(opengeode.EdgedCurve3D, filename)
    check(reloaded.nb_vertices() == 55, "pl reload: vertices")
    check(reloaded.nb_edges() == 54, "pl reload: edges")

    solid = load(opengeode.TetrahedralSolid3D, "bmsh_342.og_tso3d")
    filename = os.path.join(directory, "solid.fem")
    check(opengeode.is_saveable(solid, filename), "fem should be saveable")
    saved = opengeode.save(solid, filename)
    check(saved[0] == filename and os.path.getsize(filename) > 0, "Wrong fem files")
    # The header of the FEM file gives the number of nodes and of elements
    with open(filename) as file:
        lines = file.read().splitlines()
    dimens = lines[lines.index("DIMENS") + 1].split()
    check(
        (int(dimens[0]), int(dimens[1])) == (solid.nb_vertices(), solid.nb_polyhedra()),
        f"Wrong FEM sizes {dimens[:2]}",
    )
    varnode = lines.index("VARNODE")
    check(
        int(lines[varnode + 1].split()[0]) == solid.nb_polyhedra(),
        "Wrong number of FEM elements",
    )

    # Formats that are only read
    for obj, extension in [(surface, "omf"), (points, "omf"), (curve, "wl")]:
        name = os.path.join(directory, "none." + extension)
        check(not opengeode.is_saveable(obj, name), f"{extension} should not be saveable")
        expect_data_error(opengeode.save, obj, name)


def test_errors(directory):
    garbage = os.path.join(directory, "garbage.ts")
    with open(garbage, "w") as file:
        file.write("not a tsurf\n")
    expect_data_error(opengeode.load, opengeode.TriangulatedSurface3D, garbage)
    check(
        opengeode.is_loadable(opengeode.TriangulatedSurface3D, garbage).value() == 0,
        "garbage.ts is not loadable",
    )
    expect_data_error(
        opengeode.load,
        opengeode.TriangulatedSurface3D,
        os.path.join(directory, "missing.ts"),
    )
    error = expect_data_error(
        opengeode.load, opengeode.PointSet3D, os.path.join(DATA, "surf2d.ts")
    )
    check("Unknown extension" in error.message, "Unknown extension expected")


if __name__ == "__main__":
    test_extensions()
    test_is_loadable()
    test_loading()
    test_load_all()
    with tempfile.TemporaryDirectory() as directory:
        test_saving(directory)
        test_errors(directory)
