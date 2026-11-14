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

"""The models of OpenGeode-GeosciencesIO (Gocad ML and LSO, Shapefile, Skua
horizons stack, FEM) are read and written through the generic functions of
the opengeode package once opengeode_geosciencesio is imported, for the types
registered by OpenGeode-Geosciences, and the errors are the OpenGeode ones."""

import os
import tempfile

import geode_test_env

import opengeode
import opengeode_geosciences as geosciences
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


def sizes(model):
    return (
        model.nb_corners(),
        model.nb_lines(),
        model.nb_surfaces(),
        model.nb_blocks(),
        model.nb_horizons(),
    )


def test_extensions():
    structural = geosciences.StructuralModel
    for extension in ["ml", "lso", "og_strm"]:
        check(
            extension in opengeode.input_extensions(structural),
            f"StructuralModel should read {extension}",
        )
        check(
            extension in opengeode.output_extensions(structural),
            f"StructuralModel should write {extension}",
        )
    check(
        "ml" not in opengeode.input_extensions(geosciences.CrossSection),
        "CrossSection reads no ml",
    )
    for extension in ["shp", "shz", "svg", "og_sctn"]:
        check(
            extension in opengeode.input_extensions(opengeode.Section),
            f"Section should read {extension}",
        )
    check("shp" not in opengeode.output_extensions(opengeode.Section), "Section writes no shp")
    for extension in ["ml", "fem", "vtm", "og_brep"]:
        check(
            extension in opengeode.output_extensions(opengeode.BRep),
            f"BRep should write {extension}",
        )
    check("ml" not in opengeode.input_extensions(opengeode.BRep), "BRep reads no ml")
    for cls in [geosciences.HorizonsStack2D, geosciences.HorizonsStack3D]:
        check("xml" in opengeode.input_extensions(cls), f"{cls.__name__} should read xml")
    check(
        "geos" in geosciences.SimulationCheckerFactory.list_creators(),
        "The GEOS simulation checker should be registered",
    )


def test_structural_model():
    path = os.path.join(DATA, "test.lso")
    check(
        opengeode.is_loadable(geosciences.StructuralModel, path).value() == 1,
        "test.lso should be loadable",
    )
    model = opengeode.load(geosciences.StructuralModel, path)
    check(isinstance(model, geosciences.StructuralModel), "test.lso is a StructuralModel")
    check(sizes(model) == (22, 39, 23, 4, 4), f"Wrong test.lso: {sizes(model)}")
    check(
        opengeode.is_loadable(
            geosciences.StructuralModel, os.path.join(DATA, "2models.ml")
        ).value()
        > 0,
        "2models.ml should be loadable",
    )
    models = opengeode.load_all(
        geosciences.StructuralModel, os.path.join(DATA, "2models.ml")
    )
    check(
        [(m.nb_corners(), m.nb_lines(), m.nb_surfaces()) for m in models]
        == [(2, 2, 1), (8, 12, 8)],
        "Wrong models of 2models.ml",
    )
    expect_data_error(
        opengeode.load, geosciences.StructuralModel, os.path.join(DATA, "2models.ml")
    )
    models = opengeode.load_all(
        geosciences.StructuralModel, os.path.join(DATA, "2lightsolids.lso")
    )
    check(len(models) == 2, "2lightsolids.lso holds two models")


def test_section_and_horizons_stack():
    section = opengeode.load(opengeode.Section, os.path.join(DATA, "test.shp"))
    check(section.nb_surfaces() == 1, "test.shp holds one surface")
    check(
        opengeode.is_loadable(
            opengeode.Section, os.path.join(DATA, "test.shp")
        ).value()
        == 1,
        "test.shp should be loadable",
    )
    zipped = opengeode.load(opengeode.Section, os.path.join(DATA, "toulouse.shz"))
    check(
        (zipped.nb_corners(), zipped.nb_lines(), zipped.nb_surfaces())
        == (51718, 10, 12),
        "Wrong toulouse.shz",
    )
    path = os.path.join(DATA, "test_skua_horizons_stack.xml")
    stack = opengeode.load(geosciences.HorizonsStack3D, path)
    check(isinstance(stack, geosciences.HorizonsStack3D), "xml is a HorizonsStack3D")
    check(stack.nb_horizons() == 4, "Wrong horizons of the Skua stack")
    check(
        stack.nb_stratigraphic_units() == 5,
        "Wrong stratigraphic units of the Skua stack",
    )
    check(
        opengeode.is_loadable(geosciences.HorizonsStack3D, path).value() == 1,
        "The Skua stack should be loadable",
    )


def test_saving(directory):
    model = opengeode.load(
        geosciences.StructuralModel, os.path.join(DATA, "test.lso")
    )
    for extension in ["ml", "lso"]:
        filename = os.path.join(directory, "model." + extension)
        check(
            opengeode.is_saveable(model, filename), f"{extension} should be saveable"
        )
        check(opengeode.save(model, filename) == [filename], f"Wrong {extension} files")
        reloaded = opengeode.load(geosciences.StructuralModel, filename)
        check(
            (reloaded.nb_blocks(), reloaded.nb_surfaces()) == (4, 23),
            f"{extension} reload: {sizes(reloaded)}",
        )
    # The extension is case insensitive
    upper = os.path.join(directory, "UPPER.ML")
    opengeode.save(model, upper)
    check(
        opengeode.load(geosciences.StructuralModel, upper).nb_blocks() == 4,
        "UPPER.ML reload",
    )

    brep = opengeode.load(opengeode.BRep, os.path.join(DATA, "picasso.og_brep"))
    for extension in ["ml", "fem"]:
        filename = os.path.join(directory, "brep." + extension)
        check(opengeode.is_saveable(brep, filename), f"BRep saveable as {extension}")
        check(opengeode.save(brep, filename) == [filename], f"Wrong {extension} files")
        check(os.path.getsize(filename) > 0, f"{filename} is empty")
    # Only the structural model is read from a ml
    expect_data_error(
        opengeode.load, opengeode.BRep, os.path.join(directory, "brep.ml")
    )
    # A cross section has no lso format
    cross_section = geosciences.CrossSection()
    check(
        not opengeode.is_saveable(cross_section, os.path.join(directory, "x.lso")),
        "A CrossSection is not saveable as lso",
    )


def test_errors(directory):
    garbage = os.path.join(directory, "garbage.lso")
    with open(garbage, "w") as file:
        file.write("not a lightsolid\n")
    expect_data_error(opengeode.load, geosciences.StructuralModel, garbage)
    check(
        opengeode.is_loadable(geosciences.StructuralModel, garbage).value() == 0,
        "garbage.lso is not loadable",
    )
    expect_data_error(
        opengeode.load,
        geosciences.StructuralModel,
        os.path.join(directory, "missing.ml"),
    )
    error = expect_data_error(
        opengeode.load, geosciences.StructuralModel, os.path.join(DATA, "test.shp")
    )
    check("Unknown extension" in error.message, "Unknown extension expected")
    expect_data_error(
        opengeode.load, opengeode.Section, os.path.join(directory, "missing.shp")
    )
    expect_data_error(
        opengeode.load,
        geosciences.HorizonsStack3D,
        os.path.join(directory, "missing.xml"),
    )


if __name__ == "__main__":
    test_extensions()
    test_structural_model()
    test_section_and_horizons_stack()
    with tempfile.TemporaryDirectory() as directory:
        test_saving(directory)
        test_errors(directory)
