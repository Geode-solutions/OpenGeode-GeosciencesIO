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

"""The GEOS exporters of OpenGeode-GeosciencesIO: keyword arguments, lifetime
of the model and of the problem, errors of the OpenGeode hierarchy, written
files, and the GEOS simulation checker."""

import gc
import os
import tempfile

import geode_test_env

import opengeode
import opengeode_geosciences as geosciences
import opengeode_geosciencesio as geosciences_io

DATA = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../../../../tests/data")
)

EXPORTERS = [
    (geosciences_io.BRepGeosExporter, "BRepGeosExporter"),
    (geosciences_io.StructuralModelGeosExporter, "StructuralModelGeosExporter"),
]


def check(condition, message):
    if not condition:
        raise ValueError("[Test] " + message)


def expect_data_error(function, *args, **kwargs):
    try:
        function(*args, **kwargs)
    except opengeode.OpenGeodeDataError as error:
        check(isinstance(error, RuntimeError), "Data errors are RuntimeError")
        return error
    raise ValueError("[Test] OpenGeodeDataError expected")


def load_brep():
    return opengeode.load(opengeode.BRep, os.path.join(DATA, "picasso.og_brep"))


def load_structural_model():
    return opengeode.load(
        geosciences.StructuralModel, os.path.join(DATA, "test.lso")
    )


def create_well():
    well = opengeode.PointSet3D.create()
    opengeode.PointSetBuilder3D.create(well).create_point(
        opengeode.Point3D([1.0, 1.0, 1.0])
    )
    return well


def test_documentation():
    for exporter, name in EXPORTERS:
        check(exporter.__doc__, f"No docstring for {name}")
        check("files_directory" in exporter.__init__.__doc__, f"{name}: init arguments")
        for method in ["add_well_perforations", "set_physical_problem", "run"]:
            check(getattr(exporter, method).__doc__, f"No docstring for {name}.{method}")
        check(
            "well_perforations" in exporter.add_well_perforations.__doc__
            and "problem" in exporter.set_physical_problem.__doc__,
            f"{name}: argument names",
        )
        try:
            exporter()
        except TypeError:
            pass
        else:
            raise ValueError(f"[Test] {name} needs a model")


def test_keyword_arguments(directory):
    exporter = geosciences_io.BRepGeosExporter(
        model=load_brep(), files_directory=os.path.join(directory, "kw")
    )
    exporter.add_well_perforations(well_perforations=create_well(), name="well")
    exporter.set_physical_problem(problem=geosciences.PhysicalProblem())
    # An empty problem cannot be simulated: the checks of the problem fail
    error = expect_data_error(exporter.run)
    check("Cannot write the GEOS simulation deck" in error.message, "Deck error")
    check("no flow model nor solid mechanics" in error.message, "Problem error")


def test_directory_errors(directory):
    for exporter, name in EXPORTERS:
        model = load_brep() if name.startswith("BRep") else load_structural_model()
        for bad in ["", "/"]:
            error = expect_data_error(exporter, model, bad)
            check("must not be empty" in error.message, f"{name}: message for {bad!r}")
        directory_name = os.path.join(directory, "created_" + name)
        exporter(model, directory_name)
        check(os.path.isdir(directory_name), f"{name} should create its directory")


def test_well_errors(directory):
    for exporter_class, name in EXPORTERS:
        model = load_brep() if name.startswith("BRep") else load_structural_model()
        exporter = exporter_class(model, os.path.join(directory, "wells_" + name))
        well = create_well()
        for bad in ["", "a b", "a,b", "a{b"]:
            error = expect_data_error(exporter.add_well_perforations, well, bad)
            check("no space, comma nor brace" in error.message, f"{name}: well {bad!r}")
        exporter.add_well_perforations(well, "well")
        error = expect_data_error(exporter.add_well_perforations, well, "well")
        check("already added" in error.message, f"{name}: duplicated well")


def test_lifetime(directory):
    def make(directory_name):
        exporter = geosciences_io.BRepGeosExporter(
            load_brep(), os.path.join(directory, directory_name)
        )
        problem = geosciences.PhysicalProblem()
        geosciences.PhysicalProblemBuilder(problem).set_name("lifetime")
        exporter.set_physical_problem(problem)
        return exporter

    exporter = make("lifetime")
    gc.collect()
    # The model and the problem are kept alive by the exporter: the problem
    # is read, and reported as incomplete, instead of a dangling reference
    error = expect_data_error(exporter.run)
    check("no flow model nor solid mechanics" in error.message, "Problem error")
    without_problem = geosciences_io.BRepGeosExporter(
        load_brep(), os.path.join(directory, "lifetime_brep")
    )
    gc.collect()
    without_problem.run()
    check(
        os.path.exists(os.path.join(directory, "lifetime_brep", "lifetime_brep.vtu")),
        "The mesh of the exporter should be written",
    )
    structural = geosciences_io.StructuralModelGeosExporter(
        load_structural_model(), os.path.join(directory, "lifetime_structural")
    )
    gc.collect()
    structural.run()
    check(
        os.path.exists(
            os.path.join(directory, "lifetime_structural", "lifetime_structural.vtu")
        ),
        "The mesh of the structural exporter should be written",
    )


def test_run(directory):
    for exporter_class, name in EXPORTERS:
        model = load_brep() if name.startswith("BRep") else load_structural_model()
        folder = os.path.join(directory, "run_" + name)
        exporter = exporter_class(model, folder)
        exporter.add_well_perforations(create_well(), "w1")
        check(exporter.run() is None, f"{name}.run returns nothing")
        prefix = os.path.basename(folder)
        files = sorted(os.listdir(folder))
        check(prefix + ".vtu" in files, f"{name}: mesh missing in {files}")
        check(prefix + "_simulation.xml" in files, f"{name}: deck missing in {files}")
        check(prefix + "_well0.vtp" in files, f"{name}: well missing in {files}")
        # run can be called again
        exporter.run()
        check(sorted(os.listdir(folder)) == files, f"{name}: files after second run")


def test_gil_release(directory):
    # The callables repeated by assert_gil_released contain only the call
    # under test: the exporters and the checker are built before
    for exporter_class, name in EXPORTERS:
        model = load_brep() if name.startswith("BRep") else load_structural_model()
        folder = os.path.join(directory, "gil_" + name)
        exporter = exporter_class(model, folder)
        exporter.add_well_perforations(create_well(), "w1")
        geode_test_env.assert_gil_released(
            exporter.run, deadline=60.0, name=name + ".run"
        )
        check(
            os.path.basename(folder) + ".vtu" in os.listdir(folder),
            f"{name}: mesh missing after the released run",
        )
    checker = geosciences.SimulationCheckerFactory.create(
        "geos", geosciences.PhysicalProblem(), load_brep()
    )
    report = geode_test_env.assert_gil_released(
        checker.check, deadline=60.0, name="GEOS SimulationChecker.check"
    )
    check(not report.is_simulable(), "The released check returns its report")


def test_simulation_checker():
    problem = geosciences.PhysicalProblem()
    brep = load_brep()
    check(
        geosciences.SimulationCheckerFactory.has_creator("geos"),
        "The GEOS checker should be registered",
    )
    report = geosciences.SimulationCheckerFactory.create("geos", problem, brep).check()
    check(isinstance(report, geosciences.SimulationReport), "check returns a report")
    check(not report.is_simulable(), "An empty problem is not simulable with GEOS")
    check(report.nb_errors() > 0, "An empty problem has errors")


def test_initialize():
    for library in [
        geosciences_io.OpenGeodeGeosciencesIOMeshLibrary,
        geosciences_io.OpenGeodeGeosciencesIOModelLibrary,
    ]:
        name = library.__name__
        check(library.__doc__, f"No docstring for {name}")
        check(library.initialize.__doc__, f"No docstring for {name}.initialize")
        check(
            isinstance(library.__dict__["initialize"], staticmethod),
            f"{name}.initialize should be static",
        )
        check(
            "Initialize self" not in (library.__init__.__doc__ or ""),
            f"{name} should document its missing constructor",
        )
        try:
            library()
        except TypeError:
            pass
        else:
            raise ValueError(f"[Test] {name} should not be constructible")
        library.initialize()
        library.initialize()


if __name__ == "__main__":
    test_documentation()
    test_initialize()
    test_simulation_checker()
    with tempfile.TemporaryDirectory() as directory:
        test_keyword_arguments(directory)
        test_directory_errors(directory)
        test_well_errors(directory)
        test_lifetime(directory)
        test_run(directory)
        test_gil_release(directory)
