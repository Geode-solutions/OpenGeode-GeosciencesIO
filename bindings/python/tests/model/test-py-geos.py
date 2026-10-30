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

import os

import geode_test_env

import opengeode
import opengeode_io
import opengeode_geosciences as geosciences
import opengeode_geosciencesio_py_model as geosciences_io


# The GEOS results come from a picasso simulation run outside the test
def test_brep_export(data_dir):
    model = opengeode.load_brep(os.path.join(data_dir, "picasso.og_brep"))
    problem = geosciences.PhysicalProblem()
    builder = geosciences.PhysicalProblemBuilder(problem)
    builder.set_name("picasso")
    builder.set_model_id(model.id())
    fluid = geosciences.SinglePhaseFluid(
        [geosciences.PHASE.water],
        [geosciences.FluidComponent("water", 0.018)],
        reference_density=1000.0,
        reference_viscosity=1e-3,
        compressibility=1e-9,
        viscosibility=0.0,
        reference_pressure=1e7,
    )
    flow = geosciences.SinglePhaseFlow()
    flow.set_fluid(fluid)
    builder.set_flow(flow)
    rock = geosciences.Material()
    rock.name = "rock"
    rock.pore_compressibility = 1e-9
    rock.reference_pressure = 1e7
    rock.porosity = 0.2
    rock.permeability = [1e-13, 1e-13, 1e-14]
    builder.add_material(rock)
    builder.add_initial_condition(geosciences.UniformInitialCondition([], 1e7))
    builder.set_schedule(geosciences.Schedule(0.0, 1e8, 1e7))
    output = geosciences.OutputSettings()
    output.frequency = 5e7
    builder.set_output(output)
    exporter = geosciences_io.BRepGeosExporter(model, "python_geos")
    exporter.set_physical_problem(problem)
    exporter.run()
    opengeode.load_brep_time_series(
        model, os.path.join(data_dir, "picasso_geos", "picasso_results.pvd")
    )


def test_structural_model_export(data_dir):
    model = geosciences.load_structural_model(os.path.join(data_dir, "test.lso"))
    geosciences_io.BRepGeosExporter(model, "brep_to_geos")
    geosciences_io.StructuralModelGeosExporter(model, "SM_to_geos")


if __name__ == "__main__":
    geosciences_io.OpenGeodeGeosciencesIOModelLibrary.initialize()
    test_dir = os.path.dirname(__file__)
    data_dir = os.path.abspath(os.path.join(test_dir, "../../../../tests/data"))
    test_brep_export(data_dir)
    test_structural_model_export(data_dir)
