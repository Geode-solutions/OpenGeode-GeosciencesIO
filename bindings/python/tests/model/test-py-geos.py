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
import sys
import platform

if sys.version_info >= (3, 8, 0) and platform.system() == "Windows":
    for path in [x.strip() for x in os.environ["PATH"].split(";") if x]:
        os.add_dll_directory(path)

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


def add_block_property(model, property, name, value):
    attribute_id = opengeode.uuid()
    values = opengeode.AttributeValuesDouble()
    values.default_value = value
    for block in model.blocks():
        manager = block.mesh().polyhedron_attribute_manager()
        manager.create_attribute_variable_double(
            name, attribute_id, values, opengeode.AttributeProperties()
        )
    opengeode.BRepBuilder(model).set_physical_property(
        property, opengeode.Block3D.component_type_static(), attribute_id
    )


def surfaces_on_face(model, axis, side):
    return [
        surface
        for surface in model.surfaces()
        if abs(surface.mesh().bounding_box().min().value(axis) - side) < 1e-3
        and abs(surface.mesh().bounding_box().max().value(axis) - side) < 1e-3
    ]


# Gravity on a cube whose sides slide and whose bottom is fixed vertically
def test_mechanics_export(data_dir):
    model = opengeode.load_brep(
        os.path.join(data_dir, "adaptive_brep_perm_and_poro.og_brep")
    )
    for name, value in [
        ("bulk_modulus", 5.5556e9),
        ("shear_modulus", 4.16667e9),
        ("density", 2650.0),
    ]:
        add_block_property(
            model, getattr(opengeode.PHYSICAL_PROPERTY_NAME, name), name, value
        )
    free = float("nan")
    displacement_id = opengeode.uuid()
    for axis, side, displacement in [
        (0, 0.0, [0.0, free, free]),
        (0, 40.0, [0.0, free, free]),
        (1, 0.0, [free, 0.0, free]),
        (1, 40.0, [free, 0.0, free]),
        (2, 0.0, [free, free, 0.0]),
    ]:
        values = opengeode.AttributeValuesArrayDouble3()
        values.default_value = displacement
        for surface in surfaces_on_face(model, axis, side):
            manager = surface.mesh().polygon_attribute_manager()
            manager.create_attribute_constant_arraydouble3(
                "displacement",
                displacement_id,
                values,
                opengeode.AttributeProperties(),
            )
    opengeode.BRepBuilder(model).set_physical_property(
        opengeode.PHYSICAL_PROPERTY_NAME.boundary_displacement,
        opengeode.Surface3D.component_type_static(),
        displacement_id,
    )
    stress_id = opengeode.uuid()
    stress = opengeode.AttributeValuesArrayDouble3()
    stress.default_value = [11.25e6, 9e6, 15e6]
    for block in model.blocks():
        manager = block.mesh().polyhedron_attribute_manager()
        manager.create_attribute_constant_arraydouble3(
            "initial_stress", stress_id, stress, opengeode.AttributeProperties()
        )
    opengeode.BRepBuilder(model).set_physical_property(
        opengeode.PHYSICAL_PROPERTY_NAME.initial_stress,
        opengeode.Block3D.component_type_static(),
        stress_id,
    )
    problem = geosciences.PhysicalProblem()
    builder = geosciences.PhysicalProblemBuilder(problem)
    builder.set_name("gravity")
    builder.set_model_id(model.id())
    builder.set_solid_mechanics(geosciences.SolidMechanics())
    builder.set_schedule(geosciences.Schedule(0.0, 1.0, 1.0))
    exporter = geosciences_io.BRepGeosExporter(model, "python_gravity")
    exporter.set_physical_problem(problem)
    exporter.run()
    with open("python_gravity/python_gravity_simulation.xml") as deck:
        if "solid_stress" not in deck.read():
            raise ValueError("[Test] Mechanics deck misses the initial stress")


def test_structural_model_export(data_dir):
    model = geosciences.load_structural_model(os.path.join(data_dir, "test.lso"))
    geosciences_io.BRepGeosExporter(model, "brep_to_geos")
    geosciences_io.StructuralModelGeosExporter(model, "SM_to_geos")


if __name__ == "__main__":
    geosciences_io.OpenGeodeGeosciencesIOModelLibrary.initialize()
    test_dir = os.path.dirname(__file__)
    data_dir = os.path.abspath(os.path.join(test_dir, "../../../../tests/data"))
    test_brep_export(data_dir)
    test_mechanics_export(data_dir)
    test_structural_model_export(data_dir)
