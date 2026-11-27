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

"""Export of picasso with a single-phase problem to OpenGeoSys: the project
file names the LIQUID_FLOW process and the written meshes; with thermal
parameters, the HT process and its temperature variable; with solid
mechanics and heat sources, the THERMO_HYDRO_MECHANICS process and its
Volumetric and Nodal temperature source terms."""

import os
import tempfile
import xml.etree.ElementTree as ElementTree

import geode_test_env

import opengeode
import opengeode_io
import opengeode_geosciences as geosciences
import opengeode_geosciencesio_py_model as geosciences_io


def check(condition, message):
    if not condition:
        raise ValueError("[Test] " + message)


def picasso_problem(model):
    problem = geosciences.PhysicalProblem()
    builder = geosciences.PhysicalProblemBuilder(problem)
    builder.set_name("picasso")
    builder.set_model_id(model.id())
    fluid = geosciences.SinglePhaseFluid(
        [geosciences.Phase.water],
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
    return problem


def thermal_problem(model):
    """Thermal single-phase problem (OGS HT): the values of the C++
    test_ogs_ht_project"""
    problem = geosciences.PhysicalProblem()
    builder = geosciences.PhysicalProblemBuilder(problem)
    builder.set_name("picasso_thermal")
    builder.set_model_id(model.id())
    fluid_thermal = geosciences.FluidThermalProperties()
    fluid_thermal.specific_heat_capacity = 4182.0
    fluid_thermal.thermal_expansion = 0.0
    fluid_thermal.viscosity_temperature_coefficient = 0.0
    fluid = geosciences.SinglePhaseFluid(
        [geosciences.Phase.water],
        [geosciences.FluidComponent("water", 0.018)],
        reference_density=1000.0,
        reference_viscosity=1e-3,
        compressibility=0.0,
        viscosibility=0.0,
        reference_pressure=1e5,
        thermal=fluid_thermal,
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
    rock_thermal = geosciences.RockThermalProperties()
    rock_thermal.volumetric_heat_capacity = 2e6
    rock_thermal.thermal_conductivity = [2.0, 2.0, 3.0]
    rock.thermal = rock_thermal
    builder.add_material(rock)
    builder.add_initial_condition(
        geosciences.UniformInitialCondition([], 1e7, temperature=300.0)
    )
    builder.set_thermal(geosciences.ThermalParameters(300.0))
    builder.set_schedule(geosciences.Schedule(0.0, 1e5, 1e3, 1e4))
    output = geosciences.OutputSettings()
    output.frequency = 5e4
    builder.set_output(output)
    return problem


def add_east_pressure(model):
    side = model.bounding_box().max().value(0)
    attribute_id = opengeode.uuid()
    values = opengeode.AttributeValuesDouble()
    values.default_value = 1e7
    for surface in model.surfaces():
        box = surface.mesh().bounding_box()
        if abs(box.min().value(0) - side) < 1e-3 and abs(box.max().value(0) - side) < 1e-3:
            surface.mesh().polygon_attribute_manager().create_attribute_constant_double(
                "boundary_pressure", attribute_id, values, opengeode.AttributeProperties()
            )
    opengeode.BRepBuilder(model).set_physical_property(
        opengeode.PhysicalPropertyName.boundary_pressure,
        opengeode.Surface3D.component_type_static(),
        attribute_id,
    )


def test_brep_export(data_dir, directory):
    model_filename = os.path.join(data_dir, "picasso.og_brep")
    model = opengeode.load_one(opengeode.load_brep(model_filename), model_filename)
    add_east_pressure(model)
    problem = picasso_problem(model)
    folder = os.path.join(directory, "python_ogs")
    exporter = geosciences_io.BRepOgsExporter(model, folder)
    exporter.set_physical_problem(problem)
    exporter.run()
    project = ElementTree.parse(os.path.join(folder, "python_ogs.prj")).getroot()
    check(
        project.findtext("processes/process/type") == "LIQUID_FLOW",
        "The project should use LIQUID_FLOW",
    )
    meshes = [mesh.text for mesh in project.findall("meshes/mesh")]
    check(meshes[0] == "python_ogs.vtu", "The bulk mesh comes first")
    check(len(meshes) == 2, f"One boundary mesh expected: {meshes}")
    for mesh in meshes:
        check(os.path.exists(os.path.join(folder, mesh)), f"{mesh} should be written")


def test_thermal_export(data_dir, directory):
    model_filename = os.path.join(data_dir, "picasso.og_brep")
    model = opengeode.load_one(opengeode.load_brep(model_filename), model_filename)
    add_east_pressure(model)
    folder = os.path.join(directory, "python_ogs_ht")
    exporter = geosciences_io.BRepOgsExporter(model, folder)
    exporter.set_physical_problem(thermal_problem(model))
    exporter.run()
    project = ElementTree.parse(os.path.join(folder, "python_ogs_ht.prj")).getroot()
    process_type = project.findtext("processes/process/type")
    check(process_type == "HT", f"The project should use HT, not {process_type}")
    temperature = project.find(
        "process_variables/process_variable[name='temperature']"
    )
    check(temperature is not None, "The temperature process variable is missing")
    variables = project.find("processes/process/process_variables")
    check(
        [child.tag for child in variables] == ["temperature", "pressure"],
        "The HT process variables are temperature, then pressure",
    )


FREE = float("nan")


def add_block_property(model, name, value):
    attribute_id = opengeode.uuid()
    values = opengeode.AttributeValuesDouble()
    values.default_value = value
    for block in model.blocks():
        block.mesh().polyhedron_attribute_manager().create_attribute_constant_double(
            name, attribute_id, values, opengeode.AttributeProperties()
        )
    opengeode.BRepBuilder(model).set_physical_property(
        getattr(opengeode.PhysicalPropertyName, name),
        opengeode.Block3D.component_type_static(),
        attribute_id,
    )


def add_rollers(model):
    """Lateral faces slide, the bottom face is fixed vertically"""
    box = model.bounding_box()
    attribute_id = opengeode.uuid()
    sides = [(axis, side) for axis in (0, 1) for side in (box.min(), box.max())]
    sides.append((2, box.min()))
    for surface in model.surfaces():
        surface_box = surface.mesh().bounding_box()
        for axis, side in sides:
            value = side.value(axis)
            if (
                abs(surface_box.min().value(axis) - value) < 1e-3
                and abs(surface_box.max().value(axis) - value) < 1e-3
            ):
                displacement = [FREE, FREE, FREE]
                displacement[axis] = 0.0
                values = opengeode.AttributeValuesArrayDouble3()
                values.default_value = displacement
                surface.mesh().polygon_attribute_manager().create_attribute_constant_arraydouble3(
                    "boundary_displacement", attribute_id, values, opengeode.AttributeProperties()
                )
    opengeode.BRepBuilder(model).set_physical_property(
        opengeode.PhysicalPropertyName.boundary_displacement,
        opengeode.Surface3D.component_type_static(),
        attribute_id,
    )


def thermo_hydro_mechanics_problem(model, block_id):
    """The HT problem with solid mechanics (the Biot coefficient gives the
    pore compressibility), a Block heat source and a point one"""
    problem = thermal_problem(model)
    builder = geosciences.PhysicalProblemBuilder(problem)
    builder.set_solid_mechanics(geosciences.SolidMechanics())
    rock = problem.materials()[0]
    rock.pore_compressibility = 0.0
    rock.reference_pressure = 0.0
    builder.set_material(material_id=0, material=rock)
    builder.add_heat_source(
        geosciences.HeatSource(name="canister", block=block_id, times=[], powers=[1000.0])
    )
    builder.add_heat_source(
        geosciences.HeatSource(name="probe", block=None, times=[0.0, 1e4], powers=[150.0, 100.0])
    )
    return problem


def test_thermo_hydro_mechanics_export(data_dir, directory):
    model_filename = os.path.join(data_dir, "picasso.og_brep")
    model = opengeode.load_one(opengeode.load_brep(model_filename), model_filename)
    add_east_pressure(model)
    for name, value in [
        ("bulk_modulus", 5e9),
        ("shear_modulus", 3e9),
        ("grain_density", 2650.0),
        ("biot_coefficient", 0.8),
        ("thermal_expansion_coefficient", 1e-5),
    ]:
        add_block_property(model, name, value)
    add_rollers(model)
    block_id = next(iter(model.blocks())).id()
    problem = thermo_hydro_mechanics_problem(model, block_id)
    folder = os.path.join(directory, "python_ogs_thm")
    exporter = geosciences_io.BRepOgsExporter(model, folder)
    probe = opengeode.PointSet3D.create()
    opengeode.PointSetBuilder3D.create(probe).create_point(
        model.bounding_box().center()
    )
    exporter.add_heat_source_points(points=probe, name="probe")
    try:
        exporter.add_heat_source_points(probe, "probe")
        raise ValueError("[Test] Heat source points added twice should raise")
    except opengeode.OpenGeodeDataError:
        pass
    exporter.set_physical_problem(problem)
    exporter.run()
    project = ElementTree.parse(os.path.join(folder, "python_ogs_thm.prj")).getroot()
    process_type = project.findtext("processes/process/type")
    check(
        process_type == "THERMO_HYDRO_MECHANICS",
        f"The project should use THERMO_HYDRO_MECHANICS, not {process_type}",
    )
    sources = {
        term.findtext("mesh"): term.findtext("type")
        for term in project.findall(
            "process_variables/process_variable[name='temperature']/source_terms/source_term"
        )
    }
    check(
        sources == {"python_ogs_thm_heat_canister": "Volumetric", "python_ogs_thm_heat_probe": "Nodal"},
        f"One Volumetric and one Nodal temperature source term: {sources}",
    )
    check(
        project.find("curves/curve[name='heat_probe']") is not None,
        "The power table of the probe is a curve",
    )
    for mesh in sources:
        check(os.path.exists(os.path.join(folder, mesh + ".vtu")), f"{mesh} should be written")


if __name__ == "__main__":
    geosciences_io.OpenGeodeGeosciencesIOModelLibrary.initialize()
    data = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../../../tests/data")
    )
    with tempfile.TemporaryDirectory() as directory:
        test_brep_export(data, directory)
        test_thermal_export(data, directory)
        test_thermo_hydro_mechanics_export(data, directory)
