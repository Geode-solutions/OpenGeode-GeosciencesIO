/*
 * Copyright (c) 2019 - 2024 Geode-solutions
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include <string_view>

#include <geode/mesh/core/point_set.h>

#include <geode/geosciences/explicit/representation/core/structural_model.h>
#include <geode/geosciences_io/model/helpers/brep_geos_export.h>
#include <geode/geosciences_io/model/helpers/structural_model_geos_export.h>

#include <pybind11/pybind11.h>

PYBIND11_MODULE( opengeode_geosciencesio_py_model, module )
{
    module.doc() = "OpenGeode-GeosciencesIO Python binding for model";
    pybind11::class_< geode::GeosciencesIOModelLibrary >(
        module, "GeosciencesIOModelLibrary" )
        .def( "initialize", &geode::GeosciencesIOModelLibrary::initialize );

    pybind11::class_< geode::BRepGeosExporter >( module, "BRepGeosExporter" )
        .def( pybind11::init< const geode::BRep&, std::string_view >() )
        .def( "add_well_perforations",
            &geode::BRepGeosExporter::add_well_perforations )
        .def( "add_cell_property_1d",
            &geode::BRepGeosExporter::add_cell_property_1d )
        .def( "add_cell_property_2d",
            &geode::BRepGeosExporter::add_cell_property_2d )
        .def( "add_cell_property_3d",
            &geode::BRepGeosExporter::add_cell_property_3d )
        .def( "run", &geode::BRepGeosExporter::run );

    pybind11::class_< geode::StructuralModelGeosExporter >(
        module, "StructuralModelGeosExporter" )
        .def( pybind11::init< const geode::StructuralModel&,
            std::string_view >() )
        .def( "add_well_perforations",
            &geode::StructuralModelGeosExporter::add_well_perforations )
        .def( "add_cell_property_1d",
            &geode::StructuralModelGeosExporter::add_cell_property_1d )
        .def( "add_cell_property_2d",
            &geode::StructuralModelGeosExporter::add_cell_property_2d )
        .def( "add_cell_property_3d",
            &geode::StructuralModelGeosExporter::add_cell_property_3d )
        .def( "run", &geode::StructuralModelGeosExporter::run );
}