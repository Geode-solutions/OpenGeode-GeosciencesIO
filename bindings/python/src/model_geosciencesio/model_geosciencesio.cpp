/*
 * Copyright (c) 2019 - 2026 Geode-solutions
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

#include <geode/bindings/python/common.hpp>

#include <geode/mesh/core/point_set.hpp>

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/geosciences/simulation/core/physical_problem.hpp>

#include <geode/geosciences_io/model/helpers/brep_geos_export.hpp>
#include <geode/geosciences_io/model/helpers/structural_model_geos_export.hpp>

#include <pybind11/pybind11.h>

namespace
{
    template < typename Exporter, typename Model >
    void define_geos_exporter(
        pybind11::module_& module, const char* name, const char* class_doc )
    {
        pybind11::class_< Exporter >( module, name, class_doc )
            .def( pybind11::init< const Model&, std::string_view >(),
                pybind11::arg( "model" ), pybind11::arg( "files_directory" ),
                pybind11::keep_alive< 1, 2 >(),
                "Prepare the export of the model in the directory "
                "files_directory, which is created if needed. The directory "
                "name (without its parent) is the prefix of the written files "
                "and must not be empty (OpenGeodeDataError otherwise). The "
                "exporter keeps the model alive." )
            .def( "add_well_perforations", &Exporter::add_well_perforations,
                pybind11::arg( "well_perforations" ), pybind11::arg( "name" ),
                "Add the points of a well perforations (copied at once) under "
                "the given name, which must be non empty, contain no space, "
                "comma nor brace and be different from the previous ones "
                "(OpenGeodeDataError otherwise)." )
            .def( "set_physical_problem", &Exporter::set_physical_problem,
                pybind11::arg( "problem" ), pybind11::keep_alive< 1, 2 >(),
                "Write a full simulation deck (solvers, fluid, rock, initial "
                "conditions, events) described by the PhysicalProblem; "
                "without it, only the mesh, wells and boundary conditions "
                "are written. The exporter keeps the problem alive." )
            .def( "run", &Exporter::run,
                pybind11::call_guard< pybind11::gil_scoped_release >(),
                "Write the GEOS files in the directory; raises "
                "OpenGeodeDataError if the model or the problem cannot be "
                "exported. The GIL is released during the export." );
    }
} // namespace

PYBIND11_MODULE( opengeode_geosciencesio_py_model, module )
{
    module.doc() =
        "OpenGeode-GeosciencesIO model formats: registers in the OpenGeode "
        "factories the readers and writers of StructuralModel (Gocad ML and "
        "LSO), the Shapefile reader of Section, the Gocad ML and FEM writers "
        "of BRep, the Skua reader of HorizonsStack and the GEOS simulation "
        "checker, and binds the GEOS exporters of BRep and StructuralModel";
    geode::import_binding_dependency( module, "opengeode" );
    geode::import_binding_dependency( module, "opengeode_io" );
    geode::import_binding_dependency( module, "opengeode_geosciences" );
    pybind11::class_< geode::OpenGeodeGeosciencesIOModelLibrary >( module,
        "OpenGeodeGeosciencesIOModelLibrary",
        "Initialization of the OpenGeode-GeosciencesIO model library" )
        .def_static( "initialize",
            &geode::OpenGeodeGeosciencesIOModelLibrary::initialize,
            "Register the model readers and writers, and the GEOS "
            "simulation checker, in the OpenGeode factories; also "
            "initializes the OpenGeode-Geosciences explicit, implicit and "
            "simulation libraries and the OpenGeode-GeosciencesIO and "
            "OpenGeode-IO mesh libraries; calling it again has no effect" );
    define_geos_exporter< geode::BRepGeosExporter, geode::BRep >( module,
        "BRepGeosExporter",
        "Export of a BRep (its blocks) as a GEOS simulation input" );
    define_geos_exporter< geode::StructuralModelGeosExporter,
        geode::StructuralModel >( module, "StructuralModelGeosExporter",
        "Export of a StructuralModel (its blocks) as a GEOS simulation "
        "input" );
    geode::document_missing_constructors( module );
}
