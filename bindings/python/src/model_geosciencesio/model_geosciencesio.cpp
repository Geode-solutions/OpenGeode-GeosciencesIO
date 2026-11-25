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
#include <geode/geosciences_io/model/helpers/brep_ogs_export.hpp>
#include <geode/geosciences_io/model/helpers/structural_model_geos_export.hpp>
#include <geode/geosciences_io/model/helpers/structural_model_ogs_export.hpp>

#include <pybind11/pybind11.h>

namespace
{
    // Docstrings of the parts of an exporter that depend on its simulator
    struct SimulationExporterDocs
    {
        const char* class_doc;
        const char* well_name_rule;
        const char* problem_doc;
        const char* run_doc;
    };

    template < typename Exporter, typename Model >
    void define_simulation_exporter( pybind11::module_& module,
        const char* name,
        const SimulationExporterDocs& docs )
    {
        pybind11::class_< Exporter >( module, name, docs.class_doc )
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
                docs.well_name_rule )
            .def( "set_physical_problem", &Exporter::set_physical_problem,
                pybind11::arg( "problem" ), pybind11::keep_alive< 1, 2 >(),
                docs.problem_doc )
            .def( "run", &Exporter::run,
                pybind11::call_guard< pybind11::gil_scoped_release >(),
                docs.run_doc );
    }

    constexpr SimulationExporterDocs GEOS_DOCS{ nullptr,
        "Add the points of a well perforations (copied at once) under "
        "the given name, which must be non empty, contain no space, "
        "comma nor brace and be different from the previous ones "
        "(OpenGeodeDataError otherwise).",
        "Write a full simulation deck (solvers, fluid, rock, initial "
        "conditions, events) described by the PhysicalProblem; "
        "without it, only the mesh, wells and boundary conditions "
        "are written. The exporter keeps the problem alive.",
        "Write the GEOS files in the directory; raises "
        "OpenGeodeDataError if the model or the problem cannot be "
        "exported. The GIL is released during the export." };

    constexpr SimulationExporterDocs OGS_DOCS{ nullptr,
        "Add the points of a well perforations (copied at once) under "
        "the given name, which must be non empty, contain only letters, "
        "digits, '_' and '-' and be different from the previous ones "
        "(OpenGeodeDataError otherwise). Each point gives the vertex "
        "nearest to it of the cell containing it (points outside the "
        "model are ignored with a warning); a source term rate on the "
        "well is split evenly over these vertices.",
        "Write the OGS project file (process, media, parameters, initial "
        "and boundary conditions, time loop, solvers) described by the "
        "PhysicalProblem: LIQUID_FLOW for single-phase flow, "
        "TWOPHASE_FLOW_PP for water-oil DeadOilFluid flow. Without it, "
        "only the meshes are written. The exporter keeps the problem "
        "alive.",
        "Write the OGS files in the directory: <prefix>.vtu (bulk mesh "
        "with MaterialIDs), one mesh per boundary condition group, one "
        "point mesh per well and, with a problem, <prefix>.prj; raises "
        "OpenGeodeDataError if the model or the problem cannot be "
        "exported. The GIL is released during the export." };

    SimulationExporterDocs with_class_doc(
        SimulationExporterDocs docs, const char* class_doc )
    {
        docs.class_doc = class_doc;
        return docs;
    }
} // namespace

PYBIND11_MODULE( opengeode_geosciencesio_py_model, module )
{
    module.doc() =
        "OpenGeode-GeosciencesIO model formats: registers in the OpenGeode "
        "factories the readers and writers of StructuralModel (Gocad ML and "
        "LSO), the Shapefile reader of Section, the Gocad ML and FEM writers "
        "of BRep, the Skua reader of HorizonsStack and the GEOS and OGS "
        "simulation checkers, and binds the GEOS and OGS (OpenGeoSys) "
        "exporters of BRep and StructuralModel";
    geode::import_binding_dependency( module, "opengeode" );
    geode::import_binding_dependency( module, "opengeode_io" );
    geode::import_binding_dependency( module, "opengeode_geosciences" );
    pybind11::class_< geode::OpenGeodeGeosciencesIOModelLibrary >( module,
        "OpenGeodeGeosciencesIOModelLibrary",
        "Initialization of the OpenGeode-GeosciencesIO model library" )
        .def_static( "initialize",
            &geode::OpenGeodeGeosciencesIOModelLibrary::initialize,
            "Register the model readers and writers, and the GEOS and OGS "
            "simulation checkers, in the OpenGeode factories; also "
            "initializes the OpenGeode-Geosciences explicit, implicit and "
            "simulation libraries and the OpenGeode-GeosciencesIO and "
            "OpenGeode-IO mesh libraries; calling it again has no effect" );
    define_simulation_exporter< geode::BRepGeosExporter, geode::BRep >( module,
        "BRepGeosExporter",
        with_class_doc( GEOS_DOCS,
            "Export of a BRep (its blocks) as a GEOS simulation input" ) );
    define_simulation_exporter< geode::StructuralModelGeosExporter,
        geode::StructuralModel >( module, "StructuralModelGeosExporter",
        with_class_doc( GEOS_DOCS,
            "Export of a StructuralModel (its blocks) as a GEOS simulation "
            "input" ) );
    define_simulation_exporter< geode::BRepOgsExporter, geode::BRep >( module,
        "BRepOgsExporter",
        with_class_doc( OGS_DOCS,
            "Export of a BRep (its blocks) as an OpenGeoSys 6 simulation "
            "input" ) );
    define_simulation_exporter< geode::StructuralModelOgsExporter,
        geode::StructuralModel >( module, "StructuralModelOgsExporter",
        with_class_doc( OGS_DOCS,
            "Export of a StructuralModel (its blocks) as an OpenGeoSys 6 "
            "simulation input" ) );
    geode::document_missing_constructors( module );
}
