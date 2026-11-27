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

#include <geode/geosciences_io/model/internal/simulation_export.hpp>

#include <filesystem>
#include <limits>

#include <absl/algorithm/container.h>

#include <geode/basic/range.hpp>

#include <geode/geometry/aabb.hpp>
#include <geode/geometry/basic_objects/tetrahedron.hpp>
#include <geode/geometry/distance.hpp>
#include <geode/geometry/information.hpp>
#include <geode/geometry/point.hpp>
#include <geode/geometry/position.hpp>

#include <geode/mesh/core/point_set.hpp>
#include <geode/mesh/core/solid_mesh.hpp>

#include <geode/geosciences/explicit/representation/core/structural_model.hpp>
#include <geode/model/representation/core/brep.hpp>

namespace
{
    // Directory without trailing separator nor dot segments
    std::filesystem::path normalized_directory( std::string_view directory )
    {
        auto path = std::filesystem::path{ geode::to_string( directory ) }
                        .lexically_normal();
        if( !path.empty() && !path.has_filename() )
        {
            path = path.parent_path();
        }
        return path;
    }

    // The polyhedron is split in tetrahedra joining its barycenter to a fan
    // triangulation of its facets
    bool is_point_in_polyhedron( const geode::SolidMesh3D& solid,
        geode::index_t polyhedron_id,
        const geode::Point3D& point )
    {
        const auto barycenter = solid.polyhedron_barycenter( polyhedron_id );
        for( const auto facet :
            geode::LRange{ solid.nb_polyhedron_facets( polyhedron_id ) } )
        {
            const auto vertices =
                solid.polyhedron_facet_vertices( { polyhedron_id, facet } );
            const auto& first = solid.point( vertices[0] );
            for( const auto vertex : geode::LRange{ 1, vertices.size() - 1 } )
            {
                const geode::Tetrahedron tetrahedron{ first,
                    solid.point( vertices[vertex] ),
                    solid.point( vertices[vertex + 1] ), barycenter };
                const auto position =
                    geode::point_tetrahedron_position( point, tetrahedron );
                if( position != geode::POSITION::outside
                    && position != geode::POSITION::parallel )
                {
                    return true;
                }
            }
        }
        return false;
    }
} // namespace

namespace geode::internal
{
    index_t containing_cell(
        const SolidMesh3D& solid, const AABBTree3D& aabb, const Point3D& point )
    {
        double distance_to_nearest_cell_center{
            std::numeric_limits< double >::max()
        };
        index_t selected_cell_id{ NO_ID };
        for( const auto cell_id : aabb.containing_boxes( point ) )
        {
            if( !is_point_in_polyhedron( solid, cell_id, point ) )
            {
                continue;
            }
            const auto distance = point_point_distance< 3 >(
                point, solid.polyhedron_barycenter( cell_id ) );
            if( distance < distance_to_nearest_cell_center )
            {
                distance_to_nearest_cell_center = distance;
                selected_cell_id = cell_id;
            }
        }
        return selected_cell_id;
    }

    template < typename Model >
    SimulationExporterImpl< Model >::SimulationExporterImpl(
        std::string_view files_directory,
        const Model& model,
        std::string_view context )
        : model_( model ),
          context_{ to_string( context ) },
          files_directory_{ normalized_directory( files_directory ).string() },
          prefix_{ normalized_directory( files_directory ).stem().string() }
    {
        OpenGeodeGeosciencesIOModelException::check_exception( !prefix_.empty(),
            nullptr, OpenGeodeException::TYPE::data, context_,
            " Cannot export in directory \"", files_directory,
            "\": its name is used as file prefix and must not be empty." );
        std::filesystem::create_directories( files_directory_ );
    }

    template < typename Model >
    SimulationExporterImpl< Model >::~SimulationExporterImpl() = default;

    template < typename Model >
    void SimulationExporterImpl< Model >::set_physical_problem(
        const PhysicalProblem& problem )
    {
        problem_ = &problem;
    }

    template < typename Model >
    void SimulationExporterImpl< Model >::add_well_perforations(
        const PointSet3D& perforations, std::string_view name )
    {
        check_well_name( name );
        OpenGeodeGeosciencesIOModelException::check_exception(
            absl::c_find( well_names_, name ) == well_names_.end(), nullptr,
            OpenGeodeException::TYPE::data, context_,
            " Well perforations named ", name, " already added." );
        well_perforations_.push_back( perforations.clone() );
        well_names_.emplace_back( to_string( name ) );
    }

    template < typename Model >
    void SimulationExporterImpl< Model >::add_heat_source_points(
        const PointSet3D& points, std::string_view name )
    {
        OpenGeodeGeosciencesIOModelException::check_exception(
            absl::c_find( heat_source_names_, name )
                == heat_source_names_.end(),
            nullptr, OpenGeodeException::TYPE::data, context_,
            " Heat source points named ", name, " already added." );
        heat_source_points_.push_back( points.clone() );
        heat_source_names_.emplace_back( to_string( name ) );
    }

    template < typename Model >
    const Model& SimulationExporterImpl< Model >::model() const
    {
        return model_;
    }

    template < typename Model >
    const PhysicalProblem* SimulationExporterImpl< Model >::problem() const
    {
        return problem_;
    }

    template < typename Model >
    std::string_view SimulationExporterImpl< Model >::files_directory() const
    {
        return files_directory_;
    }

    template < typename Model >
    std::string_view SimulationExporterImpl< Model >::prefix() const
    {
        return prefix_;
    }

    template < typename Model >
    absl::Span< const std::unique_ptr< PointSet3D > >
        SimulationExporterImpl< Model >::well_perforations() const
    {
        return well_perforations_;
    }

    template < typename Model >
    absl::Span< const std::string >
        SimulationExporterImpl< Model >::well_names() const
    {
        return well_names_;
    }

    template < typename Model >
    absl::Span< const std::unique_ptr< PointSet3D > >
        SimulationExporterImpl< Model >::heat_source_points() const
    {
        return heat_source_points_;
    }

    template < typename Model >
    absl::Span< const std::string >
        SimulationExporterImpl< Model >::heat_source_names() const
    {
        return heat_source_names_;
    }

    template class opengeode_geosciencesio_model_api
        SimulationExporterImpl< BRep >;
    template class opengeode_geosciencesio_model_api
        SimulationExporterImpl< StructuralModel >;
} // namespace geode::internal
