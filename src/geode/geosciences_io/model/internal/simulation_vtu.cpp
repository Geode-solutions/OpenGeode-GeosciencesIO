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

#include <geode/geosciences_io/model/internal/simulation_vtu.hpp>

#include <geode/basic/attribute_manager.hpp>

#include <geode/mesh/core/hybrid_solid.hpp>
#include <geode/mesh/core/polyhedral_solid.hpp>
#include <geode/mesh/core/tetrahedral_solid.hpp>

#include <geode/io/mesh/detail/vtu_solid_output_impl.hpp>

namespace
{
    // Writes the solid polyhedra followed by the extra polygons in a single
    // UnstructuredGrid, all the cells sharing the solid vertices
    template < typename SolidOutputImpl >
    class SimulationVTUOutputImpl : public SolidOutputImpl
    {
    public:
        template < typename Solid >
        SimulationVTUOutputImpl( std::string_view filename,
            const Solid& solid,
            absl::Span< const geode::PolygonVertices > extra_polygons,
            const geode::AttributeManager& cell_attributes )
            : SolidOutputImpl{ filename, solid },
              extra_polygons_{ extra_polygons },
              cell_attributes_( cell_attributes )
        {
        }

    private:
        [[nodiscard]] geode::index_t nb_additional_polygons() const override
        {
            return extra_polygons_.size();
        }

        [[nodiscard]] absl::Span< const geode::index_t >
            additional_polygon_vertices(
                geode::index_t polygon_id ) const override
        {
            return extra_polygons_[polygon_id];
        }

        [[nodiscard]] const geode::AttributeManager&
            cell_attribute_manager() const override
        {
            return cell_attributes_;
        }

        absl::Span< const geode::PolygonVertices > extra_polygons_;
        const geode::AttributeManager& cell_attributes_;
    };

    template < typename SolidOutputImpl, typename Solid >
    void save_vtu( std::string_view filename,
        const Solid& solid,
        absl::Span< const geode::PolygonVertices > extra_polygons,
        const geode::AttributeManager& cell_attributes )
    {
        SimulationVTUOutputImpl< SolidOutputImpl > writer{ filename, solid,
            extra_polygons, cell_attributes };
        writer.write_file();
    }
} // namespace

namespace geode::internal
{
    bool save_simulation_vtu( std::string_view filename,
        const SolidMesh3D& solid,
        absl::Span< const PolygonVertices > extra_polygons,
        const AttributeManager& cell_attributes )
    {
        if( const auto* tetra =
                dynamic_cast< const TetrahedralSolid3D* >( &solid ) )
        {
            save_vtu< detail::VTUTetrahedralOutputImpl >(
                filename, *tetra, extra_polygons, cell_attributes );
            return true;
        }
        if( const auto* hybrid =
                dynamic_cast< const HybridSolid3D* >( &solid ) )
        {
            save_vtu< detail::VTUHybridOutputImpl >(
                filename, *hybrid, extra_polygons, cell_attributes );
            return true;
        }
        if( const auto* poly =
                dynamic_cast< const PolyhedralSolid3D* >( &solid ) )
        {
            save_vtu< detail::VTUPolyhedralOutputImpl >(
                filename, *poly, extra_polygons, cell_attributes );
            return true;
        }
        return false;
    }
} // namespace geode::internal
