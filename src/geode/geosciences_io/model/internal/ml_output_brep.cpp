/*
 * Copyright (c) 2019 - 2025 Geode-solutions
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

#include <geode/geosciences_io/model/internal/ml_output_brep.hpp>

#include <string>
#include <vector>

#include <geode/mesh/core/surface_mesh.hpp>

#include <geode/model/representation/builder/brep_builder.hpp>

#include <geode/geosciences_io/model/internal/gocad_common.hpp>
#include <geode/geosciences_io/model/internal/ml_output_impl.hpp>

namespace
{
    class MLOutputImplBRep : public geode::internal::MLOutputImpl< geode::BRep >
    {
    public:
        MLOutputImplBRep( std::string_view filename, const geode::BRep& model )
            : geode::internal::MLOutputImpl< geode::BRep >( filename, model ),
              model_( model )
        {
        }

    private:
        void write_geological_tfaces() override {}

        void write_geological_tsurfs() override {}

        std::vector< geode::uuid > unclassified_tsurfs() const override
        {
            std::vector< geode::uuid > result;
            for( const auto& surface : model_.surfaces() )
            {
                bool is_part_of_model_boundary{ false };
                for( const auto& collection :
                    model_.collections( surface.id() ) )
                {
                    if( collection.type()
                        == geode::ModelBoundary3D::component_type_static() )
                    {
                        is_part_of_model_boundary = true;
                    }
                }
                if( !is_part_of_model_boundary )
                {
                    result.push_back( surface.id() );
                }
            }
            return result;
        }

        void write_geological_regions() override {}

        void write_geological_model_surfaces() override {}

    private:
        const geode::BRep& model_;
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::vector< std::string > MLOutputBRep::write( const BRep& brep ) const
        {
            OPENGEODE_EXCEPTION( check_brep_polygons( brep ),
                "[MLOutput::write] Can not export into .ml a "
                "BRep with non triangular surface polygons." );
            if( brep.nb_model_boundaries() > 0 )
            {
                MLOutputImplBRep impl{ filename(), brep };
                impl.write_file();
            }
            else
            {
                const auto new_brep = clone_with_model_boundaries( brep );
                MLOutputImplBRep impl{ filename(), new_brep };
                impl.write_file();
            }
            return { to_string( filename() ) };
        }

        bool MLOutputBRep::is_saveable( const BRep& brep ) const
        {
            return check_brep_polygons( brep );
        }
    } // namespace internal
} // namespace geode
