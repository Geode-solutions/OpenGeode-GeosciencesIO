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

#include <geode/geosciences_io/mesh/internal/pl_output.hpp>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <geode/basic/attribute_manager.hpp>

#include <geode/basic/filename.hpp>
#include <geode/mesh/core/edged_curve.hpp>
#include <geode/mesh/core/graph.hpp>

#include <geode/geosciences_io/mesh/internal/gocad_common.hpp>

namespace
{
    class PLOutputImpl
    {
    public:
        static constexpr geode::index_t OFFSET_START{ 1 };
        static constexpr char EOL{ '\n' };
        static constexpr char SPACE{ ' ' };
        PLOutputImpl(
            std::string_view filename, const geode::EdgedCurve3D& edged_curve )
            : file_{ geode::to_string( filename ) },
              edged_curve_( edged_curve ),
              edge_done_( edged_curve.nb_edges(), false )
        {
            OPENGEODE_EXCEPTION(
                file_.good(), "Error while opening file: ", filename );
        }

        void write_prop_header()
        {
            const auto names =
                edged_curve_.vertex_attribute_manager().attribute_names();
            geode::internal::PropHeaderData prop_header;
            std::vector< geode::internal::PropClassHeaderData >
                prop_class_header;
            prop_class_header.reserve( names.size() );

            for( const auto& name : names )
            {
                VRTX_KEYWORD = "PVRTX";
                const auto attribute = edged_curve_.vertex_attribute_manager()
                                           .find_generic_attribute( name );
                if( !attribute || !attribute->is_genericable() )
                {
                    continue;
                }
                generic_att_.push_back( attribute );
                prop_header.names.emplace_back( geode::to_string( name ) );
                prop_header.prop_legal_ranges.push_back(
                    std::make_pair( "**none**", "**none**" ) );
                prop_header.no_data_values.push_back( -99999. );
                prop_header.property_classes.emplace_back(
                    geode::to_string( name ) );
                prop_header.kinds.push_back( "Real Number" );
                prop_header.property_subclass.push_back(
                    std::make_pair( "QUANTITY", "Float" ) );
                prop_header.esizes.push_back( 1 );
                prop_header.units.push_back( "unitless" );

                geode::internal::PropClassHeaderData tmp_prop_class_header;
                tmp_prop_class_header.name = geode::to_string( name );

                prop_class_header.push_back( tmp_prop_class_header );
            }
            if( !prop_header.empty() )
            {
                geode::internal::write_prop_header( file_, prop_header );
            }
            write_XYZ_prop_class_header();
            for( const auto& ch : prop_class_header )
            {
                write_property_class_header( file_, ch );
            }
        }

        void write_XYZ_prop_class_header()
        {
            geode::internal::PropClassHeaderData x_prop_header;
            x_prop_header.name = "X";
            x_prop_header.kind = "X";
            x_prop_header.unit = "m";
            x_prop_header.is_z = false;
            write_property_class_header( file_, x_prop_header );

            geode::internal::PropClassHeaderData y_prop_header;
            y_prop_header.name = "Y";
            y_prop_header.kind = "Y";
            y_prop_header.unit = "m";
            y_prop_header.is_z = false;
            write_property_class_header( file_, y_prop_header );

            geode::internal::PropClassHeaderData z_prop_header;
            z_prop_header.name = "Z";
            z_prop_header.kind = "Z";
            z_prop_header.unit = "m";
            z_prop_header.is_z = true;
            write_property_class_header( file_, z_prop_header );
        }

        void write_pvrtx(
            const geode::index_t v, const geode::index_t current_offset )
        {
            file_ << VRTX_KEYWORD << SPACE << current_offset << SPACE
                  << edged_curve_.point( v ).value( 0 ) << SPACE
                  << edged_curve_.point( v ).value( 1 ) << SPACE
                  << edged_curve_.point( v ).value( 2 );
            for( const auto& att : generic_att_ )
            {
                file_ << SPACE << att->generic_value( v );
            }
            file_ << EOL;
        }
        std::vector< geode::EdgeVertex > get_edged_vertex_on_iline(
            const geode::EdgeVertex& ev )
        {
            std::vector< geode::EdgeVertex > ev_on_iline;
            geode::EdgeVertex next_ev = ev;

            auto propagate = true;
            while( propagate )
            {
                propagate = false;
                ev_on_iline.push_back( next_ev );
                edge_done_[next_ev.edge_id] = true;
                next_ev = { next_ev.edge_id,
                    static_cast< geode::local_index_t >(
                        ( next_ev.vertex_id + 1 ) % 2 ) };
                auto edges_around = edged_curve_.edges_around_vertex(
                    edged_curve_.edge_vertex( next_ev ) );
                propagate = edges_around.size() == 2
                            && ( !edge_done_[edges_around[0].edge_id]
                                 || !edge_done_[edges_around[1].edge_id] );
                if( propagate )
                {
                    for( const auto& edge : edges_around )
                    {
                        if( edge_done_[edge.edge_id] )
                        {
                            continue;
                        }
                        next_ev = edge;
                    }
                }
                else if( edges_around.size() != 2 )
                {
                    ev_on_iline.push_back( next_ev );
                }
            }
            return ev_on_iline;
        }

        void write_ilines()
        {
            std::vector< geode::index_t > start_point;
            for( auto v : geode::Range{ edged_curve_.nb_vertices() } )
            {
                if( edged_curve_.edges_around_vertex( v ).size() != 2 )
                {
                    start_point.push_back( v );
                }
            }
            auto current_offset = OFFSET_START;
            geode::index_t nb_edges_done{ 0 };
            for( const auto& v_id : start_point )
            {
                for( const auto& edge :
                    edged_curve_.edges_around_vertex( v_id ) )
                {
                    if( edge_done_[edge.edge_id] )
                    {
                        continue;
                    }
                    write_edge_and_vertex(
                        edge, current_offset, nb_edges_done );
                }
            }
            while( nb_edges_done != edged_curve_.nb_edges() )
            {
                const auto start_edge = starting_edge();
                write_edge_and_vertex(
                    { start_edge, 0 }, current_offset, nb_edges_done );
                file_ << "SEG" << SPACE << current_offset - 1 << SPACE
                      << edged_curve_.edge_vertices( start_edge )[0] + 1 << EOL;
                nb_edges_done++;
            }
        }

        geode::index_t starting_edge()
        {
            geode::index_t starting_edge{ 0 };
            for( const auto& edge : geode::Range( edged_curve_.nb_edges() ) )
            {
                if( !edge_done_[edge] )
                {
                    starting_edge = edge;
                    break;
                }
            }
            return starting_edge;
        }

        void write_edge_and_vertex( const geode::EdgeVertex& edge,
            geode::index_t& current_offset,
            geode::index_t& nb_edges_done )
        {
            file_ << "ILINE" << EOL;
            auto ev_on_iline = get_edged_vertex_on_iline( edge );
            geode::index_t cur_v = 0;
            for( const auto& ev : ev_on_iline )
            {
                write_pvrtx(
                    edged_curve_.edge_vertex( ev ), current_offset + cur_v );
                ++cur_v;
            }
            for( auto cur_seg : geode::Range{ ev_on_iline.size() - 1 } )
            {
                file_ << "SEG" << SPACE << current_offset + cur_seg << SPACE
                      << current_offset + cur_seg + 1 << EOL;
                nb_edges_done++;
            }
            current_offset += ev_on_iline.size();
        }

        void write_file()
        {
            geode::Logger::info( "[PLOutput::write] Writting pl file." );
            file_ << "GOCAD PLine 1" << EOL;
            geode::internal::HeaderData header;
            header.name = "edged_curve_name";
            geode::internal::write_header( file_, header );
            geode::internal::write_CRS( file_, {} );
            write_prop_header();
            write_ilines();
            file_ << "END" << EOL;
            return;
        }

    private:
        std::ofstream file_;
        const geode::EdgedCurve3D& edged_curve_;
        std::vector< std::shared_ptr< geode::AttributeBase > > generic_att_;
        std::vector< bool > edge_done_;
        std::string VRTX_KEYWORD{ "VRTX" };
    };
} // namespace

namespace geode
{
    namespace internal
    {
        std::vector< std::string > PLOutput::write(
            const EdgedCurve3D& edged_curve ) const
        {
            PLOutputImpl impl{ filename(), edged_curve };
            impl.write_file();
            return { to_string( filename() ) };
        }
    } // namespace internal
} // namespace geode
