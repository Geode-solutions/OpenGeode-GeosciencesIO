/*
 * Copyright (c) 2019 - 2022 Geode-solutions
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

#include <geode/geosciences/private/gocad_common.h>

#include <fstream>
#include <queue>

#include <absl/strings/match.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_replace.h>
#include <absl/strings/str_split.h>

#include <geode/basic/logger.h>

#include <geode/geometry/basic_objects/tetra.h>
#include <geode/geometry/bounding_box.h>
#include <geode/geometry/signed_mensuration.h>

#include <geode/mesh/core/edged_curve.h>
#include <geode/mesh/core/surface_mesh.h>

#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/line.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/model/representation/core/brep.h>

namespace
{
    static constexpr char EOL{ '\n' };
    static constexpr char SPACE{ ' ' };

    using PairedSigns =
        absl::flat_hash_map< std::pair< geode::uuid, geode::uuid >, bool >;

    PairedSigns determine_paired_signs( const geode::BRep& brep )
    {
        PairedSigns paired_signs;
        paired_signs.reserve( brep.nb_lines() * 2 );
        for( const auto& line : brep.lines() )
        {
            const auto& mesh = line.mesh();
            const auto uid0 = brep.unique_vertex(
                { line.component_id(), mesh.edge_vertex( { 0, 0 } ) } );
            const auto uid1 = brep.unique_vertex(
                { line.component_id(), mesh.edge_vertex( { 0, 1 } ) } );
            const auto surface_mcvs0 = brep.mesh_component_vertices(
                uid0, geode::Surface3D::component_type_static() );
            const auto surface_mcvs1 = brep.mesh_component_vertices(
                uid1, geode::Surface3D::component_type_static() );
            absl::flat_hash_map< geode::uuid, bool > surface_direct_edges;
            surface_direct_edges.reserve(
                std::min( surface_mcvs0.size(), surface_mcvs1.size() ) );
            for( const auto& surface_mcv0 : surface_mcvs0 )
            {
                for( const auto& surface_mcv1 : surface_mcvs1 )
                {
                    if( surface_mcv1.component_id.id()
                        != surface_mcv0.component_id.id() )
                    {
                        continue;
                    }
                    const auto& surface =
                        brep.surface( surface_mcv0.component_id.id() );
                    const auto& surface_mesh = surface.mesh();
                    const auto v0v1 = surface_mesh.polygon_edge_from_vertices(
                        surface_mcv0.vertex, surface_mcv1.vertex );
                    const auto v1v0 = surface_mesh.polygon_edge_from_vertices(
                        surface_mcv1.vertex, surface_mcv0.vertex );
                    if( v0v1 && !v1v0 )
                    {
                        surface_direct_edges.emplace( surface.id(), true );
                    }
                    else if( v1v0 && !v0v1 )
                    {
                        surface_direct_edges.emplace( surface.id(), false );
                    }
                }
            }
            if( surface_direct_edges.size() < 2 )
            {
                continue;
            }
            for( const auto& s0 : surface_direct_edges )
            {
                for( const auto& s1 : surface_direct_edges )
                {
                    if( s0.first < s1.first )
                    {
                        paired_signs.emplace(
                            std::make_pair( s0.first, s1.first ),
                            s0.second != s1.second );
                    }
                }
            }
        }
        return paired_signs;
    }

    absl::FixedArray< bool > determine_relative_signs(
        absl::Span< const geode::uuid > universe_boundaries,
        const PairedSigns& paired_signs )
    {
        const auto nb_surfaces = universe_boundaries.size();
        if( nb_surfaces == 1 )
        {
            return { true };
        }
        absl::FixedArray< bool > signs( nb_surfaces, true );
        absl::FixedArray< bool > determined( nb_surfaces, false );
        determined[0] = true;
        std::queue< geode::index_t > to_process;
        to_process.push( 0 );
        while( !to_process.empty() )
        {
            const auto determined_s = to_process.front();
            const auto& determined_s_id = universe_boundaries[determined_s];
            to_process.pop();
            for( const auto s : geode::Range{ nb_surfaces } )
            {
                if( determined[s] )
                {
                    continue;
                }
                const auto& s_id = universe_boundaries[s];
                const auto itr =
                    paired_signs.find( { std::min( determined_s_id, s_id ),
                        std::max( determined_s_id, s_id ) } );
                if( itr == paired_signs.end() )
                {
                    continue;
                }
                signs[s] =
                    itr->second ? signs[determined_s] : !signs[determined_s];
                determined[s] = true;
                to_process.push( s );
            }
        }
        OPENGEODE_ASSERT( absl::c_count( determined, false ) == 0,
            "All signs should have been found" );
        return signs;
    }

    bool are_correct_sides( const geode::BRep& brep,
        absl::Span< const geode::uuid > universe_boundaries,
        absl::Span< const bool > relative_signs )
    {
        double signed_volume{ 0 };
        const auto& first_surface_mesh =
            brep.surface( universe_boundaries[0] ).mesh();
        const auto& bbox = first_surface_mesh.bounding_box();
        const auto center = ( bbox.min() + bbox.max() ) * 0.5;
        for( const auto s : geode::Indices{ universe_boundaries } )
        {
            const auto sign = relative_signs[s];
            const auto& surface_mesh =
                brep.surface( universe_boundaries[s] ).mesh();
            for( const auto t : geode::Range{ surface_mesh.nb_polygons() } )
            {
                std::array< geode::local_index_t, 3 > vertex_order{ 0, 1, 2 };
                if( !sign )
                {
                    std::swap( vertex_order[1], vertex_order[2] );
                }
                signed_volume += geode::tetra_signed_volume(
                    { surface_mesh.point( surface_mesh.polygon_vertex(
                          { t, vertex_order[0] } ) ),
                        surface_mesh.point( surface_mesh.polygon_vertex(
                            { t, vertex_order[1] } ) ),
                        surface_mesh.point( surface_mesh.polygon_vertex(
                            { t, vertex_order[2] } ) ),
                        center } );
            }
        }
        OPENGEODE_ASSERT(
            std::fabs( signed_volume ) > 0, "Null volume block is not valid" );
        return signed_volume > 0;
    }

    absl::flat_hash_map< geode::uuid, bool > determine_universe_sides(
        const geode::BRep& brep, const PairedSigns& paired_signs )
    {
        std::vector< geode::uuid > universe_boundaries;
        universe_boundaries.reserve( 2 * brep.nb_model_boundaries() );
        for( const auto& boundary : brep.model_boundaries() )
        {
            for( const auto& item : brep.model_boundary_items( boundary ) )
            {
                universe_boundaries.push_back( item.id() );
            }
        }
        const auto relative_signs =
            determine_relative_signs( universe_boundaries, paired_signs );
        const auto correct =
            are_correct_sides( brep, universe_boundaries, relative_signs );

        absl::flat_hash_map< geode::uuid, bool > sides;
        sides.reserve( universe_boundaries.size() );
        for( const auto b : geode::Indices{ universe_boundaries } )
        {
            sides.emplace( universe_boundaries[b],
                correct ? !relative_signs[b] : relative_signs[b] );
        }
        return sides;
    }

    PairedSigns determine_regions_sides(
        const geode::BRep& brep, const PairedSigns& paired_signs )
    {
        PairedSigns sides;
        for( const auto& block : brep.blocks() )
        {
            std::vector< geode::uuid > block_boundaries;
            block_boundaries.reserve( brep.nb_boundaries( block.id() ) );
            for( const auto& surface : brep.boundaries( block ) )
            {
                block_boundaries.push_back( surface.id() );
            }
            const auto relative_signs =
                determine_relative_signs( block_boundaries, paired_signs );
            const auto correct =
                are_correct_sides( brep, block_boundaries, relative_signs );

            for( const auto b : geode::Indices{ block_boundaries } )
            {
                sides.emplace(
                    std::make_pair( block.id(), block_boundaries[b] ),
                    correct ? relative_signs[b] : !relative_signs[b] );
            }
        }
        return sides;
    }

    void read_tfaces( std::ifstream& file, geode::detail::TSurfData& tsurf )
    {
        geode::detail::goto_keyword( file, "TFACE" );
        std::string line;
        while( std::getline( file, line ) )
        {
            std::istringstream iss{ line };
            std::string keyword;
            iss >> keyword;

            if( keyword == "VRTX" || keyword == "PVRTX" )
            {
                geode::index_t index;
                double x, y, z;
                iss >> index >> x >> y >> z;
                if( tsurf.points.empty() )
                {
                    tsurf.OFFSET_START = index;
                }
                tsurf.points.push_back(
                    geode::Point3D{ { x, y, tsurf.crs.z_sign * z } } );
            }
            else if( keyword == "ATOM" || keyword == "PATOM" )
            {
                geode::index_t dummy, atom_id;
                iss >> dummy >> atom_id;
                tsurf.points.push_back(
                    tsurf.points.at( atom_id - tsurf.OFFSET_START ) );
            }
            else if( keyword == "TRGL" )
            {
                geode::index_t v0, v1, v2;
                iss >> v0 >> v1 >> v2;
                tsurf.triangles.push_back( { v0 - tsurf.OFFSET_START,
                    v1 - tsurf.OFFSET_START, v2 - tsurf.OFFSET_START } );
            }
            else if( keyword == "BSTONE" )
            {
                geode::index_t bstone_id;
                iss >> bstone_id;
                tsurf.bstones.push_back( bstone_id - tsurf.OFFSET_START );
            }
            else if( keyword == "BORDER" )
            {
                geode::index_t dummy, corner, next;
                iss >> dummy >> corner >> next;
                tsurf.borders.emplace_back(
                    corner - tsurf.OFFSET_START, next - tsurf.OFFSET_START );
            }
            else if( keyword == "TFACE" )
            {
                tsurf.tface_triangles_offset.push_back(
                    tsurf.triangles.size() );
                tsurf.tface_vertices_offset.push_back( tsurf.points.size() );
            }
            else if( keyword == "END" )
            {
                tsurf.tface_triangles_offset.push_back(
                    tsurf.triangles.size() );
                tsurf.tface_vertices_offset.push_back( tsurf.points.size() );
                return;
            }
        }
        throw geode::OpenGeodeException{
            "[read_tfaces] Cannot find the end of TSurf section"
        };
    }
} // namespace

namespace geode
{
    namespace detail
    {
        bool string_starts_with(
            absl::string_view string, absl::string_view check )
        {
            return absl::StartsWith( string, check );
        }

        bool line_starts_with( std::ifstream& file, absl::string_view check )
        {
            std::string line;
            std::getline( file, line );
            return string_starts_with( line, check );
        }

        void check_keyword( std::ifstream& file, absl::string_view keyword )
        {
            OPENGEODE_EXCEPTION( line_starts_with( file, keyword ),
                absl::StrCat( "Line should starts with \"", keyword, "\"" ) );
        }

        HeaderData read_header( std::ifstream& file )
        {
            check_keyword( file, "HEADER" );
            HeaderData header;
            std::string line;
            while( std::getline( file, line ) )
            {
                if( string_starts_with( line, "}" ) )
                {
                    return header;
                }
                std::vector< absl::string_view > tokens =
                    absl::StrSplit( absl::StripAsciiWhitespace( line ), ":" );
                if( tokens.front() == "name" )
                {
                    header.name = absl::StrReplaceAll(
                        absl::StrJoin( tokens.begin() + 1, tokens.end(), " " ),
                        { { "\"", "" } } );
                }
            }
            throw geode::OpenGeodeException{
                "[read_header] Cannot find the end of \"HEADER\" section"
            };
        }

        void write_header( std::ofstream& file, const HeaderData& data )
        {
            file << "HEADER {" << EOL;
            file << "name:" << data.name << EOL;
            file << "}" << EOL;
        }

        CRSData read_CRS( std::ifstream& file )
        {
            CRSData crs;
            if( !line_starts_with( file, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) )
            {
                return crs;
            }
            std::string line;
            while( std::getline( file, line ) )
            {
                if( string_starts_with(
                        line, "END_ORIGINAL_COORDINATE_SYSTEM" ) )
                {
                    return crs;
                }
                std::istringstream iss{ line };
                std::string keyword;
                iss >> keyword;
                if( keyword == "ZPOSITIVE" )
                {
                    std::string sign;
                    iss >> sign;
                    crs.z_sign = sign == "Elevation" ? 1 : -1;
                }
            }
            throw geode::OpenGeodeException{
                "Cannot find the end of CRS section"
            };
        }

        void write_CRS( std::ofstream& file, const CRSData& data )
        {
            file << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << EOL;
            file << "NAME " << data.name << EOL;
            file << "AXIS_NAME " << data.axis_names[0] << SPACE
                 << data.axis_names[1] << SPACE << data.axis_names[2] << EOL;
            file << "AXIS_UNIT " << data.axis_units[0] << SPACE
                 << data.axis_units[1] << SPACE << data.axis_units[2] << EOL;
            file << "ZPOSITIVE " << ( data.z_sign == 1 ? "Elevation" : "Depth" )
                 << EOL;
            file << "END_ORIGINAL_COORDINATE_SYSTEM" << EOL;
        }

        void write_prop_header(
            std::ofstream& file, const PropHeaderData& data )
        {
            file << "PROPERTIES";
            for( const auto& name : data.names )
            {
                file << SPACE << name;
            }
            file << EOL;
            file << "PROP_LEGAL_RANGES";
            for( const auto& prop_range : data.prop_legal_ranges )
            {
                file << SPACE << prop_range.first << SPACE << prop_range.second;
            }
            file << EOL;
            file << "NO_DATA_VALUES";
            for( const auto prop_ndv : data.no_data_values )
            {
                file << SPACE << prop_ndv;
            }
            file << EOL;
            file << "PROPERTY_CLASSES";
            for( const auto& prop_class : data.property_classes )
            {
                file << SPACE << prop_class;
            }
            file << EOL;
            file << "PROPERTY_KINDS";
            for( const auto& kind : data.kinds )
            {
                file << SPACE << kind;
            }
            file << EOL;
            file << "PROPERTY_SUBCLASSES";
            for( const auto& prop_subclasse : data.property_subclass )
            {
                file << SPACE << prop_subclasse.first << SPACE
                     << prop_subclasse.second;
            }
            file << EOL;
            file << "ESIZES";
            for( const auto esize : data.esizes )
            {
                file << SPACE << esize;
            }
            file << EOL;
            file << "UNITS";
            for( const auto& unit : data.units )
            {
                file << SPACE << unit;
            }
            file << EOL;
        }
        void write_property_class_header(
            std::ofstream& file, const PropClassHeaderData& data )
        {
            file << "PROPERTY_CLASS_HEADER" << SPACE << data.name << SPACE
                 << "{" << EOL;
            file << "kind:" << data.kind << EOL;
            file << "unit:" << data.unit << EOL;
            if( data.is_z )
            {
                file << "is_Z: on" << EOL;
            }
            file << "}" << EOL;
        }
        std::string goto_keyword( std::ifstream& file, absl::string_view word )
        {
            std::string line;
            while( std::getline( file, line ) )
            {
                if( string_starts_with( line, word ) )
                {
                    return line;
                }
            }
            throw geode::OpenGeodeException{
                "[goto_keyword] Cannot find the requested keyword: ", word
            };
            return "";
        }

        std::string goto_keywords(
            std::ifstream& file, absl::Span< const absl::string_view > words )
        {
            std::string line;
            while( std::getline( file, line ) )
            {
                for( const auto word : words )
                {
                    if( string_starts_with( line, word ) )
                    {
                        return line;
                    }
                }
            }
            throw geode::OpenGeodeException{
                "[goto_keywords] Cannot find one of the requested keywords"
            };
            return "";
        }

        std::string read_name( std::istringstream& iss )
        {
            std::string name;
            iss >> name;
            std::string token;
            while( iss >> token )
            {
                absl::StrAppend( &name, " ", token );
            }
            return absl::StrReplaceAll( name, { { "\"", "" } } );
        }

        absl::optional< TSurfData > read_tsurf( std::ifstream& file )
        {
            if( !line_starts_with( file, "GOCAD TSurf" ) )
            {
                return absl::nullopt;
            }
            TSurfData tsurf;
            tsurf.header = read_header( file );
            tsurf.crs = read_CRS( file );
            read_tfaces( file, tsurf );
            return tsurf;
        }

        RegionSurfaceSide determine_surface_to_regions_sides( const BRep& brep )
        {
            RegionSurfaceSide result;
            const auto paired_signs = determine_paired_signs( brep );
            result.universe_surface_sides =
                determine_universe_sides( brep, paired_signs );
            result.regions_surface_sides =
                determine_regions_sides( brep, paired_signs );
            return result;
        }
    } // namespace detail
} // namespace geode
