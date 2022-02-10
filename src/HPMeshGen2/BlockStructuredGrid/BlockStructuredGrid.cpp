#include "BlockStructuredGrid.h"

#include <algorithm>
#include <array>

#include <iostream>

#include <glog/logging.h>

#include "../svg.h"

inline void quad2tri(PolyMesh& quadMesh, TriMesh& triMesh) {
	LOG_ASSERT( triMesh.n_faces() == 0 );

	// add vertices to triangle mesh
	for ( const auto& vh : quadMesh.vertices() ) {
		PolyMesh::Point p = quadMesh.point( vh );
		auto vht = triMesh.add_vertex( p );
		triMesh.data( vht ).is_feature = quadMesh.data( vh ).is_feature;
	}

	// add faces to triangle mesh
	for ( auto fh : quadMesh.faces() ) {
		// find best configuration
		std::vector<PolyMesh::VertexHandle> faceVertices;
		
		// get all vertices of the polygon
		for ( auto fv : quadMesh.fv_range( fh ) ) {
			faceVertices.push_back( fv );
		}

		// make sure polygon is a quad
		if ( faceVertices.size() != 4 ) {
			LOG( ERROR ) << "Cannot triangulate mesh as it is not a quad mesh";
		}

		std::vector<PolyMesh::VertexHandle> f1_vertices = { faceVertices[0], faceVertices[1] , faceVertices[3] };
		std::vector<PolyMesh::VertexHandle> f2_vertices = { faceVertices[2], faceVertices[3] , faceVertices[1] };

		// store faces
		triMesh.add_face( f1_vertices );
		triMesh.add_face( f2_vertices );
	}
}

BlockStructuredGrid::BlockStructuredGrid( const PolyMesh & fm, const size_t nRefinementSteps, const Contour::Shape& contour, const OceanMesh& oceanMesh, const bool positionBoundary )
	: nRefinementSteps_{ nRefinementSteps }, nGridNodes_( ( 1 << nRefinementSteps ) + 1 ), fragmentMesh{ fm }, quadMesh{ fm }, contour_{ contour }, oceanMesh_{ oceanMesh }
{
	using std::vector;
			
	edgeVertices_.resize( fragmentMesh.n_edges());
	fragmentVertices_.resize( fragmentMesh.n_faces());

	/*
					  edgeNorth
			  0     1     2     3     4
			  0 *-----------------------* 4
			  1 |   (2,0) (2,1) (2,2)   | 3
	edgeWest  2 |   (1,0) (1,1) (1,2)   | 2  edgeEast
			  3 |   (0,0) (0,1) (0,2)   | 1
			  4 *-----------------------* 0
			  0     1     2     3    4
					  edgeSouth
	*/

	// add level 1
	// vertices on edge
	for(const auto& eh : fragmentMesh.edges() ) {
		const auto& heh = eh.h0();
		edgeVertices_[eh.idx()] = {heh.from(), heh.to()};
	}
	// vertices on face
	for( const auto& fh : fragmentMesh.faces() ) {
		auto hehSouth = fh.halfedge();
		auto hehEast = hehSouth.next();
		auto hehNorth = hehEast.next();
		auto hehWest = hehNorth.next();

		auto vhSW = hehSouth.from();
		auto vhSE = hehSouth.to();
		auto vhNE = hehNorth.from();
		auto vhNW = hehNorth.to();

		// turn block s.t. triangles have best quality
		const auto& pSW = quadMesh.point( vhSW );
		const auto& pSE = quadMesh.point( vhSE );
		const auto& pNE = quadMesh.point( vhNE );
		const auto& pNW = quadMesh.point( vhNW );

		const auto angSW = HelperFunctions::calcAngle( pSE, pSW, pNW );
		const auto angSE = HelperFunctions::calcAngle( pNE, pSE, pSW );
		const auto angNE = HelperFunctions::calcAngle( pNW, pNE, pSE );
		const auto angNW = HelperFunctions::calcAngle( pSW, pNW, pNE );

		auto maxCurr = std::max( angSW, angNE );
		auto maxNew = std::max( angSE, angNW );

		bool prohibitTurn = false;
		
		// prohibit turn if new quad would be invalid
		if( maxNew < maxCurr && angSW + angSE + angNE + angNW > 360.f ) {
			LOG( WARNING ) << "A quad turned was prohibited by angle condition. fh = " << fh.idx();
			prohibitTurn = true;
		}

		if( !prohibitTurn && maxNew < maxCurr ) {
			// turn element
			const auto hehBuf = hehSouth;
			hehSouth = hehEast, hehEast = hehNorth, hehNorth = hehWest, hehWest = hehBuf;
			vhSW = hehSouth.from();
			vhSE = hehSouth.to();
			vhNE = hehNorth.from();
			vhNW = hehNorth.to();
		}
		fragmentVertices_[fh.idx()] = { {vhSW, vhSE}, {vhNW, vhNE} };
	}

	for( int i = 0; i < nRefinementSteps_; ++i ) {
		refine();

		if( positionBoundary ) {
			// adapt contour (not stable in complex regions due to nearest neighbor search)
			for( const auto& e : fragmentMesh.edges() ) {
				if( !e.is_boundary() )
					continue;

				auto h = e.halfedge( 1 );
				if( !h.is_boundary() )
					h = e.halfedge( 0 );

				const auto& vBreak = h.next().to();
				const auto& pBreak = quadMesh.point( vBreak );

				const auto& ehVec = edgeVertices_[e.idx()];
				for( int i = 1; i < ehVec.size() - 1; i += 2 ) {
					const auto& vMid = ehVec[i];
					const auto& vLeft = ehVec[i + 1];
					const auto& vRight = ehVec[i - 1];
					const auto& pMid = quadMesh.point( vMid );
					const auto& pLeft = quadMesh.point( vLeft );
					const auto& pRight = quadMesh.point( vRight );
					//const auto& pNew = contour_.computeMid( pLeft, pRight );
					const auto& pNew = contour_.computeMidWithout( pLeft, pRight, pBreak );
					quadMesh.set_point( vMid, pNew );
				}
			}
		}
		
	}
		
	quadMesh.request_face_status();
	quadMesh.request_edge_status();
	quadMesh.request_vertex_status();

	// delete old faces
	for (unsigned i = 0; i < quadMesh.n_faces(); ++i) {
		PolyMesh::FaceHandle fh = quadMesh.face_handle(i);
		quadMesh.delete_face(fh, false);
	}
	quadMesh.garbage_collection();

	// create new faces
	fragmentFaces_ = vector<vector<vector<TriMesh::FaceHandle>>>( fragmentVertices_.size(), vector<vector<TriMesh::FaceHandle>>( nGridNodes_ - 1, vector<TriMesh::FaceHandle>( nGridNodes_ - 1 ) ) );
	for (auto fragId = 0; fragId < fragmentVertices_.size(); ++fragId) {
		for (auto i = 0; i < fragmentVertices_[fragId].size() - 1; ++i) {
			for (auto j = 0; j < fragmentVertices_[fragId][i].size() - 1; ++j) {
				vector<PolyMesh::VertexHandle> fevh(4);
				fevh[0] = fragmentVertices_[fragId][i][j];
				fevh[1] = fragmentVertices_[fragId][i][j + 1];
				fevh[2] = fragmentVertices_[fragId][i + 1][j + 1];
				fevh[3] = fragmentVertices_[fragId][i + 1][j];
				fragmentFaces_[fragId][i][j] = quadMesh.add_face(fevh);
			}
		}
	}

	// create triangle mesh
	quad2tri(quadMesh, triMesh);

	// add vertex relation
	OpenMesh::VPropHandleT<OpenMesh::SmartVertexHandle> propQuadVh;
	triMesh.add_property( propQuadVh, "propQuadVh" );
	for( const auto& vhTri : triMesh.vertices() ) {
		OpenMesh::SmartVertexHandle vhQuad( vhTri.idx(), &quadMesh );
		triMesh.property( propQuadVh, vhTri ) = vhQuad;
	}
}

PolyMesh BlockStructuredGrid::get_patch(size_t patchID)
{
	LOG_ASSERT( patchID < fragmentMesh.n_faces() );
	std::vector<std::vector<PolyMesh::Point>> points(fragmentVertices_[patchID].size(), std::vector<PolyMesh::Point>(fragmentVertices_[patchID][0].size()));

	for (size_t i = 0; i < fragmentVertices_[patchID].size(); ++i) {
		for (size_t j = 0; j < fragmentVertices_[patchID][i].size(); ++j) {
			points[i][j] = quadMesh.point(fragmentVertices_[patchID][i][j]);
		}
	}

	PolyMesh patch;

	std::vector<std::vector<PolyMesh::VertexHandle>> vertices(points.size(), std::vector<PolyMesh::VertexHandle>(points[0].size()));

	// add vertices
	for (size_t i = 0; i < points.size(); ++i) {
		for (size_t j = 0; j < points[i].size(); ++j) {
			vertices[i][j] = patch.add_vertex(points[i][j]);
		}
	}

	// add faces
	for (size_t i = 0; i < vertices.size() - 1; ++i) {
		for (size_t j = 0; j < vertices[i].size() - 1; ++j) {
			std::vector<PolyMesh::VertexHandle> fevh(4);
			fevh[0] = vertices[i][j];
			fevh[1] = vertices[i][j + 1];
			fevh[2] = vertices[i + 1][j + 1];
			fevh[3] = vertices[i + 1][j];
			patch.add_face(fevh);
		}
	}

	return patch;
}

void BlockStructuredGrid::copy_positions_tri2quad()
{
	OpenMesh::VPropHandleT<OpenMesh::SmartVertexHandle> propQuadVh;
	triMesh.get_property_handle( propQuadVh, "propQuadVh" );
	
	for( const auto& vhTri : triMesh.vertices() ) {
		const auto& vhQuad = triMesh.property( propQuadVh, vhTri );
		quadMesh.set_point( vhQuad, triMesh.point( vhTri ) );
	}

	//for (uint i = 0; i < triMesh.n_vertices(); ++i) {
	//	TriMesh::VertexHandle vhTri = triMesh.vertex_handle(i);
	//	PolyMesh::VertexHandle vhPoly = quadMesh.vertex_handle(i);
	//
	//	quadMesh.set_point(vhPoly, triMesh.point(vhTri));
	//}
}

void BlockStructuredGrid::copy_positions_quad2tri()
{
	LOG( WARNING ) << " This function is deprecated and no longer supported!";
	for (uint i = 0; i < quadMesh.n_vertices(); ++i) {
		TriMesh::VertexHandle vhTri = triMesh.vertex_handle(i);
		PolyMesh::VertexHandle vhPoly = quadMesh.vertex_handle(i);

		triMesh.set_point(vhTri, quadMesh.point(vhPoly));
	}
}

void BlockStructuredGrid::exportMeshAdcirc( const std::experimental::filesystem::path& folder, const std::string& filename, OceanMesh& oceanMesh ) {
	namespace fs = std::experimental::filesystem;
	const fs::path file14 = folder / filename;
	write14file( file14 );

	fs::path file15 = file14;
	file15.replace_extension( "15" );
	fs::path ocean15 = oceanMesh.filename();
	ocean15.replace_extension( "15" );
	if( fs::exists( ocean15 ) ) {
		fs::copy_file( ocean15, file15, fs::copy_options::overwrite_existing );
	}
	
	fs::path file17 = file14;
	file17.replace_extension( "17" );
	write17file( file17, oceanMesh );
}

void BlockStructuredGrid::patch_segmentation(const size_t n_blocks)
{
	
	size_t n_frags_per_block = (fragmentMesh.n_faces() + n_blocks - 1) / n_blocks;

	blocks_.resize(n_blocks);
	for (size_t i = 0; i < n_blocks; ++i) {
		std::vector<size_t>b;
		for (size_t j = 0; j < n_frags_per_block && n_frags_per_block * i + j < fragmentMesh.n_faces(); ++j) {
			b.push_back(n_frags_per_block * i + j);
		}
		blocks_[i] = b;
	}
}

void BlockStructuredGrid::write_bsg(const std::experimental::filesystem::path& folder, const std::string & mesh_name, const bool& writeMasks )
{
	if (blocks_.size() == 0) {
		LOG( INFO ) << "no blocks defined. Using each fragment as block.";
		patch_segmentation(fragmentMesh.n_faces());
	}

	BsgName_ = mesh_name;
	BsgFolder_ = folder;

	std::experimental::filesystem::create_directories( BsgFolder_ );

	write_mesh( writeMasks );
}

void BlockStructuredGrid::print_svg(std::string filename)
{
	std::ofstream ofs(filename);

	PolyMesh mesh = quadMesh;

	// get min/max
	float x_min = std::numeric_limits<float>::max();
	float y_min = std::numeric_limits<float>::max();
	float x_max = -std::numeric_limits<float>::max();
	float y_max = -std::numeric_limits<float>::max();

	for ( auto vh : mesh.vertices() ) {
		TriMesh::Point p = mesh.point(vh);

		x_min = std::min(x_min, p[0]);
		x_max = std::max(x_max, p[0]);
		y_min = std::min(y_min, p[1]);
		y_max = std::max(y_max, p[1]);
	}

	const float x_dim = x_max - x_min;
	const float y_dim = y_max - y_min;

	const int xpx = 1000;
	const int ypx = (int)(xpx * y_dim / x_dim);

	// norm all vertices
	for ( auto vh : mesh.vertices() ) {
		TriMesh::Point p = mesh.point( vh );
		p[0] = xpx * ( p[0] - x_min ) / x_dim;
		p[1] = ypx - ypx * ( p[1] - y_min ) / y_dim;
		p[2] = 0;

		mesh.set_point( vh, p );
	}

	// header
	ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
	ofs << "<svg xmlns=\"http://www.w3.org/2000/svg\"" << std::endl;
	ofs << "version=\"1.1\" baseProfile=\"full\"" << std::endl;
	ofs << "width=\"" << xpx << "px\" height=\"" << ypx << "px\" viewBox=\"0 0 " << xpx << " " << ypx << "\">" << std::endl;
	ofs << std::endl;

	// body
	for ( const auto& fragment : fragmentVertices_ ) {
		for ( int i = 0; i < nGridNodes_ - 1; ++i ) {
			for ( int j = 0; j < nGridNodes_ - 1; ++j ) {
				// 3 4
				// 1 2
				auto p1 = mesh.point( fragment[i][j] );
				auto p2 = mesh.point( fragment[i + 1][j] );
				auto p3 = mesh.point( fragment[i][j + 1] );
				auto p4 = mesh.point( fragment[i + 1][j + 1] );
				svg::svgLine( ofs, p3, p2, "black", 1 );	// triangle edge

				if(i != 0)
					svg::svgLine( ofs, p1, p3, "black", 1 );	// quad edge
				if ( j != 0 )
					svg::svgLine( ofs, p1, p2, "black", 1 );	// quad edge
			}
		}
	}
	for ( const auto& edge : edgeVertices_ ) {
		for ( int i = 0; i < nGridNodes_ - 1; ++i ) {
			auto p1 = mesh.point( edge[i] );
			auto p2 = mesh.point( edge[i + 1] );
			svg::svgLine( ofs, p1, p2, "red", 2 );	// block edge
		}
	}

	ofs << "</svg>" << std::endl;

	ofs.close();
}

void BlockStructuredGrid::write14file( const std::experimental::filesystem::path& file14 ) {
	std::ofstream ofs( file14 );
	ofs.precision( 5 );
	ofs.setf( std::ios::fixed, std::ios::floatfield );

	ofs << file14.stem() << "\t! This file was created by HPMeshGen" << std::endl;
	ofs << triMesh.n_faces() << " " << triMesh.n_vertices() << "\t! NE,NP - NUMBER OF ELEMENTS AND NUMBER OF NODAL POINTS" << std::endl;

	// write out vertices
	for( unsigned int i = 0; i < triMesh.n_vertices(); ++i ) {
		TriMesh::Point p = triMesh.point( triMesh.vertex_handle( i ) );
		ofs << i + 1 << "   " << p[0] << "    " << p[1] << "    " << p[2] << std::endl;
	}
	// write out faces
	for( unsigned int i = 0; i < triMesh.n_faces(); ++i ) {
		TriMesh::FaceHandle fh = triMesh.face_handle( i );
		ofs << i + 1 << " " << 3;
		for( TriMesh::FaceVertexIter fv_it = triMesh.fv_begin( fh ); fv_it != triMesh.fv_end( fh ); ++fv_it ) {
			ofs << " " << fv_it->idx() + 1;
		}
		ofs << std::endl;
	}

	// boundaries
	std::vector<OpenMesh::SmartEdgeHandle> ehVec;
	for( const auto& e : triMesh.edges() ) {
		if( e.is_boundary() ) {
			ehVec.push_back( e );
		}
	}
	std::vector<bool> isStored( triMesh.n_edges(), false );					// bool checks if edge is stored in vector (avoid boundaries to appear several times)

	std::vector<std::vector<TriMesh::VertexHandle>> landBoundaries;
	std::vector<std::vector<TriMesh::VertexHandle>> islandBoundaries;
	std::vector<std::vector<TriMesh::VertexHandle>> openseaBoundaries;
	std::vector<std::vector<TriMesh::VertexHandle>> riverBoundaries;
	std::vector<std::vector<TriMesh::VertexHandle>> radiationBoundaries;

	auto getContourType = [this]( const OpenMesh::SmartHalfedgeHandle& h ) {
		return contour_.getContourType( triMesh.point( h.from() ), triMesh.point( h.to() ) );
	};

	// collect boundaries
	for( size_t i = 0; i < ehVec.size(); ++i ) {

		if( !isStored[ehVec[i].idx()] ) {
			// this is an edge which is not stored yet

			auto heh = ehVec[i].h0();
			if( !heh.is_boundary() ) {
				heh = heh.opp();
			}
			auto hehInit = heh;

			// get next heh
			auto hehNext = heh.next();

			auto contourType = getContourType( heh );
			LOG_ASSERT( contourType != Contour::ContourType::NOTHING ) << "Boundary without contour type found " << triMesh.point( heh.from() ) << "  " << triMesh.point( heh.to() );


			// check if the next heh is the same contour type. If not, the first edge of the boundary was found
			// additionally prevent endless loop in case the boundary is a closed cycle
			while( getContourType( hehNext ) == contourType && hehNext != hehInit ) {
				heh = hehNext;
				isStored[heh.edge().idx()] = true;
				// get next heh
				hehNext = heh.next();
			}

			// heh is now the first halfedge_handle of the boundary
			hehInit = heh;

			// first vertex
			auto vhBegin = heh.to();

			// find last vertex
			hehNext = heh.prev();
			while( getContourType( hehNext ) == contourType && hehNext != hehInit ) {
				heh = hehNext;
				isStored[triMesh.edge_handle( heh ).idx()] = true;
				// get prev heh
				hehNext = heh.prev();
			}

			// last vertex
			auto vhEnd = heh.from();

			// collect vertices on triMesh

			// find halfedge pointing towards vhBegin
			for( const auto& vih : vhBegin.incoming_halfedges() ) {
				if( vih.is_boundary() ) {
					heh = vih;
					break;
				}
			}


			std::vector<TriMesh::VertexHandle> boundary;
			boundary.push_back( heh.to() );
			do {
				boundary.push_back( heh.from() );
				heh = heh.prev();
			} while( heh.to() != vhEnd );

			// set boundary
			switch( contourType ) {
			case Contour::ContourType::LAND:
				landBoundaries.push_back( boundary );
				break;
			case Contour::ContourType::ISLAND:
				islandBoundaries.push_back( boundary );
				break;
			case Contour::ContourType::OPEN_SEA:
				openseaBoundaries.push_back( boundary );
				break;
			case Contour::ContourType::RIVER:
				riverBoundaries.push_back( boundary );
				break;
			case Contour::ContourType::RADIATION:
				radiationBoundaries.push_back( boundary );
				break;
			default:
				LOG( ERROR ) << "Unknown contour type in fort.14 generation. contourType = " << contourType;
				break;
			}

		}
	}

	// open up ring-boundary / delete last vertex if it is the same as the first
	for( auto& boundary : landBoundaries ) {
		if( !boundary.empty() && boundary[0].idx() == boundary.back().idx() ) {
			boundary.pop_back();
		}
	}
	for( auto& boundary : islandBoundaries ) {
		if( !boundary.empty() && boundary[0].idx() == boundary.back().idx() ) {
			boundary.pop_back();
		}
	}

	// write boundaries to file
	// Open Sea
	ofs << openseaBoundaries.size() + radiationBoundaries.size() + riverBoundaries.size() << "\t! NOPE - TOTAL NUMBER OF OPEN BOUNDARY FORCING SEGMENTS" << std::endl;
	size_t nNeumannVertices = 0;
	for( size_t i = 0; i < openseaBoundaries.size(); ++i ) {
		nNeumannVertices += openseaBoundaries[i].size();
	}
	for( size_t i = 0; i < riverBoundaries.size(); ++i ) {
		nNeumannVertices += riverBoundaries[i].size();
	}
	for( size_t i = 0; i < radiationBoundaries.size(); ++i ) {
		nNeumannVertices += radiationBoundaries[i].size();
	}
	ofs << nNeumannVertices << "\t! NETA - TOTAL NUMBER OF OPEN BOUNDARY NODES" << std::endl;
	// open sea
	for( size_t i = 0; i < openseaBoundaries.size(); ++i ) {
		ofs << openseaBoundaries[i].size();
		ofs << "\t! NUMBER OF NODES ON OPEN BOUNDERY FORCING SEGMENT NO." << ( i + 1 ) << std::endl;

		for( size_t j = 0; j < openseaBoundaries[i].size(); ++j ) {
			ofs << openseaBoundaries[i][j].idx() + 1;
			if( j == 0 )
				ofs << "\t!BEGINING NODE NUMBER ON OPEN BOUNDARY NO. " << ( i + 1 );
			else if( j == openseaBoundaries[i].size() - 1 ) ofs << "\t!ENDING NODE NUMBER ON OPEN BOUNDARY NO. " << ( i + 1 );
			ofs << std::endl;
		}
	}
	// radiation
	for( size_t i = 0; i < radiationBoundaries.size(); ++i ) {
		ofs << radiationBoundaries[i].size();
		ofs << "\t! NUMBER OF NODES ON OPEN BOUNDERY FORCING SEGMENT NO." << ( i + 1 ) << std::endl;

		for( size_t j = 0; j < radiationBoundaries[i].size(); ++j ) {
			ofs << radiationBoundaries[i][j].idx() + 1;
			if( j == 0 )
				ofs << "\t!BEGINING NODE NUMBER ON OPEN BOUNDARY NO. " << ( i + 1 );
			else if( j == radiationBoundaries[i].size() - 1 ) ofs << "\t!ENDING NODE NUMBER ON OPEN BOUNDARY NO. " << ( i + 1 );
			ofs << std::endl;
		}
	}
	// river
	for( size_t i = 0; i < riverBoundaries.size(); ++i ) {
		ofs << riverBoundaries[i].size();
		ofs << " 100 \t! NUMBER OF NODES ON OPEN BOUNDERY FORCING SEGMENT NO." << ( i + 1 ) << std::endl;

		for( size_t j = 0; j < riverBoundaries[i].size(); ++j ) {
			ofs << riverBoundaries[i][j].idx() + 1;
			if( j == 0 )
				ofs << "\t!BEGINING NODE NUMBER ON OPEN BOUNDARY NO. " << ( i + 1 );
			else if( j == riverBoundaries[i].size() - 1 ) ofs << "\t!ENDING NODE NUMBER ON OPEN BOUNDARY NO. " << ( i + 1 );
			ofs << std::endl;
		}
	}

	ofs << landBoundaries.size() + islandBoundaries.size() << "\t! NOPE - NBOU - TOTAL NUMBER OF LAND BOUNDARY SEGMENTS INCLUDING ISLANDS SEGMENTS" << std::endl;
	size_t nDirichletVertices = 0;
	for( size_t i = 0; i < landBoundaries.size(); ++i ) {
		nDirichletVertices += landBoundaries[i].size();
	}
	for( size_t i = 0; i < islandBoundaries.size(); ++i ) {
		nDirichletVertices += islandBoundaries[i].size();
	}
	// land
	ofs << nDirichletVertices << "\t! NOPE - TOTAL NUMBER OF LAND BOUNDARY NODES" << std::endl;
	for( size_t i = 0; i < landBoundaries.size(); ++i ) {
		ofs << landBoundaries[i].size();
		ofs << " 0";	// land boundary
		ofs << "\t! NETA - NUMBER OF NODES ON LAND BOUNDERY SEGMENT NO." << ( i + 1 ) << " AND BOUNDARY TYPE" << std::endl;

		for( size_t j = 0; j < landBoundaries[i].size(); ++j ) {
			ofs << landBoundaries[i][j].idx() + 1;
			if( j == 0 )
				ofs << "\t!BEGINING NODE NUMBER ON LAND BOUNDARY NO. " << ( i + 1 );
			else if( j == landBoundaries[i].size() - 1 )
				ofs << "\t!ENDING NODE NUMBER ON LAND BOUNDARY NO. " << ( i + 1 );
			ofs << std::endl;
		}
	}
	// island
	for( size_t i = 0; i < islandBoundaries.size(); ++i ) {
		ofs << islandBoundaries[i].size();
		ofs << " 1";	// island boundary
		ofs << "\t! NETA - NUMBER OF NODES ON ISLAND BOUNDERY SEGMENT NO." << ( i + 1 ) << " AND BOUNDARY TYPE" << std::endl;

		for( size_t j = 0; j < islandBoundaries[i].size(); ++j ) {
			ofs << islandBoundaries[i][j].idx() + 1;
			if( j == 0 )
				ofs << "\t!BEGINING NODE NUMBER ON ISLAND BOUNDARY NO. " << ( i + 1 );
			else if( j == islandBoundaries[i].size() - 1 )
				ofs << "\t!ENDING NODE NUMBER ON ISLAND BOUNDARY NO. " << ( i + 1 );
			ofs << std::endl;
		}
	}

	ofs.close();
}

void BlockStructuredGrid::write17file( const std::experimental::filesystem::path& file17, OceanMesh& oceanMesh ) {
	namespace fs = std::experimental::filesystem;
	auto idx = []( auto e ) {return e.idx() + 1; };

	fs::create_directories( file17.parent_path() );
	std::ofstream ofs( file17 );

	if( !ofs.is_open() ) {
		LOG( ERROR ) << "Cannot open fort17 file for generation for fort_bahamas.";
	}

	// all edges edges and incident nodes and elements
	ofs << triMesh.n_edges() << "    ! Edge connectivity list: node(1), node(2), element(1), element(2)\n";
	for( const auto& eh : triMesh.edges() ) {
		const auto& heh = !eh.h0().is_boundary() ? eh.h0() : eh.h1();
		const auto& vh1 = heh.from();
		const auto& vh2 = heh.to();
		const auto& fh1 = heh.face();
		const auto& fh2 = heh.opp().face();

		ofs << idx( eh ) << " " << idx( vh1 ) << " " << idx( vh2 ) << " " << idx( fh1 ) << " " << idx( fh2 ) << std::endl;
	}

	// all elements and incident edges
	for( auto fh : triMesh.faces() ) {
		std::vector<TriMesh::EdgeHandle> feh;
		for( auto fe_it = triMesh.fe_iter( fh ); fe_it.is_valid(); ++fe_it ) {
			feh.push_back( *fe_it );
		}

		ofs << idx( fh ) << " " << idx( feh[0] ) << " " << idx( feh[1] ) << " " << idx( feh[2] ) << std::endl;
	}

	std::vector<OpenMesh::SmartEdgeHandle> interior, land, river, radiation, sea;
	//std::vector<bool> isSea( triMesh.n_vertices(), false );
	for( const auto& eh : triMesh.edges() ) {

		auto heh = eh.h0();
		auto vh1 = heh.from();
		auto vh2 = heh.to();

		//if( !eh.is_boundary() ) {
		//	interior.push_back( eh );
		//} else if( contour_.getContourType( triMesh.point( vh1 ), triMesh.point( vh2 ) ) == Contour::ContourType::OPEN_SEA ) {
		//	sea.push_back( eh );
		//	isSea[vh1.idx()] = true;
		//	isSea[vh2.idx()] = true;
		//} else {
		//	land.push_back( eh );
		//}

		if( !eh.is_boundary() ) {
			interior.push_back( eh );
			continue;
		}
		switch( contour_.getContourType( triMesh.point( vh1 ), triMesh.point( vh2 ) ) ) {
		case Contour::ContourType::LAND:
			land.push_back( eh );
			break;
		case Contour::ContourType::ISLAND:
			land.push_back( eh );
			break;
		case Contour::ContourType::RADIATION:
			radiation.push_back( eh );
			break;
		case Contour::ContourType::RIVER:
			river.push_back( eh );
			break;
		case Contour::ContourType::OPEN_SEA:
			sea.push_back( eh );
			break;
		default:
			LOG( ERROR ) << "Unknown contour type for edge " << eh.idx() << " pos : " << triMesh.point( vh1 ) << " , " << triMesh.point( vh2 );
			break;
		}
	}

	// interior edges
	ofs << interior.size() << "    ! NUMBER OF INTERIOR EDGES" << std::endl;
	for( int i = 0; i < interior.size(); ++i ) {
		ofs << i + 1 << " " << idx( interior[i] ) << std::endl;
	}
	// land edges
	ofs << land.size() << "    ! NUMBER OF LAND EDGES" << std::endl;
	for( int i = 0; i < land.size(); ++i ) {
		ofs << i + 1 << " " << idx( land[i] ) << std::endl;
	}
	// radiation edges x = 70
	ofs << radiation.size() << "    ! NUMBER OF RADIATION EDGES" << std::endl;
	for( int i = 0; i < radiation.size(); ++i ) {
		ofs << i + 1 << " " << idx( radiation[i] ) << std::endl;
	}
	// river edges x = -20
	ofs << river.size() << "    ! NUMBER OF RIVER EDGES" << std::endl;
	for( int i = 0; i < river.size(); ++i ) {
		ofs << i + 1 << " " << idx( river[i] ) << std::endl;
	}
	// sea edges
	ofs << sea.size() << "    ! NUMBER OF SEA EDGES" << std::endl;
	for( int i = 0; i < sea.size(); ++i ) {
		ofs << i + 1 << " " << idx( sea[i] ) << std::endl;
	}

	// EMO/EFA values
	for( auto i = 0; i < oceanMesh.nEmoEfaVals; ++i ) {
		for( const auto& edge : sea ) {
			const auto& p0 = triMesh.point( edge.h0().from() );
			const auto& p1 = triMesh.point( edge.h0().to() );
			const auto& [emo, efa] = contour_.getEmoEfaVals( p0, p1 );
			ofs << std::setprecision( 10 ) << std::fixed << emo[i] << " " << efa[i] << std::endl;
		}
	}

	ofs.close();
}

///////////////////////
// private functions //

void BlockStructuredGrid::write_mesh( const bool& writeMasks )
{
	auto bsgFile = BsgFolder_ / ( BsgName_ + std::string( ".bsg" ) );
	std::ofstream ofs(bsgFile);

	if (!ofs.is_open()) {
		LOG( WARNING ) << "Could not write BSG file '" << bsgFile << "'";
		return;
	}

	ofs << "# " << BsgName_ << std::endl;
	ofs << "n_blocks " << blocks_.size() << std::endl;
	ofs << "n_fragments " << fragmentMesh.n_faces() << std::endl;
	ofs << "n_grid_nodes " << nGridNodes_ << std::endl;
	
	ofs.close();

	for (auto i = 0; i < blocks_.size(); ++i) {
		write_block(i, writeMasks );
	}
}

void BlockStructuredGrid::write_block(size_t blockID, const bool& writeMasks )
{
	std::ofstream ofs(BsgFolder_ / (std::string("b") + std::to_string(blockID) + "_" + std::to_string( nGridNodes_ ) + std::string(".block")));
	ofs.precision(5);
	ofs.setf(std::ios::fixed, std::ios::floatfield);

	ofs << "# " << BsgName_ << std::endl;
	ofs << "block_id " << blockID << std::endl;
	ofs << "n_fragments_local " << blocks_[blockID].size() << std::endl;
	ofs << "n_grid_nodes " << nGridNodes_ << std::endl;

	for (size_t i = 0; i < blocks_[blockID].size(); ++i) {
		write_fragment_head(ofs, blockID, i);
	}
	for (size_t i = 0; i < blocks_[blockID].size(); ++i) {
		write_fragment_nodes(ofs, blockID, i);
	}

	ofs.close();

	write_depth( blockID );
	write_emo_efa( blockID );
	if( writeMasks )
		write_masks( blockID );
}

void BlockStructuredGrid::write_fragment_head( std::ofstream& ofs, const size_t blockID, const size_t commID)
{
	const size_t fragId = blocks_[blockID][commID];

	ofs << "# -------------------- comm_id " << commID << " --------------------" << std::endl;
	ofs << "fragment_id " << fragId << std::endl;

	// neighbors
	TriMesh::VertexHandle vhSW = fragmentVertices_[fragId][0][0];

	// find halfedge inside current fragment pointing away from vhSW
	PolyMesh::HalfedgeHandle heh;
	for ( auto voh : fragmentMesh.voh_range( vhSW ) ) {
		if ( fragmentMesh.face_handle( voh ).idx() == fragId ) {
			heh = voh;
			break;
		}
	}

	// collect all neighbors
	std::array<int, 4> neighbors;
	for (size_t i = 0; i < neighbors.size(); ++i) {
		neighbors[i] = fragmentMesh.opposite_face_handle(heh).idx();
		heh = fragmentMesh.next_halfedge_handle(heh);
	}

	// Find blockID and commID of neighboring fragments
	std::vector<int>neighbor_blockId(4, -1);
	std::vector<int>neighbor_commId(4, -1);
	std::vector<int>neighbor_fragId(4, -1);
	for (int i = 0; i < 4; ++i) {
		if (neighbors[i] == -1)
			continue;
		for (size_t blockId = 0; blockId < blocks_.size(); ++blockId) {
			ptrdiff_t commId = std::find(blocks_[blockId].begin(), blocks_[blockId].end(), neighbors[i]) - blocks_[blockId].begin();
			if (commId < blocks_[blockId].size()) {
				neighbor_blockId[i] = blockId;
				neighbor_commId[i] = commId;
				neighbor_fragId[i] = blocks_[blockId][commId];
				break;
			}
		}
	}

	// Collect BCs and store in neighbor_fragID
	for (size_t i = 0; i < neighbors.size(); ++i) {
		if (neighbors[i] == -1) {
			PolyMesh::EdgeHandle eh = fragmentMesh.edge_handle(heh);

			///////////////////////////////
			// use contour here ...
			const auto& p1 = fragmentMesh.point( fragmentMesh.from_vertex_handle( heh ) );
			const auto& p2 = fragmentMesh.point( fragmentMesh.to_vertex_handle( heh ) );
			const auto contourType = contour_.getContourType( p1, p2 );
			

			// boundaries have internal values
			//auto bid = fragmentMesh.property( boundaryID_, eh );
			//int boundaryType;
			//if ( fragmentMesh.property( isNeumann_, eh ) ) 
			//	boundaryType = neumannBoundaryTypes_[bid];
			//else 
			//	boundaryType = dirichletBoundaryTypes_[bid];
			//if ( boundaryType < 0 ) boundaryType = 0;
			//
			//neighbor_fragId[i] = -100 - boundaryType;
			//if (fragmentMesh.property(isNeumann_, eh))
			//	neighbor_fragId[i] -= 1000;
			//
			//LOG_ASSERT( neighbor_fragId[i] == contourType );	// TODO remove assert and boundary type stuff after testing

			neighbor_fragId[i] = contourType;
		}
		heh = fragmentMesh.next_halfedge_handle(heh);
	}

	// Find boundary type and store it in neighbor_fragId
	for ( auto& n : neighbors ) {
		n = fragmentMesh.opposite_face_handle( heh ).idx();
		heh = fragmentMesh.next_halfedge_handle(heh);
	}

	// Find neighbors edge
	std::vector<std::string>neighbor_edge(4, "X");
	for (size_t i = 0; i < neighbors.size(); ++i) {
		if (neighbors[i] == -1)
			continue;

		// neighbors
		TriMesh::VertexHandle vhSW = fragmentVertices_[neighbors[i]][0][0];

		// find halfedge inside current fragment pointing away from vhSW
		PolyMesh::HalfedgeHandle nheh;
		for ( auto voh : fragmentMesh.voh_range( vhSW ) ) {
			if ( fragmentMesh.face_handle( voh ).idx() == neighbors[i] ) {
				nheh = voh;
				break;
			}
		}

		const std::string directions = "SENW";
		for (size_t j = 0; j < neighbors.size(); ++j) {
			if (fragId == fragmentMesh.opposite_face_handle(nheh).idx()) {
				neighbor_edge[i] = directions[j];
				break;
			}
			nheh = fragmentMesh.next_halfedge_handle(nheh);
		}
	}

	ofs << "# neighboring blocks: S/E/N/W  blockID  commID neighborOrientation neighborFragmentId" << std::endl;
	ofs << "W " << neighbor_blockId[3] << " " << neighbor_commId[3] << " " << neighbor_edge[3] << " " << neighbor_fragId[3] << std::endl;
	ofs << "E " << neighbor_blockId[1] << " " << neighbor_commId[1] << " " << neighbor_edge[1] << " " << neighbor_fragId[1] << std::endl;
	ofs << "S " << neighbor_blockId[0] << " " << neighbor_commId[0] << " " << neighbor_edge[0] << " " << neighbor_fragId[0] << std::endl;
	ofs << "N " << neighbor_blockId[2] << " " << neighbor_commId[2] << " " << neighbor_edge[2] << " " << neighbor_fragId[2] << std::endl;

}

void BlockStructuredGrid::write_fragment_nodes(std::ofstream & ofs, const size_t blockID, const size_t commID)
{
	const size_t fragment_id = blocks_[blockID][commID];

	ofs << "# -------------------- comm_id " << commID << " --------------------" << std::endl;

	// nodes within this fragment
	for (size_t i = 0; i < fragmentVertices_[fragment_id].size(); ++i) {
		for (size_t j = 0; j < fragmentVertices_[fragment_id][i].size(); ++j) {
			PolyMesh::Point p = quadMesh.point(fragmentVertices_[fragment_id][i][j]);
			ofs << p[0] << " " << p[1] << std::endl;
		}
	}
}

void BlockStructuredGrid::write_depth(const size_t blockID) const
{
	std::ofstream ofs(BsgFolder_ / (std::string("bath_b") + std::to_string(blockID) + "_" + std::to_string( nGridNodes_ ) + std::string(".txt")));
	ofs.precision(5);
	ofs.setf(std::ios::fixed, std::ios::floatfield);
	for (size_t commID = 0; commID < blocks_[blockID].size(); ++commID) {
		const size_t fragment_id = blocks_[blockID][commID];

		// nodes within this fragment
		for (size_t i = 0; i < fragmentVertices_[fragment_id].size(); ++i) {
			for (size_t j = 0; j < fragmentVertices_[fragment_id][i].size(); ++j) {
				auto p = quadMesh.point(fragmentVertices_[fragment_id][i][j]);
				ofs << p[2] << " ";
			}
			ofs << std::endl;
		}

		ofs << std::endl;
	}

	ofs.close();
}

void BlockStructuredGrid::write_masks( const size_t blockID ) const {
	std::ofstream ofsLower( BsgFolder_ / ( std::string( "maskLower_b" ) + std::to_string( blockID ) + "_" + std::to_string( nGridNodes_ ) + std::string( ".txt" ) ) );
	std::ofstream ofsUpper( BsgFolder_ / ( std::string( "maskUpper_b" ) + std::to_string( blockID ) + "_" + std::to_string( nGridNodes_ ) + std::string( ".txt" ) ) );

	for( size_t commID = 0; commID < blocks_[blockID].size(); ++commID ) {
		const size_t fragment_id = blocks_[blockID][commID];
		// nodes within this fragment
		for( size_t i = 0; i < fragmentVertices_[fragment_id].size() - 1; ++i ) {
			for( size_t j = 0; j < fragmentVertices_[fragment_id][i].size() - 1; ++j ) {
				const auto& vh1 = fragmentVertices_[fragment_id][i][j];
				const auto& vh2 = fragmentVertices_[fragment_id][i][j + 1];
				const auto& vh3 = fragmentVertices_[fragment_id][i + 1][j + 1];

				auto heh = OpenMesh::SmartHalfedgeHandle( quadMesh.halfedge_handle( vh1 ).idx(), &quadMesh );
				for( int k = 0; k < quadMesh.valence( vh1 ); ++k ) {
					if( heh.to() == vh2 )
						break;
					heh = heh.opp();
					heh = heh.next();
				}

				LOG_ASSERT( heh.next().to() == vh3 ) << " could not find correct quad in mask writing";

				const auto& maskL = maskLower_[heh.face().idx()];
				ofsLower << maskL << " ";
				
				const auto& maskU = maskUpper_[heh.face().idx()];
				ofsUpper << maskU << " ";
			}
			ofsLower << std::endl;
			ofsUpper << std::endl;
		}
		ofsLower << std::endl;
		ofsUpper << std::endl;
	}
	ofsLower.close();
	ofsUpper.close();

}

void BlockStructuredGrid::write_emo_efa( const size_t blockID ) const {
	namespace fs = std::experimental::filesystem;

	if( oceanMesh_.nEmoEfaVals < 1 )
		return;

	std::vector<std::ofstream> ofsEMOs( oceanMesh_.nEmoEfaVals );
	decltype( ofsEMOs ) ofsEFAs( ofsEMOs.size() );
	
	for( int i = 0; i < ofsEMOs.size(); ++i )
		ofsEMOs[i].open( BsgFolder_ / fs::path( std::string( "emo_" ) + std::to_string( i ) + std::string( "_b" ) + std::to_string( blockID ) + "_" + std::to_string( nGridNodes_ ) + std::string( ".txt" ) ) );
	for( int i = 0; i < ofsEFAs.size(); ++i )
		ofsEFAs[i].open( BsgFolder_ / fs::path( std::string( "efa_" ) + std::to_string( i ) + std::string( "_b" ) + std::to_string( blockID ) + "_" + std::to_string( nGridNodes_ ) + std::string( ".txt" ) ) );
	
	auto writeEmoEfa = [&ofsEMOs, &ofsEFAs]( auto s ) {
		for( auto& o : ofsEMOs ) o << s;
		for( auto& o : ofsEFAs ) o << s;
	};
	
	for( size_t commId = 0; commId < blocks_[blockID].size(); ++commId ) {
		auto fragmentId = blocks_[blockID][commId];
	
		// iterate over fragment
		for( size_t i = 0; i < fragmentVertices_[fragmentId].size(); ++i ) {
			for( size_t j = 0; j < fragmentVertices_[fragmentId][i].size(); ++j ) {
				const auto& vh = fragmentVertices_[fragmentId][i][j];
				const auto& p = quadMesh.point( vh );
				if( !quadMesh.is_boundary( vh ) ) {
					writeEmoEfa( 0 );
				} else {
					// get neighbors
					OpenMesh::SmartHalfedgeHandle heh;
					for( const auto& h : quadMesh.voh_range( vh ) ) {
						if( h.is_boundary() ) {
							heh = h;
							break;
						}
					}
					const auto& pNext = quadMesh.point( heh.to() );
					const auto& pPrev = quadMesh.point( heh.prev().from() );
					const auto& ct0 = contour_.getContourType( p, pNext );
					const auto& ct1 = contour_.getContourType( p, pPrev );

					if( ct0 == Contour::ContourType::OPEN_SEA || ct1 == Contour::ContourType::OPEN_SEA ) {
						const auto& [emo, efa] = contour_.getEmoEfaVals( p, p );
						for( int i = 0; i < ofsEMOs.size(); ++i ) 
							ofsEMOs[i] << emo[i];
						for( int i = 0; i < ofsEMOs.size(); ++i ) 
							ofsEFAs[i] << efa[i];
					} else {
						writeEmoEfa( 0 );
					}
				}
				writeEmoEfa( " " );
			}
			writeEmoEfa( "\n" );
		}
		writeEmoEfa( "\n" );
	}
	
	for( auto& o : ofsEMOs ) o.close();
	for( auto& o : ofsEFAs ) o.close();
}

void BlockStructuredGrid::refine() {
	decltype( edgeVertices_ ) edgeVerticesNew(edgeVertices_.size());
	decltype( fragmentVertices_ ) fragmentVerticesNew( fragmentVertices_.size() );
	//decltype( fragmentFaces_ ) fragmentFacesNew( fragmentFaces_.size() );

	// vertices on edges
	for( const auto& e : fragmentMesh.edges() ) {
		auto& edgeVhs = edgeVertices_[e.idx()];
		auto& edgeVhsNew = edgeVerticesNew[e.idx()];

		for( int i = 0; i < edgeVhs.size() - 1; ++i ) {
			// add old and new vertex
			const auto& vh1 = edgeVhs[i];
			const auto& vh2 = edgeVhs[i + 1];
			const auto& p1 = quadMesh.point( vh1 );
			const auto& p2 = quadMesh.point( vh2 );
			const auto pMid = 0.5f * ( p1 + p2 );
			auto vhMid = quadMesh.add_vertex( pMid );
			edgeVhsNew.push_back( vh1 );
			edgeVhsNew.push_back( vhMid );
		}
		// add last old vertex
		edgeVhsNew.push_back( edgeVhs[edgeVhs.size() - 1] );
	}

	// vertices on faces
	for( const auto& f : fragmentMesh.faces() ) {
		auto& faceVhs = fragmentVertices_[f.idx()];
		auto& faceVhsNew = fragmentVerticesNew[f.idx()];

		// get edge vhs
		auto hehSouth = f.halfedge();
		while( !( ( hehSouth.from() == faceVhs[0][0] && hehSouth.to() == faceVhs[0][faceVhs[0].size() - 1] ) ||
			   ( hehSouth.from() == faceVhs[0][faceVhs[0].size() - 1] && hehSouth.to() == faceVhs[0][0] ) ) )
			hehSouth = hehSouth.next();		// fragment was turned

		auto hehEast = hehSouth.next();
		auto hehNorth = hehEast.next();
		auto hehWest = hehNorth.next();

		auto vhSW = hehSouth.from();
		auto vhSE = hehSouth.to();
		auto vhNE = hehNorth.from();
		auto vhNW = hehNorth.to();

		const auto& ehSouth = hehSouth.edge();
		const auto& ehEast = hehEast.edge();
		const auto& ehNorth = hehNorth.edge();
		const auto& ehWest = hehWest.edge();

		auto edgeVhSouth = edgeVerticesNew[ehSouth.idx()];
		auto edgeVhEast = edgeVerticesNew[ehEast.idx()];
		auto edgeVhNorth = edgeVerticesNew[ehNorth.idx()];
		auto edgeVhWest = edgeVerticesNew[ehWest.idx()];

		// orientate edges according to face
		// S
		if( edgeVhSouth[0].idx() != vhSW.idx() ) {
			std::reverse( edgeVhSouth.begin(), edgeVhSouth.end() );
		}
		// E
		if( edgeVhEast[0].idx() != vhSE.idx() && edgeVhEast[0].idx() == vhNE.idx() ) {
			std::reverse( edgeVhEast.begin(), edgeVhEast.end() );
		}
		// N
		if( edgeVhNorth[0].idx() != vhNW.idx() && edgeVhNorth[0].idx() == vhNE.idx() ) {
			std::reverse( edgeVhNorth.begin(), edgeVhNorth.end() );
		}
		// W
		if( edgeVhWest[0].idx() != vhSW.idx() && edgeVhWest[0].idx() == vhNW.idx() ) {
			std::reverse( edgeVhWest.begin(), edgeVhWest.end() );
		}

		const auto n = edgeVhSouth.size();

		faceVhsNew.resize( n );
		for( auto& fvh : faceVhsNew ) {
			fvh.resize( n );
		}

		faceVhsNew[0][0] = vhSW;
		faceVhsNew[0][faceVhsNew.size() - 1] = vhSE;
		faceVhsNew[faceVhsNew.size() - 1][faceVhsNew.size() - 1] = vhNE;
		faceVhsNew[faceVhsNew.size() - 1][0] = vhNW;
		for( auto i = 1; i < n - 1; ++i ) {
			faceVhsNew[0][i] = edgeVhSouth[i];
			faceVhsNew[i][faceVhsNew.size() - 1] = edgeVhEast[i];
			faceVhsNew[faceVhsNew.size() - 1][i] = edgeVhNorth[i];
			faceVhsNew[i][0] = edgeVhWest[i];
		}
		
		for( auto i = 1; i < faceVhs.size() - 1; ++i ) {
			for( auto j = 1; j < faceVhs.size() - 1; ++j ) {
				faceVhsNew[2 * i][2 * j] = faceVhs[i][j];
			}
		}

		// set vertices on horizontal, verical, and diagonal lines
		for( auto i = 0; i < faceVhs.size() - 1; ++i ) {
			for( auto j = 0; j < faceVhs.size() - 1; ++j ) {
				const auto& vhsw = faceVhs[i][j];
				const auto& vhse = faceVhs[i][j + 1];
				const auto& vhnw = faceVhs[i + 1][j];
				const auto& vhne = faceVhs[i + 1][j + 1];
				auto psw = quadMesh.point( vhsw );
				auto pse = quadMesh.point( vhse );
				auto pnw = quadMesh.point( vhnw );
				auto pne = quadMesh.point( vhne );
				//// south
				//if( i != 0 ) {
				//	auto ps = 0.5f * ( psw + pse );
				//	faceVhsNew[2 * i][2 * j + 1] = quadMesh.add_vertex( ps );
				//}
				//// west
				//if( j != 0 ) {
				//	auto pw = 0.5f * ( psw + pnw );
				//	faceVhsNew[2 * i + 1][2 * j] = quadMesh.add_vertex( pw );
				//}
				// north
				if( i != faceVhs.size() - 2 ) {
					auto pn = 0.5f * ( pnw + pne );
					faceVhsNew[2 * i + 2][2 * j + 1] = quadMesh.add_vertex( pn );
				}
				// east
				if( j != faceVhs.size() - 2 ) {
					auto pe = 0.5f * ( pse + pne );
					faceVhsNew[2 * i + 1][2 * j + 2] = quadMesh.add_vertex( pe );
				}
				// diagonal
				auto pMid = 0.5f * ( pnw + pse );
				faceVhsNew[2 * i + 1][2 * j + 1] = quadMesh.add_vertex( pMid );
			}
		}

	}
	edgeVertices_ = edgeVerticesNew;
	fragmentVertices_ = fragmentVerticesNew;
}

