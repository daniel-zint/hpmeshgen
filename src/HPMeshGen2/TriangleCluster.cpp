#include "TriangleCluster.h"

#include <glog/logging.h>

#include <set>

#include "HelperFunctions.h"

auto trianglePoints( TriMesh& mesh, const TriMesh::FaceHandle fh ) {
	std::vector<TriMesh::Point> points;
	for ( auto fvh : mesh.fv_range( fh ) ) {
		points.push_back( mesh.point( fvh ) );
	}

	return std::make_tuple( points[0], points[1], points[2] );
	
}


TriangleCluster::TriangleCluster( TriMesh & fragMesh, TriMesh & fineMesh ) {

	LOG_ASSERT( fragMesh.n_vertices() < fineMesh.n_vertices() ) << "fineMesh and fragMesh might be swapped";

	OpenMesh::VPropHandleT<int> vPropFragAffiliation;
	fineMesh.add_property( vPropFragAffiliation, "vPropFragAffiliation" );

	auto getAffil = [ &fineMesh, &vPropFragAffiliation ]( TriMesh::VertexHandle vh ) -> int& { return fineMesh.property( vPropFragAffiliation, vh ); };

	for ( auto vh : fineMesh.vertices() ) {
		getAffil( vh ) = -1;
	}

	// iterate over fine vertices and find surrounding fragment
	for ( auto vh : fineMesh.vertices() ) {
		for ( auto fragFh : fragMesh.faces() ) {
			auto[p1, p2, p3] = trianglePoints( fragMesh, fragFh );
			if ( HelperFunctions::isInsideTriangle( p1, p2, p3, fineMesh.point( vh ) ) ) {
				getAffil( vh ) = fragFh.idx();
				break;
			}
		}
	}

	// fill vertices which have only neighbors of one fragment
	bool somethingChanged = false;
	do {
		somethingChanged = false;
		for ( auto vh : fineMesh.vertices() ) {
			if ( getAffil( vh ) != -1 )
				continue;

			std::set<int> neighFragAffils;
			for ( auto vv : fineMesh.vv_range( vh ) ) {
				if ( getAffil( vv ) != -1 ) {
					neighFragAffils.insert( getAffil( vv ) );
				}
			}
			if ( neighFragAffils.size() != 1 ) continue;
			getAffil( vh ) = *( neighFragAffils.begin() );
			somethingChanged = true;
		}
	} while ( somethingChanged );

	// fill vertices with several neighbor fragments
	for ( auto vh : fineMesh.vertices() ) {
		if ( getAffil( vh ) != -1 )
			continue;

		std::set<int> neighFragAffils;
		for ( auto vv : fineMesh.vv_range( vh ) ) {
			if ( getAffil( vv ) != -1 ) {
				neighFragAffils.insert( getAffil( vv ) );
			}
		}
		LOG_IF( ERROR, neighFragAffils.size() < 2 ) << "Vertex " << vh.idx() << " has 0 or 1 neighboring fragments.";
		
		auto p = fineMesh.point( vh );
		std::vector<std::tuple<int, std::array<float, 3>>> fragBarys;
		for ( auto affil : neighFragAffils ) {
			// calculate barycentric coordinates for vh
			auto[p1, p2, p3] = trianglePoints( fragMesh, fragMesh.face_handle( affil ) );
			auto[b1, b2, b3] = HelperFunctions::barycentricCoordinates( { p1,p2,p3 }, p );
			std::array<float, 3> barr{ b1,b2,b3 };
			std::sort( barr.begin(), barr.end() );
			fragBarys.push_back( std::make_tuple( affil, barr ) );
		}

		fragBarys.erase( std::remove_if(
			fragBarys.begin(), fragBarys.end(),
			[]( decltype( fragBarys[0] ) a ) {
				return std::get<1>( a )[1] < 0;
			} ), fragBarys.end() );
		
		std::sort( fragBarys.begin(), fragBarys.end(), []( decltype( fragBarys[0] ) a, decltype( fragBarys[0] ) b ) { return std::get<1>(a)[0] > std::get<1>( b )[0]; } );
		
		getAffil( vh ) = std::get<0>( fragBarys[0] );
	}


	// print results

	fineMesh.request_vertex_colors();
	fineMesh.request_face_colors();
	for ( auto vh : fineMesh.vertices() ) {
		auto affil = fineMesh.property( vPropFragAffiliation, vh );
		if ( affil == -1 ) {
			fineMesh.set_color( vh, { 255,255,255 } );
		} else {
			OpenMesh::Vec3uc color = { (uint)( 255 * affil / (float)( fragMesh.n_faces() - 1 ) ), 0, 0 };
			fineMesh.set_color( vh, color );
		}
	}
	for ( auto fh : fineMesh.faces() ) {
		std::set<int> affils;
		for ( auto fv : fineMesh.fv_range( fh ) ) {
			if ( getAffil( fv ) != -1 )
				affils.insert( getAffil( fv ) );
		}
		
		if ( affils.size() == 0 )
			fineMesh.set_color( fh, { 255,255,255 } );
		else if ( affils.size() == 1 ) {
			OpenMesh::Vec3uc color = { (uint)( 255 * *(affils.begin()) / (float)( fragMesh.n_faces() - 1 ) ), 0, 0 };
			fineMesh.set_color( fh, color );
		} else {
			fineMesh.set_color( fh, { 0,255,0 } );
		}

	}
	
	////for ( auto i = 0; i < vertexAffiliation.size(); ++i ) {
	////	auto vh = fineMesh.vertex_handle( i );
	////
	////	if ( vertexAffiliation[i] == -1 ) continue;
	////
	////	OpenMesh::Vec3uc color = { (uint)( 255 * vertexAffiliation[i] / (float)( fragmentElements.size() - 1 ) ), 0, 0 };
	////	fineMesh.set_color( vh, color );
	////
	////}
	//	
	//for ( auto fragId = 0; fragId < fragmentElements.size(); ++fragId ) {
	//	OpenMesh::Vec3uc color = { (uint)( 255 * fragId / (float)( fragmentElements.size() - 1 ) ), 0, 0 };
	//
	//	for ( auto fh : fragmentElements[fragId] ) {
	//		if ( fineMesh.color( fh )[0] != 255 ) {
	//			fineMesh.set_color( fh, { 0,255,0 } );
	//		} else {
	//			fineMesh.set_color( fh, color );
	//		}
	//	}
	//	//for ( auto vh : fragmentVertices[fragId] ) {
	//	//	fineMesh.set_color( vh, color );
	//	//}
	//}
	OpenMesh::IO::write_mesh( fineMesh, "C:/Users/Daniel/Documents/HPMeshGen2/output/coloredFrag.off", OpenMesh::IO::Options::FaceColor | OpenMesh::IO::Options::VertexColor );

	

	//OpenMesh::IO::write_mesh( mesh, "resources/meshes/coloredOutput.ply", OpenMesh::IO::Options::VertexColor );
}
