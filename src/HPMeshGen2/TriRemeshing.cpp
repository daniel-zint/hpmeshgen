#include "TriRemeshing.h"

#include <glog/logging.h>

#include "MeshFunctions.h"
#include "Simplify.h"

inline float calcBoundaryAngle( TriMesh& mesh, const TriMesh::VertexHandle& vh ) {
	LOG_ASSERT( mesh.is_boundary( vh ) );
	TriMesh::HalfedgeHandle hehBound;
	for( auto voh : mesh.voh_range( vh ) ) {
		if( mesh.is_boundary( voh ) ) {
			hehBound = voh;
			break;
		}
	}
	LOG_ASSERT( hehBound.is_valid() );
	auto p = mesh.point( vh );
	auto b = mesh.point( mesh.to_vertex_handle( hehBound ) ) - p;
	b.normalize();
	auto a = mesh.point( mesh.from_vertex_handle( mesh.prev_halfedge_handle( hehBound ) ) ) - p;
	a.normalize();
	decltype( p ) n = { 0,0,1 };

	auto dot = a | b;
	auto det = ( a % b ) | n;
	auto angle = std::atan2( det, dot ) * 180.f / (float)M_PI;
	if( angle < 0 )
		angle = 360 + angle;

	return angle;
}

TriRemeshing::TriRemeshing( const TriMesh& mesh ) : mesh_{ mesh } {
	// compute lAvg_
	lAvg_ = 0;
	for( auto eh : mesh_.edges() ) {
		auto l = mesh_.calc_edge_length( eh );
		lAvg_ += l;
	}
	lAvg_ /= mesh_.n_edges();

	// lMin, lMax
	lMin_ = ( 3.f / 4.f ) * lAvg_;
	lMax_ = ( 4.f / 3.f ) * lAvg_;
}

TriMesh TriRemeshing::remesh( const int nFragments ) {
	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();

	auto nFaces = mesh_.n_faces();
	for( auto i = 0; i < 100; ++i ) {
		split();
		//OpenMesh::IO::write_mesh( mesh_, "output/split.off" );
		//mesh_.garbage_collection();
		collapse();
		//mesh_.garbage_collection();
		//OpenMesh::IO::write_mesh( mesh_, "output/collapse.off" );
		flip();
		//OpenMesh::IO::write_mesh( mesh_, "output/flip.off" );
		smooth();
		//OpenMesh::IO::write_mesh( mesh_, "output/smooth.off" );
		auto nFacesNew = mesh_.n_faces();
		
		int l1 = 0, l2 = 0, l3 = 0;
		for( auto eh : mesh_.edges() ) {
			auto heh = mesh_.halfedge_handle( eh, 0 );
			auto vh1 = mesh_.from_vertex_handle( heh );
			auto vh2 = mesh_.to_vertex_handle( heh );
			auto p1 = mesh_.point( vh1 );
			auto p2 = mesh_.point( vh2 );
			auto l = mesh_.calc_edge_length( eh ) / bg_->getScalar( { p1[0],p1[1] }, { p2[0],p2[1] } );
			if( l < lMin_ )
				++l1;
			else if( l > lMax_ )
				++l3;
			else
				++l2;
		}
		float r1 = (float)l1 / (float)mesh_.n_edges();
		float r2 = (float)l2 / (float)mesh_.n_edges();
		float r3 = (float)l3 / (float)mesh_.n_edges();
		std::cout << i << ":\t " << r1 << " \t|\t " << r2 << " \t|\t " << r3 << " \t|\tnVert: " << mesh_.n_vertices() << std::endl;


		if( nFacesNew != nFaces )
			nFaces = nFacesNew;
		else
			break;
	}
	if( nFragments != -1 ) {
		split( nFragments );
		lMin_ = lMax_;
		collapse( nFragments );
		flip( 0.1f );
		smooth();
	}
		
	// set points back to contour
	if( hasContour ) {
		//for( auto vh : mesh_.vertices() ) {
		//	if( !mesh_.is_boundary( vh ) )
		//		continue;
		//	auto p = mesh_.point( vh );
		//	p = contourOld_->findNearestNeighbor( p );
		//	mesh_.set_point( vh, p );
		//}
		//smooth();
		//contourOld_->findNearestNeighbors( mesh_ );
		for( const auto& vh : mesh_.vertices() ) {
			if( !vh.is_boundary() )
				continue;

			auto p = mesh_.point( vh );
			auto pNew = contour_->mapToContour( p );
			mesh_.set_point( vh, pNew );
		}

		// smooth contour
		for( int i = 0; i < 5; ++i ) {
			MeshFunctions::boundarySmoothing( mesh_, *contour_ );
			discreteMeshOptimization( mesh_, 9, 0.5f, 1 );
		}
	}


	mesh_.garbage_collection();

	return mesh_;
}

bool TriRemeshing::split() {
	bool didSplit = false;

	//OpenMesh::EPropHandleT<bool> isNeumann;
	//mesh_.get_property_handle( isNeumann, "isNeumann" );
	//OpenMesh::EPropHandleT<int> boundaryID;
	//mesh_.get_property_handle( boundaryID, "boundaryID" );
	//OpenMesh::VPropHandleT<int> contourIdx;
	//mesh_.get_property_handle( contourIdx, "contourIdx" );
	OpenMesh::VPropHandleT<float> feature_size;
	mesh_.get_property_handle( feature_size, "feature_size" );
	
	for( auto eh : mesh_.edges() ) {
		//if( mesh_.is_boundary( eh ) )
		//	continue;
		auto heh = mesh_.halfedge_handle( eh, 0 );
		auto vh1 = mesh_.from_vertex_handle( heh );
		auto vh2 = mesh_.to_vertex_handle( heh );
		auto p1 = mesh_.point( vh1 );
		auto p2 = mesh_.point( vh2 );

		float l;
		if( hasBackgroundGrid_ ) {
			
			l = mesh_.calc_edge_length( eh ) / MeshFunctions::getNeighMinSize( *bg_, mesh_, eh );
		} else {
			l = mesh_.calc_edge_length( eh );
		}

		if( l < lMax_ )
			continue;

		//auto neumann = mesh_.property( isNeumann, eh );
		//auto bId = mesh_.property( boundaryID, eh );

		auto pNew = 0.5f * ( p1 + p2 );

		auto vhNew = mesh_.split( eh, pNew );
		mesh_.property( feature_size, vhNew ) = 180.f;
		//if( mesh_.is_boundary( vhNew ) ) {
		//	mesh_.property( contourIdx, vhNew ) = mesh_.property( contourIdx, vh1 );
		//	for( auto ve : mesh_.ve_range( vhNew ) ) {
		//		if( mesh_.is_boundary( ve ) ) {
		//			mesh_.property( isNeumann, ve ) = neumann;
		//			mesh_.property( boundaryID, ve ) = bId;
		//		}
		//	}
		//}
		
		didSplit = true;
	}
	
	return didSplit;
}

bool TriRemeshing::split( const int nFaces ) {

	if( mesh_.n_faces() >= nFaces - 1 )
		return false;
	
	//OpenMesh::EPropHandleT<bool> isNeumann;
	//mesh_.get_property_handle( isNeumann, "isNeumann" );
	//OpenMesh::EPropHandleT<int> boundaryID;
	//mesh_.get_property_handle( boundaryID, "boundaryID" );
	//OpenMesh::VPropHandleT<int> contourIdx;
	//mesh_.get_property_handle( contourIdx, "contourIdx" );
	OpenMesh::VPropHandleT<float> feature_size;
	mesh_.get_property_handle( feature_size, "feature_size" );

	std::multimap<float, TriMesh::EdgeHandle, std::greater<float>> edgeMap;
	for( auto eh : mesh_.edges() ) {
		
		float l;
		if( hasBackgroundGrid_ ) {

			l = mesh_.calc_edge_length( eh ) / MeshFunctions::getNeighMinSize( *bg_, mesh_, eh );
		} else {
			l = mesh_.calc_edge_length( eh );
		}

		edgeMap.insert( std::make_pair( l, eh ) );
	}

	for( auto entry : edgeMap ) {
		auto eh = entry.second;

		auto heh = mesh_.halfedge_handle( eh, 0 );
		auto vh1 = mesh_.from_vertex_handle( heh );
		auto vh2 = mesh_.to_vertex_handle( heh );
		auto p1 = mesh_.point( vh1 );
		auto p2 = mesh_.point( vh2 );

		//auto neumann = mesh_.property( isNeumann, eh );
		//auto bId = mesh_.property( boundaryID, eh );
		
		auto pNew = 0.5f * ( p1 + p2 );
		
		auto vhNew = mesh_.split( eh, pNew );
		mesh_.property( feature_size, vhNew ) = 180.f;
		//if( mesh_.is_boundary( vhNew ) ) {
		//	mesh_.property( contourIdx, vhNew ) = mesh_.property( contourIdx, vh1 );
		//	for( auto ve : mesh_.ve_range( vhNew ) ) {
		//		if( mesh_.is_boundary( ve ) ) {
		//			mesh_.property( isNeumann, ve ) = neumann;
		//			mesh_.property( boundaryID, ve ) = bId;
		//		}
		//	}
		//}
		
		if( mesh_.n_faces() >= nFaces - 1 )
			return true;
	}

	return true;
}

bool TriRemeshing::collapse( const int nFaces ) {

	mesh_.request_vertex_status();
	OpenMesh::Decimater::OceanDecimaterT<TriMesh> decimater( mesh_ );
	OpenMesh::Decimater::ModEdgeLengthT<TriMesh>::Handle edgeLengthMod;
	decimater.add( edgeLengthMod );
	decimater.module( edgeLengthMod ).set_min_meanratio( 0.1 );
	decimater.module( edgeLengthMod ).set_max_length( lMin_ );
	decimater.module( edgeLengthMod ).set_max_feature_size( maxFeatureSize );
	if( convHullRemeshing ) {
		decimater.module( edgeLengthMod ).use_convex_hull_decimation();
	}
	if( hasBackgroundGrid_ ) {
		decimater.module( edgeLengthMod ).set_backgroundgrid( bg_ );
	}
	decimater.initialize();
	if( nFaces == -1 )
		decimater.decimate();
	else
		decimater.decimate_to_faces( 0, nFaces );
	mesh_.garbage_collection();

	return true;
}

bool TriRemeshing::flip( const float qMin ) {
	bool didFlip = false;

	for( auto eh : mesh_.edges() ) {
		if( mesh_.is_boundary( eh ) )
			continue;

		auto heh = mesh_.halfedge_handle( eh, 0 );
		auto vh1 = mesh_.from_vertex_handle( heh );
		auto vh2 = mesh_.to_vertex_handle( heh );
		auto vh3 = mesh_.to_vertex_handle( mesh_.next_halfedge_handle( heh ) );
		auto vh4 = mesh_.to_vertex_handle( mesh_.next_halfedge_handle( mesh_.opposite_halfedge_handle( heh ) ) );
		int val1 = mesh_.valence( vh1 );
		int val2 = mesh_.valence( vh2 );
		int val3 = mesh_.valence( vh3 );
		int val4 = mesh_.valence( vh4 );
		if( mesh_.is_boundary( vh1 ) ) {
			auto a = calcBoundaryAngle( mesh_, vh1 );
			int addValence = (int)(360.f - a) / 60;
			val1 += addValence;
		}
		if( mesh_.is_boundary( vh2 ) ) {
			auto a = calcBoundaryAngle( mesh_, vh2 );
			int addValence = (int)( 360.f - a ) / 60;
			val2 += addValence;
		}
		if( mesh_.is_boundary( vh3 ) ) {
			auto a = calcBoundaryAngle( mesh_, vh3 );
			int addValence = (int)( 360.f - a ) / 60;
			val3 += addValence;
		}
		if( mesh_.is_boundary( vh4 ) ) {
			auto a = calcBoundaryAngle( mesh_, vh4 );
			int addValence = (int)( 360.f - a ) / 60;
			val4 += addValence;
		}

		//auto valOld = std::abs( val1 - 6 ) + std::abs( val2 - 6 ) + std::abs( val3 - 6 ) + std::abs( val4 - 6 );
		//auto valNew = std::abs( val1 - 7 ) + std::abs( val2 - 7 ) + std::abs( val3 - 5 ) + std::abs( val4 - 5 );

		auto valOld = std::max( std::abs( val1 - 6 ), std::abs( val2 - 6 ) ) + std::max( std::abs( val3 - 6 ), std::abs( val4 - 6 ) );
		auto valNew = std::max( std::abs( val1 - 7 ), std::abs( val2 - 7 ) ) + std::max( std::abs( val3 - 5 ), std::abs( val4 - 5 ) );

		// make sure edge flip is possible
		auto p1 = mesh_.point( vh1 );
		auto p2 = mesh_.point( vh2 );
		auto p3 = mesh_.point( vh3 );
		auto p4 = mesh_.point( vh4 );
		if( MeshFunctions::isInsideTriangle( p1, p2, p3, p4 ) || MeshFunctions::isInsideTriangle( p2, p3, p4, p1 ) ||
			MeshFunctions::isInsideTriangle( p3, p4, p1, p2 ) || MeshFunctions::isInsideTriangle( p4, p1, p2, p3 ) || !mesh_.is_flip_ok( eh ) ) {
			continue;
		}
		if( MeshFunctions::calcAngle( p1, p3, p2 ) + MeshFunctions::calcAngle( p3, p2, p4 ) + MeshFunctions::calcAngle( p2, p4, p1 ) + MeshFunctions::calcAngle( p4, p1, p3 ) < 359.9999 ) {
			continue;
		}

		float q1 = QualityMetrics::meanRatioMetric( { p1,p2,p3 } );
		float q2 = QualityMetrics::meanRatioMetric( { p1,p4,p2 } );
		float qMinOld = std::min( q1, q2 );
		float q3 = QualityMetrics::meanRatioMetric( { p1,p4,p3 } );
		float q4 = QualityMetrics::meanRatioMetric( { p4,p2,p3 } );
		float qMinNew = std::min( q3, q4 );
		bool qFine = ( qMinNew > qMin || qMinNew > qMinOld );
		if( !qFine ) {
			continue;
		}

		if( valNew < valOld && mesh_.is_flip_ok( eh ) ) {
			mesh_.flip( eh );
			didFlip = true;
		}
	}
	return didFlip;
}

void TriRemeshing::smooth() {
	discreteMeshOptimization( mesh_, 9, 0.5f, 100 );
	//for( auto i = 0; i < 10; ++i ) {
	//	PointSmoothing::smoothConvexDomain( mesh_ );
	//}
}
