#include "RatRace.h"

#include "../HPMeshGen2/MeshSmoothing/QualityMetrics.h"
#include "../HPMeshGen2/MeshFunctions.h"
#include "../HPMeshGen2/blossom_quad.h"

// Calculate angle around point p2
inline float calcAngle( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& p3 ) {
	const float dot = ( p1[0] - p2[0] ) * ( p3[0] - p2[0] ) + ( p1[1] - p2[1] ) * ( p3[1] - p2[1] );
	const float det = ( p1[0] - p2[0] ) * ( p3[1] - p2[1] ) - ( p1[1] - p2[1] ) * ( p3[0] - p2[0] );
	return ( 180.f / (float)M_PI ) * std::fabs( atan2( det, dot ) );
}

inline float calcBoundaryAngle( PolyMesh& mesh, const PolyMesh::VertexHandle& vh ) {
	LOG_ASSERT( mesh.is_boundary( vh ) );
	PolyMesh::HalfedgeHandle hehBound;
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
	decltype(p) n = { 0,0,1 };

	auto dot = a | b;
	auto det = ( a % b ) | n;
	auto angle = std::atan2( det, dot ) * 180.f / (float)M_PI;
	if( angle < 0 )
		angle = 360.f + angle;

	return angle;
}

inline void tri2hybrid( PolyMesh& pmesh ) {
	
	// map mit weight & edge
	std::map<float, PolyMesh::EdgeHandle> prioQueue;
	std::vector<std::vector<PolyMesh::VertexHandle>> quadVertices( pmesh.n_edges() );
	// iteriere ueber edges
	// fasse beide dreiecke zu viereck zusammen und berechne min/max angle criterion
	// --> schreibe weight in map

	// store boundary properties
	for ( auto e_it = pmesh.edges_begin(); e_it != pmesh.edges_end(); ++e_it ) {
		if ( !pmesh.is_boundary( *e_it ) )
			continue;
		TriMesh::HalfedgeHandle heh = pmesh.halfedge_handle( *e_it, 0 );
		if ( !pmesh.is_boundary( heh ) )
			heh = pmesh.opposite_halfedge_handle( heh );

		TriMesh::VertexHandle vh = pmesh.from_vertex_handle( heh );
	}

	for ( uint edgeID = 0; edgeID < pmesh.n_edges(); ++edgeID ) {
		PolyMesh::EdgeHandle eh = pmesh.edge_handle( edgeID );

		if ( pmesh.is_boundary( eh ) ) {
			continue;
		}

		PolyMesh::HalfedgeHandle heh1 = pmesh.halfedge_handle( eh, 0 );
		PolyMesh::HalfedgeHandle heh2 = pmesh.halfedge_handle( eh, 1 );

		std::vector<PolyMesh::VertexHandle>vVec( 4 );
		std::vector<PolyMesh::Point>pVec( 4 );

		/*
		1-----0
		|   / |
		| /   |
		2-----3
		*/

		vVec[0] = pmesh.to_vertex_handle( heh1 );
		vVec[1] = pmesh.to_vertex_handle( pmesh.next_halfedge_handle( heh1 ) );
		vVec[2] = pmesh.to_vertex_handle( heh2 );
		vVec[3] = pmesh.to_vertex_handle( pmesh.next_halfedge_handle( heh2 ) );

		pVec[0] = pmesh.point( vVec[0] );
		pVec[1] = pmesh.point( vVec[1] );
		pVec[2] = pmesh.point( vVec[2] );
		pVec[3] = pmesh.point( vVec[3] );

		// do not delete edge if quad is non-convex
		if ( calcAngle( pVec[3], pVec[0], pVec[2] ) + calcAngle( pVec[2], pVec[0], pVec[1] ) > 180. ) {
			continue;
		}
		if ( calcAngle( pVec[1], pVec[2], pVec[0] ) + calcAngle( pVec[0], pVec[2], pVec[3] ) > 180. ) {
			continue;
		}

		std::vector<float>angleVec( 4 );
		angleVec[0] = (float)calcAngle( pVec[3], pVec[0], pVec[1] );
		angleVec[1] = (float)calcAngle( pVec[0], pVec[1], pVec[2] );
		angleVec[2] = (float)calcAngle( pVec[1], pVec[2], pVec[3] );
		angleVec[3] = (float)calcAngle( pVec[2], pVec[3], pVec[0] );

		float minAngle = std::numeric_limits<float>::max();
		float maxAngle = -std::numeric_limits<float>::max();
		for ( size_t i = 0; i < angleVec.size(); ++i ) {
			minAngle = std::min( minAngle, angleVec[i] );
			maxAngle = std::max( maxAngle, angleVec[i] );
		}

		// do not add edge if its deletion would generate a degenerate quad
		bool degenerateQuad = false;
		for ( size_t i = 0; i < angleVec.size(); ++i ) {
			if ( angleVec[i] > 165 ) {
				degenerateQuad = true;
				break;
			}
		}
		if ( degenerateQuad )
			continue;

		float weight = maxAngle - minAngle;

		prioQueue[weight] = eh;
		quadVertices[edgeID] = vVec;
	}


	// delete edges (tracke welche faces zusammengefasst wurden)

	std::vector<bool>faceDeleted( pmesh.n_faces(), false );
	//std::cout << "Number of faces: " << mesh.n_faces() << std::endl;

	pmesh.request_face_status();
	pmesh.request_edge_status();
	pmesh.request_vertex_status();

	std::vector<PolyMesh::EdgeHandle> edgesToDelete;

	for ( auto pq_it = prioQueue.begin(); pq_it != prioQueue.end(); ++pq_it ) {
		PolyMesh::EdgeHandle eh = pq_it->second;
		//uint edgeID = eh.idx();
		PolyMesh::FaceHandle fh1 = pmesh.face_handle( pmesh.halfedge_handle( eh, 0 ) );
		PolyMesh::FaceHandle fh2 = pmesh.face_handle( pmesh.halfedge_handle( eh, 1 ) );

		//std::cout << "fh1 = " << fh1.idx() << "  |  fh2 = " << fh2.idx() << std::endl;

		// check if one of the faces is already part of a quad
		if ( faceDeleted[fh1.idx()] || faceDeleted[fh2.idx()] )
			continue;

		faceDeleted[fh1.idx()] = true;
		faceDeleted[fh2.idx()] = true;

		edgesToDelete.push_back( eh );
	}

	for ( size_t i = 0; i < edgesToDelete.size(); ++i ) {
		pmesh.delete_edge( edgesToDelete[i], false );
		pmesh.add_face( quadVertices[edgesToDelete[i].idx()] );
	}

	pmesh.garbage_collection();

}

RatRace::RatRace( TriMesh mesh, bool useBlossomQuad ) : tmesh_{ mesh } {
	
	LOG_ASSERT( tmesh_.n_vertices() > 0 );
	LOG_ASSERT( tmesh_.n_faces() % 2 == 0 ) << "perfect matching can only be done for an even number of faces";

	// get initial hybrid mesh
	hmesh_ = (PolyMesh)tmesh_;
	if( useBlossomQuad ) {
		LOG( INFO ) << "Use blossom quad for hybrid mesh generation";
		blossom_quad bloss;
		hmesh_ = bloss.do_blossom_algo( hmesh_ );
	} else {
		LOG( INFO ) << "Use fast hybrid mesh generation";
		tri2hybrid( hmesh_ );
	}
	
	
}

void RatRace::run( bool doPostprocessing ) {
	DualHybridGraph dhg( hmesh_ );

	auto nTriangles = dhg.nTriangles();
	auto nTrianglesOld = nTriangles;

	LOG_ASSERT( nTriangles % 2 == 0 );

	while( nTriangles > 0 ) {
		dhg.connectTriangles();
		dhg.garbageCollection();

		nTriangles = dhg.nTriangles();
		if( nTriangles == nTrianglesOld ) {
			LOG( ERROR ) << "Cannot reduce number of triangles. nTriangles = " << nTriangles;
		} else {
			nTrianglesOld = nTriangles;
		}
	}

	hmesh_ = dhg.getMesh();
	
	if( doPostprocessing ) {
		swapValence2Vertices();
		untangle();
		optimize();
	}
}

void RatRace::untangle( const int maxIter ) {
	auto[q, qMean] = QualityMetrics::minCondition( hmesh_ );
	auto qOld = q;
	auto qMeanOld = qMean;

	auto i = 0;
	for(; q < 0 && i < maxIter / 100; ++i) {
		discreteMeshOptimization( hmesh_, Q_CONDITION, 0.5f, 100 );
		std::tie( q, qMean ) = QualityMetrics::minCondition( hmesh_ );
		if( q == qOld && qMean == qMeanOld ) {
			LOG( INFO ) << "Cannot further improve quality. Min = " << q << ". Mean = " << qMean;
			break;
		} else {
			qOld = q;
			qMeanOld = qMean;
		}
	}

	if( q < 0 && i == maxIter / 10 ) {
		LOG( INFO ) << "Reached maximum number of iterations. Optimization was stopped even though it was still converging. More iterations would give better results.";
	}
}

void RatRace::optimize( const int maxIter, const int iterPerSwap ) {
	auto [q, qMean] = QualityMetrics::minCondition( hmesh_ );
	auto qOld = q;
	auto qMeanOld = qMean;

	auto i = 0;
	for(; q < 0 && i < maxIter; ++i ) {
		discreteMeshOptimization( hmesh_, Q_CONDITION, 0.5f, iterPerSwap );
		std::tie(q, qMean) = QualityMetrics::minCondition( hmesh_ );
		if( q == qOld && qMean == qMeanOld ) {
			LOG( INFO ) << "Cannot further improve quality. Min = " << q << ". Mean = " << qMean;
			break;
		} else {
			qOld = q;
			qMeanOld = qMean;
		}
	}

	if( q < 0 && i == maxIter ) {
		LOG( INFO ) << "Reached maximum number of iterations. Optimization was stopped even though it was still converging. More iterations would give better results.";
	}
}

void RatRace::swapValence2Vertices() {
	auto valence = [ this ]( PolyMesh::VertexHandle vh ) { return std::distance( hmesh_.vv_begin( vh ), hmesh_.vv_end( vh ) ); };

	hmesh_.request_edge_status();
	hmesh_.request_face_status();

	bool flippedSomething = false;
	do {
		flippedSomething = false;

		for( auto eh : hmesh_.edges() ) {
			if( hmesh_.is_boundary( eh ) ) continue;
			auto heh = hmesh_.halfedge_handle( eh, 0 );
			if( valence( hmesh_.from_vertex_handle( heh ) ) == 2 ) continue;
			if( valence( hmesh_.to_vertex_handle( heh ) ) == 2 ) continue;


			//	1 --- 0 --- 5
			//	|     ^     |
			//	|     |     |
			//	2 --- 3 --- 4

			auto hehInit = heh;
			std::vector<PolyMesh::VertexHandle> vhs;
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );
			heh = hmesh_.opposite_halfedge_handle( hehInit );
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );

			std::vector<int> vals;
			vals.reserve( vhs.size() );
			for( auto vh : vhs ) {
				vals.push_back( static_cast<size_t>( valence( vh ) ));
			}

			if( vals[0] == 3 || vals[3] == 3 ) continue;

			if( vals[2] == 2 || vals[5] == 2 ) {
				// swap to 2,5
				MeshFunctions::edgeSwapPrev( hmesh_, eh );
				flippedSomething = true;
				break;
			} else if( vals[1] == 2 || vals[4] == 2 ) {
				// swap to 1,4
				MeshFunctions::edgeSwapNext( hmesh_, eh );
				flippedSomething = true;
				break;
			}

		}
	} while( flippedSomething );

	hmesh_.garbage_collection();
}

void RatRace::swapValence2VerticesBFS() {
	hmesh_.request_edge_status();
	hmesh_.request_face_status();

	std::set<PolyMesh::VertexHandle> vhReworked;

	while( true ) {
		// find bad element
		PolyMesh::FaceHandle fhBegin;
		PolyMesh::VertexHandle vhBegin;
		for( auto f : hmesh_.faces() ) {
			std::vector<OpenMesh::Vec3f> points;
			points.reserve( 4 );
			for( auto fv : hmesh_.fv_range( f ) ) {
				points.push_back( hmesh_.point( fv ) );
			}
			auto q = QualityMetrics::condition( points );

			if( q < 0.1 ) {
				for( auto fheh : hmesh_.fh_range( f ) ) {
					auto fv = hmesh_.from_vertex_handle( fheh );
					auto p = hmesh_.point( fv );
					auto pNext = hmesh_.point( hmesh_.to_vertex_handle( fheh ) );
					auto pPrev = hmesh_.point( hmesh_.from_vertex_handle( hmesh_.prev_halfedge_handle( fheh ) ) );
					auto qv = QualityMetrics::conditionSurfVertex( { pPrev, p, pNext }, { 0,0,1 } );

					if( hmesh_.valence( fv ) == 2 && qv < 0.1 ) {
						if( !hmesh_.is_boundary( fv ) )
							continue;
						if( vhReworked.find( fv ) != vhReworked.end() ) {
							LOG( INFO ) << "Could not resolve valence 2 for vertex: " << fv.idx();
							continue;
						}
						vhBegin = fv;
						fhBegin = f;
						break;
					}
				}
				if( fhBegin.is_valid() )
					break;
			}
		}

		if( !fhBegin.is_valid() )
			break;

		if( vhReworked.find( vhBegin ) != vhReworked.end() ) {
			LOG( INFO ) << "Could not resolve all valence 2 vertices";
			break;
		}
		vhReworked.insert( vhBegin );

		bool notStuckYet;
		//int iterCount = 0;
		do {
			notStuckYet = false;

			hmesh_.request_vertex_status();
			hmesh_.request_edge_status();
			hmesh_.request_face_status();

			std::vector<int> visited( hmesh_.n_vertices(), false );
			std::vector<PolyMesh::HalfedgeHandle> hehRoot( hmesh_.n_vertices() );
			std::queue<PolyMesh::VertexHandle> queue;
			visited[vhBegin.idx()] = true;
			for( auto voh : hmesh_.voh_range( vhBegin ) ) {
				auto vh = hmesh_.to_vertex_handle( voh );
				queue.push( vh );
				visited[vh.idx()] = true;
				hehRoot[vh.idx()] = voh;
			}

			while( !queue.empty() ) {
				auto vh = queue.front();
				queue.pop();

				auto val = hmesh_.valence( vh );

				if( val == 2 )
					continue;


				// get candidates
				auto heh = hehRoot[vh.idx()];

				if( !hmesh_.is_boundary( heh ) ) {
					auto hehCand = hmesh_.next_halfedge_handle( heh );
					auto vhNext = hmesh_.to_vertex_handle( hehCand );

					if( !visited[vhNext.idx()] ) {
						auto valNext = hmesh_.valence( vhNext );

						if( val > 3 && valNext > 3 && !hmesh_.is_boundary( hmesh_.edge_handle( hehCand ) ) ) {
							// found
							MeshFunctions::edgeSwapPrev( hmesh_, hmesh_.edge_handle( hehCand ) );
							notStuckYet = true;
							break;
						}

						queue.push( vhNext );
						visited[vhNext.idx()] = true;
						hehRoot[vhNext.idx()] = hehCand;
					}
				}
				auto hehOpp = hmesh_.opposite_halfedge_handle( heh );
				if( !hmesh_.is_boundary( hehOpp ) ) {
					auto hehCand = hmesh_.prev_halfedge_handle( hehOpp );
					auto vhPrev = hmesh_.from_vertex_handle( hehCand );

					if( !visited[vhPrev.idx()] ) {
						auto valPrev = hmesh_.valence( vhPrev );

						if( val > 3 && valPrev > 3 && !hmesh_.is_boundary( hmesh_.edge_handle( hehCand ) ) ) {
							// found
							MeshFunctions::edgeSwapNext( hmesh_, hmesh_.edge_handle( hehCand ) );
							notStuckYet = true;
							break;
						}

						queue.push( vhPrev );
						visited[vhPrev.idx()] = true;
						hehRoot[vhPrev.idx()] = hmesh_.opposite_halfedge_handle( hehCand );
					}
				}

			}

			hmesh_.garbage_collection();
			//OpenMesh::IO::write_mesh( hmesh_, std::string("C:/Users/Daniel/Documents/RatRace/output/edgeSwapDMO_vh") + std::to_string(vhBegin.idx()) + "_i" + std::to_string(iterCount++) + ".off" );
		} while( hmesh_.valence( vhBegin ) == 2 && notStuckYet );
	}
}

void RatRace::swapValence() {
	hmesh_.request_edge_status();
	hmesh_.request_face_status();

	bool flippedSomething = false;
	do {
		flippedSomething = false;

		for( auto eh : hmesh_.edges() ) {
			if( hmesh_.is_boundary( eh ) ) continue;
			auto heh = hmesh_.halfedge_handle( eh, 0 );
			if( hmesh_.valence( hmesh_.from_vertex_handle( heh ) ) == 2 ) continue;
			if( hmesh_.valence( hmesh_.to_vertex_handle( heh ) ) == 2 ) continue;


			//	1 --- 0 --- 5
			//	|     ^     |
			//	|     |     |
			//	2 --- 3 --- 4

			auto hehInit = heh;
			std::vector<PolyMesh::VertexHandle> vhs;
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );	// 0
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );	// 1
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );	// 2
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );	// 3
			heh = hmesh_.opposite_halfedge_handle( hehInit );
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );	// 4
			heh = hmesh_.next_halfedge_handle( heh );
			vhs.push_back( hmesh_.to_vertex_handle( heh ) );	// 5
			heh = hmesh_.next_halfedge_handle( heh );

			std::vector<int> vals;
			vals.reserve( vhs.size() );
			for( auto vh : vhs ) {
				auto val = hmesh_.valence( vh );
				if( hmesh_.is_boundary( vh ) ) {
					++val;
				}
				vals.push_back( val );
			}

			if( hmesh_.is_boundary( vhs[0] ) ) {
				auto angle = calcBoundaryAngle( hmesh_, vhs[0] );
				if( angle > 300 && hmesh_.valence( vhs[0] ) <= 4 )
					continue;
			}
			if( hmesh_.is_boundary( vhs[3] ) ) {
				auto angle = calcBoundaryAngle( hmesh_, vhs[3] );
				if( angle > 300 && hmesh_.valence( vhs[3] ) <= 4 )
					continue;
			}

			if( vals[0] == 3 + hmesh_.is_boundary( vhs[0] ) || vals[3] == 3 + hmesh_.is_boundary( vhs[3] ) ) continue;

			auto valOld = vals[0] + vals[3];
			auto valRight = vals[5] + vals[2];
			auto valLeft = vals[1] + vals[4];

			if( valOld - valRight >= valOld - valLeft && valOld - valRight >= 3 ) {
				// swap to 2,5
				MeshFunctions::edgeSwapPrev( hmesh_, eh );
				flippedSomething = true;
				//break;
			} else if( valOld - valLeft > valOld - valRight && valOld - valLeft >= 3 ) {
				// swap to 1,4
				MeshFunctions::edgeSwapNext( hmesh_, eh );
				flippedSomething = true;
				//break;
			}

		}
	} while( flippedSomething );

	hmesh_.garbage_collection();
}

void RatRace::valence2SplitCollapse() {
	//int nSplits = 0;

	std::set<PolyMesh::VertexHandle> vhReworked;

	hmesh_.request_face_status();
	hmesh_.request_edge_status();
	hmesh_.request_vertex_status();

	while( true ) {
		// find bad vertex
		PolyMesh::VertexHandle vhBegin;
		for( auto f : hmesh_.faces() ) {
			std::vector<OpenMesh::Vec3f> points;
			points.reserve( 4 );
			for( auto fv : hmesh_.fv_range( f ) ) {
				points.push_back( hmesh_.point( fv ) );

			}
			auto q = QualityMetrics::condition( points );

			if( q < 0.1 ) {
				for( auto fheh : hmesh_.fh_range( f ) ) {
					auto fv = hmesh_.from_vertex_handle( fheh );

					auto p = hmesh_.point( fv );
					auto pNext = hmesh_.point( hmesh_.to_vertex_handle( fheh ) );
					auto pPrev = hmesh_.point( hmesh_.from_vertex_handle( hmesh_.prev_halfedge_handle( fheh ) ) );
					auto qv = QualityMetrics::conditionSurfVertex( { pPrev, p, pNext }, { 0,0,1 } );

					if( hmesh_.valence( fv ) == 2 && qv < 0.1 ) {
						if( !hmesh_.is_boundary( fv ) )
							continue;
						if( vhReworked.find( fv ) != vhReworked.end() ) {
							continue;
						}
						vhBegin = fv;
						break;
					}
				}
				if( vhBegin.is_valid() )
					break;
			}
		}
		if( !vhBegin.is_valid() ) {
			for( auto vh : hmesh_.vertices() ) {
				if( hmesh_.valence( vh ) == 2 && !hmesh_.is_boundary( vh ) ) {
					if( vhReworked.find( vh ) != vhReworked.end() ) {
						continue;
					}
					vhBegin = vh;
					break;
				}
			}
		}

		if( !vhBegin.is_valid() )
			break;

		vhReworked.insert( vhBegin );

		// illegal splits: boundary, create val2
		// preferable splits: keep valence at other vertex low

		auto isLegalSplit = [this]( PolyMesh::HalfedgeHandle heh ) {
			if( hmesh_.is_boundary( hmesh_.edge_handle( heh ) ) ) {
				return false;
			}

			auto vh = hmesh_.from_vertex_handle( heh );
			auto val = hmesh_.valence( vh );

			if( val >= 4 )
				return true;

			if( hmesh_.is_boundary( vh ) ) {

				float angle = calcBoundaryAngle( hmesh_, vh );

				if( angle < 130 ) {
					return true;
				} else {
					return false;
				}
			} else {
				return false;
			}

			return false;
		};

		// find all possible splits
		std::vector<PolyMesh::HalfedgeHandle> hehCand;
		for( auto voh : hmesh_.voh_range( vhBegin ) ) {
			if( hmesh_.is_boundary( voh ) )
				continue;
			auto next = hmesh_.next_halfedge_handle( voh );
			//if( hmesh_.is_boundary( hmesh_.edge_handle( next ) ) )
			//	continue;
			//if( hmesh_.valence( hmesh_.from_vertex_handle( next ) ) <= 3 ) 
			//	continue;
			if( !isLegalSplit( next ) ) {
				continue;
			}
			hehCand.push_back( next );
		}
		for( auto vih : hmesh_.vih_range( vhBegin ) ) {
			if( hmesh_.is_boundary( vih ) )
				continue;
			auto next = hmesh_.prev_halfedge_handle( vih );
			if( hmesh_.is_boundary( hmesh_.edge_handle( next ) ) )
				continue;
			if( hmesh_.valence( hmesh_.to_vertex_handle( next ) ) <= 3 )
				continue;
			hehCand.push_back( hmesh_.opposite_halfedge_handle( next ) );
		}

		if( hehCand.size() == 0 ) {
			//LOG( INFO ) << "No split candidates found. Cannot resolve valence 2 vertices. vh = " << vhBegin.idx() << "  |  p = " << hmesh_.point(vhBegin);
			vhReworked.insert( vhBegin );
			continue;
		}

		std::map<int, PolyMesh::HalfedgeHandle> prioSplit;
		for( auto heh : hehCand ) {
			auto vh1 = hmesh_.to_vertex_handle( hmesh_.next_halfedge_handle( hmesh_.opposite_halfedge_handle( heh ) ) );
			auto vh2 = hmesh_.from_vertex_handle( hmesh_.prev_halfedge_handle( heh ) );
			auto val = hmesh_.valence( vh1 ) + hmesh_.valence( vh2 );
			prioSplit[val] = heh;
		}

		auto hehSplit = prioSplit.begin()->second;
		MeshFunctions::vertexSplit( hmesh_, hehSplit );

		// find possible collapses in one ring
		// if no collapse is found, increase radius

		auto valCollapse = [this]( PolyMesh::HalfedgeHandle heh ) {
			if( hmesh_.is_boundary( hmesh_.edge_handle( heh ) ) )
				return -1;
			auto vh = hmesh_.from_vertex_handle( heh );
			if( hmesh_.is_boundary( vh ) )
				return -1;
			auto vhNext = hmesh_.to_vertex_handle( heh );
			auto vhPrev = hmesh_.from_vertex_handle( hmesh_.prev_halfedge_handle( heh ) );
			if( hmesh_.valence( vhNext ) <= 3 || hmesh_.valence( vhPrev ) <= 3 )
				return -1;
			auto vhOpp = hmesh_.to_vertex_handle( hmesh_.next_halfedge_handle( heh ) );
			return (int)hmesh_.valence( vh ) + (int)hmesh_.valence( vhOpp );
		};

		std::set<PolyMesh::VertexHandle> vhCand;
		vhCand.insert( vhBegin );
		std::map<int, PolyMesh::HalfedgeHandle> prioCollapse;
		while( prioCollapse.size() == 0 ) {
			for( auto vh : vhCand ) {
				for( auto vf : hmesh_.vf_range( vh ) ) {
					for( auto h : hmesh_.fh_range( vf ) ) {
						auto val = valCollapse( h );
						if( val >= 0 ) {
							prioCollapse[val] = h;
						}
					}
				}
			}

			decltype( vhCand ) vhCandNew;
			for( auto vh : vhCand ) {
				vhCandNew.insert( vh );
				for( auto vv : hmesh_.vv_range( vh ) ) {
					vhCandNew.insert( vv );
				}
			}
			vhCand = vhCandNew;
		}

		auto hehCollapse = prioCollapse.begin()->second;
		MeshFunctions::diagonalCollapse( hmesh_, hehCollapse );
	}

	hmesh_.garbage_collection();
}