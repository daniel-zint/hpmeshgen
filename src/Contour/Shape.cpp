#include "Shape.h"

#include <vector>
#include <queue>

namespace Contour
{
	Shape::Shape( OceanMesh& mesh ) {
		OpenMesh::VPropHandleT<int> propContourIdx;
		mesh.add_property( propContourIdx );

		// add points
		for( auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it ) {
			if( !mesh.is_boundary( *v_it ) ) {
				mesh.property( propContourIdx, *v_it ) = -1;
				continue;
			}

			const int idx = static_cast<int>(points_.size());
			TriMesh::Point p = mesh.point( *v_it );
			p[2] = 0;

			auto& point = points_.add();
			point.p0() = p;

			mesh.property( propContourIdx, *v_it ) = point.idx();
		}

		// add left and right index	
		for( const auto& heh : mesh.halfedges() ) {
			if( !heh.is_boundary() )
				continue;

			const auto& idxLeft = mesh.property( propContourIdx, heh.from() );
			const auto& idxRight = mesh.property( propContourIdx, heh.to() );

			auto& edge = edges_.add();
			edge.from() = idxLeft;
			edge.to() = idxRight;

			edge.contourType() = mesh.property( mesh.propContourType, heh.edge() );
			if( edge.contourType() == ContourType::OPEN_SEA ) {
				emo_[edge.idx()] = mesh.property( mesh.propEmo, heh.edge() );
				efa_[edge.idx()] = mesh.property( mesh.propEfa, heh.edge() );
			}
		}

		mesh.remove_property( propContourIdx );

		// set next and prev edges
		std::vector<std::array<PointIdxT, 2>> prevnext( points_.size() );
		for( const auto& e : edges_ ) {
			prevnext[e.from()][1] = e.idx();
			prevnext[e.to()][0] = e.idx();

			points_[e.from()].edgeIdx() = e.idx();
		}
		for( const auto& pn : prevnext ) {
			edges_[pn[0]].next() = pn[1];
			edges_[pn[1]].prev() = pn[0];
		}

		// add local surface estimation
		for( auto& contourPoint : points_ ) {
			// local coords
			const auto& edgeNextIdx = contourPoint.edgeIdx();
			const auto& pNext = points_[edges_[edgeNextIdx].to()].p0();		// new pLeft
			const auto& edgePrevIdx = edges_[edgeNextIdx].prev();
			const auto& pPrev = points_[edges_[edgePrevIdx].from()].p0();	// new pRight

			auto t = pNext - pPrev;
			decltype( t ) n{ t[1], -t[0], t[2] };
			t = t.normalized();
			n = n.normalized();

			contourPoint.normal() = n;
			contourPoint.tangential() = t;

			// polyParam
			float ATA = 0;
			float ATb = 0;
			std::vector<TriMesh::Point> points = { pNext, pPrev };
			for( int i = 0; i < points.size(); ++i ) {
				TriMesh::Point loc = contourPoint.map_glob2loc( points[i] );
				float u = loc[0];
				float v = loc[1];
				ATb += u * u * v;
				ATA += u * u * u * u;
			}
			contourPoint.surfParam() = ATb / ATA;
		}

		// validity check
		for( const auto& e : edges_ ) {
			LOG_ASSERT( e.idx() != INVALID_EDGE );
			LOG_ASSERT( e.next() != INVALID_EDGE );
			LOG_ASSERT( e.prev() != INVALID_EDGE );
			LOG_ASSERT( e.from() != INVALID_POINT );
			LOG_ASSERT( e.to() != INVALID_POINT );
			LOG_ASSERT( e.contourType() != ContourType::NOTHING );
		}
		for( const auto& p : points_ ) {
			LOG_ASSERT( p.idx() != INVALID_POINT );
			LOG_ASSERT( p.edgeIdx() != INVALID_EDGE );
		}

		// init kd tree
		kdTree_ = std::make_unique<kd_tree_t>( 3, *this );
		kdTree_->buildIndex();
    }

	Shape::Shape( TriMesh& mesh ) {
		OpenMesh::VPropHandleT<int> propContourIdx;
		mesh.add_property( propContourIdx );

		// add points
		for( auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it ) {
			if( !mesh.is_boundary( *v_it ) ) {
				mesh.property( propContourIdx, *v_it ) = -1;
				continue;
			}

			const int idx = points_.size();
			TriMesh::Point p = mesh.point( *v_it );
			p[2] = 0;

			auto& point = points_.add();
			point.p0() = p;

			mesh.property( propContourIdx, *v_it ) = point.idx();
		}

		// add left and right index	
		for( const auto& heh : mesh.halfedges() ) {
			if( !heh.is_boundary() )
				continue;

			const auto& idxLeft = mesh.property( propContourIdx, heh.from() );
			const auto& idxRight = mesh.property( propContourIdx, heh.to() );

			auto& edge = edges_.add();
			edge.from() = idxLeft;
			edge.to() = idxRight;

			edge.contourType() = ContourType::NOTHING;
		}

		mesh.remove_property( propContourIdx );

		// set next and prev edges
		std::vector<std::array<PointIdxT, 2>> prevnext( points_.size() );
		for( const auto& e : edges_ ) {
			prevnext[e.from()][1] = e.idx();
			prevnext[e.to()][0] = e.idx();

			points_[e.from()].edgeIdx() = e.idx();
		}
		for( const auto& pn : prevnext ) {
			edges_[pn[0]].next() = pn[1];
			edges_[pn[1]].prev() = pn[0];
		}

		// add local surface estimation
		for( auto& contourPoint : points_ ) {
			// local coords
			const auto& p0 = contourPoint.p0();
			const auto& edgeNextIdx = contourPoint.edgeIdx();
			const auto& pNext = points_[edges_[edgeNextIdx].to()].p0();		// new pLeft
			const auto& edgePrevIdx = edges_[edgeNextIdx].prev();
			const auto& pPrev = points_[edges_[edgePrevIdx].from()].p0();	// new pRight

			auto t = pNext - pPrev;
			decltype( t ) n{ t[1], -t[0], t[2] };
			t = t.normalized();
			n = n.normalized();

			contourPoint.normal() = n;
			contourPoint.tangential() = t;

			// polyParam
			float ATA = 0;
			float ATb = 0;
			std::vector<TriMesh::Point> points = { pNext, pPrev };
			for( int i = 0; i < points.size(); ++i ) {
				TriMesh::Point loc = contourPoint.map_glob2loc( points[i] );
				float u = loc[0];
				float v = loc[1];
				ATb += u * u * v;
				ATA += u * u * u * u;
			}
			contourPoint.surfParam() = ATb / ATA;
		}

		// validity check
		for( const auto& e : edges_ ) {
			LOG_ASSERT( e.idx() != INVALID_EDGE );
			LOG_ASSERT( e.next() != INVALID_EDGE );
			LOG_ASSERT( e.prev() != INVALID_EDGE );
			LOG_ASSERT( e.from() != INVALID_POINT );
			LOG_ASSERT( e.to() != INVALID_POINT );
		}
		for( const auto& p : points_ ) {
			LOG_ASSERT( p.idx() != INVALID_POINT );
			LOG_ASSERT( p.edgeIdx() != INVALID_EDGE );
		}

		// init kd tree
		kdTree_ = std::make_unique<kd_tree_t>( 3, *this );
		kdTree_->buildIndex();
	}

	Point Shape::findNearestNeighbor( const TriMesh::Point& p ) const {
		constexpr int k = 1;
		std::vector<size_t> ret_indexes( k );
		std::vector<float> out_dists_sqr( k );
		nanoflann::KNNResultSet<float> resultSet( k );
		resultSet.init( &ret_indexes[0], &out_dists_sqr[0] );
		kdTree_->findNeighbors( resultSet, &p[0], nanoflann::SearchParams() );

		return points_[ret_indexes[0]];
	}

	std::vector<Point> Shape::findNearestNeighborNormalPruned( const TriMesh::Point& p, const TriMesh::Point& n ) const {
		constexpr int k = 10;
		std::vector<size_t> ret_indexes( k );
		std::vector<float> out_dists_sqr( k );
		nanoflann::KNNResultSet<float> resultSet( k );
		resultSet.init( &ret_indexes[0], &out_dists_sqr[0] );
		kdTree_->findNeighbors( resultSet, &p[0], nanoflann::SearchParams() );

		std::vector<Point> nnVec;
		for( const auto& idx : ret_indexes ) {
			const auto& pn = points_[idx];
			const auto& nn = pn.normal();
			// check if normals have the same direction
			if( ( nn | n ) < 0.3 )
				continue;
			else
				nnVec.push_back( pn );
		}

		return nnVec;
	}

	TriMesh::Point Shape::mapToContour( const TriMesh::Point& p ) const {
		const auto nn = findNearestNeighbor( p );
		const auto& e = edges_[nn.edgeIdx()];
		const auto& ePrev = edges_[e.prev()];
		const auto& nnNext = points_[e.to()];
		const auto& nnPrev = points_[ePrev.from()];

		const auto& p0 = nn.p0();
		const auto& pNext = nnNext.p0();
		const auto& pPrev = nnPrev.p0();

		const auto e1 = ( pPrev - p0 );
		const auto l1 = e1.normalized() | ( p - p0 );
		LOG_ASSERT( 1.9f * l1 <= e1.length() ) << "    " << 2 * l1 << " <= " << e1.length();		// TODO remove later, just for checking that the computation is correct
		auto p1 = p0 + l1 * e1.normalized();
		
		const auto e2 = ( pNext - p0 );
		const auto l2 = e2.normalized() | ( p - p0 );
		LOG_ASSERT( 1.9f * l2 <= e2.length() );		// TODO remove later
		auto p2 = p0 + l2 * e2.normalized();
		
		if( l1 <= 0 && l2 <= 0 ) {
			return p0;
		}
		if( l1 > 0 && l2 <= 0 ) {
			return p1;
		}
		if( l1 <= 0 && l2 > 0 ) {
			return p2;
		}

		auto p1_ = p0 + l1 * e1.normalized();
		auto p2_ = p0 + l2 * e2.normalized();
		if( ( p - p1_ ).length() < ( p - p2_ ).length() ) {
			return p1_;
		} else {
			return p2_;
		}

		return p0;
	}

	std::optional<TriMesh::Point> Shape::mapToContourWithNormal( const TriMesh::Point& p, const TriMesh::Point& n ) const {
		const auto nnVec = findNearestNeighborNormalPruned( p, n );
		if( nnVec.empty() ) {
			return {};
		}

		std::vector<TriMesh::Point> projPoints;

		for( const auto& nn : nnVec ) {
			//const auto& nn = nnVec[0];
			const auto& e = edges_[nn.edgeIdx()];
			const auto& ePrev = edges_[e.prev()];
			const auto& nnNext = points_[e.to()];
			const auto& nnPrev = points_[ePrev.from()];

			const auto& p0 = nn.p0();
			const auto& pNext = nnNext.p0();
			const auto& pPrev = nnPrev.p0();

			const auto e1 = ( pPrev - p0 );
			const auto l1 = e1.normalized() | ( p - p0 );
			if( 1.9f * l1 > e1.length() ) {
				continue;
			}
			auto p1 = p0 + l1 * e1.normalized();

			const auto e2 = ( pNext - p0 );
			const auto l2 = e2.normalized() | ( p - p0 );
			if( 1.9f * l2 > e2.length() ) {
				continue;
			}
			auto p2 = p0 + l2 * e2.normalized();

			if( l1 <= 0 && l2 <= 0 ) {
				projPoints.push_back( p0 );
			} else if( l1 > 0 && l2 <= 0 ) {
				projPoints.push_back( p1 );
			} else if( l1 <= 0 && l2 > 0 ) {
				projPoints.push_back( p2 );
			} else {
				auto p1_ = p0 + l1 * e1.normalized();
				auto p2_ = p0 + l2 * e2.normalized();
				if( ( p - p1_ ).length() < ( p - p2_ ).length() ) {
					projPoints.push_back( p1_ );
				} else {
					projPoints.push_back( p2_ );
				}
			}
		}

		if( projPoints.empty() ) {
			return {};
		}

		TriMesh::Point projBest;
		float dist = std::numeric_limits<float>::max();
		for( const auto& pp : projPoints ) {
			float d = ( pp - p ).length();
			if( d < dist ) {
				dist = d;
				projBest = pp;
			}
		}

		return projBest;
	}

	std::optional<TriMesh::Point> Shape::mapToContourWithNormal( const OpenMesh::SmartVertexHandle& v, const TriMesh& m ) const {
		if( !v.is_boundary() )
			return {};

		const auto h = v.halfedge();
		LOG_ASSERT( h.is_boundary() );

		const auto& p = m.point( v );
		const auto& pPrev = m.point( h.prev().from() );
		const auto& pNext = m.point( h.to() );
		const auto e = ( pNext - pPrev );
		OpenMesh::Vec3f n{ e[1], -e[0], 0 };
		n.normalize();

		const auto pNew = mapToContourWithNormal( p, n );

		return pNew;
	}

	float Shape::signedDistanceToNearestNeighbor( const TriMesh::Point& p ) const {
		const auto nn = findNearestNeighbor( p );
		const auto& e = edges_[nn.edgeIdx()];
		const auto& ePrev = edges_[e.prev()];
		const auto& nnNext = points_[e.to()];
		const auto& nnPrev = points_[ePrev.from()];

		const auto& p0 = nn.p0();
		const auto& pNext = nnNext.p0();
		const auto& pPrev = nnPrev.p0();

		const auto e1 = ( pPrev - p0 );
		const auto l1 = e1.normalized() | ( p - p0 );
		auto p1 = p0 + l1 * e1.normalized();

		const auto e2 = ( pNext - p0 );
		const auto l2 = e2.normalized() | ( p - p0 );
		auto p2 = p0 + l2 * e2.normalized();

		OpenMesh::Vec3f n1{ e1[1], -e1[0], 0 };
		n1.normalize();
		OpenMesh::Vec3f n2{ -e2[1], e2[0], 0 };
		n2.normalize();
		auto n0 = ( n1 + n2 ).normalized();

		if( l1 <= 0 && l2 <= 0 ) {
			if( ( ( p - p0 ) | n0 ) < 0 ) {
				return -( p - p0 ).length();
			} else {
				return ( p - p0 ).length();
			}
		}
		if( l1 > 0 && l2 <= 0 ) {
			if( ( ( p - p0 ) | n1 ) < 0 ) {
				return -( p - p1 ).length();
			} else {
				return ( p - p1 ).length();
			}
		}
		if( l1 <= 0 && l2 > 0 ) {
			if( ( ( p - p0 ) | n2 ) < 0 ) {
				return -( p - p2 ).length();
			} else {
				return ( p - p2 ).length();
			}
		}

		auto p1_ = p0 + l1 * e1.normalized();
		auto p2_ = p0 + l2 * e2.normalized();
		if( ( p - p1_ ).length() < ( p - p2_ ).length() ) {
			if( ( ( p - p0 ) | n1 ) < 0 ) {
				return -( p - p1 ).length();
			} else {
				return ( p - p1 ).length();
			}
		} else {
			if( ( ( p - p0 ) | n2 ) < 0 ) {
				return -( p - p2 ).length();
			} else {
				return ( p - p2 ).length();
			}
		}
	}

	float Shape::signedDistanceExact( const TriMesh::Point& p ) const {
		float signedDist = std::numeric_limits<float>::max();
		for( const auto& e : edges_ ) {
			float s = signedDistanceToEdge( e, p );
			if( std::abs( s ) < std::abs( signedDist ) ) {
				signedDist = s;
			}
		}

		return signedDist;
	}

	TriMesh::Point Shape::computeEdgeNormal( const Edge& e ) const {
		const auto& v0 = points_[e.from()];
		const auto& v1 = points_[e.to()];
		const auto& p0 = v0.p0();
		const auto& p1 = v1.p0();

		// compute signed distance to edge
		const auto eVec = ( p1 - p0 );

		OpenMesh::Vec3f normalEdge{ -eVec[1], eVec[0], 0 };
		normalEdge.normalize();

		return normalEdge;
	}

	float Shape::distance( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const {
		const auto nn1 = findNearestNeighbor( p1 );
		const auto nn2 = findNearestNeighbor( p2 );

		// find connection between nn1 and nn2 with BFS
		std::vector<PointIdxT> root( points_.size(), INVALID_POINT );
		std::vector<bool> visited( points_.size(), false );
		std::queue<long> q;
		q.push( nn1.idx() );
		visited[nn1.idx()] = true;
		while( !q.empty() ) {
			const auto idx = q.front();
			q.pop();
			if( idx == nn2.idx() ) {
				break;
			}

			const auto& e = edges_[points_[idx].edgeIdx()];
			const auto pNextIdx = e.to();
			if( !visited[pNextIdx] ) {
				q.push( pNextIdx );
				root[pNextIdx] = idx;
				visited[pNextIdx] = true;
			}
			const auto& ePrev = edges_[e.prev()];
			const auto pPrevIdx = ePrev.from();
			if( !visited[pPrevIdx] ) {
				q.push( pPrevIdx );
				root[pPrevIdx] = idx;
				visited[pPrevIdx] = true;
			}
		}

		// use backtracking to construct path
		std::vector<PointIdxT> path;
		auto pathIter = nn2.idx();
		while( pathIter != nn1.idx() ) {
			path.push_back( pathIter );
			pathIter = root[pathIter];
			if( pathIter == INVALID_POINT ) {
				path.clear();
				break;
			}
		}
		path.push_back( nn1.idx() );

		if( path.size() == 1 ) {
			return ( p1 - p2 ).length();
		}

		float dist = (p2 - points_[path[1]].p0()).length();
		dist += ( p1 - points_[path[path.size() - 2]].p0() ).length();

		for( int i = 1; i < path.size() - 2; ++i ) {
			dist += ( points_[path[i]].p0() - points_[path[i + 1]].p0() ).length();
		}
		

		return dist;
	}

	float Shape::pointDistance( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const {
		const auto nn1 = findNearestNeighbor( p1 );
		const auto nn2 = findNearestNeighbor( p2 );

		if( nn1.idx() == nn2.idx() ) {
			const auto& e = edges_[nn1.edgeIdx()];
			const auto& pNext = points_[e.to()];
			const auto d = ( p1 - p2 ).length() / ( nn1.p0() - pNext.p0() ).length();
			return d;
		}

		// find connection between nn1 and nn2 with BFS
		std::vector<PointIdxT> root( points_.size(), INVALID_POINT );
		std::vector<bool> visited( points_.size(), false );
		std::queue<long> q;
		q.push( nn1.idx() );
		visited[nn1.idx()] = true;
		while( !q.empty() ) {
			const auto idx = q.front();
			q.pop();
			if( idx == nn2.idx() ) {
				break;
			}

			const auto& e = edges_[points_[idx].edgeIdx()];
			const auto pNextIdx = e.to();
			if( !visited[pNextIdx] ) {
				q.push( pNextIdx );
				root[pNextIdx] = idx;
				visited[pNextIdx] = true;
			}
			const auto& ePrev = edges_[e.prev()];
			const auto pPrevIdx = ePrev.from();
			if( !visited[pPrevIdx] ) {
				q.push( pPrevIdx );
				root[pPrevIdx] = idx;
				visited[pPrevIdx] = true;
			}
		}

		// use backtracking to construct path
		std::vector<PointIdxT> path;
		auto pathIter = nn2.idx();
		while( pathIter != nn1.idx() ) {
			pathIter = root[pathIter];
			if( pathIter == nn1.idx() )
				break;
			else if( pathIter == INVALID_POINT ) {
				path.clear();
				break;
			}
			path.push_back( pathIter );
		}

		if( path.empty() ) {
			return ( p1 - p2 ).length() / ( nn1.p0() - nn2.p0() ).length();
		}

		// first distance
		auto d1 = ( p1 - points_[path[path.size() - 1]].p0() ).length();
		auto d1r = ( nn1.p0() - points_[path[path.size() - 1]].p0() ).length();
		auto r1 = d1 / d1r;
		// last distance
		auto d2 = ( p2 - points_[path[0]].p0() ).length();
		auto d2r = ( nn2.p0() - points_[path[0]].p0() ).length();
		auto r2 = d2 / d2r;

		float dist = r1 + r2 + path.size() - 1;

		return dist;
	}

	TriMesh::Point Shape::computeMid( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& pOld ) const {
		const auto nn1 = findNearestNeighbor( p1 );
		const auto nn2 = findNearestNeighbor( p2 );
		const auto nnOld = findNearestNeighbor( pOld );

		

		std::vector<PointIdxT> path;
		float r1 = 0, r2 = 0;

		float dist = 0;

		if( nn1.idx() == nn2.idx() == nnOld.idx() ) {
			const auto& e = edges_[nn1.edgeIdx()];
			const auto& pNext = points_[e.to()];
			dist = ( p1 - p2 ).length() / ( nn1.p0() - pNext.p0() ).length();
		} else {
			// find connection between nn1 and nn2 with BFS
			std::vector<PointIdxT> root( points_.size(), INVALID_POINT );
			std::vector<bool> visited( points_.size(), false );
			std::queue<long> q1;
			q1.push( nn1.idx() );
			visited[nn1.idx()] = true;
			while( !q1.empty() ) {
				const auto idx = q1.front();
				q1.pop();
				if( idx == nnOld.idx() ) {
					break;
				}
				//if( idx == nn2.idx() ) {
				//	root[idx] = INVALID_POINT;
				//	visited[idx] = false;
				//	continue;
				//}

				const auto& e = edges_[points_[idx].edgeIdx()];
				const auto pNextIdx = e.to();
				if( !visited[pNextIdx] && pNextIdx != nn2.idx() ) {
					q1.push( pNextIdx );
					root[pNextIdx] = idx;
					visited[pNextIdx] = true;
				}
				const auto& ePrev = edges_[e.prev()];
				const auto pPrevIdx = ePrev.from();
				if( !visited[pPrevIdx] && pPrevIdx != nn2.idx() ) {
					q1.push( pPrevIdx );
					root[pPrevIdx] = idx;
					visited[pPrevIdx] = true;
				}
			}
			
			std::queue<long> q2;
			q2.push( nnOld.idx() );
			while( !q2.empty() ) {
				const auto idx = q2.front();
				q2.pop();
				if( idx == nn2.idx() ) {
					break;
				}

				const auto& e = edges_[points_[idx].edgeIdx()];
				const auto pNextIdx = e.to();
				if( !visited[pNextIdx] ) {
					q2.push( pNextIdx );
					root[pNextIdx] = idx;
					visited[pNextIdx] = true;
				}
				const auto& ePrev = edges_[e.prev()];
				const auto pPrevIdx = ePrev.from();
				if( !visited[pPrevIdx] ) {
					q2.push( pPrevIdx );
					root[pPrevIdx] = idx;
					visited[pPrevIdx] = true;
				}
			}

			// use backtracking to construct path
			auto pathIter = nn2.idx();
			while( pathIter != nn1.idx() ) {
				pathIter = root[pathIter];
				if( pathIter == nn1.idx() )
					break;
				else if( pathIter == INVALID_POINT ) {
					path.clear();
					break;
				}
				path.push_back( pathIter );
			}

			if( path.empty() ) {
				dist = ( p1 - p2 ).length() / ( nn1.p0() - nn2.p0() ).length();
			} else {
				// first distance
				auto d1 = ( p1 - points_[path[path.size() - 1]].p0() ).length();
				auto d1r = ( nn1.p0() - points_[path[path.size() - 1]].p0() ).length();
				r1 = d1 / d1r;
				// last distance
				auto d2 = ( p2 - points_[path[0]].p0() ).length();
				auto d2r = ( nn2.p0() - points_[path[0]].p0() ).length();
				r2 = d2 / d2r;

				dist = r1 + r2 + path.size() - 1;
			}
		}

		auto p = p1;
		
		if( path.empty() ) {
			// TODO do not just take the average but use weights based on the edge lengths at the nearest neighbors
			p = 0.5 * ( p1 + p2 );
			p = mapToContour( p );
			return p;
		}

		// move along path until mid point is reached
		float distCovered = r1;
		p = points_[path[path.size() - 1]].p0();
		if( distCovered > 0.5 * dist ) {
			p = 0.5 * ( p1 + p2 );
			p = mapToContour( p );
			return p;
		}
		int i;
		for( i = static_cast<int>(path.size()) - 2; i >= 0 && (0.5 * dist - distCovered) > 1; --i ) {
			distCovered += 1;
		}
		p = points_[path[i + 1]].p0();
		
		float delta = 0.5f * dist - distCovered;
		if( i > 0 ) {
			auto e = ( points_[path[i]].p0() - p );
			p += delta * e;
		}
		
		return p;
	}

	TriMesh::Point Shape::computeMid( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const {
		const auto nn1 = findNearestNeighbor( p1 );
		const auto nn2 = findNearestNeighbor( p2 );
		std::vector<PointIdxT> path;
		float r1 = 0, r2 = 0;

		float dist = 0;

		if( nn1.idx() == nn2.idx() ) {
			const auto& e = edges_[nn1.edgeIdx()];
			const auto& pNext = points_[e.to()];
			dist = ( p1 - p2 ).length() / ( nn1.p0() - pNext.p0() ).length();
		} else {
			// find connection between nn1 and nn2 with BFS
			std::vector<PointIdxT> root( points_.size(), INVALID_POINT );
			std::vector<bool> visited( points_.size(), false );
			std::queue<long> q1;
			q1.push( nn1.idx() );
			visited[nn1.idx()] = true;
			while( !q1.empty() ) {
				const auto idx = q1.front();
				q1.pop();
				if( idx == nn2.idx() ) {
					break;
				}
				//if( idx == nn2.idx() ) {
				//	root[idx] = INVALID_POINT;
				//	visited[idx] = false;
				//	continue;
				//}

				const auto& e = edges_[points_[idx].edgeIdx()];
				const auto pNextIdx = e.to();
				if( !visited[pNextIdx] ) {
					q1.push( pNextIdx );
					root[pNextIdx] = idx;
					visited[pNextIdx] = true;
				}
				const auto& ePrev = edges_[e.prev()];
				const auto pPrevIdx = ePrev.from();
				if( !visited[pPrevIdx] ) {
					q1.push( pPrevIdx );
					root[pPrevIdx] = idx;
					visited[pPrevIdx] = true;
				}
			}

			// use backtracking to construct path
			auto pathIter = nn2.idx();
			while( pathIter != nn1.idx() ) {
				pathIter = root[pathIter];
				if( pathIter == nn1.idx() )
					break;
				else if( pathIter == INVALID_POINT ) {
					path.clear();
					break;
				}
				path.push_back( pathIter );
			}

			if( path.empty() ) {
				dist = ( p1 - p2 ).length() / ( nn1.p0() - nn2.p0() ).length();
			} else {
				// first distance
				auto d1 = ( p1 - points_[path[path.size() - 1]].p0() ).length();
				auto d1r = ( nn1.p0() - points_[path[path.size() - 1]].p0() ).length();
				r1 = d1 / d1r;
				// last distance
				auto d2 = ( p2 - points_[path[0]].p0() ).length();
				auto d2r = ( nn2.p0() - points_[path[0]].p0() ).length();
				r2 = d2 / d2r;

				dist = r1 + r2 + path.size() - 1;
			}
		}

		auto p = p1;

		if( path.empty() ) {
			// TODO do not just take the average but use weights based on the edge lengths at the nearest neighbors
			p = 0.5 * ( p1 + p2 );
			p = mapToContour( p );
			return p;
		}

		// move along path until mid point is reached
		float distCovered = r1;
		p = points_[path[path.size() - 1]].p0();
		if( distCovered > 0.5 * dist ) {
			p = 0.5 * ( p1 + p2 );
			p = mapToContour( p );
			return p;
		}
		int i;
		for( i = static_cast<int>(path.size()) - 2; i >= 0 && ( 0.5 * dist - distCovered ) > 1; --i ) {
			distCovered += 1;
		}
		p = points_[path[i + 1]].p0();

		float delta = 0.5f * dist - distCovered;
		if( i > 0 ) {
			auto e = ( points_[path[i]].p0() - p );
			p += delta * e;
		}

		return p;
	}

	TriMesh::Point Shape::computeMidWithout( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& pBreak ) const {
		const auto nn1 = findNearestNeighbor( p1 );
		const auto nn2 = findNearestNeighbor( p2 );
		const auto nnBreak = findNearestNeighbor( pBreak );
		LOG_ASSERT( nn1.idx() != nnBreak.idx() );
		LOG_ASSERT( nn2.idx() != nnBreak.idx() );

		std::vector<PointIdxT> path;
		float r1 = 0, r2 = 0;

		float dist = 0;

		if( nn1.idx() == nn2.idx() ) {
			const auto& e = edges_[nn1.edgeIdx()];
			const auto& pNext = points_[e.to()];
			dist = ( p1 - p2 ).length() / ( nn1.p0() - pNext.p0() ).length();
		} else {
			// find connection between nn1 and nn2 with BFS
			std::vector<PointIdxT> root( points_.size(), INVALID_POINT );
			std::vector<bool> visited( points_.size(), false );
			// set break point to visited
			visited[nnBreak.idx()] = true;

			std::queue<long> q1;
			q1.push( nn1.idx() );
			visited[nn1.idx()] = true;
			while( !q1.empty() ) {
				const auto idx = q1.front();
				q1.pop();
				if( idx == nn2.idx() ) {
					break;
				}

				const auto& e = edges_[points_[idx].edgeIdx()];
				const auto pNextIdx = e.to();
				if( !visited[pNextIdx] ) {
					q1.push( pNextIdx );
					root[pNextIdx] = idx;
					visited[pNextIdx] = true;
				}
				const auto& ePrev = edges_[e.prev()];
				const auto pPrevIdx = ePrev.from();
				if( !visited[pPrevIdx] ) {
					q1.push( pPrevIdx );
					root[pPrevIdx] = idx;
					visited[pPrevIdx] = true;
				}
			}

			// use backtracking to construct path
			auto pathIter = nn2.idx();
			while( pathIter != nn1.idx() ) {
				pathIter = root[pathIter];
				if( pathIter == nn1.idx() )
					break;
				else if( pathIter == INVALID_POINT ) {
					path.clear();
					break;
				}
				path.push_back( pathIter );
			}

			if( path.empty() ) {
				dist = ( p1 - p2 ).length() / ( nn1.p0() - nn2.p0() ).length();
			} else {
				// first distance
				auto d1 = ( p1 - points_[path[path.size() - 1]].p0() ).length();
				auto d1r = ( nn1.p0() - points_[path[path.size() - 1]].p0() ).length();
				r1 = d1 / d1r;
				// last distance
				auto d2 = ( p2 - points_[path[0]].p0() ).length();
				auto d2r = ( nn2.p0() - points_[path[0]].p0() ).length();
				r2 = d2 / d2r;

				dist = r1 + r2 + path.size() - 1;
			}
		}

		auto p = p1;

		if( path.empty() ) {
			// TODO do not just take the average but use weights based on the edge lengths at the nearest neighbors
			p = 0.5 * ( p1 + p2 );
			p = mapToContour( p );
			return p;
		}

		// move along path until mid point is reached
		float distCovered = r1;
		p = points_[path[path.size() - 1]].p0();
		if( distCovered > 0.5 * dist ) {
			p = 0.5 * ( p1 + p2 );
			p = mapToContour( p );
			return p;
		}
		int i;
		for( i = static_cast<int>( path.size() ) - 2; i >= 0 && ( 0.5 * dist - distCovered ) > 1; --i ) {
			distCovered += 1;
		}
		p = points_[path[i + 1]].p0();

		float delta = 0.5f * dist - distCovered;
		if( i > 0 ) {
			auto e = ( points_[path[i]].p0() - p );
			p += delta * e;
		}

		return p;
	}

	ContourTypeT Shape::getContourType( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const {
		auto pMid = 0.5 * ( p1 + p2 );
		pMid[2] = 0;
		const auto nn = findNearestNeighbor( pMid );
		const auto& e = edges_[nn.edgeIdx()];
		const auto& ePrev = edges_[e.prev()];
		const auto& nnNext = points_[e.to()];
		const auto& nnPrev = points_[ePrev.from()];

		const auto& p0 = nn.p0();
		const auto& pNext = nnNext.p0();
		const auto& pPrev = nnPrev.p0();

		if( e.contourType() == ePrev.contourType() )
			return e.contourType();

		// ePrev
		const auto e1 = ( pPrev - p0 );
		const auto l1 = e1.normalized() | ( pMid - p0 );

		// e
		const auto e2 = ( pNext - p0 );
		const auto l2 = e2.normalized() | ( pMid - p0 );

		if( l1 < 0 && l2 < 0 ) {
			TriMesh::Point tEdge = ( p2 - p1 );
			TriMesh::Point nEdge = { tEdge[1], -tEdge[0], 0 };
			nEdge.normalize();
			TriMesh::Point ne1 = { e1[1], -e1[0], 0 };
			ne1.normalize();
			TriMesh::Point ne2 = { e2[1], -e2[0], 0 };
			ne2.normalize();
			float s1 = std::abs( nEdge | ne1 );
			float s2 = std::abs( nEdge | ne2 );
			if( s1 > s2 ) {
				return ePrev.contourType();
			} else {
				return e.contourType();
			}
		}
		if( l1 > 0 && l2 <= 0 ) {
			return ePrev.contourType();
		}
		if( l1 <= 0 && l2 > 0 ) {
			return e.contourType();
		}

		auto p1_ = p0 + l1 * e1.normalized();
		auto p2_ = p0 + l2 * e2.normalized();
		if( ( pMid - p1_ ).length() < ( pMid - p2_ ).length() ) {
			return ePrev.contourType();
		} else {
			return e.contourType();
		}
		
	}

	std::tuple<std::vector<float>, std::vector<float>> Shape::getEmoEfaVals( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const {
		auto pMid = 0.5 * ( p1 + p2 );
		pMid[2] = 0;
		const auto nn = findNearestNeighbor( pMid );
		const auto& e1 = edges_[nn.edgeIdx()];
		const auto& e0 = edges_[e1.prev()];
		const auto& nnNext = points_[e1.to()];
		const auto& nnPrev = points_[e0.from()];

		const auto emo0it = emo_.find( e0.idx() );
		const auto emo1it = emo_.find( e1.idx() );
		const auto efa0it = efa_.find( e0.idx() );
		const auto efa1it = efa_.find( e1.idx() );

		const bool e0exists = emo0it != emo_.end();
		const bool e1exists = emo1it != emo_.end();

		if( !e0exists && !e1exists ) {
			return std::tuple<std::vector<float>, std::vector<float>>();
		}
		if( !e0exists ) {
			return std::make_tuple( emo1it->second, efa1it->second );
		}
		if( !e1exists ) {
			return std::make_tuple( emo0it->second, efa0it->second );
		}

		const auto& emo0 = emo0it->second;
		const auto& emo1 = emo1it->second;
		const auto& efa0 = efa0it->second;
		const auto& efa1 = efa1it->second;

		const auto& p0 = nn.p0();
		const auto& pNext = nnNext.p0();
		const auto& pPrev = nnPrev.p0();
		
		const auto pe0 = 0.5f * ( pPrev + p0 );
		const auto pe1 = 0.5f * ( p0 + pNext );

		const auto u = ( pMid - pe0 ).length() / ( pe1 - pe0 ).length();

		std::vector<float> emoMid( emo0.size() );
		std::vector<float> efaMid( efa0.size() );

		for( auto i = 0; i < emo0.size(); ++i ) {
			emoMid[i] = ( 1 - u ) * emo0[i] + u * emo1[i];
			efaMid[i] = ( 1 - u ) * efa0[i] + u * efa1[i];
		}

		return std::make_tuple( emoMid, efaMid );
	}
}