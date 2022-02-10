#pragma once

#include <vector>
#include <map>

#include "HPMeshGen2/OceanMesh/OceanMesh.h"
#include <nanoflann.hpp>

#include "Point.h"
#include "Edge.h"

namespace Contour
{
	class Shape
	{
		using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, Shape >, Shape, 3>;
		std::unique_ptr<kd_tree_t> kdTree_;

		class : public std::vector<Point>
		{
		public:
			Point& add() {
				push_back( Point( (PointIdxT)this->size() ) );
				return this->back();
			}
		} points_;

		class : public std::vector<Edge>
		{
		public:
			Edge& add() {
				push_back( Edge( (EdgeIdxT)this->size() ) );
				return this->back();
			}
		} edges_;

		std::map<EdgeIdxT, std::vector<float>> emo_;
		std::map<EdgeIdxT, std::vector<float>> efa_;

	public:
		Shape() = default;
		Shape( OceanMesh& mesh );
		Shape( TriMesh& mesh );
		Shape( const Shape& s) {
			points_ = s.points_;
			edges_ = s.edges_;
			emo_ = s.emo_;
			efa_ = s.efa_;

			// init kd tree
			kdTree_ = std::make_unique<kd_tree_t>( 3, *this );
			kdTree_->buildIndex();
		};
		Shape& operator= ( const Shape& s ) {
			points_ = s.points_;
			edges_ = s.edges_;
			emo_ = s.emo_;
			efa_ = s.efa_;

			// init kd tree
			kdTree_ = std::make_unique<kd_tree_t>( 3, *this );
			kdTree_->buildIndex();

			return *this;
		};


		const auto& points() { return points_; }
		const auto& edges() { return edges_; }

		Point findNearestNeighbor( const TriMesh::Point& p ) const;

		/// <summary>
		/// Find nearest neighbor with approximately the same normal
		/// </summary>
		/// <param name="p">point</param>
		/// <param name="n">normal</param>
		/// <returns>optional point</returns>
		std::vector<Point> findNearestNeighborNormalPruned( const TriMesh::Point& p, const TriMesh::Point& n ) const;

		TriMesh::Point mapToContour( const TriMesh::Point& p ) const;

		std::optional<TriMesh::Point> mapToContourWithNormal( const TriMesh::Point& p, const TriMesh::Point& n ) const;

		std::optional<TriMesh::Point> mapToContourWithNormal( const OpenMesh::SmartVertexHandle& v, const TriMesh& m) const;

		float signedDistanceToNearestNeighbor( const TriMesh::Point& p ) const;

		/// <summary>
		/// Computes signed distance of a point to all edges of shape and takes the absolute minimum
		/// </summary>
		/// <param name="p">point</param>
		/// <returns></returns>
		float signedDistanceExact( const TriMesh::Point& p ) const;

		/// <summary>
		/// Signed distance of a point to an edge of shape.
		/// </summary>
		/// <param name="e">edge</param>
		/// <param name="p">point</param>
		/// <returns></returns>
		float signedDistanceToEdge( const Edge& e, const TriMesh::Point& p ) const;

		TriMesh::Point computeEdgeNormal( const Edge& e ) const;

		float distance( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const;

		float pointDistance( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const;

		TriMesh::Point computeMid( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& pOld ) const;

		TriMesh::Point computeMid( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const;

		TriMesh::Point computeMidWithout( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& pBreak ) const;

		ContourTypeT getContourType( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const;

		std::tuple< std::vector<float>, std::vector<float>> getEmoEfaVals( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const;

		// next(pointIdx)
		// prev(pointIdx)

		/*----- nanoflann functions -----*/
		// Must return the number of data points
		inline size_t kdtree_get_point_count() const { return points_.size(); }
		// Returns the dim'th component of the idx'th point in the class
		inline float kdtree_get_pt( const size_t idx, const size_t dim ) const {
			return points_[idx].p0()[dim];
		}
		template <class BBOX>
		bool kdtree_get_bbox( BBOX& ) const { return false; }
	};


	inline float Shape::signedDistanceToEdge( const Edge& e, const TriMesh::Point& p ) const {
		const auto& v0 = points_[e.from()];
		const auto& v1 = points_[e.to()];
		const auto& p0 = v0.p0();
		const auto& p1 = v1.p0();

		// compute signed distance to edge
		const auto eVec = ( p1 - p0 );
		const auto alpha = eVec.normalized() | ( p - p0 );
		auto pProj = p0 + alpha * eVec.normalized();	// projection of p onto line defined by e

		auto normalEdge = computeEdgeNormal( e );
		auto normalPrev = computeEdgeNormal( edges_[e.prev()] );
		auto normalNext = computeEdgeNormal( edges_[e.next()] );
		auto n0 = ( normalEdge + normalPrev ).normalized();
		auto n1 = ( normalEdge + normalNext ).normalized();


		if( alpha < 0 ) {
			float s0 = ( ( p - p0 ) | n0 ) < 0 ? -1 : 1;
			return s0 * ( p - p0 ).length();
		} else if( alpha > eVec.length() ) {
			float s1 = ( ( p - p0 ) | n1 ) < 0 ? -1 : 1;
			return s1 * ( p - p1 ).length();
		} else {
			float s = ( ( p - p0 ) | normalEdge ) < 0 ? -1 : 1;
			return s * ( p - pProj ).length();
		}
	}
}