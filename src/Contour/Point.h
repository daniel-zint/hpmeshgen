#pragma once

#include "../HPMeshGen2/OceanMesh/OceanMesh.h"

#include "Types.h"

namespace Contour
{
	class Point
	{
		PointIdxT idx_ = INVALID_POINT;
		TriMesh::Point p0_;
		TriMesh::Point normal_;
		TriMesh::Point tangential_;
		float surfParam_ = 0;

		EdgeIdxT edgeIdx_ = INVALID_EDGE;

	public:
		Point(int idx, TriMesh::Point p0, TriMesh::Point normal, TriMesh::Point tangential, float surfParam )
			: idx_{ idx }, p0_{ p0 }, normal_{ normal }, tangential_{ tangential }, surfParam_{ surfParam } {}
		Point() = default;
		Point( PointIdxT idx ) : idx_ { idx } {}

		const auto& idx() const { return idx_; }
		auto& p0() { return p0_; }
		const auto& p0() const { return p0_; }
		auto& normal() { return normal_; }
		const auto& normal() const { return normal_; }
		auto& tangential() { return tangential_; }
		const auto& tangential() const { return tangential_; }
		auto& surfParam() { return surfParam_; }
		const auto& surfParam() const { return surfParam_; }
		auto& edgeIdx() { return edgeIdx_; }
		const auto& edgeIdx() const { return edgeIdx_; }

		TriMesh::Point map_loc2glob( float u, float v ) const {
			auto glob = p0_ + v * normal_ + u * tangential_;
			return { glob[0], glob[1], 0 };
		}

		TriMesh::Point map_glob2loc( TriMesh::Point p ) const {
			float u = ( p - p0_ ) | ( tangential_ );
			float v = ( p - p0_ ) | ( normal_ );
			return { u, v, 0 };
		}

		float locSurf( const float u ) const {
			return surfParam_ * u * u;
		}
	};
}