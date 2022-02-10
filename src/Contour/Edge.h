#pragma once

#include "Types.h"
#include "ContourType.h"

namespace Contour
{
	class Edge
	{
		EdgeIdxT idx_ = INVALID_EDGE;
		EdgeIdxT nextEdgeIdx_ = INVALID_EDGE;
		EdgeIdxT prevEdgeIdx_ = INVALID_EDGE;
		PointIdxT fromPointIdx_ = INVALID_POINT;
		PointIdxT toPointIdx_ = INVALID_POINT;
		ContourTypeT contourType_ = ContourType::NOTHING;

	public:
		Edge( EdgeIdxT idx, PointIdxT fromIdx, PointIdxT toIdx, ContourTypeT contourId ) : idx_{ idx }, fromPointIdx_ { fromIdx }, toPointIdx_{ toIdx }, contourType_{ contourId } {}
		Edge() = default;
		Edge( EdgeIdxT idx ) : idx_{ idx } {}

		const auto& idx() const { return idx_; }
		auto& from() { return fromPointIdx_; }
		const auto& from() const { return fromPointIdx_; }
		auto& to() { return toPointIdx_; }
		const auto& to() const { return toPointIdx_; }
		auto& next() { return nextEdgeIdx_; }
		const auto& next() const { return nextEdgeIdx_; }
		auto& prev() { return prevEdgeIdx_; }
		const auto& prev() const { return prevEdgeIdx_; }
		auto& contourType() { return contourType_; }
		const auto& contourType() const { return contourType_; }
	};
}