#pragma once

#include "cuda_runtime.h"

#include "UniformGrid.h"
#include "DmoParams.h"
#include "../HPMeshGen2/SignedDistanceFunction.h"

namespace DMO
{
	namespace Metrics
	{
		struct NoMetric
		{
			__device__ __forceinline__ float vertexQuality( const float2* oneRing, const int oneRingSize, const float2& pCurr, const float2& p ) const {
				return -FLT_MAX;
			}
		};

		struct MeanRatioTriangle
		{
			static constexpr int elementSize = 3;

			__device__ float triangleQuality( const float2 p[3] ) const {
				float2 e[3];
				float e_length_squared[3];
			
				for( int i = 0; i < 3; ++i ) {
					int j = ( i + 1 ) % 3;
					e[i].x = p[j].x - p[i].x;
					e[i].y = p[j].y - p[i].y;
			
					e_length_squared[i] = e[i].x * e[i].x + e[i].y * e[i].y;
				}
			
				float l = e_length_squared[0] + e_length_squared[1] + e_length_squared[2];
				float area = e[0].x * e[1].y - e[0].y * e[1].x;
			
				if( area < 0 )
					return area;
				else
					return 2.f * sqrtf( 3.f ) * area / l;
			}
			
			__device__ __forceinline__ float vertexQuality( const float2* oneRing, const int oneRingSize, const float2& pCurr, const float2& p ) const {
				float q = FLT_MAX;
			
				for( int k = 0; k < oneRingSize - 1; ++k ) {
					float2 t[3] = { { p.x, p.y },{ oneRing[k].x, oneRing[k].y },{ oneRing[k + 1].x, oneRing[k + 1].y } };
					q = fminf( q, triangleQuality( t ) );
				}
			
				return q;
			}
		};

		struct DensityTriangle
		{
			static constexpr int elementSize = 3;

			UniformGrid grid_d_;
			
			DensityTriangle( UniformGrid grid_d ) : grid_d_{ grid_d } {}
			
			__device__ __forceinline__ float vertexQuality( const float2* oneRing, const int oneRingSize, const float2& pCurr, const float2& p ) const {
				//float q = FLT_MAX;
				float qShort = FLT_MAX;
				float qLong = 0;
			
				float hVals[MAX_ONE_RING_SIZE];
				for( int k = 0; k < oneRingSize - 1; ++k ) {
					float2 v = oneRing[k];
					float h = grid_d_.uniformGridValue( v );
					//float h = 1;
					hVals[k] = h;
				}
			
				float hp = grid_d_.uniformGridValue( p );
				for( int k = 0; k < oneRingSize - 1; ++k ) {
					float2 v = oneRing[k];
					int kLeft = ( k + oneRingSize - 1 - 1 ) % ( oneRingSize - 1 );
					int kRight = ( k + 1 ) % ( oneRingSize - 1 );
					float hLeft = hVals[kLeft];
					float hRight = hVals[kRight];
					float h = ( hVals[k] + hLeft + hRight + hp ) / 4.f;
					//float h = 1;
					float e = sqrtf( ( p.x - v.x ) * ( p.x - v.x ) + ( p.y - v.y ) * ( p.y - v.y ) );
					
					//float qEdge;
					//if( e < h ) qEdge = e / h;
					//else qEdge = h / e;
					//q = fminf( q, qEdge );
			
					float relLength = e / h;
					qShort = fminf( relLength, qShort );
					qLong = fmaxf( relLength, qLong );
				}
			
				//return q;
				return 1.f / ( qLong - qShort + 0.0001f );
			}

		};

		struct DensityWithMeanRatioTriangle
		{
			static constexpr int elementSize = 3;
			MeanRatioTriangle meanRatioTriangle_;
			DensityTriangle densityTriangle_;
			
			DensityWithMeanRatioTriangle( MeanRatioTriangle meanRatioTriangle, DensityTriangle densityTriangle )
				: meanRatioTriangle_{ meanRatioTriangle }, densityTriangle_{ densityTriangle } {}
			
			__device__ __forceinline__ float vertexQuality( const float2* oneRing, const int oneRingSize, const float2& pCurr, const float2& p ) const {
			
				float qMeanRatio = meanRatioTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
			
				const float minShapeQuality = 0.3f;
				if( qMeanRatio > minShapeQuality ) {
					float qDens = densityTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
					//return minShapeQuality + ( 1.0f * qDens + 0.0f * qMeanRatio );
					return minShapeQuality + qDens;
				} else {
					return qMeanRatio;
				}
			}
		};

		
		struct SdfTriangle
		{
			static constexpr int elementSize = 3;

			UniformGrid sdf_d_;

			SdfTriangle( UniformGrid sdf_d ) : sdf_d_{ sdf_d } {}

			__device__ __forceinline__ float vertexQuality( const float2* oneRing, const int oneRingSize, const float2& pCurr, const float2& p ) const {
				float v = sdf_d_.uniformGridValue( p );
				return v;
			}
		};

		struct SdfConstrainedMeanRatioTriangle
		{
			static constexpr int elementSize = 3;
			MeanRatioTriangle meanRatioTriangle_;
			SdfTriangle sdfTriangle_;

			SdfConstrainedMeanRatioTriangle( MeanRatioTriangle meanRatioTriangle, SdfTriangle sdfTriangle )
				: meanRatioTriangle_{ meanRatioTriangle }, sdfTriangle_{ sdfTriangle } {}

			__device__ __forceinline__ float vertexQuality( const float2* oneRing, const int oneRingSize, const float2& pCurr, const float2& p ) const {

				float qSdfOld = sdfTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, pCurr );
				float qSdf = sdfTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
				float qMeanRatio = meanRatioTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
				if( qSdf < 0 ) {
					return -FLT_MAX;
				} else if( qSdf - qSdfOld > 0 ) {
					return -FLT_MAX;
				} else {
					return qMeanRatio;
				}
			}
		};

		struct SdfMinConstrainedMeanRatioTriangle
		{
			static constexpr int elementSize = 3;
			MeanRatioTriangle meanRatioTriangle_;
			SdfTriangle sdfTriangle_;

			SdfMinConstrainedMeanRatioTriangle( MeanRatioTriangle meanRatioTriangle, SdfTriangle sdfTriangle )
				: meanRatioTriangle_{ meanRatioTriangle }, sdfTriangle_{ sdfTriangle } {}

			__device__ __forceinline__ float vertexQuality( const float2* oneRing, const int oneRingSize, const float2& pCurr, const float2& p ) const {

				float qSdfOld = sdfTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, pCurr );
				float qSdf = sdfTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
				float qMeanRatio = meanRatioTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
				if( qSdf < 0 ) {
					return -FLT_MAX;
				} else {
					return qMeanRatio;
				}
			}
		};

		struct SdfWithMeanRatioTriangle
		{
			static constexpr int elementSize = 3;
			MeanRatioTriangle meanRatioTriangle_;
			SdfTriangle sdfTriangle_;

			SdfWithMeanRatioTriangle( MeanRatioTriangle meanRatioTriangle, SdfTriangle sdfTriangle )
				: meanRatioTriangle_{ meanRatioTriangle }, sdfTriangle_{ sdfTriangle } {}

			__device__ __forceinline__ float vertexQuality( const float2* oneRing, const int oneRingSize, const float2& pCurr, const float2& p ) const {

				//float qSdfOld = sdfTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, pCurr );
				//float qSdf = sdfTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
				//float qMeanRatio = meanRatioTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
				//if( qSdf < 0 ) {
				//	return -FLT_MAX;
				//} else if( qSdf - qSdfOld > 0 ) {
				//	return -FLT_MAX;
				//} else {
				//	return qMeanRatio;
				//}

				float qMeanRatio = meanRatioTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );

				const float minShapeQuality = 0.3f;
				if( qMeanRatio > minShapeQuality ) {
					float qSdf = sdfTriangle_.vertexQuality( oneRing, oneRingSize, pCurr, p );
					return minShapeQuality + 1.f / ( fabsf( qSdf ) + 0.0001 );
				} else {
					return qMeanRatio;
				}
			}
		};
	} // namespace Metrics
} // namespace DMO

#ifndef REGISTER_METRIC
#define REGISTER_METRIC(metric)
#endif

REGISTER_METRIC( Metrics::NoMetric );
REGISTER_METRIC( Metrics::MeanRatioTriangle );
REGISTER_METRIC( Metrics::DensityTriangle );
REGISTER_METRIC( Metrics::DensityWithMeanRatioTriangle );
REGISTER_METRIC( Metrics::SdfConstrainedMeanRatioTriangle );
REGISTER_METRIC( Metrics::SdfMinConstrainedMeanRatioTriangle );
REGISTER_METRIC( Metrics::SdfWithMeanRatioTriangle );
