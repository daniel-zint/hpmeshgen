#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cstdint>
#include <stdio.h>
#include <cfloat>

#include <type_traits>

#include "math_constants.h"

#include "Vertex.h"
#include "DmoParams.h"
#include "DmoMesh.h"
#include "gpuErrchk.h"
#include "DMO_PUBLIC.h"

namespace DMO
{

	typedef union
	{
		float floats[2];             // floats[0] = lowest
		std::int32_t ints[2];        // ints[1] = lowIdx
		unsigned long long ulong;    // for atomic update
	} floatIntUnion;

	__device__ unsigned long long floatIntArgMax( unsigned long long* address, float val, std::int32_t idx ) {
		floatIntUnion loc, newValue;
		loc.floats[0] = val;
		loc.ints[1] = idx;
		newValue.ulong = *address;
		while( newValue.floats[0] < val )
			newValue.ulong = atomicCAS( address, newValue.ulong, loc.ulong );

		return newValue.ulong;
	}

	template<typename MetricT>
	__global__ void optimizeHierarchical( float2* points, const int* coloredVertexIDs, const int cOff, const DmoVertex* vertices, const int* oneRingVec, const MetricT metric ) {

		const uint2 idx1 = {
			threadIdx.x / NQ,
			threadIdx.x % NQ
		};
		const uint2 idx2 = {
			( threadIdx.x + NQ * NQ / 2 ) / NQ,
			( threadIdx.x + NQ * NQ / 2 ) % NQ
		};

		const DmoVertex& v = vertices[coloredVertexIDs[cOff + blockIdx.x]];

		float q = -FLT_MAX;

		__shared__ float2 currPos;	// current optimal position
		__shared__ float2 maxDist;

		__shared__ floatIntUnion argMaxVal;
		argMaxVal.floats[0] = -FLT_MAX;
		argMaxVal.ints[1] = NQ * NQ;

		__shared__ float2 oneRing[MAX_ONE_RING_SIZE];

		// min/max search + loading oneRing
		if( threadIdx.x == 0 ) {
			maxDist.x = -FLT_MAX;
			maxDist.y = -FLT_MAX;

			for( int k = 0; k < v.oneRingSize - 1; ++k ) {
				float2 vo = points[oneRingVec[v.oneRingID + k]];
				oneRing[k] = vo;

				float2 dist = {
					abs( points[v.idx].x - vo.x ),
					abs( points[v.idx].y - vo.y )
				};

				maxDist.x = fmaxf( maxDist.x, dist.x );
				maxDist.y = fmaxf( maxDist.y, dist.y );
			}

			oneRing[v.oneRingSize - 1] = points[oneRingVec[v.oneRingID + v.oneRingSize - 1]];

			currPos = points[v.idx];
		}

		// start depth iteration
		float depth_scale = GRID_SCALE;
		for( int depth = 0; depth < DEPTH; ++depth ) {

			float2 aabbMin, aabbMax;	// axis aligned bounding box
			aabbMin.x = currPos.x - depth_scale * maxDist.x;
			aabbMin.y = currPos.y - depth_scale * maxDist.y;
			aabbMax.x = currPos.x + depth_scale * maxDist.x;
			aabbMax.y = currPos.y + depth_scale * maxDist.y;

			float2 p1 = {
				AFFINE_FACTOR * ( idx1.x * aabbMin.x + ( NQ - 1 - idx1.x ) * aabbMax.x ),
				AFFINE_FACTOR * ( idx1.y * aabbMin.y + ( NQ - 1 - idx1.y ) * aabbMax.y )
			};
			float2 p2 = {
				AFFINE_FACTOR * ( idx2.x * aabbMin.x + ( NQ - 1 - idx2.x ) * aabbMax.x ),
				AFFINE_FACTOR * ( idx2.y * aabbMin.y + ( NQ - 1 - idx2.y ) * aabbMax.y )
			};

			float q1 = metric.vertexQuality( oneRing, v.oneRingSize, points[v.idx], p1 );
			float q2 = metric.vertexQuality( oneRing, v.oneRingSize, points[v.idx], p2 );
			//float q1 = metric.vertexQuality( oneRing, v.oneRingSize, p1 );
			//float q2 = metric.vertexQuality( oneRing, v.oneRingSize, p2 );

			float argMax = 0;
			if( q1 > q2 ) {
				q = q1;
				argMax = 1;
			} else {
				q = q2;
				argMax = 2;
			}

			__syncwarp();
			//__syncthreads();
			floatIntArgMax( ( unsigned long long* ) & ( argMaxVal.ulong ), q, idx1.x * NQ + idx1.y );

			float qOld = metric.vertexQuality( oneRing, v.oneRingSize, points[v.idx], currPos );
			//float qOld = metric.vertexQuality( oneRing, v.oneRingSize, currPos );
			if( idx1.x * NQ + idx1.y == argMaxVal.ints[1] && qOld < q ) {
				if( argMax == 1 ) {
					currPos.x = p1.x;
					currPos.y = p1.y;
				} else {
					currPos.x = p2.x;
					currPos.y = p2.y;
				}
			}

			//rescale candidate grid to the size of two cells
			depth_scale *= 2 * AFFINE_FACTOR;
		}

		// set new position if it is better than the old one
		float qOld = metric.vertexQuality( oneRing, v.oneRingSize, points[v.idx], points[v.idx] );
		//float qOld = metric.vertexQuality( oneRing, v.oneRingSize, points[v.idx] );
		if( idx1.x * NQ + idx1.y == argMaxVal.ints[1] && qOld < q ) {
			points[v.idx].x = currPos.x;
			points[v.idx].y = currPos.y;
		}
	}

	template<typename MetricT>
	void DMO_PUBLIC dmoGPU( thrust::device_vector<float2>& points_d, DmoMesh& dmoMesh1_, const MetricT& metric ) {
		for( auto colorIt = dmoMesh1_.colorOffset_.begin(); colorIt != dmoMesh1_.colorOffset_.end() - 1; ++colorIt ) {
			const int nBlocks = *( colorIt + 1 ) - *colorIt;
			optimizeHierarchical << <nBlocks, N_THREADS >> > (
				points_d.data().get(),
				dmoMesh1_.coloredVertexIDs_d.data().get(),
				*colorIt,
				dmoMesh1_.vertices_d.data().get(),
				dmoMesh1_.oneRingVec_d.data().get(),
				metric
				);
			gpuErrchk( cudaDeviceSynchronize() );
		}
	}
}

#define REGISTER_METRIC(metric) \
namespace DMO { template void DMO_PUBLIC dmoGPU<>( thrust::device_vector<float2>&, DmoMesh&, const metric& ); }

#include "Metrics.h"