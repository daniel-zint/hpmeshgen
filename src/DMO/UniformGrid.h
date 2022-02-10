#pragma once

#include "cuda_runtime.h"

namespace DMO
{
	struct UniformGrid
	{
		int2 n;
		float2 h;
		float2 aabbMin, aabbMax;
		float* vals;

		__device__  __forceinline__ float uniformGridValue( float2 p ) const  {
	
			if( p.x < aabbMin.x ) p.x = aabbMin.x;
			else if( p.x > aabbMax.x ) p.x = aabbMax.x;
			if( p.y < aabbMin.y ) p.y = aabbMin.y;
			else if( p.y > aabbMax.y ) p.y = aabbMax.y;

			int i1, i2, j1, j2;

			i1 = (int)( ( p.x - aabbMin.x ) / h.x );
			j1 = (int)( ( p.y - aabbMin.y ) / h.y );
			i2 = i1 + 1;
			j2 = j1 + 1;
			if( i2 == n.x ) {
				i1--;
				i2--;
			}
			if( j2 == n.y ) {
				j1--;
				j2--;
			}

	
			float x1 = aabbMin.x + h.x * i1;
			float y1 = aabbMin.y + h.y * j1;
			float x2 = aabbMin.x + h.x * i2;
			float y2 = aabbMin.y + h.y * j2;

			float pref = 1.f / ( ( x2 - x1 ) * ( y2 - y1 ) );
			float v1 = x2 - p.x;
			float v2 = p.x - x1;
			float v3 = vals[i1 * n.y + j1];
			float v4 = vals[i1 * n.y + j2];
			float v5 = vals[i2 * n.y + j1];
			float v6 = vals[i2 * n.y + j2];
			float v7 = y2 - p.y;
			float v8 = p.y - y1;

			return pref * ( v1 * ( v3 * v7 + v4 * v8 ) + v2 * ( v5 * v7 + v6 * v8 ) );
		}

		/// <summary>
		/// WARNING: UNTESTED CODE!!!! Compute gradient of bilinear interpolant at point p
		/// </summary>
		/// <param name="p">point</param>
		/// <returns></returns>
		__device__  __forceinline__ float2 cellGradient( float2 p ) const {
			if( p.x < aabbMin.x ) p.x = aabbMin.x;
			else if( p.x > aabbMax.x ) p.x = aabbMax.x;
			if( p.y < aabbMin.y ) p.y = aabbMin.y;
			else if( p.y > aabbMax.y ) p.y = aabbMax.y;

			int i1, i2, j1, j2;

			i1 = (int)( ( p.x - aabbMin.x ) / h.x );
			j1 = (int)( ( p.y - aabbMin.y ) / h.y );
			i2 = i1 + 1;
			j2 = j1 + 1;
			if( i2 == n.x ) {
				i1--;
				i2--;
			}
			if( j2 == n.y ) {
				j1--;
				j2--;
			}


			float x1 = aabbMin.x + h.x * i1;
			float y1 = aabbMin.y + h.y * j1;
			float x2 = aabbMin.x + h.x * i2;
			float y2 = aabbMin.y + h.y * j2;

			float a = vals[i1 * n.y + j1];
			float b = vals[i1 * n.y + j2];
			float c = vals[i2 * n.y + j1];
			float d = vals[i2 * n.y + j2];

			float u = ( p.x - x1 ) / ( x2 - x1 );
			float v = ( p.y - y1 ) / ( y2 - y1 );

			float2 g;
			g.x = ( 1 - v ) * ( b - a ) + v * ( d - c );
			g.y = ( 1 - u ) * ( c - a ) + u * ( d - b );

			return g;
		}
	};
} // namespace DMO
