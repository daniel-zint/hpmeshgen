#pragma once

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm> 
#include <tuple>
#include <array>
#include <cmath>


// Open Mesh
#include "MeshHeader.h"

namespace HelperFunctions
{
// Calculate area of triangle formed by (x1, y1), (x2, y2) and (x3, y3)
	inline float triangleArea( float x1, float y1, float x2, float y2, float x3, float y3 ) {
		return std::fabs( ( x1*( y2 - y3 ) + x2*( y3 - y1 ) + x3*( y1 - y2 ) ) * 0.5f );
	}

	inline float signedTriangleArea( float x1, float y1, float x2, float y2, float x3, float y3 ) {
		return ( x1*( y2 - y3 ) + x2*( y3 - y1 ) + x3*( y1 - y2 ) ) * 0.5f;
	}

	inline float triangleArea( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& p3 ) {
		return triangleArea( p1[0], p1[1], p2[0], p2[1], p3[0], p3[1] );
	}

	inline auto barycentricCoordinates( const std::array<TriMesh::Point, 3> &triangle, const TriMesh::Point& p ) {
		const float A = signedTriangleArea( triangle[0][0], triangle[0][1], triangle[1][0], triangle[1][1], triangle[2][0], triangle[2][1] );
		const float Ainv = 1.f / A;
		const float A0 = signedTriangleArea( p[0], p[1], triangle[1][0], triangle[1][1], triangle[2][0], triangle[2][1] );
		const float A1 = signedTriangleArea( triangle[0][0], triangle[0][1], p[0], p[1], triangle[2][0], triangle[2][1] );
		const float A2 = signedTriangleArea( triangle[0][0], triangle[0][1], triangle[1][0], triangle[1][1], p[0], p[1] );

		return std::make_tuple( A0 * Ainv, A1 * Ainv, A2 * Ainv );
	}

	inline float barycentricInterpolation( const std::vector<TriMesh::Point>& triangle, const float& a, const float& b, const float& c ) {
		return a * triangle[0][2] + b * triangle[1][2] + c * triangle[2][2];
	}

	// Check whether point p4 lies within the triangle (p1,p2,p3)
	// DEPRECATED !!!!!!
	inline bool isInsideTriangle( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& p3, const TriMesh::Point& p4 ) {
		const float A = triangleArea( p1, p2, p3 );
		const float A1 = triangleArea( p4, p2, p3 );
		const float A2 = triangleArea( p1, p4, p3 );
		const float A3 = triangleArea( p1, p2, p4 );

		// Check if sum of A1, A2 and A3 is same as A
		const float A_sum = A1 + A2 + A3;
		const float eps = 0.0001f * A;
		return ( A <= A_sum + eps ) && ( A >= A_sum - eps );
	}

	// Check if barycentric coordinates are inside the triangle, i.e. very close to zero or negative.
	inline bool isInsideTriangle( const float& a, const float& b, const float& c ) {
		const float eps = static_cast<float>(1e-5);
		if ( a < eps || b < eps || c < eps ) {
			return false;
		} else {
			return true;
		}
	}

	// Calculate angle around point p2
	inline float calcAngle( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& p3 ) {
		const float dot = ( p1[0] - p2[0] ) * ( p3[0] - p2[0] ) + ( p1[1] - p2[1] ) * ( p3[1] - p2[1] );
		const float det = ( p1[0] - p2[0] ) * ( p3[1] - p2[1] ) - ( p1[1] - p2[1] ) * ( p3[0] - p2[0] );
		return ( 180.f / (float)M_PI ) * std::fabs( atan2( det, dot ) );
	}

	// Calculate angle around point p2
	inline float calcAngle2( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& p3 ) {
		const float dot = ( p1[0] - p2[0] ) * ( p3[0] - p2[0] ) + ( p1[1] - p2[1] ) * ( p3[1] - p2[1] );
		const float det = ( p1[0] - p2[0] ) * ( p3[1] - p2[1] ) - ( p1[1] - p2[1] ) * ( p3[0] - p2[0] );
		const auto ang = ( 180.f / (float)M_PI ) * atan2( det, dot );
		if( ang > 0 )
			return ang;
		else
			return 360 + ang;
	}

	// Calculate cross product of the vectors (p1,p2) and (p1,p3)
	inline float calcCrossProduct( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& p3 ) {
		return ( p2[0] - p1[0] ) * ( p3[1] - p1[1] ) - ( p2[1] - p1[1] ) * ( p3[0] - p1[0] );
	}


	inline float calcSquaredDistance( const TriMesh::Point& p1, const TriMesh::Point& p2 ) {
		return ( p1[0] - p2[0] ) * ( p1[0] - p2[0] ) + ( p1[1] - p2[1] ) * ( p1[1] - p2[1] );
	}
	// Calculate distance between two points
	inline float calcDistance( const TriMesh::Point& p1, const TriMesh::Point& p2 ) {
		return sqrt( calcSquaredDistance( p1, p2 ) );
	}

	inline void bezierReduction( const std::vector<PolyMesh::Point>& cp, std::vector<PolyMesh::Point>& cpNew ) {
		size_t n = cp.size() - 1;

		cpNew = std::vector<PolyMesh::Point>( n, { 0,0,0 } );

		cpNew[0][0] = cp[0][0];
		cpNew[0][1] = cp[0][1];
		cpNew[n - 1][0] = cp[n][0];
		cpNew[n - 1][1] = cp[n][1];

		for ( size_t i = 0; i < n * n; ++i ) {
			for ( size_t j = 1; j < n - 1; ++j ) {
				float a = (float)( ( n - j ) * j );
				float b = (float)( ( n - j ) * ( n - j ) + ( j + 1 ) * ( j + 1 ) );
				float c = (float)( ( n - j - 1 ) * ( j + 1 ) );

				cpNew[j][0] = ( n * ( n - j ) * cp[j][0] + n * ( j + 1 ) * cp[j + 1][0]
								- a * cpNew[j - 1][0]
								- c * cpNew[j + 1][0] ) / b;
				cpNew[j][1] = ( n * ( n - j ) * cp[j][1] + n * ( j + 1 ) * cp[j + 1][1]
								- a * cpNew[j - 1][1]
								- c * cpNew[j + 1][1] ) / b;
			}
		}
	}

	inline void bezierElevation( const std::vector<PolyMesh::Point>& cp, std::vector<PolyMesh::Point>& cpNew ) {
		size_t n = cp.size();

		cpNew = std::vector<PolyMesh::Point>( n + 1, { 0,0,0 } );
		cpNew[0][0] = cp[0][0];
		cpNew[0][1] = cp[0][1];
		cpNew[n][0] = cp[n - 1][0];
		cpNew[n][1] = cp[n - 1][1];

		for ( size_t i = 1; i < n; ++i ) {
			cpNew[i][0] = (float)i / (float)n * cp[i - 1][0] + ( 1.f - (float)i / (float)n ) * cp[i][0];
			cpNew[i][1] = (float)i / (float)n * cp[i - 1][1] + ( 1.f - (float)i / (float)n ) * cp[i][1];
		}

	}

	inline void project4areaConsistency( const std::vector<PolyMesh::Point>& cp, PolyMesh::Point& p ) {

		// compute bezier reduction
		p[0] = 0.25f * ( 3.f * cp[1][0] + 3.f * cp[2][0] - cp[0][0] - cp[3][0] );
		p[1] = 0.25f * ( 3.f * cp[1][1] + 3.f * cp[2][1] - cp[0][1] - cp[3][1] );

		// get base
		float base[2] = { cp[3][0] - cp[0][0], cp[3][1] - cp[0][1] };
		float orth[2] = { -base[1], base[0] };

		float area = ( cp[2][0] - cp[0][0] ) * ( cp[1][1] - cp[3][1] ) - ( cp[2][1] - cp[0][1] ) * ( cp[1][0] - cp[3][0] );

		//float lambda = (area + base[0] * (cp[0][1] - p[1]) - base[1] * (cp[0][0] - p[0])) / (base[0] * orth[1] - base[1] * orth[0]);
		float lambda = ( area + orth[0] * ( cp[0][0] - p[0] ) + orth[1] * ( cp[0][1] - p[1] ) ) / ( orth[0] * orth[0] + orth[1] * orth[1] );

		// project point onto parallel line to base s.t. the triangle has the same area as the (unreduced) quad.
		p[0] = p[0] + lambda * orth[0];
		p[1] = p[1] + lambda * orth[1];
	}

}