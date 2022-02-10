#pragma once

#include <vector>
#include "cuda_runtime.h" 

#include "../../MeshHeader.h"	// TODO reduce to OpenMesh::VecXf s.t. functions can be used also on GPU

namespace QualityMetrics
{
	// mean ratio metric for a triangle
	inline float meanRatioMetric( const OpenMesh::Vec3f (&triangle)[3] ) {
		OpenMesh::Vec2f p[3];
		for ( int i = 0; i < 3; ++i ) {
			p[i][0] = triangle[i][0];
			p[i][1] = triangle[i][1];
		}
		
		auto e12 = p[0] - p[1];
		auto e23 = p[1] - p[2];
		auto e13 = p[0] - p[2];

		float l = e12.sqrnorm() + e23.sqrnorm() + e13.sqrnorm();
		float area = e23[1] * e13[0] - e23[0] * e13[1];

		return 2.f * std::sqrt( 3.f ) * area / l;
	}
	// mean ratio metric for a triangle
	inline float meanRatioMetric( TriMesh& mesh, const TriMesh::FaceHandle& fh ) {
		// From the paper "R. Rangarajan, A.J. Lew: Directional Vertex Relaxation"
		std::vector<TriMesh::Point> meshPoints;
		for ( auto fv_it = mesh.fv_iter( fh ); fv_it.is_valid(); ++fv_it ) {
			meshPoints.push_back( mesh.point( *fv_it ) );
		}

		std::vector<OpenMesh::Vec2f> p( 3 );
		for ( size_t i = 0; i < 3; ++i ) {
			p[i][0] = meshPoints[i][0];
			p[i][1] = meshPoints[i][1];
		}

		auto e12 = p[0] - p[1];
		auto e23 = p[1] - p[2]; 
		auto e13 = p[0] - p[2]; 

		float l = e12.sqrnorm() + e23.sqrnorm() + e13.sqrnorm();
		float area = e23[1] * e13[0] - e23[0] * e13[1];

		return 2.f * std::sqrt( 3.f ) * area / l;
	}
	// mean ratio for a point in a triangle mesh
	inline float meanRatioMetric( const std::vector<OpenMesh::Vec3f>& oneRing, const OpenMesh::Vec3f& p ) {
		float q = FLT_MAX;

		for ( int i = 0; i < oneRing.size() - 1; ++i ) {
			int j = i + 1;

			float qBuf = meanRatioMetric( { oneRing[i], oneRing[j], p } );

			q = fminf( q, qBuf );
		}

		return q;
	}
	
	// mean ratio metric for quads
	inline float meanRatioMetric( PolyMesh& mesh, const TriMesh::FaceHandle& fh ) {
		// From the paper "R. Rangarajan, A.J. Lew: Directional Vertex Relaxation"
		std::vector<TriMesh::Point> points;
		for ( auto fv_it = mesh.fv_iter( fh ); fv_it.is_valid(); ++fv_it ) {
			points.push_back( mesh.point( *fv_it ) );
		}

		float e[4][2];
		float e_length_squared[4];
		for ( size_t i = 0; i < 4; ++i ) {
			int j = ( i + 1 ) % 4;
			e[i][0] = points[j][0] - points[i][0];
			e[i][1] = points[j][1] - points[i][1];

			e_length_squared[i] = e[i][0] * e[i][0] + e[i][1] * e[i][1];
		}

		float l = e_length_squared[0] + e_length_squared[1] + e_length_squared[2] + e_length_squared[3];
		float area1 = e[0][0] * e[1][1] - e[0][1] * e[1][0];
		float area2 = e[2][0] * e[3][1] - e[2][1] * e[3][0];

		return 2.f * ( area1 + area2 ) / l;

	}

	// condition number for vertex '1'
	inline float conditionSurfVertex( const std::vector<OpenMesh::Vec3f>& points, const OpenMesh::Vec3f normal ) {
		LOG_ASSERT( points.size() == 3 );

		//			 2
		// 			 |
		//			 e1
		//	 		 |
		// 0---e2----1

		OpenMesh::Vec3f e1 = points[2] - points[1];
		OpenMesh::Vec3f e2 = points[0] - points[1];

		auto l1 = e1.sqrnorm();
		auto l2 = e2.sqrnorm();

		auto detJ = ( e1 % e2 ) | normal;

		if( detJ < 0 ) {
			return detJ;
		} else {
			auto c = 2 * detJ / ( l1 + l2 );
			return c;
		}
	}

	inline float condition( const std::vector<OpenMesh::Vec3f>& points ) {
		const auto ps = points.size();

		std::vector<OpenMesh::Vec3f> edges( ps, { 0,0,0 } );
		std::vector<float> edgeSqrNorm( ps, -1 );

		for( auto i = 0; i < ps; ++i ) {
			auto j = ( i + 1 ) % ps;
			edges[i] = points[j] - points[i];
			edges[i][2] = 0;
			edgeSqrNorm[i] = edges[i].sqrnorm();
		}

		auto det = []( const OpenMesh::Vec3f& p1, const OpenMesh::Vec3f& p2 ) {
			return ( p1 % p2 )[2];
		};

		std::vector<float> detJ( ps, -1 );
		for( auto i = 0; i < ps; ++i ) {
			auto il = ( i - 1 + ps ) % ps;
			detJ[i] = det( edges[i], -edges[il] );
		}


		std::vector<float> c( ps, -1 );

		for( auto i = 0; i < ps; ++i ) {
			auto il = ( i - 1 + ps ) % ps;
			c[i] = detJ[i] / ( edgeSqrNorm[i] + edgeSqrNorm[il] );
		}

		//auto cMin = std::min_element( c.begin(), c.end() );
		//return 2 * *cMin;

		auto detMin = std::min_element( detJ.begin(), detJ.end() );

		if( *detMin > 0 ) {
			auto cMin = std::min_element( c.begin(), c.end() );
			return 2 * *cMin;
		} else {
			return *detMin;
		}
	}

	inline float condition( PolyMesh& mesh, const PolyMesh::FaceHandle& fh ) {
		std::vector<PolyMesh::Point> p;
		p.reserve( 4 );
		for( auto fv : mesh.fv_range( fh ) ) {
			p.push_back( mesh.point( fv ) );
		}

		return condition( p );
	}

	inline auto minCondition( PolyMesh& mesh ) {
		float qMin = FLT_MAX;
		float qMean = 0;

		for( auto fh : mesh.faces() ) {
			auto q = condition( mesh, fh );
			qMin = std::min( q, qMin );
			qMean += q;

		}
		qMean /= mesh.n_faces();
		return std::make_tuple( qMin, qMean );
	}
	
	inline float minAngle( const std::vector<OpenMesh::Vec3f>& points ) {
		const auto ps = points.size();
		
		auto calcAngle = []( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& p3 ) {
			const float dot = ( p1[0] - p2[0] ) * ( p3[0] - p2[0] ) + ( p1[1] - p2[1] ) * ( p3[1] - p2[1] );
			const float det = ( p1[0] - p2[0] ) * ( p3[1] - p2[1] ) - ( p1[1] - p2[1] ) * ( p3[0] - p2[0] );
			return ( 180.f / (float)M_PI ) * std::fabs( atan2( det, dot ) );
		};

		float minAngle = FLT_MAX;
		for( auto i = 0; i < ps; ++i ) {
			auto ir = ( i + 1 ) % ps;
			auto il = ( i - 1 + ps ) % ps;

			float a = calcAngle( points[ir], points[i], points[il] );
			minAngle = std::min( minAngle, a );
		}

		return minAngle;
	}
}
