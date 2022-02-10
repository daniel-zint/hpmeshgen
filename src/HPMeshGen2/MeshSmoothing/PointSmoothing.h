#pragma once

#include <vector>

#include "../HelperFunctions.h"

// Open Mesh
#include "../../MeshHeader.h"
 
namespace PointSmoothing {
	inline void getConvexDomain(TriMesh& mesh, const TriMesh::VertexHandle& vh, std::vector<TriMesh::Point>& convexDomainPoints) {
		assert(convexDomainPoints.size() == 0);

		// find all lines
		std::vector<std::vector<TriMesh::Point>> lines;
		for (auto voh_it = mesh.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
			TriMesh::HalfedgeHandle heh = mesh.next_halfedge_handle(*voh_it);
			TriMesh::Point p1 = mesh.point(mesh.from_vertex_handle(heh));
			TriMesh::Point p2 = mesh.point(mesh.to_vertex_handle(heh));

			lines.push_back({ p1, p2 });
		}

		// find all intersection points
		std::vector<TriMesh::Point> points;
		for (size_t i = 0; i < lines.size(); ++i) {
			const OpenMesh::Vec2f alpha_i = OpenMesh::Vec2f(lines[i][0][0], lines[i][0][1]);
			const OpenMesh::Vec2f beta_i= OpenMesh::Vec2f(lines[i][1][0] - lines[i][0][0], lines[i][1][1] - lines[i][0][1]);

			for (size_t j = i + 1; j < lines.size(); ++j) {
				const OpenMesh::Vec2f alpha_j= OpenMesh::Vec2f(lines[j][0][0], lines[j][0][1]);
				const OpenMesh::Vec2f beta_j= OpenMesh::Vec2f(lines[j][1][0] - lines[j][0][0], lines[j][1][1] - lines[j][0][1]);

				// Check if vectors are linear dependent (parallel)
				if (std::abs(beta_i[0] * beta_j[1] - beta_i[1] * beta_j[0]) < 0.0001) {
					//std::cout << "Linear dependency found " << std::abs(beta_i[0] * beta_j[1] - beta_i[1] * beta_j[0]) << std::endl;
					//std::cout << "Linear dependency found. Angle = " << angle << std::endl;
					//std::cout << beta_j[0] * beta_i[1] / (beta_i[0] * beta_j[1]) << std::endl;
					continue;
				}

				float zaehler = beta_i[0] * (alpha_i[1] - alpha_j[1]) + beta_i[1] * (alpha_j[0] - alpha_i[0]);
				float nenner = beta_i[0] * beta_j[1] - beta_j[0] * beta_i[1];
				float lambda_j = zaehler / nenner;

				float lambda_i = (alpha_j[0] + beta_j[0] * lambda_j - alpha_i[0]) / beta_i[0];	// just needed for checking if the calculated point is inside the fan

				const OpenMesh::Vec2f v_intersect = alpha_j + lambda_j * beta_j;

				// catch all the points that can't be part of the convex domain because they are for sure outside of the fan
				if ((lambda_i > 1 || lambda_i < 0) && (lambda_j > 1 || lambda_j < 0))
					continue;

				TriMesh::Point p = { v_intersect[0], v_intersect[1], 0.f };
				points.push_back(p);
			}
		}

		// delete all points that are not part of the convex environment
		std::vector<bool> points_flag(points.size(), true);
		for (size_t i = 0; i < lines.size(); ++i) {
			const OpenMesh::Vec2f v1 = OpenMesh::Vec2f(lines[i][0][0], lines[i][0][1]);
			const OpenMesh::Vec2f v2 = OpenMesh::Vec2f(lines[i][1][0], lines[i][1][1]);

			const  OpenMesh::Vec2f base = v2 - v1;

			// orthonormal vector to the edge (v1,v2)
			const float norm_base = base.length();
			const  OpenMesh::Vec2f n = OpenMesh::Vec2f(-base[1]/norm_base, base[0]/norm_base);

			for (size_t j = 0; j < points.size(); ++j) {
				float distance = n[0] * points[j][0] + n[1] * points[j][1] - n[0] * v1[0] - n[1] * v1[1];
				//float length_base = base.norm();
				//already calculated in line 77 (norm_base)

				if (distance < -0.001 * norm_base) {
					// point is not part of the convex environment
					points_flag[j] = false;
				}
			}
		}

		std::vector<TriMesh::Point> points_buf;
		for (size_t i = 0; i < points.size(); ++i) {
			if (points_flag[i]) {
				convexDomainPoints.push_back(points[i]);
			}
		}
	}

	inline void smoothConvexDomain(TriMesh& mesh) {

		for( auto vh : mesh.vertices() ) {
			TriMesh::Point p_center = mesh.point( vh );

			// leave boundary as it is
			if( mesh.is_boundary( vh ) )
				continue;

			// convex Domain
			std::vector<TriMesh::Point> convexDomainPoints;
			getConvexDomain( mesh, vh, convexDomainPoints );

			if( convexDomainPoints.size() == 0 ) {
				continue;
			}

			convexDomainPoints.push_back( p_center );

			TriMesh::Point p_center_new = { 0.f,0.f,0.f };
			for( size_t i = 0; i < convexDomainPoints.size(); ++i ) {
				p_center_new += convexDomainPoints[i];
			}

			p_center_new /= convexDomainPoints.size();

			mesh.set_point( vh, p_center_new );
		}
		
	}

}
