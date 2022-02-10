#pragma once

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm> 
#include <experimental/filesystem>
#include <queue>

// Open Mesh
#include "MeshHeader.h"
#include "OceanMesh/OceanMesh.h"

// Size Grid
#include "ScalarField/ScalarField.h"

// Point Smoothing
#include "MeshSmoothing/PointSmoothing.h"

// Helpers
#include "HelperFunctions.h"
#include "MeshSmoothing/QualityMetrics.h"

#include "Contour/Shape.h"

// Block Structured Grid
#include "BlockStructuredGrid/BlockStructuredGrid.h"
#include "BlockStructuredGrid/BlockStructuredTriangleGrid.h"

#include "svg.h"


////////////////////
// cuda functions //
#define Q_MEANRATIO 0
#define Q_JACOBIAN 4
#define Q_MINANGLE 5
#define Q_RADIUSRATIO 6
#define Q_MAXANGLE 7
#define Q_CONDITION 8

//void discreteMeshOptimization(TriMesh& mesh, const int q_crit = Q_MEANRATIO, const float grid_scale = 0.5f, int n_iter = 100);
//void discreteMeshOptimization(PolyMesh& mesh, const int q_crit = Q_MEANRATIO, const float grid_scale = 0.5f, int n_iter = 100);
template <typename T> void discreteMeshOptimization( T& mesh, const int q_crit = Q_MEANRATIO, const float grid_scale = 0.5f, int n_iter = 100 );
template <typename T> inline void discreteMeshOptimizationTargets( T& mesh, float* targets, const float grid_scale = 0.5f, int n_iter = 100 );
//void dmoMeshAdaption(PolyMesh& mesh, const SizeGrid& sizegrid, float* targets, const float grid_scale = 0.5f, int n_iter = 100);
//template <typename T> void dmoMeshAdaption( T& mesh, const BackgroundGrid& sizegrid, float* targets, const float grid_scale = 0.5f, int n_iter = 100 );
//template <typename T> inline void discreteMeshOptimizationDensity( T& mesh, const BackgroundGrid& sizegrid, const float grid_scale = 0.5f, int n_iter = 100 );

extern void discreteMeshOptimizationCPU(TriMesh& mesh, const int q_crit = Q_MEANRATIO, const float grid_scale = 0.5f, int n_iter = 100);
// -------------- //
////////////////////

namespace MeshFunctions {
	using namespace HelperFunctions;

	/////////////////////////////////////////////////
	//// Input/Output Handling | Mesh Conversion ////
	//---------------------------------------------//

	// Write a mesh to ADCIRC-format
	inline void exportMeshAdcirc(const std::string& folder, const std::string& filename, TriMesh& mesh, 
				std::vector<int>& neumannBoundaryTypes, std::vector<int>& dirichletBoundaryTypes) {

		OpenMesh::EPropHandleT<bool> isNeumann;
		mesh.get_property_handle(isNeumann, "isNeumann");

		OpenMesh::EPropHandleT<int> boundaryID;
		mesh.get_property_handle(boundaryID, "boundaryID");

		std::ofstream ofs(folder + filename);

		ofs << filename << "\t! This file was created by HPMeshGen" << std::endl;
		ofs << mesh.n_faces() << " " << mesh.n_vertices() << "\t! NE,NP - NUMBER OF ELEMENTS AND NUMBER OF NODAL POINTS" << std::endl;

		// write out vertices
		for (unsigned int i = 0; i < mesh.n_vertices(); ++i) {
			TriMesh::Point p = mesh.point(mesh.vertex_handle(i));
			ofs << i + 1 << "   " << p[0] << "    " << p[1] << "    " << p[2] << std::endl;
		}
		// write out faces
		for (unsigned int i = 0; i < mesh.n_faces(); ++i) {
			TriMesh::FaceHandle fh = mesh.face_handle(i);
			ofs << i + 1 << " " << 3;
			for (TriMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); ++fv_it) {
				ofs << " " << fv_it->idx() + 1;
			}
			ofs << std::endl;
		}

		if (neumannBoundaryTypes.size() + dirichletBoundaryTypes.size() == 0) {
			ofs.close();
			return;
		}

		// boundaries
		std::vector<TriMesh::EdgeHandle> ehVec;
		for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
			if (mesh.is_boundary(*e_it)) {
				ehVec.push_back(*e_it);
			}
		}
		std::vector<bool> isStored(mesh.n_edges(), false);					// bool checks if edge is stored in vector (avoid boundaries to appear several times)

		std::vector<std::vector<TriMesh::VertexHandle>> neumannBoundaries(neumannBoundaryTypes.size());
		std::vector<std::vector<TriMesh::VertexHandle>> dirichletBoundaries(dirichletBoundaryTypes.size());

		// collect boundaries
		for (size_t i = 0; i < ehVec.size(); ++i) {
			
			if (!isStored[ehVec[i].idx()]) {
				// this is an edge which is not stored yet
				bool boundaryType = mesh.property(isNeumann, ehVec[i]);
				
				TriMesh::HalfedgeHandle heh = mesh.halfedge_handle(ehVec[i], 0);
				if (!mesh.is_boundary(heh)) {
					heh = mesh.opposite_halfedge_handle(heh);
				}
				TriMesh::HalfedgeHandle hehInit = heh;

				// get next heh
				TriMesh::HalfedgeHandle hehNext = mesh.next_halfedge_handle(heh);

				// check if the next heh is also Neumann. If not, the first edge of the boundary was found
				// additionally prevent endless loop in case the boundary is a closed cycle
				while(mesh.property(isNeumann, mesh.edge_handle(hehNext)) == boundaryType && hehNext != hehInit) {
					heh = hehNext;
					// get next heh
					hehNext = mesh.next_halfedge_handle(heh);
				}

				// heh is now the first halfedge_handle of the boundary
				hehInit = heh;

				std::vector<TriMesh::VertexHandle> boundary;
				
				// walk backwards along the boundary to collect all vertices
				boundary.push_back(mesh.to_vertex_handle(heh));
				boundary.push_back(mesh.from_vertex_handle(heh));
				isStored[mesh.edge_handle(heh).idx()] = true;
				// get prev heh
				heh = mesh.prev_halfedge_handle(heh);
				while (mesh.property(isNeumann, mesh.edge_handle(heh)) == boundaryType && heh != hehInit) {
					boundary.push_back(mesh.from_vertex_handle(heh));
					isStored[mesh.edge_handle(heh).idx()] = true;
					// get prev heh
					heh = mesh.prev_halfedge_handle(heh);
				}

				if (heh == hehInit) {
					boundary.pop_back();
				}

				if(boundaryType)
					neumannBoundaries[mesh.property(boundaryID, mesh.edge_handle(hehInit))] = boundary;
				else
					dirichletBoundaries[mesh.property(boundaryID, mesh.edge_handle(hehInit))] = boundary;

			}
		}

		// write boundaries to file
		// Neumann
		ofs << neumannBoundaries.size() << "\t! NOPE - TOTAL NUMBER OF OPEN BOUNDARY FORCING SEGMENTS" << std::endl;
		size_t nNeumannVertices = 0;
		for (size_t i = 0; i < neumannBoundaries.size(); ++i) {
			nNeumannVertices += neumannBoundaries[i].size();
		}
		ofs << nNeumannVertices << "\t! NOPE - TOTAL NUMBER OF OPEN BOUNDARY NODES" << std::endl;
		for (size_t i = 0; i < neumannBoundaries.size(); ++i) {
			ofs << neumannBoundaries[i].size(); 
			if (neumannBoundaryTypes[i] >= 0) ofs << " " << neumannBoundaryTypes[i];
			ofs << "\t! NETA - NUMBER OF NODES ON OPEN BOUNDERY FORCING SEGMENT NO." << (i + 1) << std::endl;

			for (size_t j = 0; j < neumannBoundaries[i].size(); ++j) {
				ofs << neumannBoundaries[i][j].idx() + 1;
				if (j == 0) ofs << "\t!BEGINING NODE NUMBER ON OPEN BOUNDARY NO. " << (i + 1);
				else if(j == neumannBoundaries[i].size() - 1) ofs << "\t!ENDING NODE NUMBER ON OPEN BOUNDARY NO. " << (i + 1);
				ofs << std::endl;
			}
		}

		// Dirichlet
		ofs << dirichletBoundaries.size() << "\t! NOPE - NBOU - TOTAL NUMBER OF LAND BOUNDARY SEGMENTS INCLUDING ISLANDS SEGMENTS" << std::endl;
		size_t nDirichletVertices = 0;
		for (size_t i = 0; i < dirichletBoundaries.size(); ++i) {
			nDirichletVertices += dirichletBoundaries[i].size();
		}
		ofs << nDirichletVertices << "\t! NOPE - TOTAL NUMBER OF LAND BOUNDARY NODES" << std::endl;
		for (size_t i = 0; i < dirichletBoundaries.size(); ++i) {
			ofs << dirichletBoundaries[i].size(); 
			if (dirichletBoundaryTypes[i] >= 0) ofs << " " << dirichletBoundaryTypes[i];
			ofs << "\t! NETA - NUMBER OF NODES ON LAND BOUNDERY SEGMENT NO." << (i + 1) << " AND BOUNDARY TYPE" << std::endl;

			for (size_t j = 0; j < dirichletBoundaries[i].size(); ++j) {
				ofs << dirichletBoundaries[i][j].idx() + 1;
				if (j == 0) ofs << "\t!BEGINING NODE NUMBER ON LAND BOUNDARY NO. " << (i + 1);
				else if (j == dirichletBoundaries[i].size() - 1) ofs << "\t!ENDING NODE NUMBER ON LAND BOUNDARY NO. " << (i + 1);
				ofs << std::endl;
			}
		}

		ofs.close();
	}
	
	// Write mesh in an OpenMesh-format
	inline void printMesh(const std::experimental::filesystem::path& filename, PolyMesh& mesh) {
		LOG( INFO ) << "Write '" << filename.string() << "'";
		OpenMesh::IO::write_mesh(mesh, filename.string());
	}
	inline void printMesh(const std::experimental::filesystem::path& filename, TriMesh& mesh) {
		LOG( INFO ) << "Write '" << filename.string() << "'";
		OpenMesh::IO::write_mesh(mesh, filename.string());
	}
	inline void printMeshCartesian( const std::experimental::filesystem::path& filename, TriMesh& mesh, OceanMesh& oceanMesh ) {
		TriMesh m = mesh;
		for( const auto& v : m.vertices() ) {
			auto p = m.point( v );
			p = oceanMesh.sphericalToCartesian( p );
			m.set_point( v, p );
		}
		LOG( INFO ) << "Write '" << filename.string() << "'";
		OpenMesh::IO::write_mesh( m, filename.string() );
	}
	
	// Convert a triangular mesh to a poly mesh
	inline PolyMesh tri2poly(TriMesh& mesh) {
		PolyMesh pmesh;
		OpenMesh::EPropHandleT<bool> polyIsNeumann;
		OpenMesh::EPropHandleT<int> polyBoundaryID;
		OpenMesh::VPropHandleT<int> polyContourIdx;
		pmesh.add_property( polyIsNeumann, "polyIsNeumann" );
		pmesh.add_property( polyBoundaryID, "polyBoundaryID" );
		pmesh.add_property( polyContourIdx, "polyContourIdx" );

		OpenMesh::EPropHandleT<bool> isNeumann;
		mesh.get_property_handle(isNeumann, "isNeumann");

		OpenMesh::EPropHandleT<int> boundaryID;
		mesh.get_property_handle(boundaryID, "boundaryID");

		OpenMesh::VPropHandleT<int> contourIdx;
		mesh.get_property_handle(contourIdx, "contourIdx");

		std::vector<PolyMesh::VertexHandle> vhandle_poly(mesh.n_vertices());
		for (unsigned int i = 0; i < vhandle_poly.size(); ++i) {
			vhandle_poly[i] = pmesh.add_vertex(mesh.point(mesh.vertex_handle(i)));
		}

		for (unsigned int i = 0; i < mesh.n_faces(); ++i) {
			std::vector<PolyMesh::VertexHandle>  face_vhandles;
			auto fh = mesh.face_handle(i);
			for (TriMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); ++fv_it) {
				face_vhandles.push_back(vhandle_poly[fv_it->idx()]);
			}

			pmesh.add_face(face_vhandles);
		}

		// transport information about boundaries
		for (uint i = 0; i < mesh.n_edges(); ++i) {
			TriMesh::EdgeHandle eh = mesh.edge_handle(i);
			
			if (!mesh.is_boundary(eh))
				continue;

			// find corresponding poly-edge
			/* Note:
				OpenMesh keeps the indices of the vertices but not the indices of the edges!
				That is why the search for the corresponding edge in the polygonal mesh is necessary.
			*/
			PolyMesh::EdgeHandle peh; //= pmesh.edge_handle(i);
			TriMesh::VertexHandle vh1 = mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0));
			TriMesh::VertexHandle vh2 = mesh.from_vertex_handle(mesh.halfedge_handle(eh, 1));

			PolyMesh::VertexHandle pvh1 = pmesh.vertex_handle(vh1.idx());
			for (auto voh_it = pmesh.voh_iter(pvh1); voh_it.is_valid(); ++voh_it) {
				if (pmesh.to_vertex_handle(*voh_it).idx() == vh2.idx()) {
					peh = pmesh.edge_handle(*voh_it);
					break;
				}
			}

			pmesh.property(polyIsNeumann , peh) = mesh.property(isNeumann, eh);
			pmesh.property(polyBoundaryID, peh) = mesh.property(boundaryID, eh);
			pmesh.property(polyContourIdx, vh1) = mesh.property(contourIdx, vh1);
			pmesh.property(polyContourIdx, vh2) = mesh.property(contourIdx, vh2);
		}

		return pmesh;
	}
	
	inline void printMeshSVG(const std::string& filename, TriMesh& meshInput) {
		std::ofstream ofs(filename);

		TriMesh mesh = meshInput;

		// get min/max
		float x_min = std::numeric_limits<float>::max();
		float y_min = std::numeric_limits<float>::max();
		float x_max = -std::numeric_limits<float>::max();
		float y_max = -std::numeric_limits<float>::max();

		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
			TriMesh::Point p = mesh.point(*v_it);

			x_min = std::min(x_min, p[0]);
			x_max = std::max(x_max, p[0]);
			y_min = std::min(y_min, p[1]);
			y_max = std::max(y_max, p[1]);
		}

		//float x_dim = x_max - x_min;
		//float y_dim = y_max - y_min;

		ofs << "<?xml version=\"1.0\" encoding=\"UTF - 8\"?>" << std::endl;

		ofs << "<svg xmlns=\"http://www.w3.org/2000/svg\"" << std::endl;

		ofs << "version=\"1.1\" baseProfile=\"full\"" << std::endl;
		ofs << "width=\"700px\" height=\"400px\" viewBox=\"0 0 700 400\">" << std::endl;
		ofs << "" << std::endl;

		//ofs << "" << std::endl;

		ofs << "</svg>" << std::endl;

		ofs.close();
	}

	//---------------------------------------------//
	//// Input/Output Handling | Mesh Conversion ////
	/////////////////////////////////////////////////

	//////////////////////////////////
	//// Mesh Reduction Functions ////
	//------------------------------//

	// get one ring. Works also for boundary points
	inline void getOneRing(TriMesh& mesh, const TriMesh::VertexHandle vh, std::vector<TriMesh::VertexHandle> &one_ring) {
		assert(one_ring.size() == 0);

		TriMesh::HalfedgeHandle heh = mesh.halfedge_handle(vh);

		if (mesh.is_boundary(vh)) {
			while (!mesh.is_boundary(heh)) {
				heh = mesh.prev_halfedge_handle(heh);
				heh = mesh.opposite_halfedge_handle(heh);
			}
		}

		TriMesh::HalfedgeHandle heh_init = heh;

		do {
			one_ring.push_back(mesh.to_vertex_handle(heh));
			heh = mesh.opposite_halfedge_handle(heh);
			heh = mesh.next_halfedge_handle(heh);
		} while (heh != heh_init && !mesh.is_boundary(mesh.opposite_halfedge_handle(heh)));

		one_ring.push_back(mesh.to_vertex_handle(heh));
	}
	// Check if the new point would flip a triangle or create thin triangles
	inline bool isPointAllowed(TriMesh& mesh, const TriMesh::VertexHandle& vh1, const TriMesh::VertexHandle& vh2, const TriMesh::Point& p_new) {
		const float min_q = 0.5f;

		std::vector<TriMesh::VertexHandle> vh1fan;
		getOneRing(mesh, vh1, vh1fan);

		for (size_t i = 0; i < vh1fan.size() - 1; ++i) {
			size_t j = i + 1;
			if (vh1fan[i] == vh2 || vh1fan[j] == vh2) {
				continue;
			}
			TriMesh::Point p1 = mesh.point(vh1fan[i]);
			TriMesh::Point p2 = mesh.point(vh1fan[j]);

			const float quality = QualityMetrics::meanRatioMetric({ p_new, p2, p1 });

			if (quality < 0.5f)
				return false;

		}

		std::vector<TriMesh::VertexHandle> vh2fan;
		getOneRing(mesh, vh2, vh2fan);

		for (size_t i = 0; i < vh2fan.size() - 1; ++i) {
			size_t j = i + 1;
			if (vh2fan[i] == vh1 || vh2fan[j] == vh1) {
				continue;
			}
			TriMesh::Point p1 = mesh.point(vh2fan[i]);
			TriMesh::Point p2 = mesh.point(vh2fan[j]);

			const float quality = QualityMetrics::meanRatioMetric({ p_new, p2, p1 });

			if (quality < 0.5f)
				return false;

		}

		return true;
	}

	// create Delauney triangulation
	inline void edgeFlip(TriMesh& mesh) {
		//std::cout << "Edge flip ...";
		while (true) {
			bool nothingWasFlipped = true;
			for (unsigned int edgeID = 0; edgeID < mesh.n_edges(); ++edgeID) {
				const TriMesh::EdgeHandle eh = mesh.edge_handle(edgeID);

				if (mesh.is_boundary(eh)) {
					continue;
				}
				const TriMesh::HalfedgeHandle heh1 = mesh.halfedge_handle(eh, 0);
				const TriMesh::HalfedgeHandle heh2 = mesh.halfedge_handle(eh, 1);

				// get all 4 points
				const TriMesh::Point p1 = mesh.point(mesh.from_vertex_handle(heh1));
				const TriMesh::Point p2 = mesh.point(mesh.from_vertex_handle(heh2));
				const TriMesh::Point p3 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh1)));
				const TriMesh::Point p4 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh2)));

				// make sure edge flip is possible
				if (isInsideTriangle(p1, p2, p3, p4) || isInsideTriangle(p2, p3, p4, p1) ||
					isInsideTriangle(p3, p4, p1, p2) || isInsideTriangle(p4, p1, p2, p3) || !mesh.is_flip_ok(eh)) {
					continue;
				}
				if (calcAngle(p1, p3, p2) + calcAngle(p3, p2, p4) + calcAngle(p2, p4, p1) + calcAngle(p4, p1, p3) < 359.9999) {
					continue;
				}

				// calculate minimal angles
				// initial triangles: t1 = (1,2,3) ; t2 = (1,2,4)
				const float alpha_t1 = calcAngle(p1, p2, p3);
				const float beta_t1 = calcAngle(p2, p3, p1);
				const float gamma_t1 = calcAngle(p3, p1, p2);
				const float alpha_t2 = calcAngle(p1, p2, p4);
				const float beta_t2 = calcAngle(p2, p4, p1);
				const float gamma_t2 = calcAngle(p4, p1, p2);

				float min_initial = std::min(alpha_t1, beta_t1);
				min_initial = std::min(min_initial, gamma_t1);
				min_initial = std::min(min_initial, alpha_t2);
				min_initial = std::min(min_initial, beta_t2);
				min_initial = std::min(min_initial, gamma_t2);

				const float alpha_t3 = calcAngle(p3, p4, p1);
				const float beta_t3 = calcAngle(p4, p1, p3);
				const float gamma_t3 = calcAngle(p1, p3, p4);
				const float alpha_t4 = calcAngle(p3, p4, p2);
				const float beta_t4 = calcAngle(p4, p2, p3);
				const float gamma_t4 = calcAngle(p2, p3, p4);
				float min_flipped = std::min(alpha_t3, beta_t3);
				min_flipped = std::min(min_flipped, gamma_t3);
				min_flipped = std::min(min_flipped, alpha_t4);
				min_flipped = std::min(min_flipped, beta_t4);
				min_flipped = std::min(min_flipped, gamma_t4);

				if ((min_flipped > min_initial)) {
					mesh.flip(eh);
					nothingWasFlipped = false;
				}

			}
			if (nothingWasFlipped)
				break;
		}
	}

	inline void edgeFlip(TriMesh& mesh, std::set<uint>& markedEdges, const std::vector<TriMesh::HalfedgeHandle>& hehVec, OpenMesh::EPropHandleT<uint>& origID) {
		//std::cout << "Edge flip ...";

		std::set<uint> markedEdgesBuf = markedEdges;

		while (true) {
			bool nothingWasFlipped = true;

			for (auto it = markedEdges.begin(); it != markedEdges.end(); ++it) {

				const TriMesh::HalfedgeHandle heh1 = hehVec[*it];
				const TriMesh::HalfedgeHandle heh2 = mesh.opposite_halfedge_handle(heh1);

				if (heh1.idx() < 0) {
					std::cout << "Error Found!" << std::endl;
					OpenMesh::IO::write_mesh(mesh, "../../HPMeshGenOutput/debugEdgeFlip.off");
				}

				/*
					Vertex positions:
					3 ---- 1
					|    / |
					| /    |
					2 ---- 4

					Diagonal edge is the one that should be flipped (heh1 & heh2)

					Helfedge hehXY points from vertex X to vertex Y
				*/

				const TriMesh::HalfedgeHandle heh13 = mesh.next_halfedge_handle(heh1);
				const TriMesh::HalfedgeHandle heh32 = mesh.next_halfedge_handle(heh13);
				const TriMesh::HalfedgeHandle heh24 = mesh.next_halfedge_handle(heh2);
				const TriMesh::HalfedgeHandle heh41 = mesh.next_halfedge_handle(heh24);

				const TriMesh::EdgeHandle eh = mesh.edge_handle(heh1);
				const uint edgeID = eh.idx();

				if (mesh.is_boundary(eh)) {
					continue;
				}

				// get all 4 points
				const TriMesh::Point p1 = mesh.point(mesh.from_vertex_handle(heh1));
				const TriMesh::Point p2 = mesh.point(mesh.from_vertex_handle(heh2));
				const TriMesh::Point p3 = mesh.point(mesh.to_vertex_handle(heh13));
				const TriMesh::Point p4 = mesh.point(mesh.to_vertex_handle(heh24));

				const float beta_t1 = calcAngle(p2, p3, p1);
				const float beta_t2 = calcAngle(p2, p4, p1);
				const float beta_t3 = calcAngle(p4, p1, p3);
				const float beta_t4 = calcAngle(p4, p2, p3);

				// make sure edge flip is possible
				if (isInsideTriangle(p1, p2, p3, p4) || isInsideTriangle(p2, p3, p4, p1) ||
					isInsideTriangle(p3, p4, p1, p2) || isInsideTriangle(p4, p1, p2, p3) || !mesh.is_flip_ok(eh)) {
					continue;
				}
				if (beta_t1 + beta_t2 + beta_t3 + beta_t4 < 359.9999) {
					continue;
				}

				// calculate minimal angles
				const float alpha_t1 = calcAngle(p1, p2, p3);
				const float gamma_t1 = calcAngle(p3, p1, p2);
				const float alpha_t2 = calcAngle(p1, p2, p4);
				const float gamma_t2 = calcAngle(p4, p1, p2);

				float min_initial = std::min(alpha_t1, beta_t1);
				min_initial = std::min(min_initial, gamma_t1);
				min_initial = std::min(min_initial, alpha_t2);
				min_initial = std::min(min_initial, beta_t2);
				min_initial = std::min(min_initial, gamma_t2);

				// flipped triangles: t3 = (3,4,1) ; t4 = (3,4,2)
				const float alpha_t3 = calcAngle(p3, p4, p1);
				const float gamma_t3 = calcAngle(p1, p3, p4);
				const float alpha_t4 = calcAngle(p3, p4, p2);
				const float gamma_t4 = calcAngle(p2, p3, p4);
				float min_flipped = std::min(alpha_t3, beta_t3);
				min_flipped = std::min(min_flipped, gamma_t3);
				min_flipped = std::min(min_flipped, alpha_t4);
				min_flipped = std::min(min_flipped, beta_t4);
				min_flipped = std::min(min_flipped, gamma_t4);

				if ((min_flipped > min_initial)) {
					
					markedEdges.insert(mesh.property(origID, eh));
					
					// mark surrounding edges
					markedEdges.insert(mesh.property(origID, mesh.edge_handle(heh13)));
					markedEdges.insert(mesh.property(origID, mesh.edge_handle(heh32)));
					markedEdges.insert(mesh.property(origID, mesh.edge_handle(heh24)));
					markedEdges.insert(mesh.property(origID, mesh.edge_handle(heh41)));

					mesh.flip(eh);

					nothingWasFlipped = false;
				}
			}

			std::swap(markedEdges, markedEdgesBuf);

			if (nothingWasFlipped)
				break;
		}
		//std::cout << "\b\b\bdone" << std::endl;
	}

	//------------------------------//
	//// Mesh Reduction Functions ////
	//////////////////////////////////

	///////////////////
	//// Remeshing ////
	//---------------//


	// smooth boundary. Vertices are not bound to initial vertex positions.
	inline void boundarySmoothing(TriMesh& mesh, const Contour::Shape& contour, float featureAngle = 100) {

		OpenMesh::VPropHandleT<float> feature_size;
		mesh.get_property_handle( feature_size, "feature_size" );
		
		for( const auto& vh : mesh.vertices() ) {
			const auto p = mesh.point( vh );

			if( !vh.is_boundary() ) 
				continue;
			if( mesh.data( vh ).is_feature )
				continue;

			if( feature_size.is_valid() && mesh.property( feature_size, vh ) < featureAngle )
				continue;

			std::vector<TriMesh::VertexHandle> oneRing;
			getOneRing( mesh, vh, oneRing );
			std::vector<TriMesh::Point> oneRingP;
			std::transform( oneRing.rbegin(), oneRing.rend(), std::back_inserter(oneRingP), [&mesh]( TriMesh::VertexHandle vh ) { return mesh.point( vh ); } );

			float qOld = QualityMetrics::meanRatioMetric( oneRingP, p );

			const auto heh = vh.out(); // boundary halfedge

			auto vhLeft = heh.prev().from();
			auto vhRight = heh.to();

			TriMesh::Point pLeft = mesh.point( vhLeft );
			TriMesh::Point pRight = mesh.point( vhRight );


			auto mid = contour.computeMid( pLeft, pRight, p );


			// compute mean ratio metric for new point
			float qNew = QualityMetrics::meanRatioMetric( oneRingP, mid );
			if( qNew < 0.3 && qNew < qOld ) {
				mid = 0.5 * ( mid + p );
				qNew = QualityMetrics::meanRatioMetric( oneRingP, mid );
				if( qNew < 0.3 && qNew < qOld )
					continue;
			}

			mesh.set_point( vh, mid );

		}
	}


	inline void tri2quad(PolyMesh& pmesh) {

		OpenMesh::EPropHandleT<bool> polyIsNeumann;
		pmesh.get_property_handle(polyIsNeumann, "polyIsNeumann");
		
		OpenMesh::EPropHandleT<int> polyBoundaryID;
		pmesh.get_property_handle(polyBoundaryID, "polyBoundaryID");

		// map mit weight & edge
		std::map<float, PolyMesh::EdgeHandle> prioQueue;
		std::vector<std::vector<PolyMesh::VertexHandle>> quadVertices(pmesh.n_edges());
		// iteriere ueber edges
		// fasse beide dreiecke zu viereck zusammen und berechne min/max angle criterion
		// --> schreibe weight in map

		// store boundary properties
		std::vector<bool> isNeumannBuf(pmesh.n_vertices(), false);
		std::vector<int> boundaryIdBuf(pmesh.n_vertices(), -INT_MAX);
		for (auto e_it = pmesh.edges_begin(); e_it != pmesh.edges_end(); ++e_it) {
			if (!pmesh.is_boundary(*e_it))
				continue;
			TriMesh::HalfedgeHandle heh = pmesh.halfedge_handle(*e_it, 0);
			if (!pmesh.is_boundary(heh))
				heh = pmesh.opposite_halfedge_handle(heh);

			TriMesh::VertexHandle vh = pmesh.from_vertex_handle(heh);

			isNeumannBuf[vh.idx()] = pmesh.property(polyIsNeumann, *e_it);
			boundaryIdBuf[vh.idx()] = pmesh.property(polyBoundaryID, *e_it);
		}

		for (uint edgeID = 0; edgeID < pmesh.n_edges(); ++edgeID) {
			PolyMesh::EdgeHandle eh = pmesh.edge_handle(edgeID);

			if (pmesh.is_boundary(eh)) {
				continue;
			}

			PolyMesh::HalfedgeHandle heh1 = pmesh.halfedge_handle(eh, 0);
			PolyMesh::HalfedgeHandle heh2 = pmesh.halfedge_handle(eh, 1);
			
			std::vector<PolyMesh::VertexHandle>vVec(4);
			std::vector<PolyMesh::Point>pVec(4);

			/*
				1-----0
				|   / |
				| /   |
				2-----3
			*/

			vVec[0] = pmesh.to_vertex_handle(heh1);
			vVec[1] = pmesh.to_vertex_handle(pmesh.next_halfedge_handle(heh1));
			vVec[2] = pmesh.to_vertex_handle(heh2);
			vVec[3] = pmesh.to_vertex_handle(pmesh.next_halfedge_handle(heh2));

			pVec[0] = pmesh.point(vVec[0]);
			pVec[1] = pmesh.point(vVec[1]);
			pVec[2] = pmesh.point(vVec[2]);
			pVec[3] = pmesh.point(vVec[3]);

			// do not delete edge if quad is non-convex
			if (HelperFunctions::calcAngle(pVec[3], pVec[0], pVec[2]) + HelperFunctions::calcAngle(pVec[2], pVec[0], pVec[1]) > 180.) {
				continue;
			}
			if (HelperFunctions::calcAngle(pVec[1], pVec[2], pVec[0]) + HelperFunctions::calcAngle(pVec[0], pVec[2], pVec[3]) > 180.) {
				continue;
			}

			std::vector<float>angleVec(4);
			angleVec[0] = (float)HelperFunctions::calcAngle(pVec[3], pVec[0], pVec[1]);
			angleVec[1] = (float)HelperFunctions::calcAngle(pVec[0], pVec[1], pVec[2]);
			angleVec[2] = (float)HelperFunctions::calcAngle(pVec[1], pVec[2], pVec[3]);
			angleVec[3] = (float)HelperFunctions::calcAngle(pVec[2], pVec[3], pVec[0]);

			float minAngle = std::numeric_limits<float>::max();
			float maxAngle = -std::numeric_limits<float>::max();
			for (size_t i = 0; i < angleVec.size(); ++i) {
				minAngle = std::min(minAngle, angleVec[i]);
				maxAngle = std::max(maxAngle, angleVec[i]);
			}

			// do not add edge if its deletion would generate a degenerate quad
			bool degenerateQuad = false;
			for (size_t i = 0; i < angleVec.size(); ++i) {
				if (angleVec[i] > 165) {
					degenerateQuad = true;
					break;
				}
			}
			if (degenerateQuad)
				continue;
			
			float weight = maxAngle - minAngle;

			prioQueue[weight] = eh;
			quadVertices[edgeID] = vVec;
		}
		

		// delete edges (tracke welche faces zusammengefasst wurden)

		std::vector<bool>faceDeleted(pmesh.n_faces(), false);
		//std::cout << "Number of faces: " << mesh.n_faces() << std::endl;

		pmesh.request_face_status();
		pmesh.request_edge_status();
		pmesh.request_vertex_status();

		std::vector<PolyMesh::EdgeHandle> edgesToDelete;

		for (auto pq_it = prioQueue.begin(); pq_it != prioQueue.end(); ++pq_it) {
			PolyMesh::EdgeHandle eh = pq_it->second;
			//uint edgeID = eh.idx();
			PolyMesh::FaceHandle fh1 = pmesh.face_handle(pmesh.halfedge_handle(eh, 0));
			PolyMesh::FaceHandle fh2 = pmesh.face_handle(pmesh.halfedge_handle(eh, 1));
			
			//std::cout << "fh1 = " << fh1.idx() << "  |  fh2 = " << fh2.idx() << std::endl;

			// check if one of the faces is already part of a quad
			if (faceDeleted[fh1.idx()] || faceDeleted[fh2.idx()])
				continue;

			faceDeleted[fh1.idx()] = true; 
			faceDeleted[fh2.idx()] = true;

			edgesToDelete.push_back(eh);
		}

		for (size_t i = 0; i < edgesToDelete.size(); ++i) {
			pmesh.delete_edge(edgesToDelete[i], false);
			pmesh.add_face(quadVertices[edgesToDelete[i].idx()]);
		}

		pmesh.garbage_collection();

		// restore boundary properties
		for (auto e_it = pmesh.edges_begin(); e_it != pmesh.edges_end(); ++e_it) {
			if (!pmesh.is_boundary(*e_it))
				continue;
			TriMesh::HalfedgeHandle heh = pmesh.halfedge_handle(*e_it, 0);
			if (!pmesh.is_boundary(heh))
				heh = pmesh.opposite_halfedge_handle(heh);

			TriMesh::VertexHandle vh = pmesh.from_vertex_handle(heh);

			pmesh.property(polyIsNeumann, *e_it) = isNeumannBuf[vh.idx()];
			pmesh.property(polyBoundaryID, *e_it) = boundaryIdBuf[vh.idx()];
		}
	}

	inline float getNeighMinSize( const ScalarField::ScalarField& sf, TriMesh& mesh, const TriMesh::EdgeHandle eh ) {
		auto heh = mesh.halfedge_handle( eh, 0 );
		auto vh1 = mesh.from_vertex_handle( heh );
		auto vh2 = mesh.to_vertex_handle( heh );
		auto p1 = mesh.point( vh1 );
		auto p2 = mesh.point( vh2 );

		auto s = sf.getScalar( { p1[0],p1[1] }, { p2[0],p2[1] } );

		if( !mesh.is_boundary( heh ) ) {
			auto vh3 = mesh.to_vertex_handle( mesh.next_halfedge_handle( heh ) );
			auto p3 = mesh.point( vh3 );
			s = std::min( s, sf.getScalar( { p1[0],p1[1] }, { p3[0],p3[1] } ) );
			s = std::min( s, sf.getScalar( { p2[0],p2[1] }, { p3[0],p3[1] } ) );
		}
		if( !mesh.is_boundary( mesh.opposite_halfedge_handle( heh ) ) ) {
			auto vh4 = mesh.to_vertex_handle( mesh.next_halfedge_handle( mesh.opposite_halfedge_handle( heh ) ) );
			auto p4 = mesh.point( vh4 );
			s = std::min( s, sf.getScalar( { p1[0],p1[1] }, { p4[0],p4[1] } ) );
			s = std::min( s, sf.getScalar( { p2[0],p2[1] }, { p4[0],p4[1] } ) );
		}

		return s;
	}
	
	//---------------//
	//// Remeshing ////
	///////////////////

	///////////////////////////////////
	//// Mesh Refinement Functions ////
	//-------------------------------//

	// Refine poly mesh by performing one Catmull Clark step
	inline void catmullClark2D(PolyMesh& pmesh) {

		OpenMesh::EPropHandleT<bool> polyIsNeumann;
		pmesh.get_property_handle(polyIsNeumann, "polyIsNeumann");

		OpenMesh::EPropHandleT<int> polyBoundaryID;
		pmesh.get_property_handle(polyBoundaryID, "polyBoundaryID");

		std::vector<PolyMesh::VertexHandle> edge2vertex(pmesh.n_edges());
		std::vector<PolyMesh::VertexHandle> face2vertex(pmesh.n_faces());

		std::vector<bool> isBoundaryVec(pmesh.n_edges());
		std::vector<bool> isNeumannVec(pmesh.n_edges());
		std::vector<int> boundaryIDVec(pmesh.n_edges());

		// add vertices
		// vertex on edge
		for (unsigned int edgeID = 0; edgeID < pmesh.n_edges(); ++edgeID) {
			PolyMesh::EdgeHandle eh = pmesh.edge_handle(edgeID);
			PolyMesh::HalfedgeHandle heh = pmesh.halfedge_handle(eh, 0);
			PolyMesh::Point p1 = pmesh.point(pmesh.from_vertex_handle(heh));
			PolyMesh::Point p2 = pmesh.point(pmesh.to_vertex_handle(heh));
			PolyMesh::Point p_mid;
			p_mid[0] = 0.5f * (p1[0] + p2[0]);
			p_mid[1] = 0.5f * (p1[1] + p2[1]);
			p_mid[2] = 0.5f * (p1[2] + p2[2]);

			edge2vertex[edgeID] = pmesh.add_vertex(p_mid);

			// store boundary information
			if (pmesh.is_boundary(eh)) {
				isNeumannVec[edgeID] = pmesh.property(polyIsNeumann, eh);
				boundaryIDVec[edgeID] = pmesh.property(polyBoundaryID, eh);
			}
		}
		// vertex on face
		for (unsigned int faceID = 0; faceID < pmesh.n_faces(); ++faceID) {
			PolyMesh::FaceHandle fh = pmesh.face_handle(faceID);
			PolyMesh::Point p_mid = { 0.f, 0.f, 0.f };
			size_t vertexCount = 0;
			for (PolyMesh::FaceVertexIter fv_it = pmesh.fv_begin(fh); fv_it != pmesh.fv_end(fh); ++fv_it) {
				PolyMesh::Point p = pmesh.point(*fv_it);
				p_mid[0] += p[0];
				p_mid[1] += p[1];
				p_mid[2] += p[2];
				++vertexCount;
			}
			p_mid[0] = p_mid[0] / (float)vertexCount;
			p_mid[1] = p_mid[1] / (float)vertexCount;
			p_mid[2] = p_mid[2] / (float)vertexCount;

			face2vertex[faceID] = pmesh.add_vertex(p_mid);
		}


		// do refinement
		//	e---b---*
		//  |   |   |
		//  a---c---*
		//  |   |   |
		//  *---*---*
		// e: from_vertex_handle(heh)		edge
		// a: edge2vertex[edgeID]			after edge
		// c: face2vertex[faceID]			center face
		// b: from_vertex_handle(heh-1)		before edge

		std::vector<std::vector<PolyMesh::VertexHandle>> face_vertex_handles;

		for (unsigned int faceID = 0; faceID < face2vertex.size(); ++faceID) {
			PolyMesh::FaceHandle fh = pmesh.face_handle(faceID);

			PolyMesh::HalfedgeHandle heh_init = pmesh.halfedge_handle(fh);
			PolyMesh::HalfedgeHandle heh = heh_init;

			do {
				std::vector<PolyMesh::VertexHandle> fevh;

				PolyMesh::EdgeHandle eh = pmesh.edge_handle(heh);
				PolyMesh::HalfedgeHandle heh_prev = pmesh.prev_halfedge_handle(heh);
				PolyMesh::EdgeHandle eh_prev = pmesh.edge_handle(heh_prev);

				fevh.push_back(pmesh.from_vertex_handle(heh));		// e
				fevh.push_back(edge2vertex[eh.idx()]);				// a
				fevh.push_back(face2vertex[faceID]);				// c
				fevh.push_back(edge2vertex[eh_prev.idx()]);			// b

				face_vertex_handles.push_back(fevh);

				heh = pmesh.next_halfedge_handle(heh);

			} while (heh != heh_init);

		}

		pmesh.request_face_status();
		pmesh.request_edge_status();
		pmesh.request_vertex_status();

		// delete old faces
		for (unsigned i = 0; i < pmesh.n_faces(); ++i) {
			PolyMesh::FaceHandle fh = pmesh.face_handle(i);
			pmesh.delete_face(fh, false);
		}
		pmesh.garbage_collection();

		// create new faces
		for (size_t i = 0; i < face_vertex_handles.size(); ++i) {
			pmesh.add_face(face_vertex_handles[i]);
		}

		// add boundary information to refined mesh
		for (uint i = 0; i < edge2vertex.size(); ++i) {
			if (!pmesh.is_boundary(edge2vertex[i])) {
				continue;
			}

			// i is a boundary vertex
			for (auto voh_it = pmesh.voh_iter(edge2vertex[i]); voh_it.is_valid(); ++voh_it) {
				PolyMesh::EdgeHandle eh = pmesh.edge_handle(*voh_it);

				if (pmesh.is_boundary(eh)) {
					pmesh.property(polyIsNeumann, eh) = isNeumannVec[i];
					pmesh.property(polyBoundaryID, eh) = boundaryIDVec[i];
				}
			}
		}

		/*for (auto eh_it = pmesh.edges_begin(); eh_it != pmesh.edges_end(); ++eh_it) {
			if (!pmesh.is_boundary(*eh_it)) {
				continue;
			}
			TriMesh::HalfedgeHandle heh = pmesh.halfedge_handle(*eh_it, 0);
			std::cout << pmesh.property(isNeumann, *eh_it) << "  |  " << pmesh.from_vertex_handle(heh) << "  |  " << pmesh.to_vertex_handle(heh) << std::endl;
			std::cout << pmesh.point(pmesh.from_vertex_handle(heh))[0] << "  ,  " << pmesh.point(pmesh.from_vertex_handle(heh))[1] << std::endl;
		}*/
	}

	//-------------------------------//
	//// Mesh Refinement Functions ////
	///////////////////////////////////

	template <typename T>
	inline void displayQuality( T& mesh, const int n_cols ) {
		std::vector<float> q_vec( n_cols, 0 );

		float q_min = FLT_MAX;
		// measure quality
		for( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it ) {
			float q = QualityMetrics::meanRatioMetric( mesh, *f_it );
			q_min = fminf( q_min, q );
			int idx = std::min( std::max( 0, int( q * n_cols ) ), n_cols - 1 );
			++q_vec[idx];
		}

		for( size_t i = 0; i < q_vec.size(); ++i ) {
			std::cout << "q:    " << (float)i / (float)n_cols << " - " << (float)( i + 1 ) / (float)n_cols << " = " << q_vec[i] << std::endl;
		}
		std::cout << "q_min = " << q_min << std::endl;
	}

	// rescale the mesh to the given size (minimal value is (0,0) )
	inline void rescaleMesh(PolyMesh& mesh, const float x_size, const float y_size, const float z_size) {
		float x_max = -FLT_MAX;
		float y_max = -FLT_MAX;
		float z_max = -FLT_MAX;
		float x_min = FLT_MAX;
		float y_min = FLT_MAX;
		float z_min = FLT_MAX;

		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
			TriMesh::Point p = mesh.point(*v_it);
			x_max = std::max(x_max, p[0]);
			y_max = std::max(y_max, p[1]);
			z_max = std::max(z_max, p[2]);
			x_min = std::min(x_min, p[0]);
			y_min = std::min(y_min, p[1]);
			z_min = std::min(z_min, p[2]);
		}

		float x_length = x_max - x_min;
		float y_length = y_max - y_min;
		float z_length = z_max - z_min;

		float x_factor = x_size / x_length;
		float y_factor = y_size / y_length;
		float z_factor = z_size / z_length;

		float factor = std::min(x_factor, y_factor);
		factor = std::min(factor, z_factor);

		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
			TriMesh::Point p = mesh.point(*v_it);
			p[0] = (p[0] - x_min) * factor;
			p[1] = (p[1] - y_min) * factor;
			p[2] = (p[2] - z_min) * factor;
			mesh.set_point(*v_it, p);
		}
	}

	// print max and min edge length
	inline void printMaxMinEdgeLength(PolyMesh& mesh) {
		float max_edge = -FLT_MAX;
		float min_edge = FLT_MAX;

		for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
			float l = mesh.calc_edge_length(*e_it);
			max_edge = fmaxf(max_edge, l);
			min_edge = fminf(min_edge, l);
		}

		std::cout << "e_max = " << max_edge << "  |  e_min = " << min_edge << std::endl;
	}

	////////////////////////////
	//// Quad Mesh Topology ////
	//------------------------//
	inline void edgeSwapNext(PolyMesh& mesh, const PolyMesh::EdgeHandle& eh ) {
		LOG_ASSERT( !mesh.is_boundary( eh ) );
		auto heh = mesh.halfedge_handle( eh, 0 );
		LOG_ASSERT( mesh.valence( mesh.from_vertex_handle( heh ) ) > 2 );
		LOG_ASSERT( mesh.valence( mesh.to_vertex_handle( heh ) ) > 2 );
		
		//	1 --- 0 --- 5
		//	|     ^     |
		//	|     |     |
		//	2 --- 3 --- 4

		auto hehInit = heh;
		std::vector<PolyMesh::VertexHandle> vhs;
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 0
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 1
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 2
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 3
		heh = mesh.opposite_halfedge_handle( hehInit );
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 4
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 5
		heh = mesh.next_halfedge_handle( heh );

		mesh.delete_edge( eh, false );
		auto newFh1 = mesh.add_face( { vhs[4], vhs[5], vhs[0], vhs[1] } );
		auto newFh2 = mesh.add_face( { vhs[1], vhs[2], vhs[3], vhs[4] } );

		for( auto heh_2 : mesh.fh_range( newFh1 ) ) {
			auto eh_2 = mesh.edge_handle( heh_2 );
			auto vh = mesh.from_vertex_handle( heh_2 );
		}
		for( auto heh_3 : mesh.fh_range( newFh2 ) ) {
			auto eh_2 = mesh.edge_handle( heh_3 );
			auto vh = mesh.from_vertex_handle( heh_3 );
		}

	}

	inline void edgeSwapPrev( PolyMesh& mesh, const PolyMesh::EdgeHandle& eh ) {
		LOG_ASSERT( !mesh.is_boundary( eh ) );
		auto heh = mesh.halfedge_handle( eh, 0 );
		LOG_ASSERT( mesh.valence( mesh.from_vertex_handle( heh ) ) > 2 );
		LOG_ASSERT( mesh.valence( mesh.to_vertex_handle( heh ) ) > 2 );
		
		//	1 --- 0 --- 5
		//	|     ^     |
		//	|     |     |
		//	2 --- 3 --- 4

		auto hehInit = heh;
		std::vector<PolyMesh::VertexHandle> vhs;
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 0
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 1
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 2
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 3
		heh = mesh.opposite_halfedge_handle( hehInit );
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 4
		heh = mesh.next_halfedge_handle( heh );
		vhs.push_back( mesh.to_vertex_handle( heh ) );	// 5
		heh = mesh.next_halfedge_handle( heh );

		mesh.delete_edge( eh, false );
		auto newFh1 = mesh.add_face( { vhs[5], vhs[0], vhs[1], vhs[2] } );
		auto newFh2 = mesh.add_face( { vhs[2], vhs[3], vhs[4], vhs[5] } );

		for( auto heh_2 : mesh.fh_range( newFh1 ) ) {
			auto eh_2 = mesh.edge_handle( heh_2 );
			auto vh = mesh.from_vertex_handle( heh_2 );
		}
		for( auto heh_2 : mesh.fh_range( newFh2 ) ) {
			auto eh_2 = mesh.edge_handle( heh_2 );
			auto vh = mesh.from_vertex_handle( heh_2 );
		}
	}

	inline void vertexSplit( PolyMesh& mesh, const PolyMesh::HalfedgeHandle& heh ) {
		LOG_ASSERT( !mesh.is_boundary( mesh.edge_handle( heh ) ) );

		auto vh = mesh.from_vertex_handle( heh );
		auto fh1 = mesh.face_handle( heh );
		auto fh2 = mesh.opposite_face_handle( heh );
		auto pNew = 0.5f * ( mesh.point( vh ) + mesh.point( mesh.to_vertex_handle( heh ) ) );

		PolyMesh::VertexHandle vhNew = mesh.add_vertex( pNew );

		std::vector<PolyMesh::VertexHandle> fv1, fv2, fv3;
		for( auto fv : mesh.fv_range( fh1 ) ) {
			fv1.push_back( fv );
		}
		for( auto fv : mesh.fv_range( fh2 ) ) {
			fv2.push_back( fv );
		}
		std::replace( fv1.begin(), fv1.end(), vh, vhNew );
		std::replace( fv2.begin(), fv2.end(), vh, vhNew );

		// build new face
		fv3.push_back( vhNew );
		fv3.push_back( mesh.from_vertex_handle( mesh.prev_halfedge_handle( heh ) ) );
		fv3.push_back( vh );
		fv3.push_back( mesh.to_vertex_handle( mesh.next_halfedge_handle( mesh.opposite_halfedge_handle( heh ) ) ) );

		mesh.delete_face( fh1, false );
		mesh.delete_face( fh2, false );
		mesh.add_face( fv1 );
		mesh.add_face( fv2 );
		mesh.add_face( fv3 );
	}

	inline void diagonalCollapse( PolyMesh& mesh, const PolyMesh::HalfedgeHandle& heh ) {
		LOG_ASSERT( !mesh.is_boundary( mesh.edge_handle( heh ) ) );
		auto fh = mesh.face_handle( heh );
		auto vh1 = mesh.from_vertex_handle( heh );
		auto vh2 = mesh.to_vertex_handle( mesh.next_halfedge_handle( heh ) );

		std::set<PolyMesh::FaceHandle> neighs;
		for( auto fv : mesh.fv_range( fh ) ) {
			for( auto vf : mesh.vf_range( fv ) ) {
				neighs.insert( vf );
			}
		}
		neighs.erase( fh );

		std::vector<std::vector<PolyMesh::VertexHandle>> neighVhs;
		for( auto n : neighs ) {
			std::vector<PolyMesh::VertexHandle> vhs;
			for( auto fv : mesh.fv_range( n ) ) {
				vhs.push_back( fv );
			}
			std::replace( vhs.begin(), vhs.end(), vh1, vh2 );
			neighVhs.push_back( vhs );
		}

		mesh.delete_face( fh, false );
		for( auto n : neighs ) {
			mesh.delete_face( n, false );
		}
		mesh.delete_vertex( vh1, false );
		for( auto vhs : neighVhs ) {
			mesh.add_face( vhs );
		}

	}
}

