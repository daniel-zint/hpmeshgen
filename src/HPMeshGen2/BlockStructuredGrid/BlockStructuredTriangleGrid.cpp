#include "BlockStructuredTriangleGrid.h"

#include <set>

#include "../svg.h"

BlockStructuredTriangleGrid::BlockStructuredTriangleGrid(const TriMesh & pm, const size_t n_refinement_steps, OpenMesh::EPropHandleT<bool>& isNeumann, OpenMesh::EPropHandleT<int>& boundaryID, std::vector<int>& neumannBoundaryTypes, std::vector<int>& dirichletBoundaryTypes)
{
	this->fragmentMesh = pm;
	this->isNeumann = isNeumann;
	this->boundaryID = boundaryID;
	this->neumannBoundaryTypes = neumannBoundaryTypes;
	this->dirichletBoundaryTypes = dirichletBoundaryTypes;

	// create (block structured) triangle mesh
	triMesh = fragmentMesh;

	this->n_refinement_steps_ = n_refinement_steps;

	for (int i = 0; i < n_refinement_steps; ++i) {

		std::vector<TriMesh::VertexHandle> edge2vertex(triMesh.n_edges());

		std::vector<bool> isBoundaryVec(triMesh.n_edges());
		std::vector<bool> isNeumannVec(triMesh.n_edges());
		std::vector<int> boundaryIDVec(triMesh.n_edges());

		// add vertices
		for (unsigned int edgeID = 0; edgeID < triMesh.n_edges(); ++edgeID) {
			TriMesh::EdgeHandle eh = triMesh.edge_handle(edgeID);
			TriMesh::HalfedgeHandle heh = triMesh.halfedge_handle(eh, 0);
			TriMesh::Point p1 = triMesh.point(triMesh.from_vertex_handle(heh));
			TriMesh::Point p2 = triMesh.point(triMesh.to_vertex_handle(heh));
			TriMesh::Point p_mid;
			p_mid[0] = 0.5f * (p1[0] + p2[0]);
			p_mid[1] = 0.5f * (p1[1] + p2[1]);
			p_mid[2] = 0.5f * (p1[2] + p2[2]);

			edge2vertex[edgeID] = triMesh.add_vertex(p_mid);

			// store boundary information
			if (triMesh.is_boundary(eh)) {
				isNeumannVec[edgeID] = triMesh.property(isNeumann, eh);
				boundaryIDVec[edgeID] = triMesh.property(boundaryID, eh);
			}
		}

		std::vector<std::vector<TriMesh::VertexHandle>> face_vertex_handles;

		for (auto f_it = triMesh.faces_begin(); f_it != triMesh.faces_end(); ++f_it) {
			TriMesh::FaceHandle fh = *f_it;
			// iterate over halfedges of face
			for (auto fh_it = triMesh.fh_iter(*f_it); fh_it.is_valid(); ++fh_it) {
				TriMesh::EdgeHandle eh = triMesh.edge_handle(*fh_it);
				TriMesh::EdgeHandle eh_prev = triMesh.edge_handle(triMesh.prev_halfedge_handle(*fh_it));

				std::vector<TriMesh::VertexHandle> fevh(3);
				fevh[0] = triMesh.from_vertex_handle(*fh_it);
				fevh[1] = edge2vertex[eh.idx()];
				fevh[2] = edge2vertex[eh_prev.idx()];

				face_vertex_handles.push_back(fevh);
			}

			// iterate over halfedges of face
			std::vector<TriMesh::VertexHandle> fevh;
			for (auto fh_it = triMesh.fh_iter(*f_it); fh_it.is_valid(); ++fh_it) {
				TriMesh::EdgeHandle eh = triMesh.edge_handle(*fh_it);
				fevh.push_back(edge2vertex[eh.idx()]);
			}
			face_vertex_handles.push_back(fevh);

		}

		triMesh.request_face_status();
		triMesh.request_edge_status();
		triMesh.request_vertex_status();

		// delete old faces
		for (unsigned i = 0; i < triMesh.n_faces(); ++i) {
			PolyMesh::FaceHandle fh = triMesh.face_handle(i);
			triMesh.delete_face(fh, false);
		}
		triMesh.garbage_collection();

		// create new faces
		for (size_t i = 0; i < face_vertex_handles.size(); ++i) {
			triMesh.add_face(face_vertex_handles[i]);
		}


		// add boundary information to refined mesh
		for (uint i = 0; i < edge2vertex.size(); ++i) {
			if (!triMesh.is_boundary(edge2vertex[i])) {
				continue;
			}

			// i is a boundary vertex
			for (auto voh_it = triMesh.voh_iter(edge2vertex[i]); voh_it.is_valid(); ++voh_it) {
				PolyMesh::EdgeHandle eh = triMesh.edge_handle(*voh_it);

				if (triMesh.is_boundary(eh)) {
					triMesh.property(isNeumann, eh) = isNeumannVec[i];
					triMesh.property(boundaryID, eh) = boundaryIDVec[i];
				}
			}
		}
	}

	fineMesh = &triMesh;

}


BlockStructuredTriangleGrid::~BlockStructuredTriangleGrid()
{
}

void BlockStructuredTriangleGrid::exportMeshAdcirc(const std::string & folder, const std::string & filename)
{
	std::ofstream ofs(folder + filename);

	ofs << filename << "\t! This file was created by HPMeshGen" << std::endl;
	ofs << triMesh.n_faces() << " " << triMesh.n_vertices() << "\t! NE,NP - NUMBER OF ELEMENTS AND NUMBER OF NODAL POINTS" << std::endl;

	// write out vertices
	for (unsigned int i = 0; i < triMesh.n_vertices(); ++i) {
		TriMesh::Point p = triMesh.point(triMesh.vertex_handle(i));
		ofs << i + 1 << "   " << p[0] << "    " << p[1] << "    " << p[2] << std::endl;
	}
	// write out faces
	for (unsigned int i = 0; i < triMesh.n_faces(); ++i) {
		TriMesh::FaceHandle fh = triMesh.face_handle(i);
		ofs << i + 1 << " " << 3;
		for (TriMesh::FaceVertexIter fv_it = triMesh.fv_begin(fh); fv_it != triMesh.fv_end(fh); ++fv_it) {
			ofs << " " << fv_it->idx() + 1;
		}
		ofs << std::endl;
	}

	// boundaries
	std::vector<TriMesh::EdgeHandle> ehVec;
	for (auto e_it = fragmentMesh.edges_begin(); e_it != fragmentMesh.edges_end(); ++e_it) {
		if (fragmentMesh.is_boundary(*e_it)) {
			ehVec.push_back(*e_it);
		}
	}
	std::vector<bool> isStored(fragmentMesh.n_edges(), false);					// bool checks if edge is stored in vector (avoid boundaries to appear several times)

	std::vector<std::vector<TriMesh::VertexHandle>> neumannBoundaries(neumannBoundaryTypes.size());
	std::vector<std::vector<TriMesh::VertexHandle>> dirichletBoundaries(dirichletBoundaryTypes.size());

	// collect boundaries
	for (size_t i = 0; i < ehVec.size(); ++i) {

		if (!isStored[ehVec[i].idx()]) {
			// this is an edge which is not stored yet
			bool boundaryType = fragmentMesh.property(isNeumann, ehVec[i]);

			TriMesh::HalfedgeHandle heh = fragmentMesh.halfedge_handle(ehVec[i], 0);
			if (!fragmentMesh.is_boundary(heh)) {
				heh = fragmentMesh.opposite_halfedge_handle(heh);
			}
			TriMesh::HalfedgeHandle hehInit = heh;

			// DEBUG
			//do {
			//	heh = patchMesh.next_halfedge_handle(heh);
			//	std::cout << patchMesh.from_vertex_handle(heh).idx() << ": " << patchMesh.property(isNeumann, patchMesh.edge_handle(heh)) << std::endl;
			//} while (heh != hehInit);

			// get next heh
			TriMesh::HalfedgeHandle hehNext = fragmentMesh.next_halfedge_handle(heh);

			// check if the next heh is also Neumann. If not, the first edge of the boundary was found
			// additionally prevent endless loop in case the boundary is a closed cycle
			while (fragmentMesh.property(isNeumann, fragmentMesh.edge_handle(hehNext)) == boundaryType && hehNext != hehInit) {
				heh = hehNext;
				isStored[fragmentMesh.edge_handle(heh).idx()] = true;
				// get next heh
				hehNext = fragmentMesh.next_halfedge_handle(heh);
			}

			// heh is now the first halfedge_handle of the boundary
			hehInit = heh;

			// first vertex
			TriMesh::VertexHandle vhBegin = fragmentMesh.to_vertex_handle(heh);

			// find last vertex
			hehNext = fragmentMesh.prev_halfedge_handle(heh);
			while (fragmentMesh.property(isNeumann, fragmentMesh.edge_handle(hehNext)) == boundaryType && hehNext != hehInit) {
				heh = hehNext;
				isStored[fragmentMesh.edge_handle(heh).idx()] = true;
				// get prev heh
				hehNext = fragmentMesh.prev_halfedge_handle(heh);
			}

			// last vertex
			TriMesh::VertexHandle vhEnd;
			if (heh != hehInit) {
				vhEnd = fragmentMesh.from_vertex_handle(heh);	// check if this works
			}
			else {
				vhEnd = fragmentMesh.to_vertex_handle(heh);
			}

			// collect vertices on triMesh

			// find halfedge pointing towards vhBegin
			for (auto vih_it = triMesh.vih_iter(vhBegin); vih_it.is_valid(); ++vih_it) {
				if (triMesh.is_boundary(*vih_it)) {
					heh = *vih_it;
					break;
				}
			}


			std::vector<TriMesh::VertexHandle> boundary;
			boundary.push_back(triMesh.to_vertex_handle(heh));
			//while (triMesh.to_vertex_handle(heh) != vhEnd) {
			//	boundary.push_back(triMesh.from_vertex_handle(heh));
			//	heh = triMesh.prev_halfedge_handle(heh);
			//}
			do {
				boundary.push_back(triMesh.from_vertex_handle(heh));
				heh = triMesh.prev_halfedge_handle(heh);
			} while (triMesh.to_vertex_handle(heh) != vhEnd);

			// if last vertex is the same as the first one, delete it
			if (boundary[0] == boundary[boundary.size() - 1])
				boundary.pop_back();

			if (boundaryType)
				neumannBoundaries[fragmentMesh.property(boundaryID, fragmentMesh.edge_handle(hehInit))] = boundary;
			else
				dirichletBoundaries[fragmentMesh.property(boundaryID, fragmentMesh.edge_handle(hehInit))] = boundary;

			//// walk backwards along the boundary to collect all vertices
			//boundary.push_back(patchMesh.to_vertex_handle(heh));
			//boundary.push_back(patchMesh.from_vertex_handle(heh));
			//isStored[patchMesh.edge_handle(heh).idx()] = true;
			//// get prev heh
			//heh = patchMesh.prev_halfedge_handle(heh);
			//while (patchMesh.property(isNeumann, patchMesh.edge_handle(heh)) == boundaryType && heh != hehInit) {
			//	boundary.push_back(patchMesh.from_vertex_handle(heh));
			//	isStored[patchMesh.edge_handle(heh).idx()] = true;
			//	// get prev heh
			//	heh = patchMesh.prev_halfedge_handle(heh);
			//}
			//
			//if (heh == hehInit) {
			//	boundary.pop_back();
			//}
			//
			//if (boundaryType)
			//	neumannBoundaries[patchMesh.property(boundaryID, patchMesh.edge_handle(hehInit))] = boundary;
			//else
			//	dirichletBoundaries[patchMesh.property(boundaryID, patchMesh.edge_handle(hehInit))] = boundary;

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
			else if (j == neumannBoundaries[i].size() - 1) ofs << "\t!ENDING NODE NUMBER ON OPEN BOUNDARY NO. " << (i + 1);
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

void BlockStructuredTriangleGrid::print_svg(std::string filename, const bool printValence)
{
	std::ofstream ofs(filename);
	
	TriMesh mesh = triMesh;

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

	float x_dim = x_max - x_min;
	float y_dim = y_max - y_min;

	int xpx = 1000;
	int ypx = (int)(xpx * y_dim / x_dim);

	// norm all vertices
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		TriMesh::Point p = mesh.point(*v_it);
		p[0] = 15 + xpx * (p[0] - x_min) / x_dim;
		p[1] = 15 + ypx - ypx * (p[1] - y_min) / y_dim;
		p[2] = 0;

		mesh.set_point(*v_it, p);
	}

	// header
	ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
	ofs << "<svg xmlns=\"http://www.w3.org/2000/svg\"" << std::endl;
	ofs << "version=\"1.1\" baseProfile=\"full\"" << std::endl;
	ofs << "width=\"" << xpx + 30 << "px\" height=\"" << ypx + 30 << "px\" viewBox=\"0 0 " << xpx + 30 << " " << ypx + 30 << "\">" << std::endl;
	ofs << "" << std::endl;
	
	// body

	std::set<TriMesh::EdgeHandle> eh_block;
	// find all edges that belong to block boundaries
	for (auto v_it = fragmentMesh.vertices_begin(); v_it != fragmentMesh.vertices_end(); ++v_it) {
		TriMesh::VertexHandle vh = mesh.vertex_handle(v_it->idx());
		for (auto voh_it = mesh.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
			// collect all edges until another patch vertex is reached
			eh_block.insert(mesh.edge_handle(*voh_it));

			TriMesh::HalfedgeHandle heh = *voh_it;
			TriMesh::HalfedgeHandle heh_init = heh;
			if (mesh.is_boundary(*voh_it)) {
				while (mesh.to_vertex_handle(heh).idx() >= fragmentMesh.n_vertices()) {
					heh = mesh.next_halfedge_handle(heh);
					eh_block.insert(mesh.edge_handle(heh));
				}
			}
			else {
				while (mesh.to_vertex_handle(heh).idx() >= fragmentMesh.n_vertices()) {
					heh = mesh.next_halfedge_handle(heh);
					heh = mesh.opposite_halfedge_handle(heh);
					heh = mesh.next_halfedge_handle(heh);
					heh = mesh.opposite_halfedge_handle(heh);
					heh = mesh.next_halfedge_handle(heh);
					eh_block.insert(mesh.edge_handle(heh));
				}
			}

		}
	}

	for (auto eh_it = mesh.edges_begin(); eh_it != mesh.edges_end(); ++eh_it) {
		if (eh_block.find(*eh_it) != eh_block.end())
			continue;

		TriMesh::HalfedgeHandle heh = mesh.halfedge_handle(*eh_it, 0);
		TriMesh::Point p1 = mesh.point(mesh.from_vertex_handle(heh));
		TriMesh::Point p2 = mesh.point(mesh.to_vertex_handle(heh));

		svg::svgLine(ofs, p1, p2, "black", 1);
	}

	for (auto eh_it = eh_block.begin(); eh_it != eh_block.end(); ++eh_it) {
		TriMesh::HalfedgeHandle heh = mesh.halfedge_handle(*eh_it, 0);
		TriMesh::Point p1 = mesh.point(mesh.from_vertex_handle(heh));
		TriMesh::Point p2 = mesh.point(mesh.to_vertex_handle(heh));

		svg::svgLine(ofs, p1, p2, "red", 2);
	}

	if (printValence) {
		for (auto vh_it = mesh.vertices_begin(); vh_it != mesh.vertices_end(); ++vh_it) {
			TriMesh::Point p = mesh.point(*vh_it);

			int valence;
			if (mesh.is_boundary(*vh_it))
				valence = 4;
			else
				valence = 6;

			int n_neighbors = 0;
			for (auto vv_it = mesh.vv_iter(*vh_it); vv_it.is_valid(); ++vv_it) ++n_neighbors;

			if (n_neighbors > valence)
				svg::svgCircle(ofs, p, 10.f, "blue", "black", 1.f);
			if (n_neighbors < valence)
				svg::svgCircle(ofs, p, 10.f, "orange", "black", 1.f);
		}
	}

	ofs << "</svg>" << std::endl;

	ofs.close();
}
