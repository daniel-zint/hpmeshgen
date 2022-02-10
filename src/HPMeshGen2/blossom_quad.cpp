#include "blossom_quad.h"

#include "MeshSmoothing/QualityMetrics.h"

using namespace std;

PolyMesh blossom_quad::do_blossom_algo(PolyMesh loadedMesh ) {

	// raw triangle mesh
	workMesh = loadedMesh;
	// property for edge in connection graph
	OpenMesh::HPropHandleT<int>	hprop_graph_edge_idx;
	workMesh.add_property( hprop_graph_edge_idx, "hprop_graph_edge_idx" );
	// property, if mesh edge has alredy edge in connection graph
	OpenMesh::HPropHandleT<bool> hprop_has_graph_edge;
	workMesh.add_property( hprop_has_graph_edge, "hprop_has_graph_edge" );
	// property, if a vertex was already used for an external edge
	OpenMesh::VPropHandleT<bool> vprop_has_external_edge;
	workMesh.add_property( vprop_has_external_edge, "vprop_has_external_edge" );
	// property for external edge in connection graph
	OpenMesh::VPropHandleT<int> vprop_external_edge_idx;
	workMesh.add_property( vprop_external_edge_idx, "vprop_external_edge_idx" );

	// connectivity graph 
	con_graph graph;

	/*
	 * Blossom-quad needs an even number of faces to work
	 * if number of faces is uneven, remove one face
	 */
	int node_number = workMesh.n_faces();
	if (node_number % 2 != 0) {
		PolyMesh::HalfedgeHandle largest_border;
		float largest_size = 0.f;
		/* 
		 * go through all border edges and calculate their length
		 * safe the one edge with the largest length
		 */ 
		for( auto heh : workMesh.halfedges() ) {
			if( workMesh.is_boundary( heh ) == true ) {
				float length = workMesh.calc_edge_length( heh );
				if( length > largest_size ) {
					largest_border = heh;
					largest_size = length;
				}
			}
		}
		/*
		 * move end vertex of edge to the mid point and collapse the edge
		 * now number of faces in mesh is even
		 */
		auto vh1 = workMesh.from_vertex_handle( largest_border );
		auto vh2 = workMesh.to_vertex_handle( largest_border );
		auto p1 = workMesh.point( vh1 );
		auto p2 = workMesh.point( vh2 );
		workMesh.set_point( vh2, 0.5 * ( p1 + p2 ) );
		workMesh.request_edge_status();
		workMesh.request_face_status();
		workMesh.request_vertex_status();
		workMesh.collapse(largest_border);
		workMesh.garbage_collection();
		node_number = workMesh.n_faces();
	}

	/*
	 * go through all halfedges of the mesh and create a connectivity graph containing all faces
	 * this graph is needed for the blossom algorithm
	 */
	int edge_number = 0;
	for (auto halfedge_it = workMesh.halfedges_begin(); halfedge_it != workMesh.halfedges_end(); ++halfedge_it) {
		PolyMesh::HalfedgeHandle halfedgeHandle = *halfedge_it;
		PolyMesh::HalfedgeHandle oppHalfedgeHandle = workMesh.opposite_halfedge_handle(halfedgeHandle);
		PolyMesh::FaceHandle fhs = workMesh.face_handle(halfedgeHandle);
		PolyMesh::FaceHandle fhso = workMesh.opposite_face_handle(halfedgeHandle);

		/*
		 * check, if corresponding edge is border edge
		 * check, if corresponding edge is already in connectivity graph
		 */
		if (fhs.idx() != -1 && fhso.idx() != -1 && workMesh.property(hprop_has_graph_edge, halfedgeHandle) == false) {
			// create new graph edge in connectivity graph
			graph_edge internal_edge;
			internal_edge.face0 = fhs;
			internal_edge.face1 = fhso;

			// calculate vertices for possible quad
			auto quadVhs = connect_triangles( halfedgeHandle );
			internal_edge.quad_vertices = quadVhs;
			std::vector<OpenMesh::Vec3f> quadPoints{ workMesh.point( quadVhs[0] ), workMesh.point( quadVhs[1] ), workMesh.point( quadVhs[2] ), workMesh.point( quadVhs[3] ) };

			// calculate cost for graph edge
			float qual = 1.f - QualityMetrics::condition( quadPoints );
			int converted_qual = (int)(qual * 100);
			internal_edge.cost = converted_qual;

			// add edge to connection graph
			graph.edges.push_back(internal_edge);
			// mark mesh edge as connected
			workMesh.property(hprop_graph_edge_idx, halfedgeHandle) = graph.edges_counter;
			workMesh.property(hprop_graph_edge_idx, oppHalfedgeHandle) = graph.edges_counter;
			graph.edges_counter++;
			workMesh.property(hprop_has_graph_edge, halfedgeHandle) = true;
			workMesh.property(hprop_has_graph_edge, oppHalfedgeHandle) = true;

			edge_number++;
		}
		/*
		 * check, if halfedge is boundary
		 * if so, add an external edge to the connectivity graph
		 */
		else if (workMesh.is_boundary(halfedgeHandle)) {
			PolyMesh::VertexHandle fromVertex = workMesh.from_vertex_handle(halfedgeHandle);
			PolyMesh::VertexHandle toVertex = workMesh.to_vertex_handle(halfedgeHandle);

			// add external edge to the from vertex
			if (!workMesh.property(vprop_has_external_edge, fromVertex)) {
				PolyMesh::FaceHandle externalConnectionface;
				int numberOfOutgoingEdges = 0;
				/*
				 * count number of outgoing halfedges of the vertex
				 * find border face, connected with this vertex 
				 */
				for( auto voh : workMesh.voh_range( fromVertex ) ) {
					PolyMesh::HalfedgeHandle  vohOpp = workMesh.opposite_halfedge_handle( voh );
					PolyMesh::FaceHandle fh = workMesh.face_handle( voh );
					PolyMesh::FaceHandle fhOpp = workMesh.opposite_face_handle( voh );

					if( fh.idx() == -1 && fhOpp.idx() != fhso.idx() && fhOpp.idx() != -1 ) {
						externalConnectionface = fhOpp;
					} else if( fhOpp.idx() == -1 && fh.idx() != fhso.idx() && fh.idx() != -1 ) {
						externalConnectionface = fh;
					}

					numberOfOutgoingEdges++;
				}
				// create new external edge for connectiivity graph
				graph_edge external_edge;
				external_edge.face0 = fhso;
				external_edge.face1 = externalConnectionface;
				external_edge.connecting_vertex = fromVertex;
				external_edge.cost = EXTERNAL_EDGE_COST;
			
				/*
				 * based on the number of outgoing half edges, decide, how this external edge is used
				 * with 2 or 3 outgoing halfedges, no external edge is needed
				 */
				if (numberOfOutgoingEdges == 2 || numberOfOutgoingEdges == 3) {
					external_edge.~graph_edge();
				}
				else {
					workMesh.property(vprop_has_external_edge, fromVertex) = true;
					workMesh.property(vprop_external_edge_idx, fromVertex) = graph.external_edges_counter;
					graph.external_edges_counter++;
					// with 4 outgoing halfedges, edge swap is used
					if (numberOfOutgoingEdges == 4) {
						external_edge.use_edge_swap = true;
					}
					// with 5 or more outgoing halfedges, vertex duplication is used
					else if (numberOfOutgoingEdges >= 5) {
						external_edge.use_vertex_duplication = true;	
					}
					graph.external_edges.push_back(external_edge);
					edge_number++;
				}
			}

			// add external edge to the to vertex
			if (workMesh.property(vprop_has_external_edge, toVertex)) {
				PolyMesh::FaceHandle externalConnectionface;
				int numberOfOutgoingEdges = 0;
				/*
				* count number of outgoing halfedges of the vertex
				* find border face, connected with this vertex
				*/
				for (auto vertexEdgeIdx = workMesh.voh_iter(toVertex); vertexEdgeIdx.is_valid(); ++vertexEdgeIdx) {
					PolyMesh::HalfedgeHandle vertexHalfEdge = *vertexEdgeIdx;
					PolyMesh::HalfedgeHandle  vertexOppositeHalfEdge = workMesh.opposite_halfedge_handle(vertexHalfEdge);
					PolyMesh::FaceHandle vertexFace = workMesh.face_handle(vertexHalfEdge);
					PolyMesh::FaceHandle vertexOppositeFace = workMesh.opposite_face_handle(vertexHalfEdge);

					if (vertexFace.idx() == -1 && vertexOppositeFace.idx() != fhso.idx() && vertexOppositeFace.idx() != -1) {
						externalConnectionface = vertexOppositeFace;
					}
					else if (vertexOppositeFace.idx() == -1 && vertexFace.idx() != fhso.idx() && vertexFace.idx() != -1) {
						externalConnectionface = vertexFace;
					}

					numberOfOutgoingEdges++;
				}
				// create new external edge for connectiivity graph
				graph_edge external_edge;
				external_edge.face0 = externalConnectionface;
				external_edge.face1 = fhso;
				external_edge.connecting_vertex = toVertex;
				external_edge.cost = EXTERNAL_EDGE_COST;
				
				/*
				* based on the number of outgoing half edges, decide, how this external edge is used
				* with 2 or 3 outgoing halfedges, no external edge is needed
				*/
				if (numberOfOutgoingEdges == 2 || numberOfOutgoingEdges == 3) {
					external_edge.~graph_edge();
				}
				else {
					workMesh.property(vprop_has_external_edge, toVertex) = true;
					workMesh.property(vprop_external_edge_idx, toVertex) = graph.external_edges_counter;
					graph.external_edges_counter++;
					// with 4 outgoing halfedges, edge swap is used
					if (numberOfOutgoingEdges == 4) {
						external_edge.use_edge_swap = true;
					}
					// with 5 or more outgoing halfedges, vertex duplication is used
					else if (numberOfOutgoingEdges >= 5) {
						external_edge.use_vertex_duplication = true;
					}
					graph.external_edges.push_back(external_edge);
					edge_number++;
				}
			}
		}
	
	}

	/* 
	 * create data structures needed for the blossom algorithm
	 * in the edges array there are the indices of the start and end index of every edge in the connectivity graph, so every index is a face index
	 * in the weights array there are the calculated weights of every edge in the connectivity graph
	 */
	std::vector<int> edges, weights;
	edges.reserve( 2 * edge_number );
	weights.reserve( edge_number );
	for( const auto& e : graph.edges ) {
		weights.push_back( e.cost );
		edges.push_back( e.face0.idx() );
		edges.push_back( e.face1.idx() );
	}
	for( const auto& e : graph.external_edges ) {
		weights.push_back( e.cost );
		edges.push_back( e.face0.idx() );
		edges.push_back( e.face1.idx() );
	}

	/*
	 * beginning of creation of the result mesh
	 * first all points are copied to the new mesh
	 */
	int point_count = 0;
	for( auto vh : workMesh.vertices() ) {
		blossomMesh_step1_vertex_vec.push_back( blossomMesh_step1.add_vertex( workMesh.point( vh ) ) );
		point_count++;
	}

	// use the blossom V algorithm to find a cost perfect minimal matching
	PerfectMatching pm( node_number, edge_number );
	struct PerfectMatching::Options options;	
	options.update_duals_before = true;
	options.fractional_jumpstart = false;
	bool check_perfect_matching = true;
	int e = 0;
	for( e = 0; e < edge_number; e++ )
		pm.AddEdge( edges[2 * e], edges[2 * e + 1], weights[e] );
	pm.options = options;
	pm.Solve();
	if (check_perfect_matching)
	{
		int res = CheckPerfectMatchingOptimality(node_number, edge_number, edges.data(), weights.data(), &pm);
		std::printf("check optimality: res=%d (%s)\n", res, (res == 0) ? "ok" : ((res == 1) ? "error" : "fatal error"));
	}
	double cost = ComputePerfectMatchingCost(node_number, edge_number, edges.data(), weights.data(), &pm);
	std::printf("cost = %.1f\n", cost);
	int firstFaceIdx;
	int secondFaceIdx;

	/*
	 * after a matching is found, use it to create a quad mesh from the connectivity graph
	 * if an edge is used, add its quad face to the match
	 * if an external edge is used, mark it for later usage
	 */
	for (firstFaceIdx = 0; firstFaceIdx < node_number; firstFaceIdx++) {
		secondFaceIdx = pm.GetMatch(firstFaceIdx);
		
		if (firstFaceIdx < secondFaceIdx) {
			PolyMesh::FaceHandle firstFace = workMesh.face_handle(firstFaceIdx);
			for( auto firstFaceHalfedge : workMesh.fh_range( firstFace ) ) {
				graph_edge edgeToUse = graph.edges[workMesh.property( hprop_graph_edge_idx, firstFaceHalfedge )];
				// add face to mesh
				if( ( edgeToUse.face0.idx() == firstFaceIdx && edgeToUse.face1.idx() == secondFaceIdx ) || ( edgeToUse.face0.idx() == secondFaceIdx && edgeToUse.face1.idx() == firstFaceIdx ) ) {
					blossomMesh_step1.add_face( edgeToUse.quad_vertices );
					break;
				}
				// mark external edge for usage
				PolyMesh::VertexHandle fromVertex = workMesh.from_vertex_handle( firstFaceHalfedge );
				if( workMesh.property( vprop_has_external_edge, fromVertex ) ) {
					graph_edge externalEdgeToMark = graph.external_edges[workMesh.property( vprop_external_edge_idx, fromVertex )];
					if( ( externalEdgeToMark.face0.idx() == firstFaceIdx && externalEdgeToMark.face1.idx() == secondFaceIdx ) || ( externalEdgeToMark.face0.idx() == secondFaceIdx && externalEdgeToMark.face1.idx() == firstFaceIdx ) ) {
						graph.external_edges[workMesh.property( vprop_external_edge_idx, fromVertex )].to_be_used = true;
						break;
					}
				}
			}
		}		
	}
	
	// add faces which are connected by external edges
	// TODO TEST THIS!!!!!!!!!!!!!!!
	for( int i = 0; i < graph.external_edges_counter; i++ ) {
		// find edges, that are marked for use of vertex duplication
		if( graph.external_edges[i].to_be_used && !graph.external_edges[i].used ) {
			graph.external_edges[i].used = true;

			//LOG( ERROR ) << "Untested code in blossom quad";

			// get the two faces to connect and the connecting vertex from the external edge
			PolyMesh::FaceHandle face0 = graph.external_edges[i].face0;
			PolyMesh::FaceHandle face1 = graph.external_edges[i].face1;

			/*
			 * create and add new quad face
			 * use vertices from the first triangle face and the duplicated vertex
			 */
			vector<PolyMesh::VertexHandle> quad_face0;
			for( auto vh_face0 : workMesh.fv_range( face0 ) ) {
				quad_face0.push_back( blossomMesh_step1_vertex_vec[vh_face0.idx()] );
			}
			blossomMesh_step1.add_face( quad_face0 );
			/*
			* create and add new quad face
			* use vertices from the second triangle face and the duplicated vertex
			*/
			vector<PolyMesh::VertexHandle> quad_face1;
			for( auto vh_face1 : workMesh.fv_range( face1 ) ) {
				quad_face1.push_back( blossomMesh_step1_vertex_vec[vh_face1.idx()] );
			}
			blossomMesh_step1.add_face( quad_face1 );
		}
	}

	// return the created hybrid mesh
	return blossomMesh_step1;
	
	/*
	 * use vertex duplication on external edges
	 * duplicate the connecting vertex and move the connecting vertex inside the mesh by a small bit
	 * the duplicated vertex stays on the position of the connecting vertex, so the overall border of the mesh is not altered
	 */
	for (int i = 0; i < graph.external_edges_counter; i++){
		// find edges, that are marked for use of vertex duplication
		if (graph.external_edges[i].use_vertex_duplication && graph.external_edges[i].to_be_used && !graph.external_edges[i].used) {
			graph.external_edges[i].used = true;

			// get the two faces to connect and the connecting vertex from the external edge
			PolyMesh::VertexHandle con_vert = graph.external_edges[i].connecting_vertex;
			PolyMesh::FaceHandle face0 = graph.external_edges[i].face0;
			PolyMesh::FaceHandle face1 = graph.external_edges[i].face1;

			// get the position of the connecting vertex
			OpenMesh::Vec3f y_vert = workMesh.point( con_vert );
			OpenMesh::Vec3f x_0_vert;
			OpenMesh::Vec3f x_1_vert;
			// get the faces indices
			int face0_idx = face0.idx();
			int face1_idx = face1.idx();
			// get needed vertices connected over border edges from the connecting vertex
			for (auto vheh_it = workMesh.voh_iter(con_vert); vheh_it.is_valid(); ++vheh_it) {
				PolyMesh::HalfedgeHandle heh = *vheh_it;
				PolyMesh::HalfedgeHandle oheh = workMesh.opposite_halfedge_handle(heh);
				PolyMesh::FaceHandle fhs = workMesh.face_handle(heh);
				PolyMesh::FaceHandle ofhs = workMesh.face_handle(oheh);
				if (fhs.idx() == -1 ) {
					x_0_vert = workMesh.point( workMesh.from_vertex_handle( heh ) );
				}
				else if (ofhs.idx() == -1) {
					x_1_vert = workMesh.point( workMesh.to_vertex_handle( oheh ) );
				}
			}
			// calculate new position for connecting vertex by moving its position in the direction opposite direction of the midpoint between its two neightbouring vertices in the faces to connect 
			OpenMesh::Vec3f edge_x1_x0 = x_1_vert - x_0_vert;
			OpenMesh::Vec3f mid_x1_x0 = x_0_vert + 0.5 * edge_x1_x0;
			OpenMesh::Vec3f edge_mid_x1_x0_y = mid_x1_x0 - y_vert;
			OpenMesh::Vec3f dup_vert = y_vert - 0.1 * edge_mid_x1_x0_y;
			OpenMesh::Vec3f new_vert = workMesh.point( con_vert );
			// add duplicated vertex to mesh
			blossomMesh_step1_vertex_vec.push_back(blossomMesh_step1.add_vertex(new_vert));
			// move connecting vertex inside the mesh
			blossomMesh_step1.set_point(con_vert, dup_vert);
			int dup_vert_idx = point_count;
			point_count++;
			
			/*
			 * create and add new quad face
			 * use vertices from the first triangle face and the duplicated vertex
			 */
			vector<PolyMesh::VertexHandle> quad_face0;
			for( auto vh_face0 : workMesh.fv_range(face0) ) {
				if (vh_face0.idx() == con_vert.idx()) {
					quad_face0.push_back(blossomMesh_step1_vertex_vec[vh_face0.idx()]);
					quad_face0.push_back(blossomMesh_step1_vertex_vec[dup_vert_idx]);
				}
				else {
					quad_face0.push_back(blossomMesh_step1_vertex_vec[vh_face0.idx()]);
				}
			}
			blossomMesh_step1.add_face(quad_face0);
			/*
			* create and add new quad face
			* use vertices from the second triangle face and the duplicated vertex
			*/
			vector<PolyMesh::VertexHandle> quad_face1;
			for( auto vh_face1 : workMesh.fv_range(face1) ) {
				if (vh_face1.idx() == con_vert.idx()) {
					quad_face1.push_back(blossomMesh_step1_vertex_vec[dup_vert_idx]);
					quad_face1.push_back(blossomMesh_step1_vertex_vec[vh_face1.idx()]);						
				}
				else {
					quad_face1.push_back(blossomMesh_step1_vertex_vec[vh_face1.idx()]);
				}
			}
			blossomMesh_step1.add_face(quad_face1);
		}
	}

	/*
	 * use edge swap on external edges
	 * create to quad faces from two triangle faces and a quad face, seperating these two
	 */
	for (int i = 0; i < graph.external_edges_counter; i++) {	
		// find edges, that are marked for use of edge swap
		if (graph.external_edges[i].use_edge_swap && graph.external_edges[i].to_be_used) {
			// get the two faces to connect and the connecting vertex from the external edge
			PolyMesh::VertexHandle con_vert = graph.external_edges[i].connecting_vertex;
			PolyMesh::FaceHandle old_face0 = graph.external_edges[i].face0;
			PolyMesh::FaceHandle old_face1= graph.external_edges[i].face1;
			PolyMesh::FaceHandle new_face;
			vector<PolyMesh::VertexHandle> needed_vertices(5);

			int edge_index = 0;
			// get the vertices that are connected to the connecting vertex but are no border vertices
			for (auto voh_it_blossom = workMesh.voh_iter(con_vert); voh_it_blossom.is_valid(); ++voh_it_blossom) {
				PolyMesh::HalfedgeHandle heh_new = *voh_it_blossom;
				if (edge_index == 1) {
					needed_vertices[2] = workMesh.to_vertex_handle(heh_new);
				}
				else if (edge_index == 2) {
					needed_vertices[3] = workMesh.to_vertex_handle(heh_new);
				}
				edge_index++;				
			}
			
			// get the quad face needed for edge swap
			int edge_index_b = 0;
			for (auto voh_it_blossom = blossomMesh_step1.voh_iter(con_vert); voh_it_blossom.is_valid(); ++voh_it_blossom) {
				//cout << "edge index b: " << edge_index_b << endl;
				PolyMesh::HalfedgeHandle heh_new = *voh_it_blossom;
				if (edge_index_b == 1) {
					new_face = blossomMesh_step1.face_handle(heh_new);
				}
				edge_index_b++;
			}
			if (edge_index_b > 2) {
				graph.external_edges[i].use_edge_swap = false;
				graph.external_edges[i].use_vertex_duplication = true;
				continue;
			}

			// get the last vertex from the needed quad face 
			for (auto fh_it = blossomMesh_step1.fh_iter(new_face); fh_it.is_valid(); fh_it++) {
				PolyMesh::HalfedgeHandle heh_new = *fh_it;
				PolyMesh::VertexHandle vh_new = blossomMesh_step1.from_vertex_handle(heh_new);
				if ((vh_new.idx() != con_vert.idx()) && (vh_new.idx() != needed_vertices[2].idx()) && (vh_new.idx() != needed_vertices[3].idx())) {
					needed_vertices[4] = vh_new;
				}
			}

			// get the border vertex from the first triangle face
			for (auto fh_it = workMesh.fh_iter(old_face0); fh_it.is_valid(); fh_it++) {
				PolyMesh::HalfedgeHandle heh_old = *fh_it;
				PolyMesh::VertexHandle vh_new = workMesh.to_vertex_handle(heh_old);
				if (vh_new.idx() != con_vert.idx() && vh_new.idx() != needed_vertices[2].idx() && vh_new.idx() != needed_vertices[3].idx()) {
					needed_vertices[0] = vh_new;
				}
			}

			// get the border vertex from the second triange face
			for (auto fh_it = workMesh.fh_iter(old_face1); fh_it.is_valid(); fh_it++) {
				PolyMesh::HalfedgeHandle heh_old = *fh_it;
				PolyMesh::VertexHandle vh_new = workMesh.to_vertex_handle(heh_old);
				if (vh_new.idx() != con_vert.idx() && vh_new.idx() != needed_vertices[2].idx() && vh_new.idx() != needed_vertices[3].idx()) {
					needed_vertices[1] = vh_new;
				}
			} 
			// delete old face
			blossomMesh_step1.request_face_status();
			blossomMesh_step1.delete_face(new_face, false);
			// add new faces
			blossomMesh_step1.request_face_status();
			blossomMesh_step1.add_face({ needed_vertices[4], needed_vertices[2], needed_vertices[0], con_vert });
			blossomMesh_step1.request_face_status();
			blossomMesh_step1.add_face({ needed_vertices[4], con_vert, needed_vertices[1], needed_vertices[3] });
			// garbage collection 
			blossomMesh_step1.garbage_collection();
		}
	}

	// return the created quad mesh
	return blossomMesh_step1;
}

// DEPRECATED
float blossom_quad::quality_function(const std::vector<PolyMesh::VertexHandle>& quad_vertices)
{
	// get all four quad points
	auto point0 = workMesh.point(quad_vertices[0]);
	auto point1 = workMesh.point(quad_vertices[1]);
	auto point2 = workMesh.point(quad_vertices[2]);
	auto point3 = workMesh.point(quad_vertices[3]);

	// calculate the four edges of the quad
	auto edge0 = point1 - point0;
	auto edge1 = point2 - point1;
	auto edge2 = point3 - point2;
	auto edge3 = point0 - point3;

	// calculate the angle in each of the four quad corners
	vector<float> corners;
	// corner 0
	float scal0 = edge0 | -edge3;
	float bet0 = edge0.length() * edge3.length();
	float temp0 = scal0 / bet0;
	float corner0 = std::acos(std::min(std::max((double)temp0, -1.0), 1.0));
	corners.push_back((M_PI / 2.f - corner0));
	
	// corner 1
	float scal1 = edge1 | -edge0;
	float bet1 = edge1.length() * edge0.length();
	float temp1 = scal1 / bet1;
	float corner1 = std::acos(std::min(std::max((double)temp1, -1.0), 1.0));
	corners.push_back(( M_PI / 2.f - corner1));

	// corner 2
	float scal2 = edge2 | -edge1;
	float bet2 = edge2.length() * edge1.length();
	float temp2 = scal2 / bet2;
	float corner2 = std::acos(std::min(std::max((double)temp2, -1.0), 1.0));
	corners.push_back(( M_PI / 2.f - corner2));

	// corner 3
	float scal3 = edge3 | -edge2;
	float bet3 = edge3.length() * edge2.length();
	float temp3 = scal3 / bet3;
	float corner3 = std::acos(std::min(std::max((double)temp3, -1.0), 1.0));
	corners.push_back((M_PI / 2.f - corner3));

	//find largest corner
	float max_corner = *max_element(begin(corners), end(corners));
	float temp = 1.f - ((2.f / M_PI) * max_corner);
	float quality = 1 - std::max(temp, 0.f);

	// return the calculated quality of the quad
	return quality;
}

std::vector<PolyMesh::VertexHandle> blossom_quad::connect_triangles( const PolyMesh::HalfedgeHandle& heh ) {
	LOG_ASSERT( !workMesh.is_boundary( workMesh.edge_handle( heh ) ) );

	std::vector<PolyMesh::VertexHandle> vhs(4);
	vhs[0] = workMesh.from_vertex_handle( heh );
	vhs[1] = workMesh.to_vertex_handle( workMesh.next_halfedge_handle( workMesh.opposite_halfedge_handle( heh ) ) );
	vhs[2] = workMesh.to_vertex_handle( heh );
	vhs[3] = workMesh.to_vertex_handle( workMesh.next_halfedge_handle( heh ) );

	return vhs;
}