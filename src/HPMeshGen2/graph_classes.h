#pragma once
#include <vector>
#include <iostream>
#include "MeshHeader.h"

struct sort_edge {
	int m_orig_position;
	float m_cost_value;

	sort_edge(int o, float v) :
		m_orig_position(o),
		m_cost_value(v)
	{}
	bool operator < (const sort_edge& other_s) const {
		return (m_cost_value < other_s.m_cost_value);
	}
};

struct graph_edge {
	
	PolyMesh::FaceHandle face0;
	PolyMesh::FaceHandle face1;

	std::vector<PolyMesh::VertexHandle> quad_vertices;
	PolyMesh::VertexHandle connecting_vertex;

	bool to_be_used = false;
	bool use_edge_swap = false;
	bool use_vertex_duplication = false;
	bool used = false;
	int cost;
};


struct con_graph {
	
	std::vector<graph_edge> edges;
	std::vector<graph_edge> external_edges;
	
	int edges_counter = 0;
	int external_edges_counter = 0;

	void print_graph() {
		for( int i = 0; i < edges_counter; i++ ) {
			std::cout << "edge between face " << edges[i].face0.idx() << " and face " << edges[i].face1.idx() << std::endl;
		}
		for( int j = 0; j < external_edges_counter; j++ ) {
			std::cout << "external edge between face " << external_edges[j].face0.idx() << " and face " << external_edges[j].face1.idx() << " over vertex " << external_edges[j].connecting_vertex.idx() << std::endl;
		}
	}
};

