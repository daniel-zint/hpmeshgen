#pragma once

#include "graph_classes.h"
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <glog/logging.h>

#include "MeshHeader.h"

#include "Blossom5/PerfectMatching.h"

constexpr float EXTERNAL_EDGE_COST = 10000.f;

class blossom_quad {
	PolyMesh workMesh;
	PolyMesh blossomMesh_step1;
	std::vector<PolyMesh::VertexHandle> blossomMesh_step1_vertex_vec;
	PolyMesh blossomMesh_refined;

public:
	PolyMesh do_blossom_algo(PolyMesh loadedMesh);

private:
	/*
	 * function to connect two triangles to a quad
	 */
	std::vector<PolyMesh::VertexHandle> connect_triangles( const PolyMesh::HalfedgeHandle& heh );

	/*
	 * function to calculate the quality of a quad, based on its angle sizes
	 */
	// DEPRECATED
	float quality_function( const std::vector<PolyMesh::VertexHandle>& quad_vertices );
};