
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cstdint>
#include <stdio.h>
#include <cfloat>

#include <type_traits>

#include "math_constants.h"

//#include "../BackgroundGrid/SizeGrid.h"
#include "Stopwatch/Stopwatch.h"
// Open Mesh
#include "MeshHeader.h"

#define Q_MEANRATIO 0
#define Q_JACOBIAN 4
#define Q_MINANGLE 5
#define Q_RADIUSRATIO 6
#define Q_MAXANGLE 7
#define Q_CONDITION 8
#define Q_COMBLAPLACEMEANRATIO 9


/* Keep NQ = 8 for two dimensional meshes! This value was chosen because it gives optimal
performance considering a warp-size of 32 because NQ = 8 results in 8 * 8 = 64 nodes
which is double the warp size. Each vertex is computed using one warp where each warp
computes two grid nodes.
Another implementation used 2 warps for one grid but it was slower as syncthreads is
too expensive.
*/
// Size of Quality Mesh
constexpr int NQ = 8;
// number of refinement steps within DMO
constexpr int DEPTH = 3;
// double the maximal number of allowed vertices on the one-ring neighborhood
constexpr int MAX_ONE_RING_SIZE = 32;

// For quality output
constexpr int N_QUALITY_COLS = 10;
// Set this value to print quality
#define PRINT_QUALITY 0


// Error output
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		//fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		fprintf(stderr, "GPUassert: %s. Line %d\n", cudaGetErrorString(code), line);
		if (abort) exit(code);
	}
}


typedef union {
	float floats[2];                 // floats[0] = lowest
	std::int32_t ints[2];                     // ints[1] = lowIdx
	unsigned long long ulong;    // for atomic update
} my_atomics;

__device__ unsigned long long my_atomicArgMax(unsigned long long* address, float val, std::int32_t idx)
{
	my_atomics loc, newValue;
	loc.floats[0] = val;
	loc.ints[1] = idx;
	newValue.ulong = *address;
	while (newValue.floats[0] < val)
		newValue.ulong = atomicCAS(address, newValue.ulong, loc.ulong);
		
	return newValue.ulong;
}


/*
Holds information of the vertex-id, the number of neighbors, and their location in the one-ring vector.
*/
struct Vertex {
	int oneRingID;
	int n_oneRing;
	int id;				// own vertex id
};

/////////////////////////////////
// quality metrics per element //
__host__ __device__ __forceinline__ float meanRatioTri(const float p[3][2]) {

	float e[3][2];
	float e_length_squared[3];

	for (int i = 0; i < 3; ++i) {
		int j = (i + 1) % 3;
		e[i][0] = p[j][0] - p[i][0];
		e[i][1] = p[j][1] - p[i][1];

		e_length_squared[i] = e[i][0] * e[i][0] + e[i][1] * e[i][1];
	}

	float l = e_length_squared[0] + e_length_squared[1] + e_length_squared[2];
	float area = e[0][0] * e[1][1] - e[0][1] * e[1][0];

	if( area < 0 )
		return area;
	return 2.f * sqrt(3.f) * area / l;
}

__host__ __device__ __forceinline__ float conditionQuad( const float p[4][2] ) {

	float e[4][2];
	float els[4];

	for( int i = 0; i < 4; ++i ) {
		int j = ( i + 1 ) % 4;
		e[i][0] = p[j][0] - p[i][0];
		e[i][1] = p[j][1] - p[i][1];
		els[i] = e[i][0] * e[i][0] + e[i][1] * e[i][1];
	}
	
	auto detJ1 = e[0][1] * e[3][0] - e[0][0] * e[3][1];
	auto l1 = els[0] + els[3];
	auto q1 = 2 * detJ1 / l1;

	auto detJ2 = e[1][1] * e[0][0] - e[1][0] * e[0][1];
	auto l2 = els[1] + els[0];
	auto q2 = 2 * detJ2 / l2;

	//auto detJ3 = e[2][1] * e[1][0] - e[2][0] * e[1][1];
	//auto l3 = e_length_squared[2] + e_length_squared[1];
	//auto q3 = 2 * detJ3 / l3;

	auto detJ4 = e[3][1] * e[2][0] - e[3][0] * e[2][1];
	auto l4 = els[3] + els[2];
	auto q4 = 2 * detJ4 / l4;

	auto detMin = fminf( detJ1, detJ2 );
	detMin = fminf( detMin, detJ4 );

	if( detMin < 0 )
		return detMin;

	auto qMin = fminf( q1, q2 );
	qMin = fminf( qMin, q4 );

	return qMin;
}
/////////////////////////////////

///////////////////////////////////////////
// distinguish between triangle and quad //

__host__ __device__ __forceinline__ float computeShapeQuality(const int n_oneRing, const float oneRing[MAX_ONE_RING_SIZE], const float p[2], const int element_size) {
	float q = FLT_MAX;

	switch (element_size)
	{
	case 3:
		for (int k = 0; k < n_oneRing - 1; ++k) {
			float v[3][2] = { { p[0], p[1] },{ oneRing[2 * k], oneRing[2 * k + 1] },{ oneRing[2 * (k + 1)], oneRing[2 * (k + 1) + 1] } };
			q = fminf(q, meanRatioTri(v));
		}
		break;
	default:
		break;
	}

	return q;
}

__host__ __device__ __forceinline__ float computeConditionQuality( const int n_oneRing, const float oneRing[MAX_ONE_RING_SIZE], const float p[2], const int element_size ) {
	float q = FLT_MAX;

	switch( element_size ) {
	case 3:
		break;
	case 4:
		for( int k = 0; k < n_oneRing - 1; k += 2 ) {
			float v[4][2] = { { p[0], p[1] },{ oneRing[2 * k], oneRing[2 * k + 1] },{ oneRing[2 * ( k + 1 )], oneRing[2 * ( k + 1 ) + 1] },{ oneRing[2 * ( k + 2 )], oneRing[2 * ( k + 2 ) + 1] } };
			q = fminf( q, conditionQuad( v ) );
		}
		break;
	default:
		break;
	}

	return q;
}

__host__ __device__ __forceinline__ float computeCombinedLaplaceMeanRatioQuality( const int n_oneRing, const float oneRing[MAX_ONE_RING_SIZE], const float p[2], const int element_size ) {
	float q = FLT_MAX;

	float lp[2] = { 0,0 };
	
	switch( element_size ) {
	case 3:
		// compute laplace point
		for( int k = 0; k < n_oneRing - 1; ++k ) {
			lp[0] += oneRing[2 * k];
			lp[1] += oneRing[2 * k + 1];
		}
		lp[0] /= ( n_oneRing - 1 );
		lp[1] /= ( n_oneRing - 1 );
		lp[0] = p[0] - lp[0];
		lp[1] = p[1] - lp[1];

		for( int k = 0; k < n_oneRing - 1; ++k ) {
			float v[3][2] = { { p[0], p[1] },{ oneRing[2 * k], oneRing[2 * k + 1] },{ oneRing[2 * ( k + 1 )], oneRing[2 * ( k + 1 ) + 1] } };
			q = fminf( q, meanRatioTri( v ) );
		}
		break;
	case 4:
	{
		printf( "This metric is not implemented for triangles yet.\n" );
		break;
	}
	default:
		break;
	}

	if( q < 0.5 )
		return q;
	else
		return 0.5f + 1.f / ( lp[0] * lp[0] + lp[1] * lp[1] + 1 );
}

__host__ __device__ __forceinline__ float quality(const int n_oneRing, const float oneRing[MAX_ONE_RING_SIZE], const float p[2], const int element_size, const int q_crit) {
	switch (q_crit)
	{
	case Q_MEANRATIO:
		return computeShapeQuality(n_oneRing, oneRing, p, element_size);
	case Q_CONDITION:
		return computeConditionQuality( n_oneRing, oneRing, p, element_size );
	case Q_COMBLAPLACEMEANRATIO:
		return computeCombinedLaplaceMeanRatioQuality( n_oneRing, oneRing, p, element_size );
	default:
		return -1;
	}
}


__global__ void printFaceQuality(const float* points, const int* faceVec, const int n_faces, const int element_size, const int q_crit) {
	static int counter = 0;


	int q_vec[N_QUALITY_COLS] = { 0 };
	float q_min = FLT_MAX;

	for (int i = 0; i < n_faces; ++i) {
		float q;

		switch (q_crit)
		{
		case Q_MEANRATIO:
			switch (element_size)
			{
			case 3:
			{
				int v_id[3] = { faceVec[3 * i],  faceVec[3 * i + 1], faceVec[3 * i + 2] };
				float p[3][2];
				for (int i = 0; i < 3; ++i) {
					p[i][0] = points[2 * v_id[i]];
					p[i][1] = points[2 * v_id[i] + 1];
				}

				q = meanRatioTri(p);
				break;
			}
			default:
				printf("Quality for this type of element is unknown. Element-size = %d\n", element_size);
				return;
			}
			break;
		default:
			printf("Quality metric unknown\n");
			break;
		}

		q_vec[int(q * N_QUALITY_COLS - 0.0001)] += 1;

		q_min = fminf(q_min, q);
	}

	printf("%3d: ", counter++);

	for (int i = 0; i < N_QUALITY_COLS; ++i) {
		if (q_vec[i] != 0)
			printf("%4d | ", q_vec[i]);
		else
			printf("     | ");
	}
	printf("| q_min = %1.6f", q_min);
	printf("\n");

}

__global__ void printFaceQuality(const float* points, const int* faceVec, const int n_faces, const int element_size, const int q_crit, float* q_min_vec, float* q_avg_vec) {
	static int counter = 0;

	float q_min = FLT_MAX;
	float q_avg = 0;

	for (int i = 0; i < n_faces; ++i) {
		float q;

		switch (q_crit)
		{
		case Q_MEANRATIO:
			switch (element_size)
			{
			case 3:
			{
				int v_id[3] = { faceVec[3 * i],  faceVec[3 * i + 1], faceVec[3 * i + 2] };
				float p[3][2];
				for (int i = 0; i < 3; ++i) {
					p[i][0] = points[2 * v_id[i]];
					p[i][1] = points[2 * v_id[i] + 1];
				}

				q = meanRatioTri(p);
				break;
			}
			default:
				printf("Quality for this type of element is unknown. Element-size = %d\n", element_size);
				return;
			}
			break;
		default:
			printf("Quality metric unknown\n");
			break;
		}

		q_min = fminf(q_min, q);
		q_avg += q;
	}
	q_avg /= n_faces;

	q_min_vec[counter] = q_min;
	q_avg_vec[counter++] = q_avg;
}


__global__ void optimizeHierarchical(int* coloredVertexIDs, const int cOff, const Vertex* vertices, float* points, int* oneRingVec, const float affineFactor_, const int element_size, const int q_crit, const float grid_scale) {
	const int i1 = threadIdx.x / NQ;
	const int j1 = threadIdx.x % NQ;

	const int i2 = (threadIdx.x + NQ * NQ / 2) / NQ;
	const int j2 = (threadIdx.x + NQ * NQ / 2) % NQ;

	const Vertex& v = vertices[coloredVertexIDs[cOff + blockIdx.x]];

	float q = -FLT_MAX;

	__shared__ float xPos, yPos;
	__shared__ float maxDistx;
	__shared__ float maxDisty;

	__shared__ my_atomics argMaxVal;
	argMaxVal.floats[0] = -FLT_MAX;
	argMaxVal.ints[1] = NQ*NQ;

	__shared__ float oneRing[MAX_ONE_RING_SIZE];

	// min/max search + loading oneRing
	if (threadIdx.x == 0) {
		maxDistx = -FLT_MAX;
		maxDisty = -FLT_MAX;

		for (int k = 0; k < v.n_oneRing - 1; ++k) {
			float oneRingX = points[2 * oneRingVec[v.oneRingID + k]];
			float oneRingY = points[2 * oneRingVec[v.oneRingID + k] + 1];
			oneRing[2 * k] = oneRingX;
			oneRing[2 * k + 1] = oneRingY;

			float xDist = abs(points[2 * v.id] - oneRingX);
			float yDist = abs(points[2 * v.id + 1] - oneRingY);

			maxDistx = fmaxf(maxDistx, xDist);
			maxDisty = fmaxf(maxDisty, yDist);
		}
		
		oneRing[2 * v.n_oneRing - 2] = points[2 * oneRingVec[v.oneRingID + v.n_oneRing - 1]];
		oneRing[2 * v.n_oneRing - 1] = points[2 * oneRingVec[v.oneRingID + v.n_oneRing - 1] + 1];

		xPos = points[2 * v.id];
		yPos = points[2 * v.id + 1];
	}

	// special case of valence 2 in a quad mesh
	if( element_size == 4 && v.n_oneRing == 5 ) {
		if( threadIdx.x == 0 ) {
			auto x1 = oneRing[2 * 0];
			auto y1 = oneRing[2 * 0 + 1];
			auto x2 = oneRing[2 * 2];
			auto y2 = oneRing[2 * 2 + 1];
			xPos = 0.5f * ( x1 + x2 );
			yPos = 0.5f * ( y1 + y2 );
			points[2 * v.id] = xPos;
			points[2 * v.id + 1] = yPos;
		}
		return;
	}

	// start depth iteration
	float depth_scale = grid_scale;
	float argMax = 0;
	for (int depth = 0; depth < DEPTH; ++depth) {

		float xMax, xMin, yMax, yMin;
		xMax = xPos + depth_scale * maxDistx;
		xMin = xPos - depth_scale * maxDistx;
		yMax = yPos + depth_scale * maxDisty;
		yMin = yPos - depth_scale * maxDisty;


		float pos_i1 = affineFactor_ * (i1 * xMin + (NQ - 1 - i1) * xMax);
		float pos_j1 = affineFactor_ * (j1 * yMin + (NQ - 1 - j1) * yMax);
		float pos_i2 = affineFactor_ * (i2 * xMin + (NQ - 1 - i2) * xMax);
		float pos_j2 = affineFactor_ * (j2 * yMin + (NQ - 1 - j2) * yMax);

		float p1[2] = { pos_i1, pos_j1 };
		float q1 = quality(v.n_oneRing, oneRing, p1, element_size, q_crit);
		float p2[2] = { pos_i2, pos_j2 };
		float q2 = quality(v.n_oneRing, oneRing, p2, element_size, q_crit);

		if (q1 > q2) {
			q = q1;
			argMax = 1;
		}
		else {
			q = q2;
			argMax = 2;
		}

		my_atomicArgMax( (unsigned long long *)&( argMaxVal.ulong ), q, i1 * NQ + j1 );

		float pCurrent[2] = { xPos, yPos };
		float qOld = quality(v.n_oneRing, oneRing, pCurrent, element_size, q_crit);
		if (i1 * NQ + j1 == argMaxVal.ints[1] && qOld < q) {
			if (argMax == 1) {
				xPos = pos_i1;
				yPos = pos_j1;
			}
			else {
				xPos = pos_i2;
				yPos = pos_j2;
			}
		}
		
		//depth dependent scaling factor
		depth_scale = depth_scale * (2.f / (NQ - 1));
	}

	// set new position if it is better than the old one
	float pOld[2] = { points[2 * v.id] , points[2 * v.id + 1] };
	float qOld = quality( v.n_oneRing, oneRing, pOld, element_size, q_crit );
	if (i1 * NQ + j1 == argMaxVal.ints[1] && qOld < q) {
		points[2 * v.id] = xPos;
		points[2 * v.id + 1] = yPos;
	}
}


struct UniformGrid
{
	int nx, ny;
	float hx, hy, xMin, yMin, xMax, yMax;
};


///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// OpenMesh //////////////////////////////////////

template<typename T> inline void copyOpenMeshData(T& mesh, float* points, Vertex* vertices, int* oneRingVec) {

	bool isPolyMesh = std::is_same<T, PolyMesh>::value;

	int interior_counter = 0;
	int oneRing_counter = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		TriMesh::Point p = mesh.point(*v_it);

		points[2 * v_it->idx()] = p[0];
		points[2 * v_it->idx() + 1] = p[1];


		if (!mesh.is_boundary(*v_it)) {
			// fill vertex struct

			Vertex& v = vertices[interior_counter];
			v.id = v_it->idx();

			v.n_oneRing = 0;
			for (auto voh_it = mesh.voh_iter(*v_it); voh_it.is_valid(); ++voh_it) {
				++v.n_oneRing;
				if (isPolyMesh && !mesh.is_boundary(*voh_it)) ++v.n_oneRing;
			}
			++v.n_oneRing;

			v.oneRingID = oneRing_counter;

			TriMesh::HalfedgeHandle heh = *(mesh.voh_iter(*v_it));
			TriMesh::HalfedgeHandle heh_init = heh;

			do {
				oneRingVec[oneRing_counter++] = mesh.to_vertex_handle(heh).idx();
				if (isPolyMesh) {
					heh = mesh.next_halfedge_handle(heh);
					oneRingVec[oneRing_counter++] = mesh.to_vertex_handle(heh).idx();
				}
				heh = mesh.next_halfedge_handle(heh);
				heh = mesh.next_halfedge_handle(heh);
				heh = mesh.opposite_halfedge_handle(heh);
			} while (heh.idx() != heh_init.idx());

			oneRingVec[oneRing_counter] = mesh.to_vertex_handle(heh).idx();
			++oneRing_counter;

			++interior_counter;
		}
	}

}

template <typename T> inline void createColoring(T& mesh, const int n_free_vertices, int** coloredVertexIDs, std::vector<int>& colorOffset_) {

	bool isPolyMesh = std::is_same<T, PolyMesh>::value;

	// create coloring scheme
	std::vector<int>colorScheme(mesh.n_vertices(), -1);
	int colorSchemeIt = 0;

	// set boundarys to a value that can be ignored
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		if (mesh.is_boundary(*v_it)) {
			colorScheme[v_it->idx()] = -2;
		}
	}

	while (std::find(colorScheme.begin(), colorScheme.end(), -1) != colorScheme.end()) {
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
	
			if (colorScheme[v_it->idx()] != -1) { continue; }		// vertex is already colored
	
			bool neighborIsCurrent = false;
			for (auto voh_it = mesh.voh_iter(*v_it); voh_it.is_valid(); ++voh_it) {
				PolyMesh::VertexHandle vh1 = mesh.to_vertex_handle(*voh_it);
				if (!isPolyMesh && colorScheme[vh1.idx()] == colorSchemeIt) {
					neighborIsCurrent = true;
					break;
				}
				else if (isPolyMesh) {
					PolyMesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(*voh_it));
					if (colorScheme[vh1.idx()] == colorSchemeIt || colorScheme[vh2.idx()] == colorSchemeIt) {
						neighborIsCurrent = true;
						break;
					}
				}
			}
			if (neighborIsCurrent) { continue; }			// a neighboring vertex is already in this color
	
			colorScheme[v_it->idx()] = colorSchemeIt;
		}
		++colorSchemeIt;
	}

	//// DEBUG one color for each vertex
	//for ( auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it ) {
	//	if ( colorScheme[v_it->idx()] == -1 ) {
	//		colorScheme[v_it->idx()] = colorSchemeIt++;
	//	}
	//}

	int n_colors = *(std::max_element(colorScheme.begin(), colorScheme.end())) + 1;

	if( n_colors == -1 )
		return;

	std::vector<int> n_color_vecs(n_colors, 0);
	for (int i = 0; i < colorScheme.size(); ++i) {
		if (colorScheme[i] > -1)
			++n_color_vecs[colorScheme[i]];
	}

	*coloredVertexIDs = new int[n_free_vertices];

	colorOffset_ = std::vector<int>(n_colors + 1, 0);
	for (int i = 1; i < n_colors; ++i) {
		colorOffset_[i] = colorOffset_[i - 1] + n_color_vecs[i - 1];
	}
	colorOffset_[n_colors] = n_free_vertices;		// mark the end of the colored-vertices vector

													// add vertex ids
	std::vector<int>colorCounter(n_colors, 0);
	int interior_counter = 0;
	for (int i = 0; i < colorScheme.size(); ++i) {
		if (colorScheme[i] < 0) { continue; }
		(*coloredVertexIDs)[colorOffset_[colorScheme[i]] + colorCounter[colorScheme[i]]++] = interior_counter++;
	}
}

// DEBUG VARIABLES //
int runCounter = 0;
float* vertexPosBuf;

template <typename T> void discreteMeshOptimization(T& mesh, const int q_crit = Q_MEANRATIO, const float grid_scale = 0.5f, int n_iter = 100) {

	constexpr bool isPolyMesh = std::is_same<T, PolyMesh>::value;
	constexpr int element_size = std::is_same<T, PolyMesh>::value ? 4 : 3;

	Stopwatch sw;
	int n_free_vertices = 0;
	int oneRingVecLength = 0;
#pragma omp parallel for reduction(+:n_free_vertices,oneRingVecLength)
	for (int i = 0; i < mesh.n_vertices(); ++i) {
		PolyMesh::VertexHandle vh = mesh.vertex_handle(i);
		if (mesh.is_boundary(vh)) { continue; }
		++n_free_vertices;

		for (auto voh_it = mesh.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
			++oneRingVecLength;
			if (isPolyMesh && !mesh.is_boundary(*voh_it)) ++oneRingVecLength;
		}
		++oneRingVecLength;		// additional count s.th. last element is again the first element
	}

	if( n_free_vertices == 0 ) {
		return;
	}

	// convert OpenMesh to a basic structure
	float* points = new float[2 * mesh.n_vertices()];
	Vertex* vertices = new Vertex[n_free_vertices];
	int* oneRingVec = new int[oneRingVecLength];

	float* points_d;
	Vertex* vertices_d;
	int* oneRingVec_d;
	int* coloredVertexIDs_d;

	int* coloredVertexIDs;
	std::vector<int> colorOffset_;


#pragma omp parallel sections num_threads(2)
	{
#pragma omp section
		{
			gpuErrchk(cudaMalloc((void**)&points_d, 2 * mesh.n_vertices() * sizeof(float)));
			gpuErrchk(cudaMalloc((void**)&vertices_d, n_free_vertices * sizeof(Vertex)));
			gpuErrchk(cudaMalloc((void**)&oneRingVec_d, oneRingVecLength * sizeof(int)));
			gpuErrchk(cudaMalloc((void**)&coloredVertexIDs_d, n_free_vertices * sizeof(int)));

			createColoring(mesh, n_free_vertices, &coloredVertexIDs, colorOffset_);

			gpuErrchk(cudaMemcpyAsync(coloredVertexIDs_d, coloredVertexIDs, n_free_vertices * sizeof(int), cudaMemcpyHostToDevice));
		}
#pragma omp section 
		{
			copyOpenMeshData(mesh, points, vertices, oneRingVec);
		}
	}

	gpuErrchk(cudaMemcpyAsync(points_d, points, 2 * mesh.n_vertices() * sizeof(float), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyAsync(vertices_d, vertices, n_free_vertices * sizeof(Vertex), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyAsync(oneRingVec_d, oneRingVec, oneRingVecLength * sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk( cudaDeviceSynchronize() );

	int n_colors = (int)colorOffset_.size() - 1;


	// face vector (only needed for quality evaluation)
#if PRINT_QUALITY
	int* faceVec = new int[mesh.n_faces() * 3];
	for (int i = 0; i < mesh.n_faces(); ++i) {
		TriMesh::FaceHandle fh = mesh.face_handle(i);
		int vertex_counter = 0;
		for (auto fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
			faceVec[3 * i + vertex_counter++] = fv_it->idx();
		}
	}
	int* faceVec_d;
	gpuErrchk(cudaMalloc((void**)&faceVec_d, 3 * mesh.n_faces() * sizeof(int)));
	gpuErrchk(cudaMemcpy(faceVec_d, faceVec, 3 * mesh.n_faces() * sizeof(int), cudaMemcpyHostToDevice));
#endif

	const float affineFactor_ = 1.f / (float)(NQ - 1);


#if PRINT_QUALITY
	// q_min_vec for printing
	float *q_min_vec, *q_avg_vec;
	cudaMallocManaged(&q_min_vec, (n_iter + 1) * sizeof(float));
	cudaMallocManaged(&q_avg_vec, (n_iter + 1) * sizeof(float));

	printf("    ");
	for (int i = 0; i < N_QUALITY_COLS; ++i) {
		printf("<%1.3f|", (float)(i + 1) / (float)N_QUALITY_COLS);
	}
	printf("\n\n");
	printFaceQuality << <1, 1 >> >(vertexPos_d, faceVec_d, mesh.n_faces(), 3, q_crit);
	printFaceQuality << <1, 1 >> >(vertexPos_d, faceVec_d, mesh.n_faces(), 3, q_crit, q_min_vec, q_avg_vec);
#endif // PRINT_QUALITY

	//cudaDeviceSynchronize();
	//sw.start();
	for (int i = 0; i < n_iter; ++i) {
		for (int cid = 0; cid < n_colors; ++cid) {
			const int nBlocks = colorOffset_[cid + 1] - colorOffset_[cid];
			const int nThreads = NQ * NQ / 2;
			optimizeHierarchical << <nBlocks, nThreads >> >(coloredVertexIDs_d, colorOffset_[cid], vertices_d, points_d, oneRingVec_d, affineFactor_, element_size, q_crit, grid_scale);
			gpuErrchk( cudaDeviceSynchronize() );
		}
#if PRINT_QUALITY
		gpuErrchk( cudaDeviceSynchronize() );
		printFaceQuality << <1, 1 >> >(vertexPos_d, faceVec_d, mesh.n_faces(), 3, q_crit);
		printFaceQuality << <1, 1 >> >(vertexPos_d, faceVec_d, mesh.n_faces(), 3, q_crit, q_min_vec, q_avg_vec);
#endif // PRINT_QUALITY
	}
	//cudaDeviceSynchronize();
	//sw.stop();
	//std::cout << "DMO runtime: " << sw.runtimeStr() << std::endl;

#if PRINT_QUALITY
	gpuErrchk( cudaDeviceSynchronize() );
	std::string ofs_name = "../output.txt";
	std::ofstream ofs(ofs_name);
	for (int i = 0; i < n_iter + 1; ++i) {
		ofs << i << " " << q_min_vec[i] << " " << q_avg_vec[i] << std::endl;
	}
	ofs.close();
#endif // PRINT_QUALITY
	gpuErrchk( cudaDeviceSynchronize() );
	cudaMemcpy(points, points_d, 2 * mesh.n_vertices() * sizeof(float), cudaMemcpyDeviceToHost);

	cudaFree(points_d);
	cudaFree(vertices_d);
	cudaFree(oneRingVec_d);
	cudaFree(coloredVertexIDs_d);

	delete[] vertices;
	delete[] oneRingVec;
	delete[] coloredVertexIDs;

#if PRINT_QUALITY
	delete[] faceVec;

	cudaFree(faceVec_d);
	cudaFree(q_min_vec);
	cudaFree(q_avg_vec);
#endif // PRINT_QUALITY

	//sw.start();
	// write vertex positions back to mesh
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		int id = v_it->idx();
		TriMesh::Point p = { points[2 * id], points[2 * id + 1], 0.f };
		mesh.set_point(*v_it, p);
	}
	//sw.stop();
	//std::cout << "Write back runtime: " << sw.runtimeStr() << std::endl;

	delete[] points;
}

/////////////////////////////////////////////////////
//////////////////// CPU Version ////////////////////
inline void optimizeHierarchical(const int vid, const Vertex* vertices, float* points, int* oneRingVec, const float affineFactor_, const int element_size, const int q_crit, const float grid_scale) {

	const Vertex& v = vertices[vid];

	float xPos, yPos;
	float maxDistx = 0, maxDisty = 0;

	float oneRing[MAX_ONE_RING_SIZE];

	// min/max search + loading oneRing
	for (int k = 0; k < v.n_oneRing - 1; ++k) {
		float oneRingX = points[2 * oneRingVec[v.oneRingID + k]];
		float oneRingY = points[2 * oneRingVec[v.oneRingID + k] + 1];
		oneRing[2 * k] = oneRingX;
		oneRing[2 * k + 1] = oneRingY;

		float xDist = abs(points[2 * v.id] - oneRingX);
		float yDist = abs(points[2 * v.id + 1] - oneRingY);

		maxDistx = fmaxf(maxDistx, xDist);
		maxDisty = fmaxf(maxDisty, yDist);
	}

	// set xmaxmin...
	maxDistx = grid_scale * maxDistx;
	maxDisty = grid_scale * maxDisty;

	oneRing[2 * v.n_oneRing - 2] = points[2 * oneRingVec[v.oneRingID + v.n_oneRing - 1]];
	oneRing[2 * v.n_oneRing - 1] = points[2 * oneRingVec[v.oneRingID + v.n_oneRing - 1] + 1];

	xPos = points[2 * v.id];
	yPos = points[2 * v.id + 1];

	float pOld[2] = { xPos, yPos };
	float q = quality(v.n_oneRing, oneRing, pOld, element_size, q_crit);

	// start depth iteration
	float depth_scale = grid_scale;
	for (int depth = 0; depth < DEPTH; ++depth) {

		float xMax, xMin, yMax, yMin;
		xMax = xPos + depth_scale * maxDistx;
		xMin = xPos - depth_scale * maxDistx;
		yMax = yPos + depth_scale * maxDisty;
		yMin = yPos - depth_scale * maxDisty;

#pragma omp parallel for
		for (int i = 0; i < NQ; ++i) {
			float pos_i = affineFactor_ * (i * xMin + (NQ - 1 - i) * xMax);

			for (int j = 0; j < NQ; ++j) {
				float pos_j = affineFactor_ * (j * yMin + (NQ - 1 - j) * yMax);

				float pCurrent[2] = { pos_i, pos_j };
				float qCurrent = quality(v.n_oneRing, oneRing, pCurrent, element_size, q_crit);

				#pragma omp critical
				if (qCurrent > q) {
					xPos = pos_i;
					yPos = pos_j;
					q = qCurrent;
				}
			}
		}

		//depth dependent scaling factor
		depth_scale = depth_scale * (2.f / (NQ - 1));
	}


	// set new position if it is better than the old one
	points[2 * v.id] = xPos;
	points[2 * v.id + 1] = yPos;
}

void discreteMeshOptimizationCPU(TriMesh& mesh, const int q_crit = Q_MEANRATIO, const float grid_scale = 0.5f, int n_iter = 100) {

	int n_free_vertices = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) { if (!mesh.is_boundary(*v_it)) ++n_free_vertices; }
	//printf("N free vertices = %d\n", n_free_vertices);

	int oneRingVecLength = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		if (mesh.is_boundary(*v_it)) { continue; }
		for (auto vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) { ++oneRingVecLength; }
		++oneRingVecLength;		// additional count s.th. last element is again the first element
	}

	// convert OpenMesh to a basic structure
	float* points = new float[2 * mesh.n_vertices()];
	Vertex* vertices = new Vertex[n_free_vertices];
	int* oneRingVec = new int[oneRingVecLength];

	Stopwatch sw;
	copyOpenMeshData(mesh, points, vertices, oneRingVec);

	const float affineFactor_ = 1.f / (float)(NQ - 1);

	//sw.start();
	for (int i = 0; i < n_iter; ++i) {
		for (int vid = 0; vid < n_free_vertices; ++vid) {
			optimizeHierarchical(vid, vertices, points, oneRingVec, affineFactor_, 3, q_crit, grid_scale);
		}
	}
	//sw.stop();
	//std::cout << "DMO runtime: " << sw.runtimeStr() << std::endl;

	// write vertex positions back to mesh
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		int id = v_it->idx();
		TriMesh::Point p = { points[2 * id], points[2 * id + 1], 0.f };
		mesh.set_point(*v_it, p);
	}

	delete[] points;
	delete[] vertices;
	delete[] oneRingVec;
}


template void discreteMeshOptimization( TriMesh&, const int, const float, int );
template void discreteMeshOptimization( PolyMesh&, const int, const float, int );