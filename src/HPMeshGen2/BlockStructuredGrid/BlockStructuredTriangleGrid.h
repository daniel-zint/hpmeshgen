#pragma once

#include "MeshHeader.h"
#include "HPMeshGen2/HelperFunctions.h"

#include <vector>

class BlockStructuredTriangleGrid
{

	int n_mpi_threads;
	int* mpi_patch_id;
	// number of refinement steps within a patch
	size_t n_refinement_steps_;


	// All boundary information is stored on the patch mesh
	// specifies if the boundary belongs to Neumann or Dirichlet boundaries
	OpenMesh::EPropHandleT<bool> isNeumann;
	// stores the boundary index. For specification of the boundary the property "isNeumann" is required
	OpenMesh::EPropHandleT<int> boundaryID;
	// for each Neumann boundary the type is specified here
	std::vector<int> neumannBoundaryTypes;
	// for each Dirichlet boundary the type is specified here
	std::vector<int> dirichletBoundaryTypes;

	// block structure (segmentation of patchMesh)
	std::vector<std::vector<size_t>> blocks;

	std::string mesh_name_;
	std::string mesh_folder_;
public:
	// Patch mesh of BSG
	TriMesh fragmentMesh;
	// BSG as TriangleMesh (do not change topology of this mesh!!!!)
	TriMesh triMesh;

	// Pointer to quadMesh. This just helps to make the code nicer.
	TriMesh* fineMesh;

	// Default constructor
	BlockStructuredTriangleGrid() = delete;
	// Requires a patch structure consisting of quads and the number of nodes within a patch in each dimension
	BlockStructuredTriangleGrid(const TriMesh& pm, const size_t n_refinement_steps, OpenMesh::EPropHandleT<bool>& isNeumann, OpenMesh::EPropHandleT<int>& boundaryID, std::vector<int>& neumannBoundaryTypes, std::vector<int>& dirichletBoundaryTypes);
	// Copy constructor
	BlockStructuredTriangleGrid(const BlockStructuredTriangleGrid&) = delete;

	// Destructor
	~BlockStructuredTriangleGrid();

	// Copy assignment
	const BlockStructuredTriangleGrid& operator=(const BlockStructuredTriangleGrid&) = delete;

	inline size_t n_patches();

	inline size_t n_refinement_steps();

	// export mesh in ADCIRC format (including boundary information)
	void exportMeshAdcirc(const std::string& folder, const std::string& filename);

	void print_svg(std::string filename, const bool printValence);

};


inline size_t BlockStructuredTriangleGrid::n_patches()
{
	return fragmentMesh.n_faces();
}

inline size_t BlockStructuredTriangleGrid::n_refinement_steps()
{
	return n_refinement_steps_;
}