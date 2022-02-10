#pragma once

#include "MeshHeader.h"
#include "HPMeshGen2/HelperFunctions.h"
#include "Contour/Shape.h"

#include <vector>
#include <experimental/filesystem>

#include <glog/logging.h>

/* Store all information of a Block Structured Grid
- Mesh containing the patch structure (using OpenMesh)
- Curvilinear grid for each patch (using a Struct)
- Information of mesh partitioning on patch-level for mpi distribution

  Patch structure contains quads exclusively
*/

class BlockStructuredGrid
{
	// number of nodes in one dimension within a patch (including first and last node)
	size_t nGridNodes_ = 0;
	size_t nRefinementSteps_ = 0;
	// vertices on the patch edges
	std::vector<std::vector<PolyMesh::VertexHandle>> edgeVertices_;
	// vertices within a patch including the boundaries/edges
	std::vector<std::vector<std::vector<PolyMesh::VertexHandle>>> fragmentVertices_;
	// faces within a patch (might not be required)
	std::vector<std::vector<std::vector<PolyMesh::FaceHandle>>> fragmentFaces_;

	std::vector<int> maskLower_;
	std::vector<int> maskUpper_;

	Contour::Shape contour_;

	OceanMesh oceanMesh_;

	// block structure (segmentation of patchMesh)
	std::vector<std::vector<size_t>> blocks_;

	std::string BsgName_;
	std::experimental::filesystem::path BsgFolder_;
public:
	// Fragment mesh of BSG
	PolyMesh fragmentMesh;
	// BSG as PolyMesh (do not change topology of this mesh!!!!)
	PolyMesh quadMesh;
	// BSG as TriangleMesh (do not change topology of this mesh!!!!)
	TriMesh triMesh;

	BlockStructuredGrid() {}
	// Requires a patch structure consisting of quads and the number of nodes within a patch in each dimension
	BlockStructuredGrid( const PolyMesh& pm, const size_t n_grid_nodes, const Contour::Shape& contour, const OceanMesh& oceanMesh, const bool positionBoundary = false );
	BlockStructuredGrid& operator=(BlockStructuredGrid&)= default;
	inline size_t nFragments();

	inline size_t nGridNodes();

	inline auto blocks() { return blocks_; }

	inline auto fragmentVertices() { return fragmentVertices_; }

	inline auto& maskLower() { return maskLower_; }

	inline auto& maskUpper() { return maskUpper_; }

	const auto& contour() const { return contour_; }

	PolyMesh get_patch(size_t patchID);

	// Adapt the positions of the quad mesh according to the tri mesh
	void copy_positions_tri2quad();
	// Adapt the positions of the tri mesh according to the quad mesh
	// DEPRECATED!!!
	void copy_positions_quad2tri();

	// export mesh in ADCIRC format (including boundary information)
	void exportMeshAdcirc( const std::experimental::filesystem::path& folder, const std::string& filename, OceanMesh& oceanMesh );

	// create blocks
	void patch_segmentation(const size_t n_blocks);

	// write block structured grid to the given folder. Within this folder another folder with the mesh_name is created. This folder stores all the files.
	void write_bsg(const std::experimental::filesystem::path& folder, const std::string& mesh_name, const bool& writeMasks = false);

	void print_svg(std::string filename);
	
private:
	void write14file( const std::experimental::filesystem::path& file14 );

	void write17file( const std::experimental::filesystem::path& file17, OceanMesh& oceanMesh );

	// write the mesh file of bsg
	void write_mesh( const bool& writeMasks );
	// write block of bsg
	void write_block(size_t blockID, const bool& writeMasks );
	// write fragment of bsg
	void write_fragment_head( std::ofstream& ofs, const size_t blockID, const size_t commID );
	void write_fragment_nodes( std::ofstream& ofs, const size_t blockID, const size_t commID );
	// write depth (z-value)
	void write_depth(const size_t blockID) const;
	// write masks
	void write_masks( const size_t blockID ) const;
	// write emo/efa values to txt file
	void write_emo_efa( const size_t blockID ) const;

	void refine();
};


inline size_t BlockStructuredGrid::nFragments()
{
	return fragmentMesh.n_faces();
}

inline size_t BlockStructuredGrid::nGridNodes()
{
	return nGridNodes_;
}
