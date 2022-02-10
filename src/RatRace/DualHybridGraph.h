#pragma once

#include <vector>
#include <array>
#include <map>
#include <experimental/filesystem>

#include <glog/logging.h>

#include "MeshHeader.h"

class DualHybridGraph
{
	
public:
	using NodeIdxT = long;
	using HalfedgeIdxT = long;

	class Halfedge
	{
		NodeIdxT fromPointIdx_ = -1;
		NodeIdxT toPointIdx_ = -1;
		HalfedgeIdxT idx_ = -1;
		bool isValid_ = false;

		HalfedgeIdxT oppositeIdx_ = -1;
		HalfedgeIdxT nextIdx_ = -1;
		HalfedgeIdxT prevIdx_ = -1;

		PolyMesh::VertexHandle vertex_;


	public:
		Halfedge() {}
		Halfedge( HalfedgeIdxT idx ) : idx_{ idx }, isValid_{ true } {}
		//Halfedge( NodeIdxT from, NodeIdxT to ) : from_{ from }, to_{ to }, isValid_{ true } {}

		auto& from() { return fromPointIdx_; }
		auto from() const { return fromPointIdx_; }
		auto& to() { return toPointIdx_; }
		auto to() const { return toPointIdx_; }
		auto idx() const { return idx_; }
		auto isValid() const { return isValid_; }
		auto& isValid() { return isValid_; }

		auto& oppositeIdx() { return oppositeIdx_; }
		auto& nextIdx() { return nextIdx_; }
		auto& prevIdx() { return prevIdx_; }

		auto& vertex() { return vertex_; }
		auto vertex() const { return vertex_; }
	};

	class : public std::vector<Halfedge>
	{
	public:
		Halfedge& add() {
			push_back( Halfedge( (HalfedgeIdxT)this->size() ) );
			return this->back();
		}

	} halfedges;

	class Node
	{
		NodeIdxT idx_ = -1;
		bool isValid_ = false;
		HalfedgeIdxT halfedgeIdx_ = -1;
		std::vector<PolyMesh::VertexHandle> vertices_;
	public:
		Node() {}
		Node( NodeIdxT idx ) : idx_{ idx }, isValid_{ true } {}
		auto idx() const { return idx_; }
		auto isValid() const { return isValid_; }
		auto& isValid() { return isValid_; }
		auto& halfedgeIdx() { return halfedgeIdx_; }
		auto halfedgeIdx() const { return halfedgeIdx_; }
		auto& vertices() { return vertices_; }
		auto vertices() const { return vertices_; }

	};

	class : public std::vector<Node>
	{
	public:
		Node& add() {
			push_back( Node((NodeIdxT)this->size()) );
			return this->back();
		}
	} nodes;

	PolyMesh mesh_;
	

	DualHybridGraph( PolyMesh& mesh );
	
	// return ALL neighbors including boundaries, e.g. invalid node handles
	std::vector<NodeIdxT> neighborIdxs( const NodeIdxT& idx );
	// return neighbors without boundaries (no invalid node handles)
	std::vector<NodeIdxT> innerNeighborIdxs( const NodeIdxT& idx );

	// outgoing halfedges
	std::vector<HalfedgeIdxT> ohIdxs( const NodeIdxT& idx );
	// ingoing halfedges
	std::vector<HalfedgeIdxT> ihIdxs( const NodeIdxT& idx );

	// return halfedge index
	HalfedgeIdxT findHalfedge( const NodeIdxT& nodeFrom, const NodeIdxT& nodeTo );

	PolyMesh getMesh();

	void printMesh( const std::experimental::filesystem::path& filename );

	void halfedgeCollapse( const HalfedgeIdxT& idx );

	/* Split nodeFrom of halfedges[idx] s. th. it contains a triangle neighboring nodeTo.
	   The split is performed s.th. the best quad is generated.
	   The new halfedge points to nodeFrom, halfedges[idx] remains valid.
	*/
	void split3to4( const HalfedgeIdxT& idx );

	void garbageCollection();

	bool isTriangle( const NodeIdxT& idx );

	bool hasTriangles();

	int nTriangles();

	// Search for the first triangle in nodes and connect it with the nearest triangle
	void connectTriangles();

	void printGridInfo();
};