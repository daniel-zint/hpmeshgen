#pragma once

#include <experimental/filesystem>

#include <glog/logging.h>

#include "MeshHeader.h"
#include "DualHybridGraph.h"

class RatRace
{
	using fsp = std::experimental::filesystem::path;

	const fsp meshFile_;
	TriMesh tmesh_;
	PolyMesh hmesh_;

public:
	RatRace( TriMesh mesh, bool useBlossomQuad = false );
	void run(bool doPostprocessing = true);
	void untangle( int maxIter = 1000 );
	void optimize( int maxIter = 1000, int iterPerSwap = 5 );
	void swapValence2Vertices();
	void swapValence2VerticesBFS();
	void swapValence();
	void valence2SplitCollapse();

	auto hmesh() { return hmesh_; }

private:

};