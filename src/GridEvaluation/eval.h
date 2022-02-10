#pragma once

#include <experimental/filesystem>
#include <fstream>
#include <vector>
#include <glog/logging.h>
#include "MeshHeader.h"
#include "HPMeshGen2/OceanMesh/OceanMesh.h"

// CFL Condition Test for numerical stability
// per node: compute sqrt(depth)
// per triangle: compute d / min(sqrt(depth))   with d: longest edge of triangle

namespace GridEvaluation
{
	void eval( TriMesh& mesh, OceanMesh& inputMesh, const std::experimental::filesystem::path& outputFile );

}
