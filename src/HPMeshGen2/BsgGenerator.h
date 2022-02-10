#pragma once

#include <glog/logging.h>

#include "Config.h"
#include "MeshHeader.h"
#include "OceanMesh/OceanMesh.h"
#include "BlockStructuredGrid/BlockStructuredGrid.h"
#include "Simplificator.h"
#include "Contour/Shape.h"

class BsgGenerator
{
	Config config_;
	OceanMesh inputMesh_;
	TriMesh simplifiedMesh_;
	PolyMesh fragmentMesh_;
	Contour::Shape contour_;

	BlockStructuredGrid bsg_;

	Simplificator simplificator_;

public:
	BsgGenerator( const Config& config ) : config_{ config }, inputMesh_( config.meshFile() ), contour_( inputMesh_ ), simplificator_( inputMesh_ ) {}

	void generate();

	auto& bsg() { return bsg_; };

	auto& inputMesh() { return inputMesh_; }
	auto& simplifiedMesh() { return simplifiedMesh_; }
	auto& fragmentMesh() { return fragmentMesh_; }

	auto& contour() { return contour_; }

	auto& simplificator() { return simplificator_; }
		
	void processInput();

	void simplify();

	void remeshSimplifiedMesh();

	/// <summary>
	/// DEPRECATED use generateFragments() instead
	/// </summary>
	void generateQuadFragments();

	void generateFragments();

	void optimizeSimplifiedMesh();

	void adaptElementsToContour();

private:

};