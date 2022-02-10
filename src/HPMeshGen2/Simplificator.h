#pragma once

#include <experimental/filesystem>
#include <glog/logging.h>

#include "Config.h"
#include "ScalarField/ScalarField.h"
#include "MeshHeader.h"
#include "HelperFunctions.h"
#include "MeshSmoothing/QualityMetrics.h"
#include "MeshSmoothing/PointSmoothing.h"
#include "OceanMesh/OceanMesh.h"

class Simplificator
{
	OceanMesh mesh_;

	// settings
	bool useCache_{ false };
	std::experimental::filesystem::path cacheFolder_;
	std::experimental::filesystem::path meshFile_;

	bool doEdgeFlips{ true };
	float minFeatureAngle{ 100.f };

	bool printSimplificationOutput_{ false };
	std::experimental::filesystem::path outputFolder_;


public:
	Simplificator( OceanMesh mesh ) : mesh_( mesh ) {
		for ( const auto& vh : mesh_.vertices() ) {
			auto p = mesh_.point( vh );
			p[2] = 0;
			mesh_.set_point( vh, p );
		}
	}

	bool setCache( decltype( cacheFolder_ ) cacheFolder, decltype(meshFile_) meshFile );

	bool simplificationOutput( decltype( outputFolder_ ) outputFolder );
	TriMesh simplify( const ScalarField::ScalarField& backgroundgrid, const size_t &nPatches, const bool convHullSimpl = false );

};