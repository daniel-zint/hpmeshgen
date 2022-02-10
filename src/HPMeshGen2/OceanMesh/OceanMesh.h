#pragma once

#include <vector>
#include <experimental/filesystem>
#include <functional>

#include "MeshHeader.h"
#include "ScalarField/ScalarField.h"
#include "HPMeshGen2/SizeField.h"
#include "HPMeshGen2/DepthField.h"

struct OceanParameters
{
	int CARTESIAN = 1;
	int SPHERICAL = 2;
	std::string run_description = "";
	std::string run_id = "";
	int hot_input = -1;
	int coordinate_system_type = CARTESIAN;
	int tidal_potential_parameter = -1;
	int ramp = -1;
	double earth_radius = -1;
	double ref_lam = std::numeric_limits<double>::max();
	double ref_phi = std::numeric_limits<double>::max();
};

/// <summary>
/// Stores all information from a fort.14 file
/// </summary>
class OceanMesh : public TriMesh
{
	std::experimental::filesystem::path filename_;

	OceanParameters parameters;

public:
	// property storing the contour type, e.g. Land, Open Sea, etc.
	OpenMesh::EPropHandleT<int> propContourType;
	// boundary angle at vertex in degrees between 0 (feature) and 180 (no feature).
	OpenMesh::VPropHandleT<float> feature_size;

	int nEmoEfaVals = -1;
	OpenMesh::EPropHandleT<std::vector<float>> propEmo;
	OpenMesh::EPropHandleT<std::vector<float>> propEfa;

	OceanMesh() : TriMesh() {}

	OceanMesh( const std::experimental::filesystem::path& filename );
	OceanMesh& operator=(OceanMesh&)= default;
	OceanMesh(const OceanMesh&)= default;
	void loadBackgroundGrids( const std::experimental::filesystem::path& cacheFolder, const std::experimental::filesystem::path& meshFile, const size_t& nx, const size_t& ny );

	// set z-values to 0
	void flatten();

	ScalarField::ScalarField& sizeField() {
		LOG_ASSERT( sizefield );
		return sizefield.value();
	}
	ScalarField::ScalarField& depthField() {
		LOG_ASSERT( depthfield );
		return depthfield.value();
	}

	template<typename MeshT>
	void applyDepth(MeshT& m);

	TriMesh::Point sphericalToCartesian( const TriMesh::Point& ps ) const;

	const auto& filename() const { return filename_; }
	
private:
	std::optional<ScalarField::ScalarField> sizefield = {};
	std::optional<ScalarField::ScalarField> depthfield = {};

	void readFortFiles();

	void read14File();

	void read15File( const std::experimental::filesystem::path& file15 );

	void read17File( const std::experimental::filesystem::path& file17 );
	
	void readOffFile();

	void readOmFile();
};

template<typename MeshT>
inline void OceanMesh::applyDepth( MeshT& m ) {
	for( const auto& v : m.vertices() ) {
		auto p = m.point( v );
		p[2] = depthField().getScalar( p[0], p[1] );
		m.set_point( v, p );
	}
}
