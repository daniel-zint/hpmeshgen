#pragma once

#include <vector>
#include <assert.h>
#include <experimental/filesystem>

#include "ScalarField/ScalarField.h"
#include "MeshHeader.h"

class SignedDistanceFunction
{
protected:
	ScalarField::ScalarField field_;

public:
	SignedDistanceFunction( TriMesh& mesh, const size_t& Nx, const size_t& Ny );

	const auto& field() const { return field_; }

	// load SDF if it already exists in cache, otherwise generate and store it.
	SignedDistanceFunction( TriMesh& mesh, const std::experimental::filesystem::path& cacheFolder, const std::experimental::filesystem::path& meshFile, const size_t& sizeGridSizeX, const size_t& sizeGridSizeY ) {
		namespace fs = std::experimental::filesystem;

		fs::path file = cacheFolder / ( meshFile.stem().string() + "_SignedDistanceFunction_" + std::to_string( sizeGridSizeX ) + "_" + std::to_string( sizeGridSizeY ) + ".bin" );

		auto dfLoad = ScalarField::load( file, sizeGridSizeX, sizeGridSizeY );

		if( dfLoad ) {
			field_ = dfLoad.value();
		} else {
			LOG( INFO ) << "Generate SignedDistanceFunction";
			SignedDistanceFunction sdf( mesh, sizeGridSizeX, sizeGridSizeY );
			LOG( INFO ) << "Store in cache file '" << file << "'";
			sdf.field().writeBinary( file.string() );
			field_ = sdf.field();
		}
	}

private:
	static ScalarField::AxisAlignedBoundingBox computeAABB( TriMesh& mesh ) {
		const auto n_vertices = mesh.n_vertices();
		auto xMax = -std::numeric_limits<float>::max();
		auto yMax = -std::numeric_limits<float>::max();
		auto xMin = std::numeric_limits<float>::max();
		auto yMin = std::numeric_limits<float>::max();
		for( unsigned int i = 0; i < n_vertices; ++i ) {
			TriMesh::Point point = mesh.point( mesh.vertex_handle( i ) );

			xMax = std::max( xMax, point[0] );
			yMax = std::max( yMax, point[1] );
			xMin = std::min( xMin, point[0] );
			yMin = std::min( yMin, point[1] );
		}

		const float x_range_old = xMax - xMin;
		const float y_range_old = yMax - yMin;
		xMax += 0.1 * x_range_old;
		xMin -= 0.1 * x_range_old;
		yMax += 0.1 * y_range_old;
		yMin -= 0.1 * y_range_old;

		return { xMin, xMax, yMin, yMax };
	}

};