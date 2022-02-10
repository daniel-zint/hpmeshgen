#include "eval.h"
#include "VtkWriter.h"
#include "HPMeshGen2/MeshSmoothing/QualityMetrics.h"

namespace GridEvaluation
{
	void eval( TriMesh& mesh, OceanMesh& inputMesh, const std::experimental::filesystem::path& outputFile ) {

		std::vector<float> edgeData( mesh.n_edges() );

		for( const auto& e : mesh.edges() ) {
			const auto h = e.halfedge( 0 );
			auto p0 = mesh.point( h.from() );
			auto p1 = mesh.point( h.to() );
			float d0 = p0[2];
			float d1 = p1[2];
			p0[2] = 0;
			p1[2] = 0;
			p0 = inputMesh.sphericalToCartesian( p0 );
			p1 = inputMesh.sphericalToCartesian( p1 );

			float l = ( p1 - p0 ).length();
			float d = 0.5f * ( d0 + d1 );
			d = std::sqrt( d );
			edgeData[e.idx()] = l / d;	// maybe the other way round
		}

		// compute factor to average cfl
		//float avgCfl = 0;
		//for( const auto& e : edgeData ) avgCfl += e;
		//avgCfl /= edgeData.size();
		//for( auto& e : edgeData ) e /= avgCfl;

		VtkWriter::FaceData<float> faceData{ "cfl_quotient" };
		faceData.data.resize( mesh.n_faces() );

		for( const auto& f : mesh.faces() ) {
			float minCFL = std::numeric_limits<float>::max();
			for( const auto& e : f.edges() ) {
				minCFL = std::min( minCFL, edgeData[e.idx()] );
			}
			faceData.data[f.idx()] = minCFL;
		}

		VtkWriter::PointData<float> pointDepth{ "depth" };
		pointDepth.data.resize( mesh.n_vertices() );
		for( const auto& v : mesh.vertices() ) {
			pointDepth.data[v.idx()] = mesh.point( v )[2];
		}

		VtkWriter::PointData<float> pointEdgeLength{ "edge_length" };
		pointEdgeLength.data.resize( mesh.n_vertices() );
		for( const auto& v : mesh.vertices() ) {
			pointEdgeLength.data[v.idx()] = 0;
			auto p = mesh.point( v );
			p[2] = 0;
			for( const auto& vv : v.vertices() ) {
				auto pp = mesh.point( vv );
				pp[2] = 0;
				pointEdgeLength.data[v.idx()] += ( p - pp ).length();
			}
			pointEdgeLength.data[v.idx()] /= v.valence();
		}

		VtkWriter::FaceData<float> faceQuality{ "mean_ratio_metric" };
		faceQuality.data.resize( mesh.n_faces() );

		for( const auto& f : mesh.faces() ) {
			faceQuality.data[f.idx()] = QualityMetrics::meanRatioMetric( mesh, f );
		}

		VtkWriter vtkWrite( mesh, outputFile, true );
		vtkWrite.addFaceData( faceData );
		vtkWrite.addPointData( pointDepth );
		vtkWrite.addFaceData( faceQuality );
		vtkWrite.addPointData( pointEdgeLength );

	}
}