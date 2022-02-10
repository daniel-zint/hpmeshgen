#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

#include <cstdint>
#include <stdio.h>
#include <cfloat>

#include <type_traits>

#include "math_constants.h"

// Open Mesh
#include "MeshHeader.h"

#include "Metrics.h"
#include "DmoParams.h"
#include "Vertex.h"
#include "DmoMesh.h"
#include "gpuErrchk.h"
#include "DMO_PUBLIC.h"

namespace DMO
{
	// CUDA kernel
	template<typename MetricT>
	void DMO_PUBLIC dmoGPU( thrust::device_vector<float2>& points_d, DmoMesh& dmoMesh1_, const MetricT& metric );

	template<typename MetricT, typename Metric2T = DMO::Metrics::NoMetric>
	class Solver
	{
		TriMesh& mesh_;
		MetricT* metric1_;
		Metric2T* metric2_;

		thrust::host_vector<float2> points_;
		thrust::device_vector<float2> points_d_;

		DmoMesh* dmoMesh1_;
		DmoMesh* dmoMesh2_;

	public:
		Solver<MetricT, Metric2T>( TriMesh& mesh, MetricT* metric1, DmoMesh* dmoMesh1, Metric2T* metric2 = nullptr, DmoMesh* dmoMesh2 = nullptr );

		void solve( int nIterations = 100 );
	};


	/*
	*		========================
	*		Implementation of Solver
	*		========================
	*/

	template<typename MetricT, typename Metric2T>
	Solver<MetricT, Metric2T>::Solver( TriMesh& mesh, MetricT* metric1, DmoMesh* dmoMesh1, Metric2T* metric2, DmoMesh* dmoMesh2 )
		: mesh_( mesh ), metric1_( metric1 ), metric2_( metric2 ), dmoMesh1_(dmoMesh1), dmoMesh2_( dmoMesh2 )
	{
		// copy points
		points_.resize( mesh_.n_vertices() );
		for( auto vh : mesh_.vertices() ) {
			auto p = mesh_.point( vh );
			points_[vh.idx()] = { p[0],p[1] };
		}
		points_d_ = points_;

		dmoMesh1_->copyHostToDevice();
		if constexpr( !std::is_same<Metric2T, Metrics::NoMetric>::value ) {
			dmoMesh2_->copyHostToDevice();
		}

		gpuErrchk( cudaDeviceSynchronize() );

	}

	template<typename MetricT, typename Metric2T>
	void Solver<MetricT, Metric2T>::solve( int nIterations ) {
		for( int i = 0; i < nIterations; ++i ) {
			dmoGPU<MetricT>( points_d_, *dmoMesh1_, *metric1_ );
			if constexpr( !std::is_same<Metric2T, Metrics::NoMetric>::value ) {
				dmoGPU<Metric2T>( points_d_, *dmoMesh2_, *metric2_ );
			}
		}

		points_ = points_d_;
		
		for( const auto& vh : mesh_.vertices() ) {
			int idx = vh.idx();
			TriMesh::Point p = { points_[idx].x, points_[idx].y, 0.f };
			mesh_.set_point( vh, p );
		}

	}

}
