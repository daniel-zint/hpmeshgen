#pragma once

#include <ScalarField/ScalarField.h>
#include "MeshHeader.h"
#include "HelperFunctions.h"
#include <nanoflann.hpp>

class SizeField
{
	struct MeshPointsNN
	{
		using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, MeshPointsNN >, MeshPointsNN, 3>;
		std::unique_ptr<kd_tree_t> kdTree_;
		std::vector<TriMesh::Point> points_;

		MeshPointsNN( TriMesh& mesh ) {
			points_.reserve( mesh.n_vertices() );
			for( const auto& vh : mesh.vertices() ) {
				points_.push_back( mesh.point( vh ) );
			}
			for( auto& p : points_ ) {
				p[2] = 0;
			}

			// init kd tree
			kdTree_ = std::make_unique<kd_tree_t>( 3, *this );
			kdTree_->buildIndex();
		}

		std::vector<size_t> findNearestNeighbor( const TriMesh::Point& p, int k ) const {
			std::vector<size_t> ret_indexes( k );
			std::vector<float> out_dists_sqr( k );
			nanoflann::KNNResultSet<float> resultSet( k );
			resultSet.init( &ret_indexes[0], &out_dists_sqr[0] );
			kdTree_->findNeighbors( resultSet, &p[0], nanoflann::SearchParams() );
			return ret_indexes;
		}

		/*----- nanoflann functions -----*/
		// Must return the number of data points
		inline size_t kdtree_get_point_count() const { return points_.size(); }
		// Returns the dim'th component of the idx'th point in the class
		inline float kdtree_get_pt( const size_t idx, const size_t dim ) const {
			return points_[idx][dim];
		}
		template <class BBOX>
		bool kdtree_get_bbox( BBOX& ) const { return false; }
	};

	ScalarField::ScalarField field_;

public:
	SizeField( const std::experimental::filesystem::path& filename ) : field_(filename) {}
	SizeField( const ScalarField::ScalarField& field ) : field_( field ) {}
	SizeField( TriMesh& mesh, const size_t& Nx, const size_t& Ny, const bool useInterpolation = true ) 
		: field_( useInterpolation ? computeAabbLecacy( mesh ) : computeAabbLecacy( mesh ), Nx, Ny )
	{
		if( useInterpolation ) {
			fillGrid( mesh );
		} else {
			fillGridFromNearestNeighbors( mesh );
		}
	}

	const auto& field() const { return field_; }

private:

	void gaussSeidelStep( ScalarField::ScalarField& sf, int i, int j, float* v, float* vNew ) {
		int n = 0;
		float val = 0;
		if( i > 0 && v[sf.lex( i - 1, j )] != -FLT_MAX ) {
			val += v[sf.lex( i - 1, j )];
			++n;
		}
		if( i < sf.Nx() - 1 && v[sf.lex( i + 1, j )] != -FLT_MAX ) {
			val += v[sf.lex( i + 1, j )];
			++n;
		}
		if( j > 0 && v[sf.lex( i, j - 1 )] != -FLT_MAX ) {
			val += v[sf.lex( i, j - 1 )];
			++n;
		}
		if( j < sf.Ny() - 1 && v[sf.lex( i, j + 1 )] != -FLT_MAX ) {
			val += v[sf.lex( i, j + 1 )];
			++n;
		}

		if( n != 0 )
			vNew[sf.lex( i, j )] = val / n;
	}

	void gaussSeidelStepInterior( ScalarField::ScalarField& sf, int i, int j, float* v, float* vNew ) {
		vNew[sf.lex( i, j )] = 0.25f * (
			v[sf.lex( i - 1, j )] +
			v[sf.lex( i + 1, j )] +
			v[sf.lex( i, j - 1 )] +
			v[sf.lex( i, j + 1 )] );
	}

	void fillGrid( TriMesh& mesh ) {
		// fill scalar field
		std::vector<float> vertexAvg( mesh.n_vertices() );
		std::vector<bool> isDomain( field_.scalars().size(), false );
		for( const auto& v : mesh.vertices() ) {
			float avg{ 0.f };
			auto nEdges{ 0 };

			for( const auto& voh : v.outgoing_halfedges() ) {
				++nEdges;
				auto p1 = mesh.point( voh.from() );
				auto p2 = mesh.point( voh.to() );
				p1[2] = 0;
				p2[2] = 0;
				avg += ( p1 - p2 ).length();
			}

			vertexAvg[v.idx()] = avg / (float)nEdges;
		}

		for( const auto& fh : mesh.faces() ) {
			std::vector<TriMesh::Point> triangle;	// points of the triangle

			for( const auto& v : fh.vertices() ) {
				TriMesh::Point p{ mesh.point( v )[0], mesh.point( v )[1], vertexAvg[v.idx()] };
				triangle.push_back( p );
			}


			// get the rectangle of grid points that surrounds the triangle
			float x_t_max = std::fmaxf( std::max( triangle[0][0], triangle[1][0] ), triangle[2][0] );
			float x_t_min = std::fminf( std::min( triangle[0][0], triangle[1][0] ), triangle[2][0] );
			float y_t_max = std::fmaxf( std::max( triangle[0][1], triangle[1][1] ), triangle[2][1] );
			float y_t_min = std::fminf( std::min( triangle[0][1], triangle[1][1] ), triangle[2][1] );

			size_t i_min = field_.x2col( x_t_min );
			size_t i_max = field_.x2col( x_t_max ) + 1;
			size_t j_min = field_.y2row( y_t_min );
			size_t j_max = field_.y2row( y_t_max ) + 1;

			// check for all inner points of rectangle if they are inside the triangle. If yes, calculate the depth
			for( size_t i = i_min + 1; i < i_max; ++i ) {
				const float xv = field_.x( i );

				for( size_t j = j_min + 1; j < j_max; ++j ) {
					const float yv = field_.y( j );
					const auto [a, b, c] = HelperFunctions::barycentricCoordinates( { triangle[0], triangle[1], triangle[2] }, { xv,yv,0 } );

					if( HelperFunctions::isInsideTriangle( a, b, c ) ) {
						// calculate average edge length
						field_( i, j ) = HelperFunctions::barycentricInterpolation( triangle, a, b, c );
						isDomain[field_.lex( i, j )] = true;
					}
				}
			}
		}

		auto v2 = field_.scalars();
		auto vDiff = decltype( v2 )( v2.size(), 0 );
		float* pv1 = field_.scalars().data();
		float* pv2 = v2.data();

		const auto& v1Size = field_.scalars().size();

		std::vector<std::array<int, 2>> warmUpRed;
		std::vector<std::array<int, 2>> warmUpBlack;

		for( int i = 0; i < field_.Nx(); ++i ) {
			for( int j = 0; j < field_.Ny(); ++j ) {
				if( isDomain[field_.lex( i, j )] ) continue;
				if( ( i + j ) % 2 == 0 )
					warmUpRed.push_back( { i,j } );
				else
					warmUpBlack.push_back( { i,j } );
			}
		}

		const auto warmUpRedSize = warmUpRed.size();
		const auto warmUpBlackSize = warmUpBlack.size();

		// warm-up phase (spread reasonable values fast)
		for( long k = 0; k < 1e6; ++k ) {

		#pragma omp parallel for
			for( int idx = 0; idx < warmUpRedSize; ++idx ) {
				const int i = warmUpRed[idx][0];
				const int j = warmUpRed[idx][1];
				gaussSeidelStep( field_, i, j, pv1, pv2 );
			}
		#pragma omp parallel for
			for( int idx = 0; idx < warmUpBlackSize; ++idx ) {
				const int i = warmUpBlack[idx][0];
				const int j = warmUpBlack[idx][1];
				gaussSeidelStep( field_, i, j, pv2, pv2 );
			}

			std::swap( pv1, pv2 );
			auto minIt = std::min_element( field_.scalars().begin(), field_.scalars().end() );

			if( *minIt != -FLT_MAX ) {
				break;
			}
		}

		std::vector<std::array<int, 2>> interiorRed;
		std::vector<std::array<int, 2>> interiorBlack;

		for( int i = 1; i < field_.Nx() - 1; ++i ) {
			for( int j = 1; j < field_.Ny() - 1; ++j ) {
				if( isDomain[field_.lex( i, j )] ) continue;
				if( ( i + j ) % 2 == 0 )
					interiorRed.push_back( { i,j } );
				else
					interiorBlack.push_back( { i,j } );
			}
		}

		const auto interiorRedSize = interiorRed.size();
		const auto interiorBlackSize = interiorBlack.size();

		// main loop (after warm-up)
		for( long k = 0; k < 1e6; ++k ) {

		#pragma omp parallel for
			for( int idx = 0; idx < interiorRedSize; ++idx ) {
				const int i = interiorRed[idx][0];
				const int j = interiorRed[idx][1];
				gaussSeidelStepInterior( field_, i, j, pv1, pv2 );
			}
			// edges red
			for( int i = 2; i < field_.Nx() - 1; i += 2 ) {
				const int j = 0;
				pv2[field_.lex( i, j )] = (
					pv1[field_.lex( i - 1, j )] +
					pv1[field_.lex( i + 1, j )] +
					pv1[field_.lex( i, j + 1 )] ) / 3;
			}
			for( int j = 2; j < field_.Ny() - 1; j += 2 ) {
				const int i = 0;
				pv2[field_.lex( i, j )] = (
					pv1[field_.lex( i + 1, j )] +
					pv1[field_.lex( i, j - 1 )] +
					pv1[field_.lex( i, j + 1 )] ) / 3;
			}
			for( int i = field_.Ny() % 2 == 0 ? 1 : 2; i < field_.Nx() - 1; i += 2 ) {
				const int j = field_.Ny() - 1;
				pv2[field_.lex( i, j )] = (
					pv1[field_.lex( i - 1, j )] +
					pv1[field_.lex( i + 1, j )] +
					pv1[field_.lex( i, j - 1 )] ) / 3;
			}
			for( int j = field_.Nx() % 2 == 0 ? 1 : 2; j < field_.Ny() - 1; j += 2 ) {
				const int i = field_.Nx() - 1;
				pv2[field_.lex( i, j )] = (
					pv1[field_.lex( i - 1, j )] +
					pv1[field_.lex( i, j - 1 )] +
					pv1[field_.lex( i, j + 1 )] ) / 3;
			}
			// corners red
			pv2[field_.lex( 0, 0 )] = ( pv1[field_.lex( 0 + 1, 0 )] + pv1[field_.lex( 0, 0 + 1 )] ) / 2;
			if( field_.Ny() % 2 == 1 )
				pv2[field_.lex( field_.Nx() - 1, 0 )] = ( pv1[field_.lex( field_.Nx() - 1 - 1, 0 )] + pv1[field_.lex( field_.Nx() - 1, 0 + 1 )] ) / 2;
			if( field_.Ny() % 2 == 1 )
				pv2[field_.lex( 0, field_.Ny() - 1 )] = ( pv1[field_.lex( 0 + 1, field_.Ny() - 1 )] + pv1[field_.lex( 0, field_.Ny() - 1 - 1 )] ) / 2;
			if( ( field_.Ny() + field_.Nx() ) % 2 == 0 )
				pv2[field_.lex( field_.Nx() - 1, field_.Ny() - 1 )] = ( pv1[field_.lex( field_.Nx() - 1 - 1, field_.Ny() - 1 )] + pv1[field_.lex( field_.Nx() - 1, field_.Ny() - 1 - 1 )] ) / 2;

		#pragma omp parallel for
			for( int idx = 0; idx < interiorBlackSize; ++idx ) {
				const int i = interiorBlack[idx][0];
				const int j = interiorBlack[idx][1];
				gaussSeidelStepInterior( field_, i, j, pv2, pv2 );
			}
			// edges black
			for( int i = 1; i < field_.Nx() - 1; i += 2 ) {
				const int j = 0;
				pv2[field_.lex( i, j )] = (
					pv2[field_.lex( i - 1, j )] +
					pv2[field_.lex( i + 1, j )] +
					pv2[field_.lex( i, j + 1 )] ) / 3;
			}
			for( int j = 1; j < field_.Ny() - 1; j += 2 ) {
				const int i = 0;
				pv2[field_.lex( i, j )] = (
					pv2[field_.lex( i + 1, j )] +
					pv2[field_.lex( i, j - 1 )] +
					pv2[field_.lex( i, j + 1 )] ) / 3;
			}
			for( int i = field_.Ny() % 2 == 1 ? 1 : 2; i < field_.Nx() - 1; i += 2 ) {
				const int j = field_.Ny() - 1;
				pv2[field_.lex( i, j )] = (
					pv2[field_.lex( i - 1, j )] +
					pv2[field_.lex( i + 1, j )] +
					pv2[field_.lex( i, j - 1 )] ) / 3;
			}
			for( int j = field_.Nx() % 2 == 1 ? 1 : 2; j < field_.Ny() - 1; j += 2 ) {
				const int i = field_.Nx() - 1;
				pv2[field_.lex( i, j )] = (
					pv2[field_.lex( i - 1, j )] +
					pv2[field_.lex( i, j - 1 )] +
					pv2[field_.lex( i, j + 1 )] ) / 3;
			}
			// corners black
			if( field_.Ny() % 2 == 0 )
				pv2[field_.lex( field_.Nx() - 1, 0 )] = ( pv2[field_.lex( field_.Nx() - 1 - 1, 0 )] + pv2[field_.lex( field_.Nx() - 1, 0 + 1 )] ) / 2;
			if( field_.Ny() % 2 == 0 )
				pv2[field_.lex( 0, field_.Ny() - 1 )] = ( pv2[field_.lex( 0 + 1, field_.Ny() - 1 )] + pv2[field_.lex( 0, field_.Ny() - 1 - 1 )] ) / 2;
			if( ( field_.Ny() + field_.Nx() ) % 2 == 1 )
				pv2[field_.lex( field_.Nx() - 1, field_.Ny() - 1 )] = ( pv2[field_.lex( field_.Nx() - 1 - 1, field_.Ny() - 1 )] + pv2[field_.lex( field_.Nx() - 1, field_.Ny() - 1 - 1 )] ) / 2;

		#pragma omp parallel for
			for( int i = 0; i < v1Size; ++i ) {
				vDiff[i] = std::fabs( pv2[i] - pv1[i] );
			}
			auto maxIt = std::max_element( vDiff.begin(), vDiff.end() );
			std::swap( pv1, pv2 );

			if( *maxIt <= 1e-5 * pv1[std::distance( vDiff.begin(), maxIt )] ) {
				break;
			}
		}
	}

	void fillGridFromNearestNeighbors( TriMesh& mesh ) {
		std::vector<float> vertexAvg( mesh.n_vertices() );
		for( const auto& v : mesh.vertices() ) {
			float avg{ 0.f };
			auto nEdges{ 0 };

			for( const auto& voh : v.outgoing_halfedges() ) {
				++nEdges;
				auto p1 = mesh.point( voh.from() );
				auto p2 = mesh.point( voh.to() );
				p1[2] = 0;
				p2[2] = 0;
				avg += ( p1 - p2 ).length();
			}

			vertexAvg[v.idx()] = avg / (float)nEdges;
		}

		MeshPointsNN knn( mesh );
		for( int i = 0; i < field_.Nx(); ++i ) {
			for( int j = 0; j < field_.Ny(); ++j ) {
				const auto& pos = field_.pos( i, j );
				TriMesh::Point p{ pos.x, pos.y, 0 };

				const auto nnIdxs = knn.findNearestNeighbor( p, 5 );

				float invDistSum = 0;
				float avg = 0;
				for( const auto& idx : nnIdxs ) {
					const auto& pn = mesh.point( mesh.vertex_handle( idx ) );
					auto dist = ( pn - p ).sqrnorm();
					auto invDist = 1.f / ( dist + 0.0001f );
					auto val = vertexAvg[idx];
					avg += val * invDist;
					invDistSum += invDist;
				}
				avg /= invDistSum;
				field_( i, j ) = avg;
			}
		}
	}

	static ScalarField::AxisAlignedBoundingBox computeAabbLecacy( TriMesh& mesh ) {
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

	static ScalarField::AxisAlignedBoundingBox computeAabb( TriMesh& mesh ) {
		auto xMax = -std::numeric_limits<float>::max();
		auto yMax = -std::numeric_limits<float>::max();
		auto xMin = std::numeric_limits<float>::max();
		auto yMin = std::numeric_limits<float>::max();
		for( const auto& v : mesh.vertices() ) {
			const auto& point = mesh.point( v );
			xMax = std::max( xMax, point[0] );
			yMax = std::max( yMax, point[1] );
			xMin = std::min( xMin, point[0] );
			yMin = std::min( yMin, point[1] );
		}
		return { xMin, xMax, yMin, yMax };
	}

public:
	// load size field if it already exists in cache, otherwise generate and store it.
	static ScalarField::ScalarField load( TriMesh& mesh, const std::experimental::filesystem::path& cacheFolder, const std::experimental::filesystem::path& meshFile, const size_t& sizeGridSizeX, const size_t& sizeGridSizeY ) {
		namespace fs = std::experimental::filesystem;

		fs::path file = cacheFolder / ( meshFile.stem().string() + "_SizeField_" + std::to_string( sizeGridSizeX ) + "_" + std::to_string( sizeGridSizeY ) + ".bin" );

		auto sfLoad = ScalarField::load( file, sizeGridSizeX, sizeGridSizeY );
		//std::optional<ScalarField::ScalarField> sfLoad = {};

		if( sfLoad ) {
			return sfLoad.value();
		} else {
			LOG( INFO ) << "Generate SizeField";
			SizeField sf( mesh, sizeGridSizeX, sizeGridSizeY, true );
			LOG( INFO ) << "Store in cache file '" << file << "'";
			sf.field().writeBinary( file.string() );
			return sf.field();
		}
	}

};