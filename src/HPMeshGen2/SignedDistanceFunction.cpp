#include "SignedDistanceFunction.h"

#include "Contour/Shape.h"

SignedDistanceFunction::SignedDistanceFunction( TriMesh& mesh, const size_t& Nx, const size_t& Ny ) : field_(computeAABB(mesh), Nx, Ny) {
	Contour::Shape shape( mesh );

	const auto nx = field_.Nx();
	const auto ny = field_.Ny();

	// initialize all values
	#pragma omp parallel for
	for( int i = 0; i < nx; ++i ) {
		for( int j = 0; j < ny; ++j ) {
			auto fieldPos = field_.pos( i, j );
			TriMesh::Point p{ fieldPos[0], fieldPos[1], 0 };
			//field_( i, j ) = shape.signedDistanceToNearestNeighbor( p );
			field_( i, j ) = shape.signedDistanceExact( p );
		}
	}


	//for( int i = 0; i < nx; ++i ) {
	//	for( int j = 0; j < ny; ++j ) {
	//		field_( i, j ) = std::numeric_limits<float>::max();
	//	}
	//}
	//
	//for( const auto& e : shape.edges() ) {
	//	#pragma omp parallel for
	//	for( int i = 0; i < nx; ++i ) {
	//		for( int j = 0; j < ny; ++j ) {
	//			auto fieldPos = field_.pos( i, j );
	//			TriMesh::Point p{ fieldPos[0], fieldPos[1], 0 };
	//			const auto s = shape.signedDistanceToEdge( e, p );
	//			if( std::abs( s ) < std::abs( field_( i, j ) ) ) {
	//				field_( i, j ) = s;
	//			}
	//		}
	//	}
	//	std::cout << e.idx() << std::endl;
	//}

}
