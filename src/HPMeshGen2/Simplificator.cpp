#include "Simplificator.h"

#define Q_MEANRATIO 0
#define Q_JACOBIAN 4
#define Q_MINANGLE 5
#define Q_RADIUSRATIO 6
#define Q_MAXANGLE 7

#include "MeshFunctions.h"

#include "Simplify.h"


template<typename KeyType, typename ValueType>
inline std::pair<KeyType, ValueType> get_min( const std::map<KeyType, ValueType>& x ) {
	using pairtype = std::pair<KeyType, ValueType>;
	return *std::max_element( x.begin(), x.end(), []( const pairtype & p1, const pairtype & p2 ) {
		return p1.second > p2.second;
							  } );
}


inline bool isPointAllowed( TriMesh& mesh, const TriMesh::VertexHandle& vh1, const TriMesh::VertexHandle& vh2, const TriMesh::Point& p_new ) {
	const float min_q = 0.1f;

	std::vector<TriMesh::VertexHandle> vh1fan;
	MeshFunctions::getOneRing( mesh, vh1, vh1fan );

	for ( size_t i = 0; i < vh1fan.size() - 1; ++i ) {
		size_t j = i + 1;
		if ( vh1fan[i] == vh2 || vh1fan[j] == vh2 ) {
			continue;
		}
		TriMesh::Point p1 = mesh.point( vh1fan[i] );
		TriMesh::Point p2 = mesh.point( vh1fan[j] );

		const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

		if ( quality < min_q )
			return false;

	}

	std::vector<TriMesh::VertexHandle> vh2fan;
	MeshFunctions::getOneRing( mesh, vh2, vh2fan );

	for ( size_t i = 0; i < vh2fan.size() - 1; ++i ) {
		size_t j = i + 1;
		if ( vh2fan[i] == vh1 || vh2fan[j] == vh1 ) {
			continue;
		}
		TriMesh::Point p1 = mesh.point( vh2fan[i] );
		TriMesh::Point p2 = mesh.point( vh2fan[j] );

		const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

		if ( quality < min_q )
			return false;

	}

	return true;
}




bool Simplificator::setCache( decltype( cacheFolder_ ) cacheFolder, decltype( meshFile_ ) meshFile ) {
	namespace fs = std::experimental::filesystem;

	cacheFolder_ = cacheFolder;
	meshFile_ = meshFile;

	if ( !fs::exists( cacheFolder_ ) ) {
		LOG( INFO ) << "Cache directory '" << cacheFolder_ << "' does not exist yet and is therefore now created";
		if ( !fs::create_directories( cacheFolder_ ) ) {
			LOG( WARNING ) << "Directory '" << cacheFolder_ << "' cannot be created!" << std::endl;
			return false;
		}
	}

	useCache_ = true;
	return useCache_;
}

bool Simplificator::simplificationOutput( decltype( outputFolder_ ) outputFolder ) {
	namespace fs = std::experimental::filesystem;

	outputFolder_ = outputFolder;

	if ( !fs::exists( outputFolder_ ) ) {
		LOG( INFO ) << "Output directory '" << outputFolder_ << "' does not exist yet and is therefore now created";
		if ( !fs::create_directories( outputFolder_ ) ) {
			LOG( WARNING ) << "Directory '" << outputFolder_ << "' cannot be created!" << std::endl;
			return false;
		}
	}

	printSimplificationOutput_ = true;
	return printSimplificationOutput_;
}


TriMesh Simplificator::simplify( const ScalarField::ScalarField& backgroundgrid, const size_t &nPatches, const bool convHullSimpl ) {
	namespace fs = std::experimental::filesystem;

	OceanMesh meshRet;

	fs::path patchStructureFile;
	if( !convHullSimpl ) {
		patchStructureFile = cacheFolder_ / ( meshFile_.stem().string() + "_PatchStructure_" + std::to_string( nPatches ) + ".om" );
	} else {
		patchStructureFile = cacheFolder_ / ( meshFile_.stem().string() + "_PatchStructureConvHull_" + std::to_string( nPatches ) + ".om" );
	}

	if ( useCache_ && fs::exists( patchStructureFile ) ) {
		LOG( INFO ) << "Read PatchStructure from cache file '" << fs::canonical( patchStructureFile ) << "'";
		OceanMesh meshRet_tmp( patchStructureFile.string() );
        meshRet = meshRet_tmp;
	} else {
		LOG( INFO ) << "Generate patch structure";
		meshRet = mesh_;
		
		//elbReductionFast( meshRet, backgroundgrid, nPatches * 1.5 );
		//elbReduction( meshRet, backgroundgrid, nPatches );
		meshRet.request_vertex_status();
		
		if( !convHullSimpl ) {
			OpenMesh::Decimater::OceanDecimaterT<TriMesh> decimater( meshRet );
			OpenMesh::Decimater::ModPointDistT<TriMesh>::Handle pointDistMod;
			decimater.add( pointDistMod );
			decimater.module( pointDistMod ).set_backgroundgrid( &backgroundgrid );
			decimater.module( pointDistMod ).set_min_meanratio( 0.1 );

			decimater.initialize();
			decimater.decimate_to_faces( 0, nPatches );
			meshRet.garbage_collection();
		} else {
			// fill small interior boundaries
			for( auto heh : meshRet.halfedges() ) {
				if( !heh.is_boundary() )
					continue;
				if( heh.prev().from() == heh.next().to() ) {
					std::vector<OpenMesh::SmartVertexHandle> vhs{ heh.prev().from(), heh.from(), heh.to() };
					std::vector<TriMesh::Point> pVec( 3 );
					std::transform( vhs.begin(), vhs.end(), pVec.begin(), [&meshRet]( const auto& v ) { return meshRet.point( v ); } );
					auto a = HelperFunctions::calcAngle2( pVec[2], pVec[1], pVec[0] );
					if( a < 180 )
						meshRet.add_face( vhs );
				}
			}
			OpenMesh::Decimater::OceanDecimater2T decimater( meshRet );
			OpenMesh::Decimater::ModConvHullPointDist2T<TriMesh>::Handle convHullPointDistMod;
			decimater.add( convHullPointDistMod );
			decimater.module( convHullPointDistMod ).set_backgroundgrid( &backgroundgrid );
			decimater.module( convHullPointDistMod ).set_min_meanratio( 0.01 );

			decimater.initialize();
			decimater.decimate_to_faces( 0, nPatches );
			meshRet.garbage_collection();
		}
		
		
		if ( useCache_ ) {
			LOG( INFO ) << "Store patch structure in cache: '" << patchStructureFile << "'";
			meshRet.property( meshRet.feature_size ).set_persistent( true );
			OpenMesh::IO::write_mesh( meshRet, patchStructureFile.string() );
		}
	}
	return meshRet;
}
