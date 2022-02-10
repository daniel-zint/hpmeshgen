#include "BsgGenerator.h"

#include <experimental/filesystem>
#include "MeshFunctions.h"
#include "RatRace/RatRace.h"
#include "TriRemeshing.h"
#include "SignedDistanceFunction.h"
#include "DMO/Solver.h"

#include "GridEvaluation/eval.h"

void BsgGenerator::generate() {
	
	processInput();

	if( config_.initialMeshOutput() ) {
		MeshFunctions::printMesh( config_.outputFolder() / "initialMesh.off", inputMesh_ );
		MeshFunctions::printMeshCartesian( config_.outputFolder() / "initialMesh_cart.off", inputMesh_, inputMesh_ );
	}

	simplify();
	remeshSimplifiedMesh();
	// fragment mesh
	generateFragments();
	
	// BSG
	BlockStructuredGrid bsg_tmp( fragmentMesh_, config_.nRefinementSteps(), contour_, inputMesh_ );
    bsg_ = bsg_tmp;
    
	LOG( INFO ) << "# quads after refinement: " << bsg_.quadMesh.n_faces() << std::endl;
}

void BsgGenerator::processInput() {
	GridEvaluation::eval( inputMesh_, inputMesh_, config_.outputFolder() / "cflInputQuotient.vtk" );
	inputMesh_.loadBackgroundGrids( config_.cacheFolder(), config_.meshFile(), config_.sizegridSizeX(), config_.sizegridSizeY() );
	inputMesh_.flatten();
}

void BsgGenerator::simplify() {

	simplificator_.setCache( config_.cacheFolder(), config_.meshFile() );
	simplifiedMesh_ = simplificator_.simplify( inputMesh_.sizeField(), 2 * config_.nPatches(), config_.convexHullDecimation() );
	//MeshFunctions::printMesh( config_.outputFolder() / "debSimpl.off", simplifiedMesh_ );

	if( !config_.convexHullDecimation() ) {
		//contourOld_.mapBoundaryToContour( simplifiedMesh_ );

		for( const auto& vh : simplifiedMesh_.vertices() ) {
			if( !vh.is_boundary() )
				continue;
		
			auto p = simplifiedMesh_.point( vh );
			auto pNew = contour_.mapToContour( p );
			simplifiedMesh_.set_point( vh, pNew );
		}
		//MeshFunctions::printMesh( config_.outputFolder() / "debSimplCont.off", simplifiedMesh_ );
	}

	LOG( INFO ) << "# coarse triangles: " << simplifiedMesh_.n_faces();
	if( !config_.convexHullDecimation() )
		optimizeSimplifiedMesh();
}

void BsgGenerator::remeshSimplifiedMesh() {
	
	namespace fs = std::experimental::filesystem;
	
	constexpr bool useCache = true;

	fs::path remeshedFile;
	remeshedFile = ( config_.cacheFolder() / config_.meshFile().stem() );
	remeshedFile += "_SimplifiedRemeshed";
	if( config_.convexHullDecimation() ) {
		remeshedFile += "ConvHull";
	}
	if( config_.forceFragmentNumber() ) {
		remeshedFile += "ForceFragments";
	}
	remeshedFile += std::string( "_" ) + std::to_string( 2 * config_.nPatches() ) + ".om";

	if( useCache && fs::exists( remeshedFile ) ) {
		LOG( INFO ) << "Read simplified and remeshed file from cache file '" << fs::canonical( remeshedFile ) << "'";
		OpenMesh::VPropHandleT<bool> is_feature;
		simplifiedMesh_.add_property( is_feature, "is_feature" );

		//meshRet = OceanMesh( patchStructureFile.string() );
		OpenMesh::IO::read_mesh( simplifiedMesh_, remeshedFile.string() );

		for( const auto& vh : simplifiedMesh_.vertices() ) {
			simplifiedMesh_.data( vh ).is_feature = simplifiedMesh_.property( is_feature, vh );
		}
		simplifiedMesh_.remove_property( is_feature );
	} else {
		LOG( INFO ) << "Remesh simplified file";
		
		TriRemeshing triRem( simplifiedMesh_ );
		triRem.set_backgroundgrid( &inputMesh_.sizeField() );
		if( !config_.convexHullDecimation() ) {
			triRem.set_contour( &contour_ );
		} else {
			triRem.use_convex_hull_remeshing();
		}
		if( config_.forceFragmentNumber() ) {
			simplifiedMesh_ = triRem.remesh( 2 * config_.nPatches() );
		} else {
			simplifiedMesh_ = triRem.remesh();
		}
		
		if( useCache ) {
			//OpenMesh::EPropHandleT<bool> isNeumann;
			//simplifiedMesh_.get_property_handle( isNeumann, "isNeumann" );
			//OpenMesh::EPropHandleT<int> boundaryID;
			//simplifiedMesh_.get_property_handle( boundaryID, "boundaryID" );
			//OpenMesh::VPropHandleT<int> contourIdx;
			//simplifiedMesh_.get_property_handle( contourIdx, "contourIdx" );
			//OpenMesh::VPropHandleT<float> feature_size;
			//simplifiedMesh_.get_property_handle( feature_size, "feature_size" );
			//
			OpenMesh::VPropHandleT<bool> is_feature;
			simplifiedMesh_.add_property( is_feature, "is_feature" );
			for( const auto& vh : simplifiedMesh_.vertices() ) {
				simplifiedMesh_.property( is_feature, vh ) = simplifiedMesh_.data( vh ).is_feature;
			}
			
			//LOG( INFO ) << "Store remeshed mesh in cache: '" << remeshedFile << "'";
			//simplifiedMesh_.property( contourIdx ).set_persistent( true );
			//simplifiedMesh_.property( isNeumann ).set_persistent( true );
			//simplifiedMesh_.property( boundaryID ).set_persistent( true );
			//simplifiedMesh_.property( feature_size ).set_persistent( true );
			simplifiedMesh_.property( is_feature ).set_persistent( true );
			OpenMesh::IO::write_mesh( simplifiedMesh_, remeshedFile.string() );
			simplifiedMesh_.remove_property( is_feature );
		}
		
	}
	
}

/// <summary>
/// DEPRECATED use generateFragments() instead
/// </summary>
void BsgGenerator::generateQuadFragments() {
	fragmentMesh_ = MeshFunctions::tri2poly( simplifiedMesh_ );
	MeshFunctions::tri2quad( fragmentMesh_ );
	MeshFunctions::catmullClark2D( fragmentMesh_ );

	//MeshFunctions::fitBoundaryToContour( simplifiedMesh_, fragmentMesh_, contourOld_ );	// removed contourOld
	
	LOG( INFO ) << "# quad patches: " << fragmentMesh_.n_faces() << std::endl;
	discreteMeshOptimization( fragmentMesh_, Q_MEANRATIO, 0.5f, 100 );
}


inline void evenTriangleNumber( BsgGenerator& bsgGenerator, const Config& config ) {
	if ( bsgGenerator.simplifiedMesh().n_faces() % 2 == 1 ) {
		auto& mesh = bsgGenerator.simplifiedMesh();
		auto& bg = bsgGenerator.inputMesh().sizeField();

		// get edge which is relatively the longest
		TriMesh::EdgeHandle ehLongest;
		float relativeLength = -FLT_MAX;
		for ( auto eh : mesh.edges() ) {
			if ( !mesh.is_boundary( eh ) ) continue;
			auto heh = mesh.halfedge_handle( eh, 0 );
			auto vh1 = mesh.from_vertex_handle( heh );
			auto vh2 = mesh.to_vertex_handle( heh );
			auto p1 = mesh.point( vh1 );
			auto p2 = mesh.point( vh2 );

			auto s = bg.getScalar( { p1[0],p1[1] }, { p2[0], p2[1] } );
			auto l = ( p2 - p1 ).length();
			auto rl = l / s;
			if ( rl > relativeLength ) {
				relativeLength = rl;
				ehLongest = eh;
			}
		}
		auto heh = mesh.halfedge_handle( ehLongest, 0 );
		auto vh1 = mesh.from_vertex_handle( heh );
		auto vh2 = mesh.to_vertex_handle( heh );
		auto p1 = mesh.point( vh1 );
		auto p2 = mesh.point( vh2 );

		// perform edge split
		auto vhNew = mesh.split( ehLongest, 0.5f* ( p2 + p1 ) );
	}
}

void BsgGenerator::generateFragments() {
	evenTriangleNumber( ( *this ), config_ );
	//optimizeSimplifiedMesh();
	//discreteMeshOptimizationDensity<TriMesh>( simplifiedMesh_, inputMesh_.sizegrid, 0.5f, 100 );
	//MeshFunctions::printMesh( config_.outputFolder() / "simplifiedMesh_.off", simplifiedMesh_ );

	RatRace rr( simplifiedMesh_, true );
	//MeshFunctions::printMesh( config_.outputFolder() / "hmesh_.off", rr.hmesh() );
	rr.run( false );
	
	fragmentMesh_ = rr.hmesh();
		
	for( const auto& vh : simplifiedMesh_.vertices() ) {
		auto fragVh = fragmentMesh_.vertex_handle( vh.idx() );
		fragmentMesh_.data( fragVh ).is_feature = simplifiedMesh_.data( vh ).is_feature;
	}

	LOG( INFO ) << "# quad patches: " << fragmentMesh_.n_faces() << std::endl;
	//MeshFunctions::printMesh( config_.outputFolder() / "fragmentMesh_.off", fragmentMesh_ );
	//MeshFunctions::printMesh( config_.outputFolder() / "simplifiedMesh_.off", simplifiedMesh_ );

	// fix bad elements on boundary (at least for convexHullDecimation)
	if( config_.convexHullDecimation() ) {
		for( const auto& v : fragmentMesh_.vertices() ) {
			if( !v.is_boundary() )
				continue;

			if( v.valence() > 2 ) 
				continue;
			std::vector<OpenMesh::SmartFaceHandle> faces;
			for( const auto& f : v.faces() ) {
				faces.push_back( f );
			}
			LOG_ASSERT( faces.size() == 1 );
			auto q = QualityMetrics::condition( fragmentMesh_, faces[0] );
			if( q > 0 )
				continue;

			// set vertex in between its neighbors
			const auto& vNext = v.halfedge().to();
			const auto& vPrev = v.halfedge().prev().from();
			const auto& pNext = fragmentMesh_.point( vNext );
			const auto& pPrev = fragmentMesh_.point( vPrev );

			if( vNext.valence() == 2 || vPrev.valence() == 2 )
				continue;

			fragmentMesh_.set_point( v, 0.5f * ( pNext + pPrev ) );
		}
	}
}

void BsgGenerator::optimizeSimplifiedMesh() {
	for ( size_t i = 0; i < 10; ++i ) {
		MeshFunctions::boundarySmoothing( simplifiedMesh_, contour_ );
		//discreteMeshOptimizationDensity<TriMesh>( simplifiedMesh_, inputMesh_.sizegrid, 0.5f, 100 );		// causes race conditions
		MeshFunctions::edgeFlip( simplifiedMesh_ );
		//MeshFunctions::printMesh( config_.outputFolder() / "deb.off", simplifiedMesh_ );
	}
}

void BsgGenerator::adaptElementsToContour() {
	constexpr bool fixContour = true;

	//MeshFunctions::printMesh( config_.outputFolder() / "deb0.off", bsg_.triMesh );
	auto triComplete = bsg_.triMesh;

	auto& sizefield = inputMesh_.sizeField();
	DMO::UniformGrid grid_d;
	grid_d.n = { (int)sizefield.Nx(), (int)sizefield.Ny() };
	grid_d.h = { sizefield.hx(), sizefield.hy() };
	grid_d.aabbMin = { sizefield.aabb().xMin, sizefield.aabb().yMin };
	grid_d.aabbMax = { sizefield.aabb().xMax, sizefield.aabb().yMax };
	thrust::host_vector<float> sgVals;
	sgVals.resize( grid_d.n.x * grid_d.n.y );
	for( int i = 0; i < grid_d.n.x * grid_d.n.y; ++i ) {
		sgVals[i] = sizefield( i );
	}
	thrust::device_vector<float> sgVals_d = sgVals;
	grid_d.vals = sgVals_d.data().get();

	
	DMO::Metrics::MeanRatioTriangle metricMeanRatio;
	DMO::Metrics::DensityTriangle metricDensity( grid_d );
	DMO::Metrics::DensityWithMeanRatioTriangle metricInner( metricMeanRatio, metricDensity );

	//ScalarField::ScalarField sdf = SignedDistanceFunction( inputMesh_, config_.sizegridSizeX(), config_.sizegridSizeY() ).field();
	ScalarField::ScalarField sdf = SignedDistanceFunction( inputMesh_, config_.cacheFolder(), config_.meshFile(), config_.sizegridSizeX(), config_.sizegridSizeY() ).field();
	sdf.toVTK( config_.outputFolder() / "sdf.vtk", true );

	DMO::UniformGrid sdf_d;
	sdf_d.n = { (int)sdf.Nx(), (int)sdf.Ny() };
	sdf_d.h = { sdf.hx(), sdf.hy() };
	sdf_d.aabbMin = { sdf.aabb().xMin, sdf.aabb().yMin };
	sdf_d.aabbMax = { sdf.aabb().xMax, sdf.aabb().yMax };
	thrust::host_vector<float> valsSdf;
	valsSdf.resize( sdf_d.n.x * sdf_d.n.y );
	for( int j = 0; j < sdf_d.n.y; ++j ) {
		for( int i = 0; i < sdf_d.n.x; ++i ) {
			valsSdf[j + sdf_d.n.y * i] = sdf( i, j );
		}
	}
	thrust::device_vector<float> valsSdf_d = valsSdf;
	sdf_d.vals = valsSdf_d.data().get();
	DMO::Metrics::SdfTriangle sdfTriangle( sdf_d );
	DMO::Metrics::SdfConstrainedMeanRatioTriangle sdfcmrm( metricMeanRatio, sdfTriangle );

	DMO::DmoMesh dmoMeshInner = DMO::DmoMesh::create<DMO::Set::Inner>( bsg_.triMesh );
	std::vector<OpenMesh::SmartVertexHandle> boundaryVertices;
	for( const auto& vh : bsg_.triMesh.vertices() ) {
		if( vh.is_boundary() && !bsg_.triMesh.data( vh ).is_feature ) {
			boundaryVertices.push_back( vh );
		}
	}
	DMO::DmoMesh dmoMeshBound( boundaryVertices );
	if constexpr( fixContour ) {
		DMO::Solver dmo( bsg_.triMesh, &metricInner, &dmoMeshInner, &sdfcmrm, &dmoMeshBound );
		dmo.solve(100);
	} else {
		DMO::Solver dmo( bsg_.triMesh, &metricInner, &dmoMeshInner );
		dmo.solve( 100 );
	}
	bsg_.copy_positions_tri2quad();

	//MeshFunctions::printMesh( config_.outputFolder() / "deb1.off", bsg_.triMesh );

	/////////////////////////////////////////////////////////////////
	
	auto& m = bsg_.triMesh;
	//{
	//	// HACK set a vertex on open sea boundary s.th. it is removed
	//	auto v = m.vertex_handle( 5413 );
	//	auto p = m.point( v );
	//	p[0] = -55;
	//	m.set_point( v, p );
	//	bsg_.copy_positions_tri2quad();
	//}

	OpenMesh::VPropHandleT<float> vsdf;
	m.add_property( vsdf );

	for( const auto& v : m.vertices() ) {
		const auto& p = m.point( v );
		m.property( vsdf, v ) = sdf.getScalar( p[0], p[1] );
	}

	//MeshFunctions::printMesh(config_.outputFolder() / "debbeforeremove.off", m);
	m.request_face_status();
	m.request_edge_status();
	m.request_vertex_status();
	m.request_halfedge_status();

	// remove triangles that have a certain percentage outside of the boundary
	for( const auto& f : m.faces() ) {
		bool del = true;
		std::vector<TriMesh::Point> points;
		points.reserve( 3 );
		for( auto v : f.vertices() ) {
			if( m.property( vsdf, v ) <= 0 || ( bsg_.triMesh.data( v ).is_feature && v.valence() == 2) ) {
				del = false;
				//break;
			}
			points.push_back( m.point( v ) );
		}

		//triangle is not compeletly outside the boundary so look at the area
		if( !del ) {
			std::vector<TriMesh::VertexHandle> triangle_vh;
			triangle_vh.reserve( 3 );
			std::vector<TriMesh::VertexHandle> triangle_inside;
			triangle_inside.reserve( 3 );

			int vertices_inside = 0; //counts vertices inside the boundary
			for (const auto& v : f.vertices()) {
				triangle_vh.push_back( v );
				if (m.property(vsdf, v) <= 0 ) {
					vertices_inside++;
					triangle_inside.push_back( v );
				}
			}
			//triangle has only one vertex inside and could potentially be removed
			if (vertices_inside == 1) {

				//pushes the vertices into triangle_inside so that the vertex inside the boundary is at position 0
				for( const auto& v : triangle_vh ) {
					if( v != triangle_inside[0] ) {
						triangle_inside.push_back( v );
					}
				}

				//linear interpolation to find the point where the edge meets the boundary
				float lerp_param_1 = m.property( vsdf, triangle_inside[0] ) / ( m.property( vsdf, triangle_inside[0] ) - m.property( vsdf, triangle_inside[1] ) );
				float lerp_param_2 = m.property( vsdf, triangle_inside[0] ) / ( m.property( vsdf, triangle_inside[0] ) - m.property( vsdf, triangle_inside[2] ) );
					
				auto intersect_1 = ( 1.0f - lerp_param_1 ) * m.point( triangle_inside[0] ) + lerp_param_1 * m.point( triangle_inside[1] );
				auto intersect_2 = ( 1.0f - lerp_param_2 ) * m.point( triangle_inside[0] ) + lerp_param_2 * m.point( triangle_inside[2] );


				////calculate the areas inside and outside the boundary
				//double area_small = 0.5 * ((m.point( triangle_inside[0] ) - value_1) % (m.point( triangle_inside[0] ) - value_2)).norm();
				//
				//double area_big = 0.5 * ((points[0] - points[1]) % (points[0] - points[2])).norm();
				//
				//double ratio = area_small / area_big;

				//caluclate the smallest percentage of edge inside the boundary
				float ratio_1 = ( ( intersect_1 - m.point( triangle_inside[0] ) ).norm() ) / ( ( m.point( triangle_inside[1] ) - m.point( triangle_inside[0] ) ).norm() );
				float ratio_2 = ( ( intersect_2 - m.point( triangle_inside[0] ) ).norm() ) / ( ( m.point( triangle_inside[2] ) - m.point( triangle_inside[0] ) ).norm() );
				
				float ratio = std::min( ratio_1, ratio_2 );
				
				const float ratio_tolerance = 0.85; //tolerance for how many percent of a triangle can lie inside the boundary

				if( ratio <= 1.0f - ratio_tolerance ) {
					del = true;
				}

			}
		}

			
		if( del ) {
			// delete triangle if its midpoint is also outside the domain
			auto pMid = ( 1.f / 3.f ) * ( points[0] + points[1] + points[2] );

			if( sdf.getScalar( pMid[0], pMid[1] ) > 0 )
				m.delete_face( f );
		}
	}

	// remove boundary vertices if their neighbors approximate the boundary better
	bool smthDeleted = true;
	while( smthDeleted ) {
		smthDeleted = false;
		for( const auto& v : m.vertices() ) {
			if( !v.is_boundary() || bsg_.triMesh.data( v ).is_feature )
				continue;
			float vs = m.property( vsdf, v );
			if( vs <= 0 )
				continue;
			float sNeighMin = 0;
			for( const auto& vv : v.vertices() ) {
				const float& vvs = m.property( vsdf, vv );
				sNeighMin = std::min( sNeighMin, vvs );
			}
			if( vs + sNeighMin <= 0 )
				continue;

			// do not 'destroy' land bridges
			int nBoundaryNeighs = 0;
			for( const auto& vv : v.vertices() ) {
				if( vv.is_boundary() ) ++nBoundaryNeighs;
			}
			if( nBoundaryNeighs > 2 )
				continue;

			// if all vertices of a triangle are on the boundary, check also for cog
			if( v.valence() == 2 ) {
				TriMesh::Point p = m.point( v );
				for( const auto& vv : v.vertices() ) {
					p += m.point( vv );
				}
				p /= 3;
				auto s = sdf.getScalar( p[0], p[1] );
				sNeighMin = std::min( sNeighMin, s );
			}
			if( vs + sNeighMin > 0 ) {
				m.delete_vertex( v );
				smthDeleted = true;
			}
		}
	}

	// remove non-manifold vertices
	bool removedNonManifoldVertex = false;
	do {
		removedNonManifoldVertex = false;
		for( const auto& v : m.vertices() ) {
			int c = 0;
			for( const auto& e : v.edges() ) {
				if( e.is_boundary() ) ++c;
			}
			if( c == 0 || c == 2 ) {
				continue;
			}

			m.delete_vertex( v );
			removedNonManifoldVertex = true;
		}
	} while( removedNonManifoldVertex );
	
	//m.garbage_collection();
	//MeshFunctions::printMesh( config_.outputFolder() / "deb2.off", m );

	m.remove_property( vsdf );

	if constexpr( fixContour ) {

		// get boundary back on contour
		dmoMeshInner = DMO::DmoMesh::create<DMO::Set::Inner>( bsg_.triMesh );
		boundaryVertices.clear();
		for( const auto& vh : bsg_.triMesh.vertices() ) {
			if( vh.is_boundary() && !bsg_.triMesh.data( vh ).is_feature ) {
				boundaryVertices.push_back( vh );
			}
		}
		
		dmoMeshBound = DMO::DmoMesh( boundaryVertices );
		DMO::Metrics::SdfWithMeanRatioTriangle metricSdfWithMrm( metricMeanRatio, sdfTriangle );
		
		
		DMO::Solver dmo2( bsg_.triMesh, &metricInner, &dmoMeshInner, &metricSdfWithMrm, &dmoMeshBound ); 
		
		dmo2.solve( 100 );
		
		bsg_.copy_positions_tri2quad();
		
		//m.garbage_collection();
		//MeshFunctions::printMesh( config_.outputFolder() / "deb3.off", m );
		
		// smooth boundary
		for( int i = 0; i < 5; ++i ) {
			//MeshFunctions::boundarySmoothingFree( bsg_.triMesh, contourOld_, inputMesh_.sizegrid );
			MeshFunctions::boundarySmoothing( bsg_.triMesh, contour_ );
			DMO::Solver( bsg_.triMesh, &metricMeanRatio, &dmoMeshInner ).solve( 5 );
		}
		bsg_.copy_positions_tri2quad();
		
		// remove boundary triangles that would give low quality if mapped on contour
		for( const auto& f : bsg_.triMesh.faces() ) {
			if( !f.is_boundary() )
				continue;
			bool onlyBoundaryVertices = true;
			std::vector<TriMesh::Point> p;
			p.reserve( 3 );
			for( const auto& fv : f.vertices() ) {
				// do not remove triangles that contain a feature vertex
				if( !fv.is_boundary() || bsg_.triMesh.data(fv).is_feature ) {
					onlyBoundaryVertices = false;
					break;
				}
				p.push_back( bsg_.triMesh.point( fv ) );
			}
			if( !onlyBoundaryVertices )
				continue;
			int nBoundaryEdges = 0;
			for( const auto& fe : f.edges() ) {
				if( fe.is_boundary() ) ++nBoundaryEdges;
			}
			if( nBoundaryEdges < 2 ) {
				continue;
			}
			for( auto& e : p ) {
				e = contour_.mapToContour( e );
			}
			auto q = QualityMetrics::meanRatioMetric( {p[0], p[1], p[2]} );
		
			if( q < 0.3 ) {
				bsg_.triMesh.delete_face( f );
			}
		}
		
		// smooth boundary again
		dmoMeshInner = DMO::DmoMesh::create<DMO::Set::Inner>( bsg_.triMesh );
		for( int i = 0; i < 5; ++i ) {
			MeshFunctions::boundarySmoothing( bsg_.triMesh, contour_ );
			DMO::Solver( bsg_.triMesh, &metricInner, &dmoMeshInner ).solve( 10 );
		}
		DMO::Solver( bsg_.triMesh, &metricMeanRatio, &dmoMeshInner ).solve( 1 );	// adapting to size function can look quite ugly. This resets the vertices a bit.
		bsg_.copy_positions_tri2quad();
	} else {
		//bsg_.triMesh.garbage_collection();
		//MeshFunctions::printMesh( config_.outputFolder() / "deb0.off", bsg_.triMesh );
		
		// set vertices to contour if possible
		OpenMesh::VPropHandleT<bool> onContour;
		m.add_property( onContour );
		for( const auto& v : m.vertices() ) {
			m.property( onContour, v ) = false;
			if( !v.is_boundary() || m.data( v ).is_feature )
				continue;

			const auto pNewOptional = contour_.mapToContourWithNormal( v, m );

			if( !pNewOptional.has_value() )
				continue;

			const auto& pOld = m.point( v );

			const auto& pNew = pNewOptional.value();

			// check validity of new point
			float qMin = std::numeric_limits<float>::max();
			for( const auto& voh : v.outgoing_halfedges() ) {
				if( voh.is_boundary() )
					continue;
				const auto& p1 = m.point( voh.to() );
				const auto& p2 = m.point( voh.next().to() );

				auto q = QualityMetrics::meanRatioMetric( { pNew, p1, p2 } );
				qMin = std::min( q, qMin );

				// do not allow vertices that are far away
				const auto dOld = ( pOld - p1 ).length();
				const auto dNew = ( pNew - p1 ).length();
				if( dNew > 2 * dOld ) {
					qMin = 0;
					break;
				}
			}

			if( qMin > 0.3 ) {
				m.set_point( v, pNew );
				m.property( onContour, v ) = true;
			}
		}

		//MeshFunctions::printMesh( config_.outputFolder() / "deb1.off", bsg_.triMesh );

		// smooth vertices if neighbors are also on contour
		for( const auto& v : m.vertices() ) {
			if( !v.is_boundary() || m.data( v ).is_feature || !m.property( onContour, v ) )
				continue;

			const auto h = v.halfedge();
			LOG_ASSERT( h.is_boundary() );

			if( !m.property( onContour, h.to() ) || !m.property( onContour, h.prev().from() ) ) {
				continue;
			}

			const auto& pMid = m.point( v );
			const auto& pPrev = m.point( h.prev().from() );
			const auto& pNext = m.point( h.to() );
			const auto ePrev = ( pNext - pMid );
			const auto eNext = ( pMid - pPrev );
			OpenMesh::Vec3f nPrev{ ePrev[1], -ePrev[0], 0 };
			nPrev.normalize();
			OpenMesh::Vec3f nNext{ eNext[1], -eNext[0], 0 };
			nNext.normalize();
			if( ( nNext | nPrev ) < 0.3 )
				continue;


			const auto e = ( pNext - pPrev );
			OpenMesh::Vec3f n{ e[1], -e[0], 0 };
			n.normalize();

			const auto& p = 0.5f * ( pPrev + pNext );
			const auto pNew = contour_.mapToContourWithNormal( p, n );
			if( pNew.has_value() ) {
				// check validity of new point
				float qMin = std::numeric_limits<float>::max();
				for( const auto& voh : v.outgoing_halfedges() ) {
					if( voh.is_boundary() )
						continue;
					const auto& p1 = m.point( voh.to() );
					const auto& p2 = m.point( voh.next().to() );

					auto q = QualityMetrics::meanRatioMetric( { pNew.value(), p1, p2 } );
					qMin = std::min( q, qMin );
				}

				if( qMin > 0.3 ) {
					m.set_point( v, pNew.value() );
					m.property( onContour, v ) = true;
				}
			}
		}

		dmoMeshInner = DMO::DmoMesh::create<DMO::Set::Inner>( bsg_.triMesh );
		DMO::Solver( bsg_.triMesh, &metricInner, &dmoMeshInner ).solve( 10 );
		DMO::Solver( bsg_.triMesh, &metricMeanRatio, &dmoMeshInner ).solve( 1 );

		// smooth vertices if neighbors are on contour
		for( const auto& v : m.vertices() ) {
			if( !v.is_boundary() || m.data( v ).is_feature || m.property( onContour, v ) )
				continue;

			const auto h = v.halfedge();
			LOG_ASSERT( h.is_boundary() );

			if( !m.property( onContour, h.to() ) || !m.property( onContour, h.prev().from() ) ) {
				continue;
			}

			const auto& pMid = m.point( v );
			const auto& pPrev = m.point( h.prev().from() );
			const auto& pNext = m.point( h.to() );

			const auto e = ( pNext - pPrev );
			OpenMesh::Vec3f n{ e[1], -e[0], 0 };
			n.normalize();

			const auto& p = 0.5f * ( pPrev + pNext );
			const auto pNew = contour_.mapToContourWithNormal( p, n );
			if( pNew.has_value() ) {
				// check validity of new point
				float qMin = std::numeric_limits<float>::max();
				for( const auto& voh : v.outgoing_halfedges() ) {
					if( voh.is_boundary() )
						continue;
					const auto& p1 = m.point( voh.to() );
					const auto& p2 = m.point( voh.next().to() );

					auto q = QualityMetrics::meanRatioMetric( { pNew.value(), p1, p2 } );
					qMin = std::min( q, qMin );
				}

				if( qMin > 0.3 ) {
					m.set_point( v, pNew.value() );
					m.property( onContour, v ) = true;
				}
			}
		}

		// remove boundary triangles that would give low quality if mapped on contour
		for( const auto& f : m.faces() ) {
			if( !f.is_boundary() )
				continue;
			bool onlyBoundaryVertices = true;
			for( const auto& fv : f.vertices() ) {
				// do not remove triangles that contain a feature vertex
				if( !fv.is_boundary() || bsg_.triMesh.data( fv ).is_feature ) {
					onlyBoundaryVertices = false;
					break;
				}
			}
			if( !onlyBoundaryVertices )
				continue;
			int nBoundaryEdges = 0;
			for( const auto& fe : f.edges() ) {
				if( fe.is_boundary() ) ++nBoundaryEdges;
			}
			if( nBoundaryEdges < 2 ) {
				continue;
			}
			std::vector<TriMesh::Point> p;
			p.reserve( 3 );
			for( const auto& fv : f.vertices() ) {
				auto mapped = contour_.mapToContourWithNormal( fv, m );
				if( mapped.has_value() ) {
					p.push_back( mapped.value() );
				} else {
					p.push_back( m.point( fv ) );
				}
			}
			auto q = QualityMetrics::meanRatioMetric( { p[0], p[1], p[2] } );

			if( q < 0.3 ) {
				bsg_.triMesh.delete_face( f );
			}
		}

		m.remove_property( onContour );

		//MeshFunctions::printMesh( config_.outputFolder() / "deb2.off", bsg_.triMesh );

		dmoMeshInner = DMO::DmoMesh::create<DMO::Set::Inner>( bsg_.triMesh );
		DMO::Solver( bsg_.triMesh, &metricInner, &dmoMeshInner ).solve( 10 );
		DMO::Solver( bsg_.triMesh, &metricMeanRatio, &dmoMeshInner ).solve( 1 );	// adapting to size function can look quite ugly. This resets the vertices a bit.
		bsg_.copy_positions_tri2quad();

		//MeshFunctions::printMesh( config_.outputFolder() / "deb3.off", bsg_.triMesh );
	}

	//bsg_.triMesh.garbage_collection();
	//MeshFunctions::printMesh( config_.outputFolder() / "deb4.off", m );

	bsg_.maskLower().resize( bsg_.quadMesh.n_faces() );
	bsg_.maskUpper().resize( bsg_.quadMesh.n_faces() );
	for( const auto& fh : bsg_.quadMesh.faces() ) {
		auto idxLower = 2 * fh.idx();
		auto idxUpper = 2 * fh.idx() + 1;
		auto fhLower = bsg_.triMesh.face_handle( idxLower );
		auto fhUpper = bsg_.triMesh.face_handle( idxUpper );
		if(bsg_.triMesh.status(fhLower).deleted() )
			bsg_.maskLower()[fh.idx()] = Contour::ContourType::NOTHING;
		else
			bsg_.maskLower()[fh.idx()] = Contour::ContourType::INSIDE;
		
		if( bsg_.triMesh.status( fhUpper ).deleted() )
			bsg_.maskUpper()[fh.idx()] = Contour::ContourType::NOTHING;
		else
			bsg_.maskUpper()[fh.idx()] = Contour::ContourType::INSIDE;
	}

	bsg_.triMesh.request_vertex_status();
	bsg_.triMesh.request_face_status();
		
	for( const auto& heh : bsg_.triMesh.halfedges() ) {

		if( !heh.is_boundary() )
			continue;

		auto v0 = heh.from();
		auto v1 = heh.to();

		const auto contourType = contour_.getContourType( bsg_.triMesh.point( v0 ), bsg_.triMesh.point( v1 ) );
		
		// find quads with these two vertices
		OpenMesh::SmartVertexHandle v0q( v0.idx(), &bsg_.quadMesh );
		for( const auto& voh : v0q.outgoing_halfedges() ) {
			if( voh.is_boundary())
				continue;

			auto quad = voh.face();
			if( voh.to().idx() == v1.idx() ) {
				// heh and voh are the same halfedge but in tri and quad

				//if( voh.is_boundary() )
				//	break;

				// find corresponding triangle
				const auto idxLower = 2 * quad.idx();
				const auto idxUpper = 2 * quad.idx() + 1;
				const OpenMesh::SmartFaceHandle fhLower( idxLower, &triComplete );
				const OpenMesh::SmartFaceHandle fhUpper( idxUpper, &triComplete );
				const auto hehLower = fhLower.halfedge();
				const auto hehUpper = fhUpper.halfedge();
				bool foundEdge = false;
				for( const auto& h : fhLower.halfedges() ) {
					if( h.from().idx() == v0.idx() && h.to().idx() == v1.idx() ) {
						bsg_.maskLower()[quad.idx()] = contourType;
						foundEdge = true;
						break;
					}
				}
				if( foundEdge )
					break;
				for( const auto& h : fhUpper.halfedges() ) {
					if( h.from().idx() == v0.idx() && h.to().idx() == v1.idx() ) {
						bsg_.maskUpper()[quad.idx()] = contourType;
						foundEdge = true;
						break;
					}
				}
				if( foundEdge )
					break;
			}
			if( voh.next().to().idx() == v1.idx() ) {
				const auto idxLower = 2 * quad.idx();
				const auto idxUpper = 2 * quad.idx() + 1;
				const OpenMesh::SmartFaceHandle fhLower( idxLower, &bsg_.triMesh );
				const OpenMesh::SmartFaceHandle fhUpper( idxUpper, &bsg_.triMesh );
				const bool deletedLower = bsg_.triMesh.status( fhLower ).deleted();
				const bool deletedUpper = bsg_.triMesh.status( fhUpper ).deleted();
				
				if( deletedLower && deletedUpper ) {
					LOG( ERROR );
				} else if( deletedLower && !deletedUpper ) {
					for( const auto& h : fhUpper.halfedges() ) {
						if( h.opp().is_boundary() ) {
							//const auto v = h.to();
							//const int& idx = m.property( contourIdx, v );
							//const auto contourPoint = contourOld_( idx );
							//const auto bid = contourPoint.boundaryIdRight;
							//LOG_ASSERT( bid == contourType );
							bsg_.maskLower()[quad.idx()] = contourType;
							break;
						}
					}
				} else if( !deletedLower && deletedUpper ) {
					for( const auto& h : fhLower.halfedges() ) {
						if( h.opp().is_boundary() ) {
							//const auto v = h.to();
							//const int& idx = m.property( contourIdx, v );
							//const auto contourPoint = contourOld_( idx );
							//const auto bid = contourPoint.boundaryIdRight;
							//LOG_ASSERT( bid == contourType );
							bsg_.maskUpper()[quad.idx()] = contourType;
							break;
						}
					}
				} else {
					LOG( ERROR );
				}
			}
		}
	}

	//// smooth masked vertices ////
	// copy to triComplete
	for( const auto& vQuad : bsg_.quadMesh.vertices() ) {
		const auto vTriComplete = OpenMesh::SmartVertexHandle( vQuad.idx(), &triComplete );
		triComplete.set_point( vTriComplete, bsg_.quadMesh.point( vQuad ) );
	}
	//OpenMesh::IO::write_mesh( triComplete, ( config_.outputFolder() / "deb.off" ).string() );
	std::vector<OpenMesh::SmartVertexHandle> outerVertices;
	for( const auto& vTriComplete : triComplete.vertices() ) {
		const auto vTri = OpenMesh::SmartVertexHandle( vTriComplete.idx(), &bsg_.triMesh );
		if( bsg_.triMesh.status( vTri ).deleted() ) {
			outerVertices.push_back( vTriComplete );
		}
	}
	DMO::DmoMesh dmoMeshOuter( outerVertices );
	//DMO::Solver dmoOuter( triComplete, &metricMeanRatio, &dmoMeshOuter );
	DMO::Metrics::SdfMinConstrainedMeanRatioTriangle sdfmincmrm( metricMeanRatio, sdfTriangle );
	DMO::Solver dmoOuter( triComplete, &sdfmincmrm, &dmoMeshOuter );
	dmoOuter.solve( 100 );
	// copy to quadMesh
	for( const auto& vQuad : bsg_.quadMesh.vertices() ) {
		const auto vTriComplete = OpenMesh::SmartVertexHandle( vQuad.idx(), &triComplete );
		bsg_.quadMesh.set_point( vQuad, triComplete.point( vTriComplete ) );
	}
	//OpenMesh::IO::write_mesh( triComplete, ( config_.outputFolder() / "deb.off" ).string() );


	m.garbage_collection();

	//// test masking
	//bsg_.quadMesh.request_face_colors();
	//for( const auto& fh : bsg_.quadMesh.faces() ) {
	//	auto l = bsg_.maskLower()[fh.idx()];
	//	if( l < 0 )
	//		l = l <= -1000 ? 255 : 127;
	//	auto u = bsg_.maskUpper()[fh.idx()];
	//	if( u < 0 )
	//		u = u <= -1000 ? 255 : 127;
	//	if( l == 5000 ) {
	//		l = 50;
	//	}
	//	if( u == 5000 ) {
	//		u = 50;
	//	}
	//	bsg_.quadMesh.set_color( fh, { l, u, 0 } );
	//}
	//OpenMesh::IO::write_mesh( bsg_.quadMesh, ( config_.outputFolder() / "deb5.off" ).string(), OpenMesh::IO::Options::FaceColor );
	//MeshFunctions::printMesh( config_.outputFolder() / "deb5tri.off", m );
}
