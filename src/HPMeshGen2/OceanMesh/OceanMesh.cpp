#include "OceanMesh.h"

#include <fstream>

#include <glog/logging.h>

#include "HPMeshGen2/HelperFunctions.h"
#include "Contour/ContourType.h"

using std::vector;
namespace fs = std::experimental::filesystem;

// Helper functions for import_mesh_adcirc
void deleteComments( std::string &str ) {
	if ( str.find( '!' ) == std::string::npos ) {
		// no comment in this string
		return;
	}
	// delete comment
	str.erase( str.find( '!' ) );
	// delete spaces at the beginning of string
	str.erase( str.begin(), std::find_if( str.begin(), str.end(), []( int ch ) { return !std::isspace( ch ); } ) );
	// delete spaces at the end of string
	str.erase( std::find_if( str.rbegin(), str.rend(), []( int ch ) { return !std::isspace( ch ); } ).base(), str.end() );
}
bool getLine( std::ifstream& ifs, std::string& line ) {
	std::getline( ifs, line );
	deleteComments( line );
	while ( line.size() == 0 ) {
		if ( !std::getline( ifs, line ) ) {
			return false;
		}

		deleteComments( line );
	}
	return true;
}

bool getLine( std::ifstream& ifs, std::istringstream& iss ) {
	std::string line;
	bool ret = getLine( ifs, line );
	iss.clear();
	iss.str( line );

	return ret;
}

void importMeshAdcirc( const std::experimental::filesystem::path& filename, std::vector<TriMesh::Point>& vertices, vector<vector<size_t>>& faces,
					   vector<vector<size_t>>& neumannBoundaries, vector<int>& neumannBoundaryTypes, vector<std::vector<size_t>>& dirichletBoundaries, vector<int>& dirichletBoundaryTypes ) {
	std::ifstream ifs;

	ifs.open( filename );
	if ( !ifs.is_open() ) {
		LOG( ERROR ) << "Unable to open file '" << filename << "'";
	}

	std::string line;

	// first line - mesh name
	getLine( ifs, line );
	LOG( INFO ) << "input-mesh name: " << line;

	// second line - ne & np
	getLine( ifs, line );
	std::istringstream iss( line );
	size_t nTriangles;
	size_t nVertices;
	iss >> nTriangles >> nVertices;

	vertices.resize( nVertices );
	faces.resize( nTriangles );

	// vertices
	for ( size_t i = 0; i < nVertices; ++i ) {
		getLine( ifs, line );

		iss.clear();
		iss.str( line );

		size_t dummy; // this is needed to ignore the vertex ID
		float x;
		float y;
		float z;

		iss >> dummy >> x >> y >> z;

		vertices[i] = { x, y, z };
	}

	// elements
	for ( size_t i = 0; i < nTriangles; ++i ) {
		getLine( ifs, line );
		iss.clear();
		iss.str( line );

		size_t idDummy;
		size_t nDummy;
		size_t id1;
		size_t id2;
		size_t id3;

		iss >> idDummy >> nDummy >> id1 >> id2 >> id3;

		faces[i] = { id1 - 1, id2 - 1, id3 - 1 };
	}

	// Neumann Boundaries
	if ( !getLine( ifs, line ) ) {
		// no boundary information given --> return
		return;
	}
	iss.clear();
	iss.str( line );
	size_t nNeumannBoundaries;
	iss >> nNeumannBoundaries;

	neumannBoundaries.resize( nNeumannBoundaries );
	neumannBoundaryTypes.resize( nNeumannBoundaries );

	getLine( ifs, line );

	for ( size_t i = 0; i < nNeumannBoundaries; ++i ) {
		getLine( ifs, line );
		iss.clear();
		iss.str( line );
		size_t nBoundaryVertices;	// number of vertices in this neumann boundary
		iss >> nBoundaryVertices;
		int boundaryType;
		if ( iss >> boundaryType ) {
			neumannBoundaryTypes[i] = boundaryType;
		} else {
			neumannBoundaryTypes[i] = 0;
		}

		std::vector<size_t> neumannBoundary( nBoundaryVertices );

		for ( size_t j = 0; j < nBoundaryVertices; ++j ) {
			getLine( ifs, line );
			iss.clear();
			iss.str( line );
			size_t id;
			iss >> id;

			neumannBoundary[j] = id - 1;
		}

		neumannBoundaries[i] = neumannBoundary;
	}

	// Dirichlet Boundaries
	getLine( ifs, line );
	iss.clear();
	iss.str( line );
	size_t nDirichletBoundaries;
	iss >> nDirichletBoundaries;

	dirichletBoundaries.resize( nDirichletBoundaries );
	dirichletBoundaryTypes.resize( nDirichletBoundaries );

	getLine( ifs, line );

	for ( size_t i = 0; i < nDirichletBoundaries; ++i ) {
		getLine( ifs, line );
		iss.clear();
		iss.str( line );
		size_t nDirichletVertices;	// number of vertices in this neumann boundary
		iss >> nDirichletVertices;
		int boundaryType;
		if ( iss >> boundaryType ) {
			dirichletBoundaryTypes[i] = boundaryType;
		} else {
			dirichletBoundaryTypes[i] = 0;
		}

		std::vector<size_t> dirichletBoundary( nDirichletVertices );

		for ( size_t j = 0; j < nDirichletVertices; ++j ) {
			getLine( ifs, line );
			iss.clear();
			iss.str( line );
			size_t id;
			iss >> id;

			dirichletBoundary[j] = id - 1;
		}

		dirichletBoundaries[i] = dirichletBoundary;
	}

	ifs.close();
}

void storeFeatureSize(OceanMesh& mesh) {
	#pragma omp parallel for
	for ( int i = 0; i < mesh.n_vertices(); ++i ) {
		TriMesh::VertexHandle vh = mesh.vertex_handle( i );
		// give all interior points a dummy value (max angle --> 180 degree)
		if ( !mesh.is_boundary( vh ) ) {
			mesh.property( mesh.feature_size, vh ) = 180.f;
			
			continue;
		}

		// compute feature size for boundary point
		std::vector<TriMesh::Point>boundaryNeighbors;
		std::vector<TriMesh::EdgeHandle> incidentEdgeHandles;
		for ( const auto& voh : mesh.voh_range( vh ) ) {
			if ( mesh.is_boundary( voh.edge() ) ) {
				boundaryNeighbors.push_back( mesh.point( voh.to() ) );
				incidentEdgeHandles.push_back( voh.edge() );
			}
		}

		LOG_ASSERT( boundaryNeighbors.size() == 2 );

		float angle = HelperFunctions::calcAngle( boundaryNeighbors[0], mesh.point( vh ), boundaryNeighbors[1] );
		mesh.property( mesh.feature_size, vh ) = angle;

		// set feature size at the vertex between Dirichlet and Neumann to 0
		if( mesh.property( mesh.propContourType, incidentEdgeHandles[0] ) != mesh.property( mesh.propContourType, incidentEdgeHandles[1] ) ) {
			mesh.property( mesh.feature_size, vh ) = 0.f;
			mesh.data( vh ).is_feature = true;
		}
	}
}

OceanMesh::OceanMesh( const fs::path& filename ) : TriMesh(), filename_{ filename } {
	add_property( feature_size, "feature_size" );
	add_property( propContourType, "propContourType" );
	add_property( propEmo, "propEmo" );
	add_property( propEfa, "propEfa" );
	
	if ( filename_.extension() == ".off" || filename_.extension() == ".obj" ) {
		readOffFile();
	} else if ( filename_.extension() == ".om" ) {
		readOmFile();
	} else if ( filename_.extension() == ".14" ) {
		readFortFiles();
	} else {
		LOG( ERROR ) << "Unknown mesh file type. Mesh file is " << filename_;
	}
	
	storeFeatureSize(*this);
}

void OceanMesh::loadBackgroundGrids( const fs::path& cacheFolder, const fs::path& meshFile, const size_t & nx, const size_t& ny ) {
	sizefield = SizeField::load( *this, cacheFolder, meshFile, nx, ny );
	depthfield = DepthField::load( *this, cacheFolder, meshFile, nx, ny );
}

void OceanMesh::flatten() {
	for ( const auto& vh : vertices() ) {
		auto p = point( vh );
		p[2] = 0.f;
		set_point( vh, p );
	}
}


TriMesh::Point OceanMesh::sphericalToCartesian( const TriMesh::Point& ps ) const {
	if( parameters.coordinate_system_type == parameters.CARTESIAN ) {
		return ps;
	}
	if( parameters.coordinate_system_type != parameters.SPHERICAL ) {
		LOG( WARNING ) << "Unknown coordinate system type! Value is " << parameters.coordinate_system_type;
		return ps;
	}

	//const auto& r = parameters.earth_radius;
	//const auto& lo = ps[0] * M_PI / 180.;	// longitude
	//const auto& la = ps[1] * M_PI / 180.;	// latitude
	//
	//TriMesh::Point pc;
	//pc[0] = r * std::cos( lo ) * std::cos( la );
	//pc[1] = r * std::sin( lo ) * std::cos( la );
	//pc[2] = r * std::sin( la );


	// transformation as it is done in UTBEST
	const auto& r = parameters.earth_radius;
	const auto& lambda0 = parameters.ref_lam;
	const auto& phi0 = parameters.ref_phi;
	const auto& lambda = ps[0] * M_PI / 180.;	// longitude
	const auto& phi = ps[1] * M_PI / 180.;		// latitude
	const auto x = r * ( lambda - lambda0 ) * std::cos( phi0 );
	const auto y = phi * r;

	TriMesh::Point pc( x, y, 0 );

	return pc;
}

void OceanMesh::readFortFiles() {
	read14File();

	fs::path file15 = filename_;
	file15.replace_extension( "15" );
	read15File( file15 );

	fs::path file17 = filename_;
	file17.replace_extension( "17" );
	read17File( file17 );
}

void OceanMesh::read14File() {
	std::vector<TriMesh::Point> points;	// TODO change this to TriMesh::Point
	std::vector<std::vector<size_t>> faces;

	std::vector<std::vector<size_t>> neumannBoundaries;
	std::vector<std::vector<size_t>> dirichletBoundaries;
	std::vector<int> neumannBoundaryTypes;
	std::vector<int> dirichletBoundaryTypes;

	importMeshAdcirc( filename_, points, faces, neumannBoundaries, neumannBoundaryTypes, dirichletBoundaries, dirichletBoundaryTypes );

	std::vector<TriMesh::VertexHandle> vhandle( points.size() );
	std::transform( points.begin(), points.end(), vhandle.begin(), [this]( TriMesh::Point& p ) { return add_vertex( p ); } );

	for( const auto& f : faces ) {
		add_face( { vhandle[f[0]], vhandle[f[1]],vhandle[f[2]] } );
	}

	// add boundary properties
	// Neumann
	for ( size_t i = 0; i < neumannBoundaries.size(); ++i ) {
		for ( size_t j = 0; j < neumannBoundaries[i].size() - 1; ++j ) {
			TriMesh::VertexHandle vh1 = vertex_handle( static_cast<uint>( neumannBoundaries[i][j] ) );
			TriMesh::VertexHandle vh2 = vertex_handle( static_cast<uint>( neumannBoundaries[i][j + 1] ) );

			for( const auto& voh : voh_range( vh1 ) ) {
				if( voh.to() == vh2 ) {
					const auto& eh = voh.edge();
					switch( neumannBoundaryTypes[i] ) {
					case 0:
						property( propContourType, eh ) = Contour::ContourType::OPEN_SEA;
						break;
					case 100:
						property( propContourType, eh ) = Contour::ContourType::RIVER;
						break;
					default:
						LOG( WARNING ) << "Unknown boundary type: " << neumannBoundaryTypes[i] << ". Set to OpenSea";
						property( propContourType, eh ) = Contour::ContourType::OPEN_SEA;
						break;
					}
				}
			}
		}
	}
	// Dirichlet
	for ( size_t i = 0; i < dirichletBoundaries.size(); ++i ) {
		for ( size_t j = 0; j < dirichletBoundaries[i].size(); ++j ) {
			TriMesh::VertexHandle vh1 = vertex_handle( (uint)dirichletBoundaries[i][j] );
			TriMesh::VertexHandle vh2 = vertex_handle( (uint)dirichletBoundaries[i][( j + 1 ) % dirichletBoundaries[i].size()] );

			for( const auto& voh : voh_range( vh1 ) ) {
				if( voh.to() == vh2 ) {
					const auto& eh = voh.edge();
					switch( dirichletBoundaryTypes[i] ) {
					case 0:
						property( propContourType, eh ) = Contour::ContourType::LAND;
						break;
					case 1:
						property( propContourType, eh ) = Contour::ContourType::ISLAND;
						break;
					default:
						LOG( WARNING ) << "Unknown boundary type: " << dirichletBoundaryTypes[i] << ". Set to Land";
						property( propContourType, eh ) = Contour::ContourType::LAND;
						break;
					}
				}
			}
		}
	}
}

void OceanMesh::read15File( const fs::path& file15 ) {
	if( !fs::exists( file15 ) ) {
		LOG( INFO ) << ".15 file for '" << filename_.filename() << "' not found.";
		return;
	}

	std::ifstream ifs;
	std::istringstream iss;

	ifs.open( file15 );
	if( !ifs.is_open() ) {
		LOG( ERROR ) << "Unable to open file '" << file15 << "'";
	}

	
	// Input run description
	getLine( ifs, iss );	// 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION
	iss >> parameters.run_description;

	// Input run id
	getLine( ifs, iss );	// 24 CHARACTER ALPANUMERIC RUN IDENTIFICATION 
	iss >> parameters.run_id;

	// Input hot input data
	getLine( ifs, iss );	// IHOT - parameter fort hot start, 0=cold, 67-read fort.67 (binary)
	iss >> parameters.hot_input;
	std::string hot_input_file;
	switch( parameters.hot_input ) {
	case 0: break;
	case 1:
		getLine( ifs, iss );
		iss >> hot_input_file;
		break;
	default:
		LOG( WARNING ) << "Cold/Hot start switch out of range! Value is " << parameters.hot_input;
		return;
	}

	// Input coordinate system type (Cartesian/Spherical)
	getLine( ifs, iss );	// ICS - COORDINATE SYSTEM SELECTION PARAMETER, 1 = cartesian, 2 = spherical
	iss >> parameters.coordinate_system_type;
	if( parameters.coordinate_system_type != parameters.CARTESIAN && parameters.coordinate_system_type != parameters.SPHERICAL ) {
		LOG( WARNING ) << "Unknown coordinate system type! Value is " << parameters.coordinate_system_type;
	}
	
	getLine( ifs, iss );	// NOLIBF - BOTTOM FRICTION TERM SELECTION PARAMETER
	getLine( ifs, iss );	// NWP - VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER

	getLine( ifs, iss );	// NCOR - VARIABLE CORIOLIS IN SPACE OPTION PARAMETER
	int coriolis_parameter;
	iss >> coriolis_parameter;
	if( coriolis_parameter != 0 && coriolis_parameter != 1 ) {
		LOG( WARNING ) << "Coriolis const/variable switch out of range! Value is " << coriolis_parameter;
	}

	getLine( ifs, iss );	// NTIP - TIDAL POTENTIAL OPTION PARAMETER
	iss >> parameters.tidal_potential_parameter;
	if( parameters.tidal_potential_parameter != 0 && parameters.tidal_potential_parameter != 1 ) {
		LOG( WARNING ) << "Tidal potential const/variable switch out of range! Value is " << parameters.tidal_potential_parameter;
	}

	getLine( ifs, iss );	// NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER
	getLine( ifs, iss );	// NRAMP - RAMP FUNCTION OPTION
	iss >> parameters.ramp;
	if( parameters.ramp != 0 && parameters.ramp != 1 ) {
		LOG( WARNING ) << "Ramp on/off switch out of range! Value is " << parameters.ramp;
	}

	getLine( ifs, iss );	// G - ACCELERATION DUE TO GRAVITY - DETERMINES UNITS 
	double gravity;
	iss >> gravity;
	double unit_conversion_factor = gravity / 9.81;
	parameters.earth_radius = 6378206.4 * unit_conversion_factor;

	getLine( ifs, iss );	// NUMBER OF QUADRATURE POINTS
	int n_quadrature_points;
	iss >> n_quadrature_points;
	for( int i = 0; i < n_quadrature_points; ++i ) {
		getLine( ifs, iss );
	}

	getLine( ifs, iss );	// NDTVAR - AUTOMATIC TIME STEP SELECTION PARAMETER, 0=NO, 1=YES
	int time_step_selection_parameter;
	iss >> time_step_selection_parameter;
	if( time_step_selection_parameter != 0 && time_step_selection_parameter != 1 ) {
		LOG( WARNING ) << "Time step const/variable switch out of range! Value is " << time_step_selection_parameter;
	}

	getLine( ifs, iss );	// DT - TIME STEP( IN SECONDS )
	getLine( ifs, iss );	// STATIM - STARTING TIME (IN DAYS)
	getLine( ifs, iss );	// REFTIM - REFERENCE TIME (IN DAYS)
	getLine( ifs, iss );	// RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS) (12.0 DAYS)
	getLine( ifs, iss );	// IRK - TIME STEPPING SCHEME: 1=EULER FORWARD, 2=2nd ORDER RUNGE-KUTTA
	getLine( ifs, iss );	// ISLOPE - SLOPE LIMITING PARAMETER: 1=YES, 0=NO.
	getLine( ifs, iss );	// ITRANS - 0=TRANSIENT, 1=STEADY STATE CALCULATION
	getLine( ifs, iss );	// CONVCR - CONV. CRITERIA FOR STEADY STATE CALCULATIONS
	getLine( ifs, iss );	// DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)
	getLine( ifs, iss );	// H0 - MIMINMUM CUTOFF DEPTH

	getLine( ifs, iss );	// SLAM0,SFEA0 - CENTER OF CPP PROJECTION IN DEGREES LONG/LAT
	iss >> parameters.ref_lam >> parameters.ref_phi;
	// Transform lon/lat to fit the correct range [0,360]x[-90,90]
	parameters.ref_lam = parameters.ref_lam < 0. ? parameters.ref_lam + 360. : parameters.ref_lam;
	double coef = parameters.ref_phi / 180., abs_coef = std::abs( coef );
	parameters.ref_phi -= abs_coef > 0.5 ? coef / abs_coef * 180. : 0.;
	// Convert lam and phi from degrees to radians
	parameters.ref_lam *= M_PI / 180.;
	parameters.ref_phi *= M_PI / 180.;

	getLine( ifs, iss );	// TAU - LINEARIZED FRICTION COEF USED WHEN NOLI = 0 OR 1
	getLine( ifs, iss );	// CF - BOTTOM FRICTION COEFFICIENT USED WHEN NOLI = 2 OR 3; IGNORED IF NWP = 1
	getLine( ifs, iss );	// NVISC - VISCOUS TERMS SELECTION PARAMETER, 0=NO, 1=YES
	getLine( ifs, iss );	// ESL - LATERAL EDDY VISCOSITY COEFFICIENT; IGNORED IF NWP =1
	getLine( ifs, iss );	// CORI - CORIOLIS PARAMETER - IGNORED IF NCOR = 1
	getLine( ifs, iss );	// NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED
	int n_tidal_potential_forcings;
	iss >> n_tidal_potential_forcings;
	for( int i = 0; i < n_tidal_potential_forcings; ++i ) {
		if( parameters.tidal_potential_parameter == 1 ) {
			getLine( ifs, iss );
			getLine( ifs, iss );
		} else 	{
			getLine( ifs, iss );
			getLine( ifs, iss );
		}
	}

	getLine( ifs, iss );	// NBFR - TOTAL NUMBER OF FORCING FREQUENCIES ON OPEN BOUNDARIES
	int n_forcing_frequencies;
	iss >> n_forcing_frequencies;
	if( n_forcing_frequencies < 0 ) {
		LOG( WARNING ) << "Negative number of boundary tidal potential components! Value is " << n_forcing_frequencies;
	}
	for( int i = 0; i < n_forcing_frequencies; ++i ) {
		getLine( ifs, iss );
		getLine( ifs, iss );
	}

	getLine( ifs, iss );	// NOUTE,TOUTSE,TOUTFE,NSPOOLE:ELEV STATION OUTPUT INFO (UNIT  61)
	getLine( ifs, iss );	// TOTAL NUMBER OF ELEVATION RECORDING STATIONS
	int n_elevation_recording_stations;
	iss >> n_elevation_recording_stations;
	for( int i = 0; i < n_elevation_recording_stations; ++i ) {
		getLine( ifs, iss );
	}
	getLine( ifs, iss );	// NOUTV,TOUTSV,TOUTFV,NSPOOLV:ELEV STATION OUTPUT INFO (UNIT  62)
	getLine( ifs, iss );	// TOTAL NUMBER OF VELOCITY RECORDING STATIONS
	int n_velocity_recording_stations;
	iss >> n_velocity_recording_stations;
	for( int i = 0; i < n_velocity_recording_stations; ++i ) {
		getLine( ifs, iss );
	}
	getLine( ifs, iss );	// NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE : GLOBAL OUTPUT INFO (UNIT  63)
	getLine( ifs, iss );	// NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV : GLOBAL OUTPUT INFO (UNIT  64)
	getLine( ifs, iss );	// hot start output control

	ifs.close();
}

void OceanMesh::read17File( const fs::path& file17 ) {
	if( !fs::exists( file17 ) ) {
		LOG( INFO ) << ".17 file for '" << filename_.filename() << "' not found.";
		return;
	}

	std::ifstream ifs;
	std::istringstream iss;

	ifs.open( file17 );
	if( !ifs.is_open() ) {
		LOG( ERROR ) << "Unable to open file '" << file17 << "'";
	}


	// first line - n_edges
	getLine( ifs, iss );
	int n_17edges;
	iss >> n_17edges;
	if( n_17edges != this->n_edges() ) {
		LOG( WARNING ) << "Number of edges in .17 file does not fit to .14 file --> .17 file will be ignored.";
		return;
	}

	// all edges
	std::vector<std::array<int, 4>> edges( n_17edges );
	auto findEdge = [&edges, this]( int eId ) {
		const auto& edgeInfo = edges[eId];
		const auto& v0id = edgeInfo[0];
		const auto& v1id = edgeInfo[1];
		const auto v0 = OpenMesh::SmartVertexHandle( v0id, this );
		const auto v1 = OpenMesh::SmartVertexHandle( v1id, this );
		bool edgeFound = false;
		for( const auto& voh : v0.outgoing_halfedges() ) {
			if( voh.to() == v1 ) {
				return voh.edge();
			}
		}
		return OpenMesh::SmartEdgeHandle();
	};

	for( auto& e : edges ) {
		getLine( ifs, iss );
		int edgeId, n1, n2, t1, t2;
		iss >> edgeId >> n1 >> n2 >> t1 >> t2;
		e[0] = n1 - 1;
		e[1] = n2 - 1;
		e[2] = t1 - 1;
		e[3] = t2 - 1;
	}

	// all elements and incident edges
	std::vector<std::array<int, 3>> faces( n_faces() );
	for( auto& f : faces ) {
		getLine( ifs, iss );
		int faceId, e1, e2, e3;
		iss >> faceId >> e1 >> e2 >> e3;
		f[0] = e1 - 1;
		f[1] = e2 - 1;
		f[2] = e3 - 1;
	}

	getLine( ifs, iss );
	int n_interior;
	iss >> n_interior;
	std::vector<int> interior( n_interior );
	for( auto& i : interior ) {
		getLine( ifs, iss );
		int buf, id;
		iss >> buf >> id;
		i = id - 1;
	}

	getLine( ifs, iss );
	int n_land;
	iss >> n_land;
	std::vector<int> land( n_land );
	for( auto& i : land ) {
		getLine( ifs, iss );
		int buf, id;
		iss >> buf >> id;
		i = id - 1;
	}

	getLine( ifs, iss );
	int n_radiation;
	iss >> n_radiation;
	std::vector<int> radiation( n_radiation );
	for( auto& i : radiation ) {
		getLine( ifs, iss );
		int buf, id;
		iss >> buf >> id;
		i = id - 1;
	}
	

	getLine( ifs, iss );
	int n_river;
	iss >> n_river;
	std::vector<int> river( n_river );
	for( auto& i : river ) {
		getLine( ifs, iss );
		int buf, id;
		iss >> buf >> id;
		i = id - 1;
	}

	getLine( ifs, iss );
	int n_sea;
	iss >> n_sea;
	std::vector<int> sea( n_sea );
	for( auto& i : sea ) {
		getLine( ifs, iss );
		int buf, id;
		iss >> buf >> id;
		i = id - 1;
	}

	for( int i = 0; i < n_land; ++i ) {
		const auto& eId = land[i];
		const auto e = findEdge( eId );
		if( !e.is_valid() ) {
			LOG( WARNING ) << "Could not find edge in mesh: fort.17 file might be corrupted.";
		} else {
			if( property( propContourType, e ) != Contour::ContourType::LAND && property( propContourType, e ) != Contour::ContourType::ISLAND ) {
				if( property( propContourType, e ) == Contour::ContourType::INSIDE )
					LOG( WARNING ) << "Contour type was not set in .14 file for edge " << eId << ". Use .17 file value";
				else
					LOG( ERROR ) << "Types in .14 and .17 file do not fit for edge " << eId << ". Use .17 file value";
				property( propContourType, e ) = Contour::ContourType::LAND;
			}
		}
	}
	for( int i = 0; i < n_radiation; ++i ) {
		const auto& eId = radiation[i];
		const auto e = findEdge( eId );
		if( !e.is_valid() ) {
			LOG( WARNING ) << "Could not find edge in mesh: fort.17 file might be corrupted.";
		} else {
			if( property( propContourType, e ) != Contour::ContourType::RADIATION ) {
				if( property( propContourType, e ) == Contour::ContourType::INSIDE )
					LOG( WARNING ) << "Contour type was not set in .14 file for edge " << eId << ". Use .17 file value";
				else
					LOG( ERROR ) << "Types in .14 and .17 file do not fit for edge " << eId << ". Use .17 file value";
				property( propContourType, e ) = Contour::ContourType::RADIATION;
			}
		}
	}
	for( int i = 0; i < n_river; ++i ) {
		const auto& eId = river[i];
		const auto e = findEdge( eId );
		if( !e.is_valid() ) {
			LOG( WARNING ) << "Could not find edge in mesh: fort.17 file might be corrupted.";
		} else {
			if( property( propContourType, e ) != Contour::ContourType::RIVER ) {
				if( property( propContourType, e ) == Contour::ContourType::INSIDE )
					LOG( WARNING ) << "Contour type was not set in .14 file for edge " << eId << ". Use .17 file value";
				else
					LOG( ERROR ) << "Types in .14 and .17 file do not fit for edge " << eId << ". Use .17 file value";
				property( propContourType, e ) = Contour::ContourType::RIVER;
			}
		}
	}
	for( int i = 0; i < n_sea; ++i ) {
		const auto& eId = sea[i];
		const auto e = findEdge( eId );
		if( !e.is_valid() ) {
			LOG( WARNING ) << "Could not find edge in mesh: fort.17 file might be corrupted.";
		} else {
			if( property( propContourType, e ) != Contour::ContourType::OPEN_SEA ) {
				if( property( propContourType, e ) == Contour::ContourType::INSIDE )
					LOG( WARNING ) << "Contour type was not set in .14 file for edge " << eId << ". Use .17 file value";
				else
					LOG( ERROR ) << "Types in .14 and .17 file do not fit for edge " << eId << ". Use .17 file value";
				property( propContourType, e ) = Contour::ContourType::OPEN_SEA;
			}
		}
	}

	// emo/efa
	std::vector<std::vector<float>> emo;
	std::vector<std::vector<float>> efa;
	bool reachedEnd = n_sea == 0;
	while( !reachedEnd ) {
		std::vector<float> emoBuf( n_sea );
		std::vector<float> efaBuf( n_sea );
		for( int i = 0; i < n_sea; ++i ) {
			reachedEnd = !getLine( ifs, iss );
			if( reachedEnd )
				break;
			iss >> emoBuf[i] >> efaBuf[i];
		}
		if( reachedEnd )
			break;
		emo.push_back( emoBuf );
		efa.push_back( efaBuf );
	}
	LOG( INFO ) << "Found " << emo.size() << " emo/eva values in .17 file.";
	nEmoEfaVals = emo.size();

	ifs.close();

	// update emo/efa properties
	for( int i = 0; i < n_sea; ++i ) {
		const auto& eId = sea[i];

		const auto e = findEdge( eId );
		if( !e.is_valid() ) {
			LOG( WARNING ) << "Could not find sea edge in mesh: fort.17 file might be corrupted.";
		} else {
			auto& emoVals = property( propEmo, e );
			for( auto& val : emo ) {
				emoVals.push_back( val[i] );
			}
			auto& efaVals = property( propEfa, e );
			for( auto& val : efa ) {
				efaVals.push_back( val[i] );
			}
		}

		//const auto& edgeInfo = edges[eId];
		//const auto& v0id = edgeInfo[0];
		//const auto& v1id = edgeInfo[1];
		//
		//const auto v0 = OpenMesh::SmartVertexHandle( v0id, this );
		//const auto v1 = OpenMesh::SmartVertexHandle( v1id, this );
		//bool edgeFound = false;
		//for( const auto& voh : v0.outgoing_halfedges() ) {
		//	if( voh.to() == v1 ) {
		//		const auto& e = voh.edge();
		//		auto& emoVals = property( propEmo, e );
		//		for( auto& val : emo ) {
		//			emoVals.push_back( val[i] );
		//		}
		//		auto& efaVals = property( propEfa, e );
		//		for( auto& val : efa ) {
		//			efaVals.push_back( val[i] );
		//		}
		//
		//		edgeFound = true;
		//		break;
		//	}
		//}
		//if( !edgeFound ) {
		//	LOG( WARNING ) << "Could not find sea edge in mesh: fort.17 file might be corrupted.";
		//}
	}

	// sanity check
	int n_sea_sanity = 0;
	for( const auto& e : this->edges() ) {
		const auto& emoVal = property( propEmo, e );
		const auto& efaVal = property( propEfa, e );
		LOG_ASSERT( emoVal.size() == efaVal.size() );
		if( e.is_boundary() ) {
			if( emoVal.size() > 0 )
				++n_sea_sanity;
		} else 	{
			LOG_ASSERT( emoVal.empty() );
		}
	}
	LOG_ASSERT( n_sea == n_sea_sanity );
}

void OceanMesh::readOffFile() {
	OpenMesh::IO::read_mesh( *this, filename_.string() );
}

void OceanMesh::readOmFile() {
	//add_property( feature_size, "feature_size" );
	OpenMesh::IO::read_mesh( *this, filename_.string() );
}
