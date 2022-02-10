#include <experimental/filesystem>
// glog
#include <glog/logging.h>

#include "Config.h"

#include "BsgGenerator.h"

// Mesh Function
#include "MeshFunctions.h"

#include "Simplificator.h"

#include "DMO/Solver.h"
#include "GridEvaluation/eval.h"

#include "baseDirectory.h"
#include "Simplify.h"
#include "ScalarField/ScalarField.h"
#include "Stopwatch/Stopwatch.h"

namespace fs = std::experimental::filesystem;

void generateBsg( const Config& config ) {

	BsgGenerator bsgGenerator( config );

	if ( config.reductionOutput() ) bsgGenerator.simplificator().simplificationOutput( config.outputFolder() );


	bsgGenerator.generate();
	//MeshFunctions::printMesh( config.outputFolder() / "fragmentMesh.off", bsgGenerator.fragmentMesh() );
	//MeshFunctions::printMesh( config.outputFolder() / "fragmentMeshRefined.off", bsgGenerator.bsg().quadMesh );
	//MeshFunctions::printMesh( config.outputFolder() / "fragmentMeshRefinedTri.off", bsgGenerator.bsg().triMesh );
	if( config.convexHullDecimation() ) {
		bsgGenerator.adaptElementsToContour();
	}

	BlockStructuredGrid& bsg = bsgGenerator.bsg();

	LOG( INFO ) << "# quads after refinement: " << bsg.quadMesh.n_faces() << std::endl;

	if ( config.sizegridOutput() ) { 
		bsgGenerator.inputMesh().sizeField().toVTK( config.outputFolder() / "sizegrid.vtk", false );
		bsgGenerator.inputMesh().depthField().toVTK( config.outputFolder() / "depthgrid.vtk", false );
	}
	if ( config.blockMeshOutput() ) {
		MeshFunctions::printMesh( config.outputFolder() / "simplifiedMesh.off", bsgGenerator.simplifiedMesh() );
		MeshFunctions::printMesh( config.outputFolder() / "fragmentMesh.off", bsgGenerator.fragmentMesh() );
		MeshFunctions::printMesh( config.outputFolder() / "fragmentMeshRefined.off", bsgGenerator.bsg().quadMesh );
	}

	//bsgGenerator.optimizeBsg();
	
	
	// evaluation
	MeshFunctions::printMesh( config.outputFolder() / "fineMeshFinal.off", bsg.triMesh );
	MeshFunctions::displayQuality( bsg.triMesh, 10 );

	// apply depth to BSG
	//bsgGenerator.inputMesh().applyDepth( bsg.triMesh );
	//bsgGenerator.inputMesh().applyDepth( bsg.quadMesh );

	///////////////////////
	// write bsg to file //

	// rescale mesh
	//MeshFunctions::rescaleMesh(*bsg.fineMesh, 1000, 1000, 1000);
	
	//bsg.patch_segmentation(nBlocks);
	//bsg.write_bsg(config.outputFolder(), "BSG_east_b" + std::to_string(nBlocks) + "_f" + std::to_string(bsg.n_patches()));
	//bsg.write_bsg( config.outputFolder(), "BSG_east_b" + std::to_string(nBlocks) + "_f" + std::to_string(bsg.n_patches()) + "_scaled");
	
	//if (config.meshFile().extension() == ".14") {
	//	bsg.exportMeshAdcirc( config.outputFolder(), config.meshFile().stem().string() + "_" + std::to_string(config.nRefinementSteps()) + ".14");
	//}

	std::string bsgName = std::string( "BSG_" ) + config.meshFile().stem().string() + "_b" + std::to_string( config.nBlocks() ) + "_f" + std::to_string( bsg.nFragments() );
	if( config.convexHullDecimation() ) {
		bsgName += "_m";
	}
	fs::path bsgFolder = config.outputFolder() / bsgName;

	fs::create_directories( bsgFolder );

	if( !config.convexHullDecimation() ) {
		const auto& sf = bsgGenerator.inputMesh().sizeField();

		DMO::UniformGrid grid_d;
		grid_d.n = { (int)sf.Nx(), (int)sf.Ny() };
		grid_d.h = { sf.hx(), sf.hy() };
		grid_d.aabbMin = { sf.aabb().xMin, sf.aabb().yMin };
		grid_d.aabbMax = { sf.aabb().xMax, sf.aabb().yMax };
		thrust::host_vector<float> sgVals;
		sgVals.resize( grid_d.n.x * grid_d.n.y );
		for( int i = 0; i < grid_d.n.x * grid_d.n.y; ++i ) {
			sgVals[i] = sf( i );
		}
		thrust::device_vector<float> sgVals_d = sgVals;
		grid_d.vals = sgVals_d.data().get();


		DMO::Metrics::MeanRatioTriangle metricMeanRatio;
		DMO::Metrics::DensityTriangle metricDensity( grid_d );
		DMO::Metrics::DensityWithMeanRatioTriangle metricInner( metricMeanRatio, metricDensity );
		DMO::DmoMesh dmoMeshInner = DMO::DmoMesh::create<DMO::Set::Inner>( bsg.triMesh );

		DMO::Solver dmo( bsg.triMesh, &metricInner, &dmoMeshInner );
		dmo.solve( 100 );
		DMO::Solver( bsg.triMesh, &metricMeanRatio, &dmoMeshInner ).solve( 1 );

		bsg.copy_positions_tri2quad();
	}

	// evaluation
	MeshFunctions::printMesh( bsgFolder / ( bsgName + "_" + std::to_string( bsg.nGridNodes() ) + ".off" ), bsg.triMesh );
	MeshFunctions::printMeshCartesian( bsgFolder / ( bsgName + "_" + std::to_string( bsg.nGridNodes() ) + "_cart.off" ), bsg.triMesh, bsgGenerator.inputMesh() );
	MeshFunctions::displayQuality( bsg.triMesh, 10 );

	// apply depth to BSG
	bsgGenerator.inputMesh().applyDepth( bsg.triMesh );
	bsgGenerator.inputMesh().applyDepth( bsg.quadMesh );

	GridEvaluation::eval( bsg.triMesh, bsgGenerator.inputMesh(), bsgFolder / ( bsgName + "_" + std::to_string( bsg.nGridNodes() ) + ".vtk" ) );
	MeshFunctions::printMesh( bsgFolder / "triBlocks.off", bsgGenerator.simplifiedMesh() );
	MeshFunctions::printMeshCartesian( bsgFolder / "triBlocks_cart.off", bsgGenerator.simplifiedMesh(), bsgGenerator.inputMesh() );

	bsg.patch_segmentation( config.nBlocks() );
	bsg.write_bsg( bsgFolder, bsgName, config.convexHullDecimation() );

	// only write BSG and fort files if input is a fort.14 file
	if( config.meshFile().extension() == ".14" ) {
		bsg.exportMeshAdcirc( bsgFolder, config.meshFile().stem().string() + "_" + std::to_string( bsg.nGridNodes() ) + ".14", bsgGenerator.inputMesh() );
	}
	
	//bsg.print_svg( config.outputFolder().string() + "bsgPrint.svg" );
}

int main( int argc, char* argv[] ) 
{
	google::InitGoogleLogging( argv[0] );
	//FLAGS_timestamp_in_logfile_name = false;
	FLAGS_logtostderr = true;
	FLAGS_colorlogtostderr = true;
	
	fs::current_path( BASE_DIRECTORY );

	// load config file
	const Config config( "HPMeshGenParameter.txt" );
	
	// set logging flags
	//FLAGS_logtostderr = false;
	//FLAGS_alsologtostderr = true;
	//FLAGS_log_dir = config.outputFolder().string();
	// print config to log
	config.log();

	// do it!
	Stopwatch sw;
	sw.start();
	generateBsg( config );
	sw.stop();
	LOG( WARNING ) << "Runtime: " << sw.runtimeStr<Stopwatch::Milliseconds>();

	LOG( INFO ) << "HPMeshGen2 finished";
}