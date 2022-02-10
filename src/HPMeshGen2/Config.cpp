#include "Config.h"

#include <string>

#include "InputFileReader/InputFileReader.h"

Config::Config( const fsp& configFile ) : configFile_( configFile ) {
	namespace fs = std::experimental::filesystem;

	if ( !fs::exists( configFile_ ) || !fs::is_regular_file( configFile_ ) )
		LOG( FATAL ) << "Config file " << configFile_ << " does not exist";
	configFile_ = fs::canonical( configFile_ );

	InputFileReader ifr( configFile_ );

	workingDir_ = ifr.getValue<std::string>( "workingDirectory" );
	if ( !fs::exists( workingDir_ ) || !fs::is_directory( workingDir_ ) )
		LOG( FATAL ) << "Working directory " << workingDir_ << " does not exist or is not a directory";
	workingDir_ = fs::canonical( workingDir_ );
	LOG( INFO ) << "Set working directory to " << workingDir_;
	fs::current_path( workingDir_ );

	cacheFolder_ = ifr.getValue<std::string>( "cacheFolder" );
	if ( !fs::exists( cacheFolder_ ) ) {
		LOG( INFO ) << "Create cache folder " << cacheFolder_;
		fs::create_directories( cacheFolder_ );
	} else if ( !fs::is_directory( cacheFolder_ ) ) {
		LOG( FATAL ) << "Cache folder " << cacheFolder_ << " is not a directory!";
	}
	cacheFolder_ = fs::canonical( cacheFolder_ );

	nBlocks_ = ifr.getValue<size_t>( "nBlocks" );

	meshFile_ = ifr.getValue<std::string>( "meshFileName" );
	if ( !fs::exists( meshFile_ ) || !fs::is_regular_file( meshFile_ ) )
		LOG( FATAL ) << "Mesh file " << meshFile_ << " does not exist or is not a file";
	meshFile_ = fs::canonical( meshFile_ );

	nRefinementSteps_ = ifr.getValue<size_t>( "nRefinementSteps" );
	nFragments_ = ifr.getValue<size_t>( "nPatches" );
	sizegridSizeX_ = ifr.getValue<size_t>( "sizegridSizeX" );
	sizegridSizeY_ = ifr.getValue<size_t>( "sizegridSizeY" );

	outputFolder_ = ifr.getValue<std::string>( "outputFolder" );
	if ( !fs::exists( outputFolder_ ) ) {
		LOG( INFO ) << "Create output folder " << outputFolder_;
		fs::create_directories( outputFolder_ );
	} else if ( !fs::is_directory( outputFolder_ ) ) {
		LOG( FATAL ) << "Output folder " << outputFolder_ << " is not a directory!";
	}


	initialMeshOutput_ = ifr.getValue<bool>( "initialMeshOutput" );
	sizegridOutput_ = ifr.getValue<bool>( "sizegridOutput" );
	reductionOutput_ = ifr.getValue<bool>( "reductionOutput" );
	blockMeshOutput_ = ifr.getValue<bool>( "blockMeshOutput" );
	meshQualityOutput_ = ifr.getValue<bool>( "meshQualityOutput" );
	forceFragmentNumber_ = ifr.getValue<bool>( "forceFragmentNumber" );
	convexHullDecimation_ = ifr.getValue<bool>( "convexHullDecimation" );
}

void Config::log() const {
	LOG( INFO ) << "Working directory: " << workingDir_;
	LOG( INFO ) << "Config: " << configFile_;
	LOG( INFO ) << "Mesh: " << meshFile_;
	LOG( INFO ) << "Cache: " << cacheFolder_;
	LOG( INFO ) << "# Refinement Steps: " << nRefinementSteps_;
	LOG( INFO ) << "# blocks: " << nBlocks_;
	LOG( INFO ) << "Input # fragments: " << nFragments_;
	LOG( INFO ) << "Input Size Grid dimensions:    X = " << sizegridSizeX_ << "  |  Y = " << sizegridSizeY_;
	google::FlushLogFiles( google::GLOG_INFO );
}