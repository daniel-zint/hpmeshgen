#pragma once

#include <experimental/filesystem>

#include <glog/logging.h>

class Config
{
	using fsp = std::experimental::filesystem::path;
	fsp configFile_;
	fsp workingDir_;
	fsp cacheFolder_;
	size_t nBlocks_;
	fsp meshFile_;
	size_t nRefinementSteps_;
	size_t nFragments_;
	size_t sizegridSizeX_;
	size_t sizegridSizeY_;
	fsp outputFolder_;
	bool initialMeshOutput_;
	bool sizegridOutput_;
	bool reductionOutput_;
	bool blockMeshOutput_;
	bool meshQualityOutput_;
	bool forceFragmentNumber_;
	bool convexHullDecimation_;

public:
	Config( const fsp& configFile );

	void log() const;

	auto configFile() const { return configFile_; }
	auto workingDir() const { return workingDir_; }
	auto cacheFolder() const { return cacheFolder_; }
	auto nBlocks() const { return nBlocks_; }
	auto meshFile() const { return meshFile_; }
	auto nRefinementSteps() const { return nRefinementSteps_; }
	auto nPatches() const { return nFragments_; }
	auto sizegridSizeX() const { return sizegridSizeX_; }
	auto sizegridSizeY() const { return sizegridSizeY_; }
	auto outputFolder() const { return outputFolder_; }
	auto initialMeshOutput() const { return initialMeshOutput_; }
	auto sizegridOutput() const { return sizegridOutput_; }
	auto reductionOutput() const { return reductionOutput_; }
	auto blockMeshOutput() const { return blockMeshOutput_; }
	auto meshQualityOutput() const { return meshQualityOutput_; }
	auto forceFragmentNumber() const { return forceFragmentNumber_; }
	auto convexHullDecimation() const { return convexHullDecimation_; }
};