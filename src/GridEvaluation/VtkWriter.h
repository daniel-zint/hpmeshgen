#pragma once

#include <experimental/filesystem>
#include <fstream>
#include <vector>
#include <glog/logging.h>
#include "../MeshHeader.h"

class VtkWriter
{
public:
	template<typename T>
	struct FaceData
	{
		std::string name = "FaceData";
		std::vector<T> data;
	};

	template<typename T>
	struct PointData
	{
		std::string name = "PointData";
		std::vector<T> data;
	};

private:
	TriMesh& mesh_;

	std::ofstream ofs_;
	bool isBinary_;
	

public:
	VtkWriter( TriMesh& mesh, const std::experimental::filesystem::path& outputFile, const bool& isBinary ) : mesh_{ mesh }, isBinary_{ isBinary } {
		
		if( isBinary_ )
			ofs_.open( outputFile, std::ios::binary | std::ios::out );
		else
			ofs_.open( outputFile );

		ofs_ << "# vtk DataFile Version 3.0" << std::endl;
		ofs_ << outputFile.stem().string() << std::endl;
		if( isBinary_ )
			ofs_ << "BINARY" << std::endl;
		else
			ofs_ << "ASCII" << std::endl;

		meshToVtk();
	}

	~VtkWriter() {
		ofs_.close();
	}

	template<typename T>
	void addFaceData( const FaceData<T>& fd );

	template<typename T>
	void addPointData( const PointData<T>& pd );
		
private:
	void meshToVtk();
	
};

template<typename T>
inline void swapEnd( T& var ) {
	char* varArray = reinterpret_cast<char*>( &var );
	for( long i = 0; i < static_cast<long>( sizeof( var ) / 2 ); i++ )
		std::swap( varArray[sizeof( var ) - 1 - i], varArray[i] );
}

template<typename T>
inline void VtkWriter::addFaceData( const FaceData<T>& fd ) {
	ofs_ << "CELL_DATA " << mesh_.n_faces() << std::endl;
	ofs_ << "SCALARS " << fd.name << " " << typeid( T ).name() << " 1" << std::endl;
	ofs_ << "LOOKUP_TABLE default" << std::endl;
	for( const auto& d : fd.data ) {
		if( isBinary_ ) {
			T val = d;
			swapEnd( val );
			ofs_.write( (char*)&val, sizeof( T ) );
		} else {
			ofs_ << d << "\n";
		}
	}
	if( isBinary_ )
		ofs_ << "\n";
}

template<typename T>
inline void VtkWriter::addPointData( const PointData<T>& pd ) {
	ofs_ << "POINT_DATA " << mesh_.n_vertices() << std::endl;
	ofs_ << "SCALARS " << pd.name << " " << typeid( T ).name() << " 1" << std::endl;
	ofs_ << "LOOKUP_TABLE default" << std::endl;
	for( const auto& d : pd.data ) {
		if( isBinary_ ) {
			T val = d;
			swapEnd( val );
			ofs_.write( (char*)&val, sizeof( T ) );
		} else {
			ofs_ << d << "\n";
		}
	}
	if( isBinary_ )
		ofs_ << "\n";
}
