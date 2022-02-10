#include "ScalarField.h"

#include <algorithm>
#include <fstream>
#include <tuple>
#include <functional>
#include <iostream>

// For reading in pictures
//#if __has_include(<CImg.h>)
//#include <CImg.h>
//#endif

namespace ScalarField
{

	ScalarField::ScalarField( const std::experimental::filesystem::path& filename ) {
		if( !std::experimental::filesystem::exists( filename ) ) {
			LOG( ERROR ) << "File with name '" << filename.string() << "' does not exist!";
		}

		if( filename.extension() == ".bin" ) {
			readBinary( filename );
		} else if( filename.extension() == ".jpg" || filename.extension() == ".bmp" ) {
			//loadPicture( filename );
			LOG( ERROR ) << "Extension of file '" << filename.string() << "' is not compatible with ScalarField!";
		} else {
			LOG( ERROR ) << "Extension of file '" << filename.string() << "' is not compatible with ScalarField!";
		}
	}

	ScalarField::ScalarField( const AxisAlignedBoundingBox& aabb, const size_t& Nx, const size_t& Ny ) : aabb_{ aabb }, Nx_{ Nx }, Ny_{ Ny } {
		if( Nx_ < 2 || Ny_ < 2 ) {
			std::cerr << "Dimensions of scalar field must be larger than 1." << std::endl;
			std::cerr << "Nx = " << Nx_ << std::endl;
			std::cerr << "Ny = " << Ny_ << std::endl;
			return;
		}

		scalars_.resize( Nx_* Ny_ );

		for( auto& e : scalars_ ) e = -std::numeric_limits<float>::max();

		// Grid spacing
		hx_ = ( aabb_.xMax - aabb_.xMin ) / ( Nx_ - 1 );
		hy_ = ( aabb_.yMax - aabb_.yMin ) / ( Ny_ - 1 );

		//fillGrid( mesh, vertexAvg );
	}

	void ScalarField::toVTK( const std::experimental::filesystem::path& name, const bool& isBinary ) const {
		using std::endl;

		auto swapEnd = []( float& var ) {
			char* varArray = reinterpret_cast<char*>( &var );
			for( long i = 0; i < static_cast<long>( sizeof( var ) / 2 ); i++ )
				std::swap( varArray[sizeof( var ) - 1 - i], varArray[i] );
		};

		std::ofstream ofs;
		if( isBinary )
			ofs.open( name, std::ios::binary | std::ios::out );
		else
			ofs.open( name );

		ofs << "# vtk DataFile Version 3.0" << endl;
		ofs << name.stem().string() << endl;
		if( isBinary )
			ofs << "BINARY" << endl;
		else
			ofs << "ASCII" << endl;
		ofs << "DATASET STRUCTURED_POINTS" << endl;
		ofs << "DIMENSIONS " << Nx_ << " " << Ny_ << " " << 1 << endl;
		ofs << "ORIGIN " << aabb_.xMin << " " << aabb_.yMin << " " << 0 << endl;
		ofs << "SPACING " << hx_ << " " << hy_ << " " << 1 << endl;

		ofs << "POINT_DATA " << scalars_.size() << endl;
		ofs << "SCALARS value float 1" << endl;
		ofs << "LOOKUP_TABLE data_table" << endl;

		for( size_t j = 0; j < Ny_; ++j ) {
			for( size_t i = 0; i < Nx_; ++i ) {
				if( isBinary ) {
					float val;
					val = ( *this )( i, j );
					swapEnd( val );
					ofs.write( (char*)&val, sizeof( val ) );
				} 			else {
					ofs << ( *this )( i, j ) << "\n";
				}
			}
		}
		if( isBinary )
			ofs << "\n";

		ofs.close();
	}

	void ScalarField::writeBinary( const std::experimental::filesystem::path& name ) const {
		std::ofstream ofs( name, std::ios::binary | std::ios::out );
		if( !ofs.is_open() ) {
			std::cerr << "Cannot write to file with name '" << name << "'" << std::endl;
			std::cerr << "Line: " << __LINE__ << "\nFile: " << __FILE__ << std::endl;
		}

		auto write = [&ofs]( const auto& x ) {
			ofs.write( (char*)&x, sizeof( x ) );
		};

		write( Nx_ );
		write( Ny_ );
		write( aabb_.xMin );
		write( aabb_.xMax );
		write( aabb_.yMin );
		write( aabb_.yMax );
		ofs.write( (char*)&scalars_[0], scalars_.size() * sizeof( scalars_[0] ) );

		ofs.close();
	}

	void ScalarField::readBinary( const std::experimental::filesystem::path& name ) {
		if( !std::experimental::filesystem::exists( name ) ) {
			std::cerr << "File with name '" << name.string() << "' does not exist!" << std::endl;
			std::cerr << "Line: " << __LINE__ << "\nFile: " << __FILE__ << std::endl;
		}

		std::ifstream ifs( name, std::ios::binary | std::ios::in );
		if( !ifs.is_open() ) {
			std::cerr << "Cannot read file with name '" << name << "'" << std::endl;
			std::cerr << "Line: " << __LINE__ << "\nFile: " << __FILE__ << std::endl;
		}

		auto read = [&ifs]( const auto& x ) {
			ifs.read( (char*)&x, sizeof( x ) );
		};

		read( Nx_ );
		read( Ny_ );
		read( aabb_.xMin );
		read( aabb_.xMax );
		read( aabb_.yMin );
		read( aabb_.yMax );

		scalars_.resize( Nx_ * Ny_ );

		hx_ = ( aabb_.xMax - aabb_.xMin ) / ( Nx_ - 1 );
		hy_ = ( aabb_.yMax - aabb_.yMin ) / ( Ny_ - 1 );

		ifs.read( (char*)&scalars_[0], scalars_.size() * sizeof( scalars_[0] ) );

		ifs.close();
	}
/*
	void ScalarField::loadPicture( const std::experimental::filesystem::path& name ) {
	#ifndef SCALARFIELD_USING_CIMG
		LOG( ERROR ) << "The optional library CImg was not found";
	#else
		using namespace cimg_library;

		CImg<unsigned char> image( name.string().c_str() );

		Nx_ = image.width();
		Ny_ = image.height();
		aabb_.xMin = 0;
		aabb_.xMax = static_cast<float>( Nx_ );
		aabb_.yMin = 0;
		aabb_.yMax = static_cast<float>( Ny_ );

		scalars_.resize( Nx_ * Ny_ );

		hx_ = ( aabb_.xMax - aabb_.xMin ) / ( Nx_ - 1 );
		hy_ = ( aabb_.yMax - aabb_.yMin ) / ( Ny_ - 1 );

		auto minVal = std::min( hx_, hy_ );
		auto maxVal = std::max( aabb_.xMax - aabb_.xMin, aabb_.yMax - aabb_.yMin );

		for( int i = 0; i < Nx_; ++i ) {
			for( int j = 0; j < Ny_; ++j ) {
				//( *this )( i, j );
				auto r = image( i, static_cast<unsigned int>( Ny_ - 1 - j ), 0 );
				auto g = image( i, static_cast<unsigned int>( Ny_ - 1 - j ), 1 );
				auto b = image( i, static_cast<unsigned int>( Ny_ - 1 - j ), 2 );
				auto greyVal = ( r + g + b ) / 3;
				//( *this )( i, j ) = 255 - greyVal;
				auto u = greyVal / 255.f;
				( *this )( i, j ) = ( 1 - u ) * minVal + u * maxVal;

				//image( i, j, 0 ) = greyVal;
				//image( i, j, 1 ) = greyVal;
				//image( i, j, 2 ) = greyVal;
			}
		}

		//image.save( ( name.parent_path() / std::experimental::filesystem::path( "grey.bmp" ) ).string().c_str() );
	#endif // !SCALARFIELD_USING_CIMG
	}
*/
}