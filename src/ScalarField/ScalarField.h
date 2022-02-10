#pragma once

#include <string>
#include <assert.h>
#include <vector>
#include <list>
#include <experimental/filesystem>
#include <optional>
#include <glog/logging.h>

#include <functional>

namespace ScalarField
{

	struct AxisAlignedBoundingBox
	{
		float xMin = 0, xMax = 0, yMin = 0, yMax = 0;
	};

	struct Point
	{
		float x = 0, y = 0;
		float& operator[]( size_t i ) {
			switch( i ) {
			case 0: return x;
			case 1: return y;
			default: LOG( ERROR ) << "Point[i] with i = " << i;
			}
		}
		float operator[]( size_t i ) const {
			switch( i ) {
			case 0: return x;
			case 1: return y;
			default: LOG( ERROR ) << "Point[i] with i = " << i;
			}
		}
	};


	class ScalarField
	{
	protected:
		AxisAlignedBoundingBox aabb_; // AABB of domain

		size_t Nx_ = 0, Ny_ = 0;	// Number of grid points
		float hx_ = 0, hy_ = 0;		// Stepsize

		std::vector<float> scalars_;

	public:
		//ScalarField() = delete;
		ScalarField() {}
		ScalarField( const std::experimental::filesystem::path& filename );
		ScalarField( const AxisAlignedBoundingBox& aabb, const size_t& Nx, const size_t& Ny );

		//ScalarField( const ScalarField& ) = delete;
		//ScalarField& operator=( const ScalarField& ) = delete;

		// Access functions
		float getScalar( float x, float y ) const;		// get size-value at a certain position
		float getScalar( const Point& p1, const Point& p2 ) const;
		float getMin( const Point& p1, const Point& p2 ) const;
		auto& scalars() { return scalars_; }
		const auto& Nx() const { return Nx_; }
		const auto& Ny() const { return Ny_; }
		const auto& hx() const { return hx_; }
		const auto& hy() const { return hy_; }
		const auto& aabb() const { return aabb_; };

		// Print functions
		void toVTK( const std::experimental::filesystem::path& name, const bool& isBinary = true ) const;

		// Access functions
		float& operator()( size_t i, size_t j ) { return scalars_[lex( i, j )]; }
		float operator()( size_t i, size_t j ) const { return scalars_[lex( i, j )]; }
		float& operator()( size_t i ) { return scalars_[i]; }
		float operator()( size_t i ) const { return scalars_[i]; }

		Point pos( size_t i, size_t j ) const { return { aabb_.xMin + i * hx_ , aabb_.yMin + j * hy_ }; }

		size_t lex( const size_t& i, const size_t& j ) const { 
			assert( i < Nx_&& j < Ny_ ); 
			return i * Ny_ + j; 
		}
		std::tuple<size_t, size_t> lex2ij( const size_t& i ) const {
			assert( i < scalars_.size() );
			return { i % Ny_, i / Ny_ };
		}

		size_t x2col( float x ) const {
			if( x > aabb_.xMax ) x = aabb_.xMax;
			if( x < aabb_.xMin ) x = aabb_.xMin;
			return (size_t)( ( x - aabb_.xMin ) / hx_ );
		}
		size_t y2row( float y ) const {
			if( y > aabb_.yMax ) y = aabb_.yMax;
			if( y < aabb_.yMin ) y = aabb_.yMin;
			return (size_t)( ( y - aabb_.yMin ) / hy_ );
		}

		void writeBinary( const std::experimental::filesystem::path& name ) const;
		void readBinary( const std::experimental::filesystem::path& name );

		float x( const size_t& i ) const { assert( i < Nx_ ); return aabb_.xMin + i * hx_; }
		float y( const size_t& j ) const { assert( j < Ny_ ); return aabb_.yMin + j * hy_; }
		
	private:

		// use jpg as backgroundgrid, set dimensions to [0,1]x[0,1]
		//void loadPicture( const std::experimental::filesystem::path& name );
	};


	inline float ScalarField::getScalar( float x, float y ) const {
		// this function was created according to the wikipedia article: "Bilinear interpolation"

		size_t i1 = this->x2col( x );
		size_t i2 = i1 + 1;
		if( i2 >= Nx_ ) {
			i1--;
			i2--;
		}
		size_t j1 = this->y2row( y );
		size_t j2 = j1 + 1;
		if( j2 >= Ny_ ) {
			j1--;
			j2--;
		}

		const float x1 = this->x( i1 );
		const float y1 = this->y( j1 );
		const float x2 = this->x( i2 );
		const float y2 = this->y( j2 );

		const float pref = 1.f / ( ( x2 - x1 ) * ( y2 - y1 ) );
		const float a = x2 - x;
		const float b = x - x1;
		const float c = ( *this )( i1, j1 );
		const float d = ( *this )( i1, j2 );
		const float e = ( *this )( i2, j1 );
		const float f = ( *this )( i2, j2 );
		const float g = y2 - y;
		const float h = y - y1;

		return pref * ( a * ( c * g + d * h ) + b * ( e * g + f * h ) );
	}
	inline float ScalarField::getScalar( const Point& p1, const Point& p2 ) const {
		float val = 0;
		const size_t n_steps = 10;
		for( size_t i = 0; i < n_steps; ++i ) {
			float x = ( ( n_steps - 1 - i ) * p1.x + i * p2.x ) / ( n_steps - 1 );
			float y = ( ( n_steps - 1 - i ) * p1.y + i * p2.y ) / ( n_steps - 1 );
			val += this->getScalar( x, y );
		}
		val /= (float)n_steps;

		return val;
	}
	inline float ScalarField::getMin( const Point& p1, const Point& p2 ) const {
		float val = std::numeric_limits<float>::max();
		const size_t n_steps = 10;
		for( size_t i = 0; i < n_steps; ++i ) {
			float x = ( ( n_steps - 1 - i ) * p1.x + i * p2.x ) / ( n_steps - 1 );
			float y = ( ( n_steps - 1 - i ) * p1.y + i * p2.y ) / ( n_steps - 1 );
			val = std::min( val, this->getScalar( x, y ) );
		}

		return val;
	}


	inline std::optional<ScalarField> load( const std::experimental::filesystem::path& file, const size_t& sizeGridSizeX, const size_t& sizeGridSizeY ) {
		namespace fs = std::experimental::filesystem;

		if( !fs::exists( file.parent_path() ) ) {
			LOG( INFO ) << "Directory '" << file.parent_path() << "' does not exist yet and is therefore now created";
			if( !fs::create_directories( file.parent_path() ) ) {
				LOG( ERROR ) << "Directory '" << file.parent_path() << "' cannot be created!";
			}
		}
		if( fs::exists( file ) ) {
			LOG( INFO ) << "Read SizeField from cache file '" << fs::canonical( file ) << "'";
			return ScalarField( file.string() );
		} else {
			return {};
		}
	}
}
