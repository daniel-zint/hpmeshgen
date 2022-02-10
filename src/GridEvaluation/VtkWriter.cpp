#include "VtkWriter.h"



void VtkWriter::meshToVtk() {
	using std::endl;

	ofs_ << "DATASET UNSTRUCTURED_GRID" << endl;
	ofs_ << "POINTS " << mesh_.n_vertices() << " float" << endl;
	for( const auto& v : mesh_.vertices() ) {
		const auto& p = mesh_.point( v );
		if( isBinary_ ) {
			float p0 = p[0];
			float p1 = p[1];
			//float p2 = p[2];
			float p2 = 0;
			swapEnd( p0 );
			swapEnd( p1 );
			swapEnd( p2 );
			ofs_.write( (char*)&p0, sizeof( float ) );
			ofs_.write( (char*)&p1, sizeof( float ) );
			ofs_.write( (char*)&p2, sizeof( float ) );
		} else {
			ofs_ << p[0] << " " << p[1] << " " << p[2] << "\n";
		}
	}
	if( isBinary_ )
		ofs_ << endl;
	ofs_ << "CELLS " << mesh_.n_faces() << " " << mesh_.n_faces() * 4 << endl;
	for( const auto& f : mesh_.faces() ) {
		std::vector<int> vids;
		for( const auto& v : f.vertices() ) {
			vids.push_back( v.idx() );
		}
		if( isBinary_ ) {
			int i0 = 3;
			int i1 = vids[0];
			int i2 = vids[1];
			int i3 = vids[2];
			swapEnd( i0 );
			swapEnd( i1 );
			swapEnd( i2 );
			swapEnd( i3 );
			ofs_.write( (char*)&i0, sizeof( int ) );
			ofs_.write( (char*)&i1, sizeof( int ) );
			ofs_.write( (char*)&i2, sizeof( int ) );
			ofs_.write( (char*)&i3, sizeof( int ) );
		} else {
			ofs_ << "3 " << vids[0] << " " << vids[1] << " " << vids[2] << "\n";
		}
	}
	if( isBinary_ )
		ofs_ << endl;
	ofs_ << "CELL_TYPES " << mesh_.n_faces() << endl;
	for( const auto& f : mesh_.faces() ) {
		if( isBinary_ ) {
			int i = 5;
			swapEnd( i );
			ofs_.write( (char*)&i, sizeof( int ) );
		} else {
			ofs_ << "5\n";
		}
	}
	if( isBinary_ )
		ofs_ << endl;
}
