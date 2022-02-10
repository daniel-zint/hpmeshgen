#pragma once

#include "MeshHeader.h"
#include "ScalarField/ScalarField.h"
//#include "ContourOld.h"
#include "Contour/Shape.h"

class TriRemeshing
{
	TriMesh mesh_;
	float lAvg_ = -1;
	float lMin_ = -1;
	float lMax_ = -1;
	bool hasBackgroundGrid_ = false;
	const ScalarField::ScalarField* bg_ = nullptr;
	bool hasContour = false;
	//const ContourOld* contourOld_ = nullptr;
	const Contour::Shape* contour_ = nullptr;
	const float maxFeatureSize = 150;
	bool convHullRemeshing = false;
public:
	TriRemeshing(const TriMesh& mesh);
	void set_backgroundgrid( const ScalarField::ScalarField* bg ) {
		bg_ = bg;
		hasBackgroundGrid_ = true;

		// compute lAvg_
		lAvg_ = 0;
		for( auto eh : mesh_.edges() ) {
			auto heh = mesh_.halfedge_handle( eh, 0 );
			auto vh1 = mesh_.from_vertex_handle( heh );
			auto vh2 = mesh_.to_vertex_handle( heh );
			auto p1 = mesh_.point( vh1 );
			auto p2 = mesh_.point( vh2 );
			auto l = mesh_.calc_edge_length( eh ) / bg_->getScalar( { p1[0],p1[1] }, { p2[0],p2[1] } );
			lAvg_ += l;
		}
		lAvg_ /= mesh_.n_edges();

		// lMin, lMax
		lMin_ = ( 3.f / 4.f ) * lAvg_;
		lMax_ = ( 4.f / 3.f ) * lAvg_;
	}
	void set_contour( const Contour::Shape* c ) {
		//contourOld_ = co;
		contour_ = c;
		hasContour = true;
	}
	void use_convex_hull_remeshing() {
		convHullRemeshing = true;
	}
	TriMesh remesh( const int nFragments = -1 );
private:
	/// Split edges longer than lMax. Iterate once through all edges.
	/// Return true if some edge was splitted, false otherwise.
	bool split();
	bool split( const int nFaces );
	bool collapse( const int nFaces = -1 );
	bool flip( const float qMin = -1 );
	void smooth();
};
