#pragma once

#include <glog/logging.h>

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

#include "MeshHeader.h"

namespace DMO
{
	enum Set { Inner, Boundary, All };

	class DmoMesh
	{		
	public:
		thrust::host_vector<DmoVertex> vertices;		// vertices considered for smoothing
		thrust::device_vector<DmoVertex> vertices_d;
		thrust::host_vector<int> oneRingVec;			// neighborhood of  vertices
		thrust::device_vector<int> oneRingVec_d;
		thrust::host_vector<int> coloredVertexIDs;		// vertex coloring
		thrust::device_vector<int> coloredVertexIDs_d;

		std::vector<int> colorOffset_;

		DmoMesh() = default;
		DmoMesh( const DmoMesh& m ) {
			vertices = m.vertices;
			oneRingVec = m.oneRingVec;
			coloredVertexIDs = m.coloredVertexIDs;
			colorOffset_ = m.colorOffset_;
		}
		DmoMesh( const std::vector<OpenMesh::SmartVertexHandle>& vhs ) {
			copyMeshData( vhs );
			createColoring( vhs );
		}

		template <int set = Set::Inner> static DmoMesh create( TriMesh& mesh )
		{
			DmoMesh m;
			m.copyMeshData<set>( mesh );
			m.createColoring<set>( mesh );
			return m;
		}

		DmoMesh& operator= ( const DmoMesh& m ) {
			vertices = m.vertices;
			oneRingVec = m.oneRingVec;
			coloredVertexIDs = m.coloredVertexIDs;
			colorOffset_ = m.colorOffset_;
			return *this;
		}

		void copyHostToDevice() {
			vertices_d = vertices;
			oneRingVec_d = oneRingVec;
			coloredVertexIDs_d = coloredVertexIDs;
		}

	private:
		
		void copyMeshData( const std::vector<OpenMesh::SmartVertexHandle>& vhs ) {

			vertices.resize( vhs.size() );
			// compute oneRing size
			int oneRingVecLength = 0;
			for( const auto& vh : vhs ) {
				oneRingVecLength += vh.valence();
				if( !vh.is_boundary() )
					oneRingVecLength++; // additional count s.th. last element is again the first element
			}
			oneRingVec.resize( oneRingVecLength );

			int freeVertexCounter = 0;
			int oneRingCounter = 0;
			for( const auto& vh : vhs ) {

				DMO::DmoVertex& v = vertices[freeVertexCounter++];
				v.idx = vh.idx();
				v.oneRingSize = vh.valence();
				if( !vh.is_boundary() )
					v.oneRingSize++;
				v.oneRingID = oneRingCounter;

				auto heh = vh.out().prev().opp();
				const auto heh_init = heh;
				do {
					oneRingVec[oneRingCounter++] = heh.to().idx();
					heh = heh.next().next().opp();
				} while( !heh.is_boundary() && heh != heh_init );
				oneRingVec[oneRingCounter++] = heh.to().idx();
			}

		}

		template<int set>
		void copyMeshData( TriMesh& mesh ) {

			int nVertices = 0;
			int oneRingVecLength = 0;

			for( const auto& vh : mesh.vertices() ) {
				if( set == Set::Inner && vh.is_boundary() )
					continue;
				else if( set == Set::Boundary && !vh.is_boundary() )
					continue;

				nVertices++;
				oneRingVecLength += vh.valence();
				if( !vh.is_boundary() )
					oneRingVecLength++; // additional count s.th. last element is again the first element
			}

			vertices.resize( nVertices );
			oneRingVec.resize( oneRingVecLength );

			int freeVertexCounter = 0;
			int oneRingCounter = 0;
			for( const auto& vh : mesh.vertices() ) {
				if( set == Set::Inner && vh.is_boundary() )
					continue;
				else if( set == Set::Boundary && !vh.is_boundary() )
					continue;

				DMO::DmoVertex& v = vertices[freeVertexCounter++];
				v.idx = vh.idx();
				v.oneRingSize = vh.valence();
				if( !vh.is_boundary() )
					v.oneRingSize++;
				v.oneRingID = oneRingCounter;

				auto heh = vh.out().prev().opp();
				const auto heh_init = heh;
				do {
					oneRingVec[oneRingCounter++] = heh.to().idx();
					heh = heh.next().next().opp();
				} while( !heh.is_boundary() && heh != heh_init );
				oneRingVec[oneRingCounter++] = heh.to().idx();
			}
		}

		void createColoring( const std::vector<OpenMesh::SmartVertexHandle>& vhs ) {

			const auto& mesh = vhs[0].mesh();

			// create coloring scheme
			std::vector<int>colorScheme( mesh->n_vertices(), -2 );

			for( const auto& vh : vhs ) {
				colorScheme[vh.idx()] = -1;
			}

			for( const auto& vh : vhs ) {			
				unsigned long colorBits = 0;
				for( auto vv : vh.vertices() ) {
					int c = colorScheme[vv.idx()];
					if( c >= 0 )
						colorBits |= 1 << c;
				}
				int color = 0;
				while( ( colorBits & ( 1 << color ) ) ) {
					++color;
				}
				colorScheme[vh.idx()] = color;
			}

			int n_colors = *( std::max_element( colorScheme.begin(), colorScheme.end() ) ) + 1;

			if( n_colors == -1 )
				return;

			std::vector<int> n_color_vecs( n_colors, 0 );
			for( int i = 0; i < colorScheme.size(); ++i ) {
				if( colorScheme[i] > -1 )
					++n_color_vecs[colorScheme[i]];
			}

			coloredVertexIDs.resize( vhs.size() );

			colorOffset_.resize( n_colors + 1, 0 );
			for( int i = 1; i < n_colors; ++i ) {
				colorOffset_[i] = colorOffset_[i - 1] + n_color_vecs[i - 1];
			}
			colorOffset_[n_colors] = vhs.size();		// mark the end of the colored-vertices vector

			// add vertex ids
			std::vector<int>colorCounter( n_colors, 0 );
			int interior_counter = 0;
			for( int i = 0; i < colorScheme.size(); ++i ) {
				if( colorScheme[i] < 0 ) { continue; }
				coloredVertexIDs[colorOffset_[colorScheme[i]] + colorCounter[colorScheme[i]]++] = interior_counter++;
			}
		}

		template<int set>
		void createColoring( TriMesh& mesh ) {
		
			// create coloring scheme
			std::vector<int>colorScheme( mesh.n_vertices(), -2 );
		
			for( const auto& vh : mesh.vertices() ) {
				if( set == Set::Inner && vh.is_boundary() )
					continue;
				else if( set == Set::Boundary && !vh.is_boundary() )
					continue;

				colorScheme[vh.idx()] = -1;
			}
		
			for( const auto& vh : mesh.vertices() ) {
				if( set == Set::Inner && vh.is_boundary() )
					continue;
				else if( set == Set::Boundary && !vh.is_boundary() )
					continue;
		
				unsigned long colorBits = 0;
				for( auto vv : vh.vertices() ) {
					int c = colorScheme[vv.idx()];
					if( c >= 0 )
						colorBits |= 1 << c;
				}
				int color = 0;
				while( ( colorBits & ( 1 << color ) ) ) {
					++color;
				}
				colorScheme[vh.idx()] = color;
			}
		
			int n_colors = *( std::max_element( colorScheme.begin(), colorScheme.end() ) ) + 1;
		
			if( n_colors == -1 )
				return;
		
			std::vector<int> n_color_vecs( n_colors, 0 );
			for( int i = 0; i < colorScheme.size(); ++i ) {
				if( colorScheme[i] > -1 )
					++n_color_vecs[colorScheme[i]];
			}
		
			coloredVertexIDs.resize( vertices.size() );
		
			colorOffset_.resize( n_colors + 1, 0 );
			for( int i = 1; i < n_colors; ++i ) {
				colorOffset_[i] = colorOffset_[i - 1] + n_color_vecs[i - 1];
			}
			colorOffset_[n_colors] = vertices.size();		// mark the end of the colored-vertices vector
		
			// add vertex ids
			std::vector<int>colorCounter( n_colors, 0 );
			int interior_counter = 0;
			for( int i = 0; i < colorScheme.size(); ++i ) {
				if( colorScheme[i] < 0 ) { continue; }
				coloredVertexIDs[colorOffset_[colorScheme[i]] + colorCounter[colorScheme[i]]++] = interior_counter++;
			}
		}
	};
}
