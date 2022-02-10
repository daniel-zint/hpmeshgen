#pragma once

#include "MeshHeader.h"
#include "MeshFunctions.h"
#include "MeshSmoothing/QualityMetrics.h"
#include "HelperFunctions.h"
#include "ScalarField/ScalarField.h"

#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Tools/Utils/HeapT.hh>
#include <OpenMesh/Tools/Decimater/BaseDecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>

#include <type_traits>
#include <algorithm>


//== NAMESPACE ================================================================

namespace OpenMesh
{ // BEGIN_NS_OPENMESH
	namespace Decimater
	{ // BEGIN_NS_DECIMATER

		struct PointDistErrorMetric
		{
			std::vector<OpenMesh::Vec3f> pVec;
			const ScalarField::ScalarField* bg = nullptr;
			bool hasBackgroundGrid = false;

			PointDistErrorMetric() {}
			PointDistErrorMetric( decltype( pVec ) pVeci ) : pVec{ pVeci } {}
			PointDistErrorMetric( decltype( pVec ) pVeci, decltype( bg ) bgi ) : pVec{ pVeci }, bg{ bgi }, hasBackgroundGrid{ true } {}
			PointDistErrorMetric( const OpenMesh::Vec3f& p ) : pVec{ p } {}
			PointDistErrorMetric( const OpenMesh::Vec3f& p, decltype( bg ) bgi ) : pVec{ p }, bg{ bgi }, hasBackgroundGrid{ true } {}

			PointDistErrorMetric operator+( const PointDistErrorMetric& pd ) const {
				auto pVecRet{ pVec };
				pVecRet.insert( pVecRet.end(), pd.pVec.begin(), pd.pVec.end() );
				if ( hasBackgroundGrid )
					return { pVecRet, bg };
				else
					return { pVecRet };
			}

			PointDistErrorMetric& operator+=( const PointDistErrorMetric& pd ) {
				pVec.insert( pVec.end(), pd.pVec.begin(), pd.pVec.end() );
				return *this;
			}

			float eval( const OpenMesh::Vec3f& x ) const {
				double q = 0;
				if ( hasBackgroundGrid ) {
					for ( const auto& p : pVec ) {
						auto s = bg->getScalar( { p[0],p[1] }, { x[0],x[1] } );
						s = s * s;
						auto l = ( p - x ).sqrnorm();
						q += l / s;
					}
				} else {
					for ( const auto& p : pVec ) {
						q += ( p - x ).sqrnorm();
					}
				}
				return (float)q;
			}

			float operator*( const OpenMesh::Vec3f& x ) const {
				return eval( x );
			}
		};

		//== CLASS DEFINITION =========================================================

		template <class MeshT>
		class ModPointDistT : public ModBaseT<MeshT>
		{
		public:
			DECIMATING_MODULE( ModPointDistT, MeshT, PointDistanceModule );
			
			using value_type = vector_traits<TriMesh::Point>::value_type;
		public:

			/// Constructor
			ModPointDistT( MeshT &_dec ) :
				Base( _dec, false ) {
				Base::mesh().add_property( vErrors_ );
				Base::mesh().add_property( ePoints_ );
			}

			/// Destructor
			~ModPointDistT() {
				Base::mesh().remove_property( vErrors_ );
				Base::mesh().remove_property( ePoints_ );
			}

		public: // inherited

			/// Initalize the module and prepare the mesh for decimation.
			virtual void initialize( void ) {
				auto& mesh = Base::mesh();

				if ( !vErrors_.is_valid() ) {
					mesh.add_property( vErrors_ );
					mesh.add_property( ePoints_ );
				}

				// init errors
				if ( hasBackgroundGrid ) {
					for ( auto vh : mesh.vertices() ) {
						mesh.property( vErrors_, vh ) = PointDistErrorMetric( mesh.point( vh ), bg_ );
					}
				} else {
					for ( auto vh : mesh.vertices() ) {
						mesh.property( vErrors_, vh ) = PointDistErrorMetric( mesh.point( vh ) );
					}
				}
				// init edge points
				for ( auto eh : mesh.edges() ) {
					auto heh = mesh.halfedge_handle( eh, 0 );
					mesh.property(ePoints_, eh) = mesh.point( mesh.to_vertex_handle( heh ) );
				}
				
				Base::set_binary( false );
			}

			float collapse_priority( const CollapseInfo& _ci ) {
				
				auto& mesh = Base::mesh();

				OpenMesh::VPropHandleT<float> feature_size;
				mesh.get_property_handle( feature_size, "feature_size" );

				float priority = Base::ILLEGAL_COLLAPSE;

				if ( true ) // continuous mode
				{
					const auto& v0 = _ci.v0;
					const auto& v1 = _ci.v1;
					const auto& eh = mesh.edge_handle( _ci.v0v1 );
										
					auto& p = mesh.property( ePoints_, eh );

					// jump over edges that connect two boundary vertices but are themselves interior (don't reduce canals)
					if ( mesh.is_boundary( v0 ) && mesh.is_boundary( v1 ) && !mesh.is_boundary( eh ) ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}

					// jump over edge that connects two features
					if( mesh.data( v0 ).is_feature && mesh.data( v1 ).is_feature ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}
					if ( mesh.property( feature_size, v0 ) == 0.f && mesh.property( feature_size, v1 ) == 0.f ) {
						LOG( ERROR ) << " Executing deprecated code!!";
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}

					// do not reduce boundarys further than 3 edges
					// DEPRECATED is checked already
					if ( mesh.is_boundary( eh ) ) {
						auto heh = _ci.v0v1;
						if ( !mesh.is_boundary( heh ) ) {
							heh = _ci.v1v0;
						}

						auto hehInit = heh;
						size_t nBoundaryEdges = 0;
						do {
							++nBoundaryEdges;
							heh = mesh.next_halfedge_handle( heh );
						} while ( heh != hehInit && nBoundaryEdges < 4 );

						if ( nBoundaryEdges <= 3 ) {
							LOG( ERROR ) << " Executing deprecated code!!";
							priority = Base::ILLEGAL_COLLAPSE;
							return priority;
						}
					}

					p = computeEdgePoint( _ci );
										
					// assure minimal quality
					if ( minMeanRatioMetric( mesh, v0, v1, p ) < min_r_ ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}
										
					auto e1 = mesh.property( vErrors_, _ci.v0 );
					auto e2 = mesh.property( vErrors_, _ci.v1 );
					auto e = e1 + e2;
					priority = e * p;

				} else // binary mode
				{
					LOG( ERROR ) << "Binary mode for point distance not implemented";
					priority = Base::ILLEGAL_COLLAPSE;
				}

				return priority;
			}
			
			void postprocess_collapse( const CollapseInfoT< MeshT >& _ci ) {
				auto& mesh = Base::mesh();
				mesh.property( vErrors_, _ci.v1 ) += mesh.property( vErrors_, _ci.v0 );
				
				auto pNew = mesh.property( ePoints_, mesh.edge_handle( _ci.v0v1 ) );
				mesh.set_point( _ci.v1, pNew );

				OpenMesh::VPropHandleT<float> feature_size;
				mesh.get_property_handle( feature_size, "feature_size" );
				
				mesh.property( feature_size, _ci.v1 ) = std::min( mesh.property( feature_size, _ci.v0 ), mesh.property( feature_size, _ci.v1 ) );
				mesh.data( _ci.v1 ).is_feature |= mesh.data( _ci.v0 ).is_feature;
			}

			void set_backgroundgrid( const ScalarField::ScalarField* bg ) {
				bg_ = bg;
				hasBackgroundGrid = true;
			}

			void set_min_meanratio( value_type _min_meanratio, bool _binary = true ) {
				assert( 0.0 <= _min_meanratio && _min_meanratio <= 1.0 );
				min_r_ = _min_meanratio;
			}

		private:
			OpenMesh::VPropHandleT<PointDistErrorMetric> vErrors_;
			OpenMesh::EPropHandleT<OpenMesh::Vec3f> ePoints_;	// store target positions
			const ScalarField::ScalarField* bg_;
			bool hasBackgroundGrid = false;

			value_type min_r_ = -1.f;

			inline OpenMesh::Vec3f computeEdgePoint( const CollapseInfo& _ci ) {

				auto& mesh = Base::mesh();

				OpenMesh::VPropHandleT<float> feature_size;
				mesh.get_property_handle( feature_size, "feature_size" );
				
				const auto& v0 = _ci.v0;
				const auto& v1 = _ci.v1;
				const auto& eh = mesh.edge_handle( _ci.v0v1 );

				const auto b0 = mesh.is_boundary( v0 );
				const auto b1 = mesh.is_boundary( v1 );

				OpenMesh::Vec3f p( 0,0,0 );

				// calculate new point
				if ( b0 && b1 ) {

					float angle1 = mesh.property( feature_size, v0 );
					float angle2 = mesh.property( feature_size, v1 );

					if ( angle1 < 100 && angle1 < angle2 ) {
						p = _ci.p0;
					} else if ( angle2 < 100 && angle2 <= angle1 ) {
						p = _ci.p1;
					} else {

						auto heh = _ci.v0v1;
						if ( !mesh.is_boundary( heh ) ) {
							heh = _ci.v1v0;
						}

						std::vector<TriMesh::VertexHandle> vhVec( 4 );
						vhVec[0] = mesh.from_vertex_handle( mesh.prev_halfedge_handle( heh ) );
						vhVec[1] = mesh.from_vertex_handle( heh );
						vhVec[2] = mesh.to_vertex_handle( heh );
						vhVec[3] = mesh.to_vertex_handle( mesh.next_halfedge_handle( heh ) );

						std::vector<TriMesh::Point> pVec;
						pVec.reserve( 4 );
						for ( auto vh : vhVec ) {
							pVec.push_back( mesh.point( vh ) );
						}

						HelperFunctions::project4areaConsistency( pVec, p );
					}

					//p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );
				} else if ( b0 && !b1 ) {
					p = _ci.p0;
				} else if ( !b0 && b1 ) {
					p = _ci.p1;
				} else {
					p = 0.5f * ( _ci.p0 + _ci.p1 );
				}

				return p;
			}

			inline float minMeanRatioMetric( MeshT& mesh, const TriMesh::VertexHandle& vh1, const TriMesh::VertexHandle& vh2, const TriMesh::Point& p_new ) {
				float minQ = FLT_MAX;
				
				std::vector<TriMesh::VertexHandle> vh1fan;
				MeshFunctions::getOneRing( mesh, vh1, vh1fan );

				for ( size_t i = 0; i < vh1fan.size() - 1; ++i ) {
					size_t j = i + 1;
					if ( vh1fan[i] == vh2 || vh1fan[j] == vh2 ) {
						continue;
					}
					TriMesh::Point p1 = mesh.point( vh1fan[i] );
					TriMesh::Point p2 = mesh.point( vh1fan[j] );

					const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

					minQ = std::min( minQ, quality );

				}

				std::vector<TriMesh::VertexHandle> vh2fan;
				MeshFunctions::getOneRing( mesh, vh2, vh2fan );

				for ( size_t i = 0; i < vh2fan.size() - 1; ++i ) {
					size_t j = i + 1;
					if ( vh2fan[i] == vh1 || vh2fan[j] == vh1 ) {
						continue;
					}
					TriMesh::Point p1 = mesh.point( vh2fan[i] );
					TriMesh::Point p2 = mesh.point( vh2fan[j] );

					const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

					minQ = std::min( minQ, quality );
				}
				
				return minQ;
			}
		};


		//=============================================================================

		//== CLASS DEFINITION =========================================================

		template <class MeshT>
		class ModMeanRatioT : public ModBaseT<MeshT>
		{
		public:
			DECIMATING_MODULE( ModMeanRatioT, MeshT, MeanRatioModule );

			using value_type = vector_traits<TriMesh::Point>::value_type;

		public:

			/// Constructor
			ModMeanRatioT( MeshT &_dec ) :
				Base( _dec, false ) {
				static_assert( std::is_same<MeshT, TriMesh>::value, "ModMeanRatio only works with type TriMesh" );
			}

			/// Destructor
			~ModMeanRatioT() {
			}

		public: // inherited

			/// Initalize the module and prepare the mesh for decimation.
			virtual void initialize( void ) {
			}

			float collapse_priority( const CollapseInfo& _ci ) {

				auto& mesh = Base::mesh();

				double priority = minMeanRatioMetric( mesh, _ci.v0, _ci.v1, _ci.p0 );

				if ( min_r_ < 0.f ) // continuous mode
				{

				} else // binary mode
				{
					priority = ( priority < min_r_ ) ? Base::ILLEGAL_COLLAPSE : Base::LEGAL_COLLAPSE;
				}

				return (float)priority;
			}

			void set_min_meanratio( value_type _min_meanratio, bool _binary = true ) {
				assert( 0.0 <= _min_meanratio && _min_meanratio <= 1.0 );
				min_r_ = _min_meanratio;
				Base::set_binary( _binary );
			}

		private:
			value_type min_r_ = -1.f;

			inline float minMeanRatioMetric( MeshT& mesh, const TriMesh::VertexHandle& vh1, const TriMesh::VertexHandle& vh2, const TriMesh::Point& p_new ) {
				float minQ = FLT_MAX;
				
				std::vector<TriMesh::VertexHandle> vh1fan;
				MeshFunctions::getOneRing( mesh, vh1, vh1fan );

				for ( size_t i = 0; i < vh1fan.size() - 1; ++i ) {
					size_t j = i + 1;
					if ( vh1fan[i] == vh2 || vh1fan[j] == vh2 ) {
						continue;
					}
					TriMesh::Point p1 = mesh.point( vh1fan[i] );
					TriMesh::Point p2 = mesh.point( vh1fan[j] );

					const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

					minQ = std::min( minQ, quality );

				}

				std::vector<TriMesh::VertexHandle> vh2fan;
				MeshFunctions::getOneRing( mesh, vh2, vh2fan );

				for ( size_t i = 0; i < vh2fan.size() - 1; ++i ) {
					size_t j = i + 1;
					if ( vh2fan[i] == vh1 || vh2fan[j] == vh1 ) {
						continue;
					}
					TriMesh::Point p1 = mesh.point( vh2fan[i] );
					TriMesh::Point p2 = mesh.point( vh2fan[j] );

					const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

					minQ = std::min( minQ, quality );
				}

				return minQ;
			}
		};


		//=============================================================================

		//== CLASS DEFINITION =========================================================

		template <class MeshT>
		class ModEdgeLengthT : public ModBaseT<MeshT>
		{
		public:
			DECIMATING_MODULE( ModEdgeLengthT, MeshT, EdgeLengthModule );

			using value_type = vector_traits<TriMesh::Point>::value_type;
		public:

			/// Constructor
			ModEdgeLengthT( MeshT& _dec ) :
				Base( _dec, false ) {
				Base::mesh().add_property( ePoints_ );
			}

			/// Destructor
			~ModEdgeLengthT() {
				Base::mesh().remove_property( ePoints_ );
			}

		public: // inherited

			/// Initalize the module and prepare the mesh for decimation.
			virtual void initialize( void ) {
				auto& mesh = Base::mesh();

				if( !ePoints_.is_valid() ) {
					mesh.add_property( ePoints_ );
				}

				// init edge points
				for( auto eh : mesh.edges() ) {
					auto heh = mesh.halfedge_handle( eh, 0 );
					mesh.property( ePoints_, eh ) = mesh.point( mesh.to_vertex_handle( heh ) );
				}
	
				Base::set_binary( false );
			}

			float collapse_priority( const CollapseInfo& _ci ) {

				auto& mesh = Base::mesh();

				OpenMesh::VPropHandleT<float> feature_size;
				mesh.get_property_handle( feature_size, "feature_size" );

				float priority = Base::ILLEGAL_COLLAPSE;

				if( true ) // continuous mode
				{
					const auto& v0 = _ci.v0;
					const auto& v1 = _ci.v1;
					const auto& eh = mesh.edge_handle( _ci.v0v1 );

					auto& p = mesh.property( ePoints_, eh );
					
					// jump over edges that connect two boundary vertices but are themselves interior (don't reduce canals)
					if( mesh.is_boundary( v0 ) && mesh.is_boundary( v1 ) && !mesh.is_boundary( eh ) ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}

					// jump over edge that connects two features
					if( mesh.data( v0 ).is_feature && mesh.data( v1 ).is_feature ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}
					if( mesh.property( feature_size, v0 ) == 0.f && mesh.property( feature_size, v1 ) == 0.f ) {
						LOG( ERROR ) << " Executing deprecated code!!";
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}
									
					p = computeEdgePoint( _ci );

					// assure minimal quality
					if( minMeanRatioMetric( mesh, v0, v1, p ) < min_r_ ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}

					value_type l;
					if( hasBackgroundGrid ) {
						l = ( mesh.point( v0 ) - mesh.point( v1 ) ).length() / MeshFunctions::getNeighMinSize( *bg_, mesh, eh );
					} else {
						l = ( mesh.point( v0 ) - mesh.point( v1 ) ).length();
					}

					// assure maximal length
					if( l > max_l_ ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}

					priority = l;

				} else // binary mode
				{
					LOG( ERROR ) << "Binary mode for point distance not implemented";
					priority = Base::ILLEGAL_COLLAPSE;
				}

				return priority;
			}

			void postprocess_collapse( const CollapseInfoT< MeshT >& _ci ) {
				auto& mesh = Base::mesh();

				auto pNew = mesh.property( ePoints_, mesh.edge_handle( _ci.v0v1 ) );
				mesh.set_point( _ci.v1, pNew );

				OpenMesh::VPropHandleT<float> feature_size;
				mesh.get_property_handle( feature_size, "feature_size" );
				
				auto& f0 = mesh.property( feature_size, _ci.v0 );
				auto& f1 = mesh.property( feature_size, _ci.v1 );

				if( mesh.is_boundary( _ci.v0 ) && !mesh.is_boundary( _ci.v1 ) ) {
					f1 = f0;
				} else if( !mesh.is_boundary( _ci.v0 ) && mesh.is_boundary( _ci.v1 ) ) {
					// do nothing. Just keep properties of v1
					//f1 = f1;
					//c1 = c1;
				} else {
					f1 = std::min( f0, f1 );
				}
				mesh.data( _ci.v1 ).is_feature |= mesh.data( _ci.v0 ).is_feature;

			}

			void set_backgroundgrid( const ScalarField::ScalarField* bg ) {
				bg_ = bg;
				hasBackgroundGrid = true;
			}

			void set_min_meanratio( value_type _min_meanratio, bool _binary = true ) {
				assert( _min_meanratio <= 1.0 );
				min_r_ = _min_meanratio;
			}

			void set_max_length( value_type _max_length, bool _binary = true ) {
				assert( 0.0 <= _max_length );
				max_l_ = _max_length;
			}

			void set_max_feature_size( value_type _max_feature_size, bool _binary = true ) {
				max_feature_size = _max_feature_size;
			}

			void use_convex_hull_decimation() {
				convHullDecimation = true;
			}

		private:
			OpenMesh::EPropHandleT<OpenMesh::Vec3f> ePoints_;	// store target positions
			const ScalarField::ScalarField* bg_;
			bool hasBackgroundGrid = false;
			bool convHullDecimation = false;

			value_type min_r_ = -FLT_MAX;
			value_type max_l_ = FLT_MAX;
			value_type max_feature_size = 0.f;

			inline OpenMesh::Vec3f computeEdgePoint( const CollapseInfo& _ci ) {

				auto& mesh = Base::mesh();

				OpenMesh::VPropHandleT<float> feature_size;
				mesh.get_property_handle( feature_size, "feature_size" );

				const auto& v0 = _ci.v0;
				const auto& v1 = _ci.v1;
				const auto& eh = mesh.edge_handle( _ci.v0v1 );

				OpenMesh::Vec3f p( 0, 0, 0 );

				const auto b0 = mesh.is_boundary( v0 );
				const auto b1 = mesh.is_boundary( v1 );

				// calculate new point
				if( b0 && b1 ) {

					if( !convHullDecimation ) {
						float angle1 = mesh.property( feature_size, v0 );
						float angle2 = mesh.property( feature_size, v1 );

						if( angle1 < max_feature_size && angle1 < angle2 ) {
							p = _ci.p0;
						} else if( angle2 < max_feature_size && angle2 <= angle1 ) {
							p = _ci.p1;
						} else {

							auto heh = _ci.v0v1;
							if( !mesh.is_boundary( heh ) ) {
								heh = _ci.v1v0;
							}

							std::vector<TriMesh::VertexHandle> vhVec( 4 );
							vhVec[0] = mesh.from_vertex_handle( mesh.prev_halfedge_handle( heh ) );
							vhVec[1] = mesh.from_vertex_handle( heh );
							vhVec[2] = mesh.to_vertex_handle( heh );
							vhVec[3] = mesh.to_vertex_handle( mesh.next_halfedge_handle( heh ) );

							std::vector<TriMesh::Point> pVec;
							pVec.reserve( 4 );
							for( auto vh : vhVec ) {
								pVec.push_back( mesh.point( vh ) );
							}

							HelperFunctions::project4areaConsistency( pVec, p );
						}
					} else {
						// compute boundary angles
						auto heh = _ci.v0v1;
						if( !mesh.is_boundary( heh ) )
							heh = _ci.v1v0;

						std::vector<TriMesh::VertexHandle> vhs( 4 );
						vhs[0] = mesh.from_vertex_handle( mesh.prev_halfedge_handle( heh ) );
						vhs[1] = mesh.from_vertex_handle( heh );
						vhs[2] = mesh.to_vertex_handle( heh );
						vhs[3] = mesh.to_vertex_handle( mesh.next_halfedge_handle( heh ) );

						std::vector<TriMesh::Point> pVec( 4 );
						std::transform( vhs.begin(), vhs.end(), pVec.begin(), [&mesh]( const auto& v ) { return mesh.point( v ); } );

						auto angle1 = 180.f - HelperFunctions::calcAngle2( pVec[0], pVec[1], pVec[2] );
						auto angle2 = 180.f - HelperFunctions::calcAngle2( pVec[1], pVec[2], pVec[3] );

						if( std::abs( angle1 + angle2 ) > 170.f ) {
							//p = { FLT_MAX,FLT_MAX,FLT_MAX };
							p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );
							return p;
						}

						//p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );

						if( angle1 > 0.1 && angle2 > 0.1 ) {
							// region is convex
							// compute intersection of neighboring edges
							auto A = pVec[0];
							auto B = pVec[1];
							auto C = pVec[2];
							auto D = pVec[3];
							// Line AB represented as a1x + b1y = c1
							auto a = B[1] - A[1];
							auto b = A[0] - B[0];
							auto c = a * A[0] + b * A[1];
							// Line CD represented as a2x + b2y = c2
							auto a1 = D[1] - C[1];
							auto b1 = C[0] - D[0];
							auto c1 = a1 * C[0] + b1 * C[1];
							auto det = a * b1 - a1 * b;
							if( det == 0 ) {
								p = { FLT_MAX,FLT_MAX,FLT_MAX };
							} else {
								auto x = ( b1 * c - b * c1 ) / det;
								auto y = ( a * c1 - a1 * c ) / det;
								p = { x,y,0 };
							}

						} else if( angle1 <= 0.1 && angle2 <= 0.1 ) {
							// region is concave --> use midpoint
							p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );
						} else if( angle1 <= 0.1 && angle2 > 0.1 ) {
							// collapse to convex point
							p = pVec[2];
						} else if( angle1 > 0.1 && angle2 <= 0.1 ) {
							// collapse to convex point
							p = pVec[1];
						} else {
							LOG( ERROR ) << "Impossible case in boundary computation of convex hull";
							LOG( ERROR ) << "angle1 = " << angle1 << "  |  angle2 = " << angle2;
							LOG( ERROR ) << pVec[0] << "  |  " << pVec[1] << "  |  " << pVec[2];
						}
					}
					

					//p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );
				} else if( b0 && !b1 ) {
					p = mesh.point( v0 );
				} else if( !b0 && b1 ) {
					p = mesh.point( v1 );
				} else {
					p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );
				}

				return p;
			}

			inline float minMeanRatioMetric( MeshT& mesh, const TriMesh::VertexHandle& vh1, const TriMesh::VertexHandle& vh2, const TriMesh::Point& p_new ) {
				float minQ = FLT_MAX;

				std::vector<TriMesh::VertexHandle> vh1fan;
				MeshFunctions::getOneRing( mesh, vh1, vh1fan );

				for( size_t i = 0; i < vh1fan.size() - 1; ++i ) {
					size_t j = i + 1;
					if( vh1fan[i] == vh2 || vh1fan[j] == vh2 ) {
						continue;
					}
					TriMesh::Point p1 = mesh.point( vh1fan[i] );
					TriMesh::Point p2 = mesh.point( vh1fan[j] );

					const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

					minQ = std::min( minQ, quality );

				}

				std::vector<TriMesh::VertexHandle> vh2fan;
				MeshFunctions::getOneRing( mesh, vh2, vh2fan );

				for( size_t i = 0; i < vh2fan.size() - 1; ++i ) {
					size_t j = i + 1;
					if( vh2fan[i] == vh1 || vh2fan[j] == vh1 ) {
						continue;
					}
					TriMesh::Point p1 = mesh.point( vh2fan[i] );
					TriMesh::Point p2 = mesh.point( vh2fan[j] );

					const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

					minQ = std::min( minQ, quality );
				}

				return minQ;
			}
		};


		//=============================================================================

		//== CLASS DEFINITION =========================================================

		template < typename MeshT >
		class OceanDecimaterT : virtual public BaseDecimaterT<MeshT> //virtual especially for the mixed decimater
		{
		public: //-------------------------------------------------------- public types

			typedef OceanDecimaterT< MeshT >           Self;
			typedef MeshT                         Mesh;
			typedef CollapseInfoT<MeshT>          CollapseInfo;
			typedef ModBaseT<MeshT>               Module;
			typedef std::vector< Module* >        ModuleList;
			typedef typename ModuleList::iterator ModuleListIterator;

		public: //------------------------------------------------------ public methods

				/// Constructor
			OceanDecimaterT( Mesh& _mesh ) :
				BaseDecimaterT<Mesh>( _mesh ),
				mesh_( _mesh ),
				#if (defined(_MSC_VER) && (_MSC_VER >= 1800)) || __cplusplus > 199711L || defined( __GXX_EXPERIMENTAL_CXX0X__ )
				heap_( nullptr )
				#else
				heap_( NULL )
				#endif

			{

				// private vertex properties
				mesh_.add_property( collapse_target_ );
				mesh_.add_property( priority_ );
				mesh_.add_property( heap_position_ );
			}

			/// Destructor
			~OceanDecimaterT() {

				// private vertex properties
				mesh_.remove_property( collapse_target_ );
				mesh_.remove_property( priority_ );
				mesh_.remove_property( heap_position_ );

			}

		public:

			/**
			* @brief Perform a number of collapses on the mesh.
			* @param _n_collapses Desired number of collapses. If zero (default), attempt
			*                     to do as many collapses as possible.
			* @return Number of collapses that were actually performed.
			* @note This operation only marks the removed mesh elements for deletion. In
			*       order to actually remove the decimated elements from the mesh, a
			*       subsequent call to ArrayKernel::garbage_collection() is required.
			*/
			size_t decimate( size_t _n_collapses = 0 ) {

				if ( !this->is_initialized() )
					return 0;

				typename Mesh::VertexIter v_it, v_end( mesh_.vertices_end() );
				typename Mesh::VertexHandle vp;
				typename Mesh::HalfedgeHandle v0v1;
				typename Mesh::VertexVertexIter vv_it;
				typename Mesh::VertexFaceIter vf_it;
				unsigned int n_collapses( 0 );

				typedef std::set<typename Mesh::VertexHandle> Support;
				typedef typename Support::iterator SupportIterator;

				Support support;
				SupportIterator s_it, s_end;

				// check _n_collapses
				if ( !_n_collapses )
					_n_collapses = mesh_.n_vertices();

				// initialize heap
				HeapInterface HI( mesh_, priority_, heap_position_ );

				#if (defined(_MSC_VER) && (_MSC_VER >= 1800)) || __cplusplus > 199711L || defined( __GXX_EXPERIMENTAL_CXX0X__ )
				heap_ = std::unique_ptr<DeciHeap>( new DeciHeap( HI ) );
				#else
				heap_ = std::auto_ptr<DeciHeap>( new DeciHeap( HI ) );
				#endif


				heap_->reserve( mesh_.n_vertices() );

				for ( v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it ) {
					heap_->reset_heap_position( *v_it );
					if ( !mesh_.status( *v_it ).deleted() )
						heap_vertex( *v_it );
				}

				const bool update_normals = mesh_.has_face_normals();

				// process heap
				while ( ( !heap_->empty() ) && ( n_collapses < _n_collapses ) ) {
					// get 1st heap entry
					vp = heap_->front();
					v0v1 = mesh_.property( collapse_target_, vp );
					heap_->pop_front();

					// setup collapse info
					CollapseInfo ci( mesh_, v0v1 );

					// check topological correctness AGAIN !
					if ( !this->is_collapse_legal( ci ) )
						continue;

					// store support
					// the two-ring neighborhood is required for updating the element quality correctly
					support.clear();
					for( auto vv : mesh_.vv_range( ci.v0 ) ) {
						support.insert( vv );
						for( auto vvv : mesh_.vv_range( vv ) ) {
							support.insert( vvv );
						}
					}
					for( auto vv : mesh_.vv_range( ci.v1 ) ) {
						support.insert( vv );
						for( auto vvv : mesh_.vv_range( vv ) ) {
							support.insert( vvv );
						}
					}
					support.erase( ci.v0 );

					// pre-processing
					this->preprocess_collapse( ci );

					// perform collapse
					mesh_.collapse( v0v1 );
					++n_collapses;

					if ( update_normals ) {
						// update triangle normals
						vf_it = mesh_.vf_iter( ci.v1 );
						for ( ; vf_it.is_valid(); ++vf_it )
							if ( !mesh_.status( *vf_it ).deleted() )
								mesh_.set_normal( *vf_it, mesh_.calc_face_normal( *vf_it ) );
					}

					// post-process collapse
					this->postprocess_collapse( ci );

					// update heap (former one ring of decimated vertex)
					for ( s_it = support.begin(), s_end = support.end(); s_it != s_end; ++s_it ) {
						assert( !mesh_.status( *s_it ).deleted() );
						heap_vertex( *s_it );
					}

					// notify observer and stop if the observer requests it
					if ( !this->notify_observer( n_collapses ) )
						return n_collapses;
				}

				// delete heap
				heap_.reset();



				// DON'T do garbage collection here! It's up to the application.
				return n_collapses;
			}

			/**
			* @brief Decimate the mesh to a desired target vertex complexity.
			* @param _n_vertices Target complexity, i.e. desired number of remaining
			*                    vertices after decimation.
			* @return Number of collapses that were actually performed.
			* @note This operation only marks the removed mesh elements for deletion. In
			*       order to actually remove the decimated elements from the mesh, a
			*       subsequent call to ArrayKernel::garbage_collection() is required.
			*/
			size_t decimate_to( size_t  _n_vertices ) {
				return ( ( _n_vertices < this->mesh().n_vertices() ) ?
						 decimate( this->mesh().n_vertices() - _n_vertices ) : 0 );
			}

			/**
			* @brief Attempts to decimate the mesh until a desired vertex or face
			*        complexity is achieved.
			* @param _n_vertices Target vertex complexity.
			* @param _n_faces Target face complexity.
			* @return Number of collapses that were actually performed.
			* @note Decimation stops as soon as either one of the two complexity bounds
			*       is satisfied.
			* @note This operation only marks the removed mesh elements for deletion. In
			*       order to actually remove the decimated elements from the mesh, a
			*       subsequent call to ArrayKernel::garbage_collection() is required.
			*/
			size_t decimate_to_faces( size_t  _nv = 0, size_t _nf = 0 ) {

				if ( !this->is_initialized() )
					return 0;

				if ( _nv >= mesh_.n_vertices() || _nf >= mesh_.n_faces() )
					return 0;

				typename Mesh::VertexIter v_it, v_end( mesh_.vertices_end() );
				typename Mesh::VertexHandle vp;
				typename Mesh::HalfedgeHandle v0v1;
				typename Mesh::VertexVertexIter vv_it;
				typename Mesh::VertexFaceIter vf_it;
				size_t nv = mesh_.n_vertices();
				size_t nf = mesh_.n_faces();
				unsigned int n_collapses = 0;

				typedef std::set<typename Mesh::VertexHandle> Support;
				typedef typename Support::iterator SupportIterator;

				Support support;

				// initialize heap
				HeapInterface HI( mesh_, priority_, heap_position_ );
				#if (defined(_MSC_VER) && (_MSC_VER >= 1800)) || __cplusplus > 199711L || defined( __GXX_EXPERIMENTAL_CXX0X__ )
				heap_ = std::unique_ptr<DeciHeap>( new DeciHeap( HI ) );
				#else
				heap_ = std::auto_ptr<DeciHeap>( new DeciHeap( HI ) );
				#endif
				heap_->reserve( mesh_.n_vertices() );

				for ( v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it ) {
					heap_->reset_heap_position( *v_it );
					if ( !mesh_.status( *v_it ).deleted() )
						heap_vertex( *v_it );
				}

				const bool update_normals = mesh_.has_face_normals();

				// process heap
				while ( ( !heap_->empty() ) && ( _nv < nv ) && ( _nf < nf ) ) {
					// get 1st heap entry
					vp = heap_->front();
					v0v1 = mesh_.property( collapse_target_, vp );
					heap_->pop_front();

					// setup collapse info
					CollapseInfo ci( mesh_, v0v1 );

					// check topological correctness AGAIN !
					if ( !this->is_collapse_legal( ci ) )
						continue;

					// store support
					// the two-ring neighborhood is required for updating the element quality correctly
					support.clear();
					for ( auto vv : mesh_.vv_range( ci.v0 ) ) {
						support.insert( vv );
						for ( auto vvv : mesh_.vv_range( vv ) ) {
							support.insert( vvv );
						}
					}
					for ( auto vv : mesh_.vv_range( ci.v1 ) ) {
						support.insert( vv );
							for ( auto vvv : mesh_.vv_range( vv ) ) {
								support.insert( vvv );
							}
					}
					support.erase( ci.v0 );

					// adjust complexity in advance (need boundary status)
					++n_collapses;
					--nv;
					if ( mesh_.is_boundary( ci.v0v1 ) || mesh_.is_boundary( ci.v1v0 ) )
						--nf;
					else
						nf -= 2;

					// pre-processing
					this->preprocess_collapse( ci );

					// perform collapse
					mesh_.collapse( v0v1 );

					// update triangle normals
					if ( update_normals ) {
						vf_it = mesh_.vf_iter( ci.v1 );
						for ( ; vf_it.is_valid(); ++vf_it )
							if ( !mesh_.status( *vf_it ).deleted() )
								mesh_.set_normal( *vf_it, mesh_.calc_face_normal( *vf_it ) );
					}

					// post-process collapse
					this->postprocess_collapse( ci );

					// update heap (former one ring of decimated vertex)
					for ( const auto& s : support ) {
						assert( !mesh_.status( s ).deleted() );
						heap_vertex( s );
					}

					// notify observer and stop if the observer requests it
					if ( !this->notify_observer( n_collapses ) )
						return n_collapses;
				}

				// delete heap
				heap_.reset();


				// DON'T do garbage collection here! It's up to the application.
				return n_collapses;
			}

		public:

			typedef typename Mesh::VertexHandle    VertexHandle;
			typedef typename Mesh::HalfedgeHandle  HalfedgeHandle;

			/// Heap interface
			class HeapInterface
			{
			public:

				HeapInterface( Mesh&               _mesh,
								VPropHandleT<float> _prio,
								VPropHandleT<int>   _pos )
					: mesh_( _mesh ), prio_( _prio ), pos_( _pos ) {
				}

				inline bool
					less( VertexHandle _vh0, VertexHandle _vh1 ) {
					return mesh_.property( prio_, _vh0 ) < mesh_.property( prio_, _vh1 );
				}

				inline bool
					greater( VertexHandle _vh0, VertexHandle _vh1 ) {
					return mesh_.property( prio_, _vh0 ) > mesh_.property( prio_, _vh1 );
				}

				inline int
					get_heap_position( VertexHandle _vh ) {
					return mesh_.property( pos_, _vh );
				}

				inline void
					set_heap_position( VertexHandle _vh, int _pos ) {
					mesh_.property( pos_, _vh ) = _pos;
				}


			private:
				Mesh&                mesh_;
				VPropHandleT<float>  prio_;
				VPropHandleT<int>    pos_;
			};

			typedef Utils::HeapT<VertexHandle, HeapInterface>  DeciHeap;


		private: //---------------------------------------------------- private methods

					/// Insert vertex in heap
			void heap_vertex( VertexHandle _vh ) {
				//   std::clog << "heap_vertex: " << _vh << std::endl;

				float prio, best_prio( FLT_MAX );
				typename Mesh::HalfedgeHandle heh, collapse_target;

				// find best target in one ring
				typename Mesh::VertexOHalfedgeIter voh_it( mesh_, _vh );
				for ( ; voh_it.is_valid(); ++voh_it ) {
					heh = *voh_it;
					CollapseInfo ci( mesh_, heh );

					if ( this->is_collapse_legal( ci ) ) {
						prio = this->collapse_priority( ci );
						if ( prio >= 0.0 && prio < best_prio ) {
							best_prio = prio;
							collapse_target = heh;
						}
					}
				}

				// target found -> put vertex on heap
				if ( collapse_target.is_valid() ) {
					//     std::clog << "  added|updated" << std::endl;
					mesh_.property( collapse_target_, _vh ) = collapse_target;
					mesh_.property( priority_, _vh ) = best_prio;

					if ( heap_->is_stored( _vh ) )
						heap_->update( _vh );
					else
						heap_->insert( _vh );
				}

				// not valid -> remove from heap
				else {
					//     std::clog << "  n/a|removed" << std::endl;
					if ( heap_->is_stored( _vh ) )
						heap_->remove( _vh );

					mesh_.property( collapse_target_, _vh ) = collapse_target;
					mesh_.property( priority_, _vh ) = -1;
				}
			}

		private: //------------------------------------------------------- private data


					// reference to mesh
			Mesh&      mesh_;

			// heap
			#if (defined(_MSC_VER) && (_MSC_VER >= 1800)) || __cplusplus > 199711L || defined( __GXX_EXPERIMENTAL_CXX0X__ )
			std::unique_ptr<DeciHeap> heap_;
			#else
			std::auto_ptr<DeciHeap> heap_;
			#endif

			// vertex properties
			VPropHandleT<HalfedgeHandle>  collapse_target_;
			VPropHandleT<float>           priority_;
			VPropHandleT<int>             heap_position_;

		};


		//=============================================================================

		class MeshPointsNN
		{
			using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, MeshPointsNN >, MeshPointsNN, 3>;
			std::unique_ptr<kd_tree_t> kdTree_;

			std::vector<TriMesh::VertexHandle> vhs_;
			TriMesh* mesh_;

		public:
			MeshPointsNN() {}
			MeshPointsNN( TriMesh* mesh ) : mesh_{ mesh } {
				for( auto vh : mesh_->vertices() ) {
					if( !vh.is_boundary() )
						continue;
					vhs_.push_back( vh );
				}

				// init kd tree
				kdTree_ = std::make_unique<kd_tree_t>( 3, *this );
				kdTree_->buildIndex();
			}

			MeshPointsNN& operator = ( const MeshPointsNN& t ) {
				mesh_ = t.mesh_;
				vhs_ = t.vhs_;

				// init kd tree
				kdTree_ = std::make_unique<kd_tree_t>( 3, *this );
				kdTree_->buildIndex();

				return *this;
			}

			void buildIndex() {
				vhs_.erase( std::remove_if(vhs_.begin(), vhs_.end(),
					[this]( const TriMesh::VertexHandle& x ) {
						return mesh_->status(x).deleted(); // put your condition here
					} ), vhs_.end() 
						);

				kdTree_->buildIndex();
			}

			TriMesh::VertexHandle findNearestNeighbor( const TriMesh::Point& p ) const {
				// the nearest neighbor should be the point itself --> search for the second nearest neighbor
				constexpr int k = 2;
				std::vector<size_t> ret_indexes( k );
				std::vector<float> out_dists_sqr( k );
				nanoflann::KNNResultSet<float> resultSet( k );
				resultSet.init( &ret_indexes[0], &out_dists_sqr[0] );
				kdTree_->findNeighbors( resultSet, &p[0], nanoflann::SearchParams() );

				return vhs_[ret_indexes[1]];
			}

			/*----- nanoflann functions -----*/
			// Must return the number of data points
			inline size_t kdtree_get_point_count() const { return vhs_.size(); }
			// Returns the dim'th component of the idx'th point in the class
			inline float kdtree_get_pt( const size_t idx, const size_t dim ) const {
				return mesh_->point( vhs_[idx] )[dim];
			}
			template <class BBOX>
			bool kdtree_get_bbox( BBOX& ) const { return false; }
		};

		template < typename MeshT >
		class MyCollapseInfoT : virtual public CollapseInfoT<MeshT>
		{
		public:
			using CollapseInfoT<MeshT>::CollapseInfoT;
			template < typename Mesh > 
			MyCollapseInfoT( Mesh& _mesh, typename Mesh::VertexHandle _v0, typename Mesh::VertexHandle _v1 )
				: CollapseInfoT<MeshT>( _mesh, _mesh.halfedge_handle(0) ) {
				CollapseInfoT<MeshT>::v0 = _v0; 
				CollapseInfoT<MeshT>::v1 = _v1; 
				CollapseInfoT<MeshT>::p0 = _mesh.point( CollapseInfoT<MeshT>::v0 );
				CollapseInfoT<MeshT>::p1 = _mesh.point( CollapseInfoT<MeshT>::v1 );
				CollapseInfoT<MeshT>::v0v1 = typename Mesh::HalfedgeHandle();
				CollapseInfoT<MeshT>::v1v0 = typename Mesh::HalfedgeHandle();

				CollapseInfoT<MeshT>::vlv1 = typename Mesh::HalfedgeHandle();
				CollapseInfoT<MeshT>::v0vl = typename Mesh::HalfedgeHandle();
				CollapseInfoT<MeshT>::vrv0 = typename Mesh::HalfedgeHandle();
				CollapseInfoT<MeshT>::v1vr = typename Mesh::HalfedgeHandle();
				CollapseInfoT<MeshT>::fl = typename Mesh::FaceHandle();
				CollapseInfoT<MeshT>::fr = typename Mesh::FaceHandle();
				CollapseInfoT<MeshT>::vl = typename Mesh::VertexHandle();
				CollapseInfoT<MeshT>::vr = typename Mesh::VertexHandle();
			}
		};

		template <class MeshT>
		class ModConvHullPointDist2T : public ModBaseT<MeshT>
		{
		public:
			DECIMATING_MODULE( ModConvHullPointDist2T, MeshT, PointDistanceModule );

			using value_type = vector_traits<TriMesh::Point>::value_type;
		public:

			/// Constructor
			ModConvHullPointDist2T( MeshT& _dec ) :
				Base( _dec, false ) {
				Base::mesh().add_property( vErrors_ );
				Base::mesh().add_property( ePoints_ );
			}

			/// Destructor
			~ModConvHullPointDist2T() {
				Base::mesh().remove_property( vErrors_ );
				Base::mesh().remove_property( ePoints_ );
			}

		public: // inherited

			/// Initalize the module and prepare the mesh for decimation.
			virtual void initialize( void ) {
				auto& mesh = Base::mesh();

				if( !vErrors_.is_valid() ) {
					mesh.add_property( vErrors_ );
					mesh.add_property( ePoints_ );
				}

				// init errors
				if( hasBackgroundGrid ) {
					for( auto vh : mesh.vertices() ) {
						mesh.property( vErrors_, vh ) = PointDistErrorMetric( mesh.point( vh ), bg_ );
					}
				} else {
					for( auto vh : mesh.vertices() ) {
						mesh.property( vErrors_, vh ) = PointDistErrorMetric( mesh.point( vh ) );
					}
				}
				// init edge points
				for( auto eh : mesh.edges() ) {
					auto heh = mesh.halfedge_handle( eh, 0 );
					mesh.property( ePoints_, eh ) = mesh.point( mesh.to_vertex_handle( heh ) );
				}

				Base::set_binary( false );
			}

			float collapse_priority( const CollapseInfo& _ci ) {

				auto& mesh = Base::mesh();

				OpenMesh::VPropHandleT<float> feature_size;
				mesh.get_property_handle( feature_size, "feature_size" );

				float priority = Base::ILLEGAL_COLLAPSE;

				if( true ) // continuous mode
				{
					const auto& v0 = _ci.v0;
					const auto& v1 = _ci.v1;

					if( !_ci.v0v1.is_valid() ) {
						auto p = 0.5 * ( _ci.p0 + _ci.p1 );
						// assure minimal quality
						if( minMeanRatioMetric( mesh, v0, v1, p ) < min_r_ ) {
							priority = Base::ILLEGAL_COLLAPSE;
							return priority;
						}
						auto e1 = mesh.property( vErrors_, _ci.v0 );
						auto e2 = mesh.property( vErrors_, _ci.v1 );
						auto e = e1 + e2;
						priority = e * p;
						return priority;
					}

					const auto& eh = mesh.edge_handle( _ci.v0v1 );

					auto& p = mesh.property( ePoints_, eh );

					// jump over edges that connect two boundary vertices but are themselves interior (don't reduce canals)
					if( mesh.is_boundary( v0 ) && mesh.is_boundary( v1 ) && !mesh.is_boundary( eh ) ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}

					p = computeEdgePoint( _ci );

					// assure minimal quality
					if( minMeanRatioMetric( mesh, v0, v1, p ) < min_r_ ) {
						priority = Base::ILLEGAL_COLLAPSE;
						return priority;
					}

					auto e1 = mesh.property( vErrors_, _ci.v0 );
					auto e2 = mesh.property( vErrors_, _ci.v1 );
					auto e = e1 + e2;
					priority = e * p;

				} else // binary mode
				{
					LOG( ERROR ) << "Binary mode for point distance not implemented";
					priority = Base::ILLEGAL_COLLAPSE;
				}

				return priority;
			}

			void postprocess_collapse( const CollapseInfoT< MeshT >& _ci ) {
				auto& mesh = Base::mesh();
				mesh.property( vErrors_, _ci.v1 ) += mesh.property( vErrors_, _ci.v0 );

				typename MeshT::Point pNew;
				if( _ci.v0v1.is_valid() )
					pNew = mesh.property( ePoints_, mesh.edge_handle( _ci.v0v1 ) );
				else
					pNew = 0.5f * ( _ci.p0 + _ci.p1 );
				mesh.set_point( _ci.v1, pNew );

				OpenMesh::VPropHandleT<float> feature_size;
				mesh.get_property_handle( feature_size, "feature_size" );

				mesh.property( feature_size, _ci.v1 ) = std::min( mesh.property( feature_size, _ci.v0 ), mesh.property( feature_size, _ci.v1 ) );
				mesh.data( _ci.v1 ).is_feature |= mesh.data( _ci.v0 ).is_feature;

				//// transfer edge properties (seems to work also without this)
				//OpenMesh::EPropHandleT<bool> isNeumann;
				//mesh.get_property_handle( isNeumann, "isNeumann" );
				//OpenMesh::EPropHandleT<int> boundaryID;
				//mesh.get_property_handle( boundaryID, "boundaryID" );
				//auto heh = _ci.v0v1;
				//if ( mesh.is_boundary( _ci.v1 ) && !mesh.is_boundary( _ci.v0 ) ) {
				//	auto prev = mesh.prev_halfedge_handle( heh );
				//	auto next = mesh.next_halfedge_handle( heh );
				//	mesh.property( isNeumann, mesh.edge_handle( next ) ) = mesh.property( isNeumann, mesh.edge_handle( prev ) );
				//	heh = _ci.v1v0;
				//	prev = mesh.prev_halfedge_handle( heh );
				//	next = mesh.next_halfedge_handle( heh );
				//	mesh.property( isNeumann, mesh.edge_handle( prev ) ) = mesh.property( isNeumann, mesh.edge_handle( next ) );
				//}

			}

			void set_backgroundgrid( const ScalarField::ScalarField* bg ) {
				bg_ = bg;
				hasBackgroundGrid = true;
			}

			void set_min_meanratio( value_type _min_meanratio, bool _binary = true ) {
				assert( 0.0 <= _min_meanratio && _min_meanratio <= 1.0 );
				min_r_ = _min_meanratio;
			}

		private:
			OpenMesh::VPropHandleT<PointDistErrorMetric> vErrors_;
			OpenMesh::EPropHandleT<OpenMesh::Vec3f> ePoints_;	// store target positions
			const ScalarField::ScalarField* bg_;
			bool hasBackgroundGrid = false;

			value_type min_r_ = -FLT_MAX;

			inline OpenMesh::Vec3f computeEdgePoint( const CollapseInfo& _ci ) {

				auto& mesh = Base::mesh();

				const auto& v0 = _ci.v0;
				const auto& v1 = _ci.v1;

				OpenMesh::Vec3f p( 0, 0, 0 );

				// calculate new point
				if( mesh.is_boundary( v0 ) && mesh.is_boundary( v1 ) ) {

					// compute boundary angles
					auto heh = _ci.v0v1;
					if( !mesh.is_boundary( heh ) )
						heh = _ci.v1v0;

					std::vector<TriMesh::VertexHandle> vhs( 4 );
					vhs[0] = mesh.from_vertex_handle( mesh.prev_halfedge_handle( heh ) );
					vhs[1] = mesh.from_vertex_handle( heh );
					vhs[2] = mesh.to_vertex_handle( heh );
					vhs[3] = mesh.to_vertex_handle( mesh.next_halfedge_handle( heh ) );

					std::vector<TriMesh::Point> pVec( 4 );
					std::transform( vhs.begin(), vhs.end(), pVec.begin(), [&mesh]( const auto& v ) { return mesh.point( v ); } );

					auto angle1 = 180.f - HelperFunctions::calcAngle2( pVec[0], pVec[1], pVec[2] );
					auto angle2 = 180.f - HelperFunctions::calcAngle2( pVec[1], pVec[2], pVec[3] );

					if( angle1 + angle2 > 170.f ) {
						p = { FLT_MAX,FLT_MAX,FLT_MAX };
						return p;
					}

					//p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );

					if( angle1 > 3 && angle2 > 3 ) {
						// region is convex
						// compute intersection of neighboring edges
						auto A = pVec[0];
						auto B = pVec[1];
						auto C = pVec[2];
						auto D = pVec[3];
						// Line AB represented as a1x + b1y = c1
						auto a = B[1] - A[1];
						auto b = A[0] - B[0];
						auto c = a * A[0] + b * A[1];
						// Line CD represented as a2x + b2y = c2
						auto a1 = D[1] - C[1];
						auto b1 = C[0] - D[0];
						auto c1 = a1 * C[0] + b1 * C[1];
						auto det = a * b1 - a1 * b;
						if( det == 0 ) {
							p = { FLT_MAX,FLT_MAX,FLT_MAX };
						} else {
							auto x = ( b1 * c - b * c1 ) / det;
							auto y = ( a * c1 - a1 * c ) / det;
							p = { x,y,0 };
						}

					} else if( angle1 <= 3 && angle2 <= 3 ) {
						// region is concave --> use midpoint
						p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );
					} else if( angle1 <= 3 && angle2 > 3 ) {
						// collapse to convex point
						p = pVec[2];
					} else if( angle1 > 3 && angle2 <= 3 ) {
						// collapse to convex point
						p = pVec[1];
					} else {
						LOG( ERROR ) << "Impossible case in boundary computation of convex hull";
					}

					//p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );
				} else if( mesh.is_boundary( v0 ) && !mesh.is_boundary( v1 ) ) {
					p = mesh.point( v0 );
				} else if( !mesh.is_boundary( v0 ) && mesh.is_boundary( v1 ) ) {
					p = mesh.point( v1 );
				} else {
					p = 0.5f * ( mesh.point( v0 ) + mesh.point( v1 ) );
				}

				return p;
			}

			inline float minMeanRatioMetric( MeshT& mesh, const TriMesh::VertexHandle& vh1, const TriMesh::VertexHandle& vh2, const TriMesh::Point& p_new ) {
				float minQ = FLT_MAX;

				std::vector<TriMesh::VertexHandle> vh1fan;
				MeshFunctions::getOneRing( mesh, vh1, vh1fan );

				for( size_t i = 0; i < vh1fan.size() - 1; ++i ) {
					size_t j = i + 1;
					if( vh1fan[i] == vh2 || vh1fan[j] == vh2 ) {
						continue;
					}
					TriMesh::Point p1 = mesh.point( vh1fan[i] );
					TriMesh::Point p2 = mesh.point( vh1fan[j] );

					const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

					minQ = std::min( minQ, quality );

				}

				std::vector<TriMesh::VertexHandle> vh2fan;
				MeshFunctions::getOneRing( mesh, vh2, vh2fan );

				for( size_t i = 0; i < vh2fan.size() - 1; ++i ) {
					size_t j = i + 1;
					if( vh2fan[i] == vh1 || vh2fan[j] == vh1 ) {
						continue;
					}
					TriMesh::Point p1 = mesh.point( vh2fan[i] );
					TriMesh::Point p2 = mesh.point( vh2fan[j] );

					const float quality = QualityMetrics::meanRatioMetric( { p_new, p2, p1 } );

					minQ = std::min( minQ, quality );
				}

				return minQ;
			}
		};

		class OceanDecimater2T : virtual public BaseDecimaterT<TriMesh> //virtual especially for the mixed decimater
		{
		public: //-------------------------------------------------------- public types

			typedef OceanDecimater2T     Self;
			typedef TriMesh                         Mesh;
			typedef MyCollapseInfoT<TriMesh>        CollapseInfo;
			typedef ModBaseT<TriMesh>               Module;
			typedef std::vector< Module* >        ModuleList;
			typedef typename ModuleList::iterator ModuleListIterator;
			MeshPointsNN meshPointsNN;

		public: //------------------------------------------------------ public methods

				/// Constructor
			OceanDecimater2T( Mesh& _mesh ) :
				BaseDecimaterT<Mesh>( _mesh ),
				mesh_( _mesh ),
			#if (defined(_MSC_VER) && (_MSC_VER >= 1800)) || __cplusplus > 199711L || defined( __GXX_EXPERIMENTAL_CXX0X__ )
				heap_( nullptr )
			#else
				heap_( NULL )
			#endif
				
			{

				// private vertex properties
				mesh_.add_property( collapse_target_ );
				mesh_.add_property( collapse_target_vh_ );
				mesh_.add_property( priority_ );
				mesh_.add_property( heap_position_ );

				meshPointsNN = MeshPointsNN( &mesh_ );
			}

			/// Destructor
			~OceanDecimater2T() {

				// private vertex properties
				mesh_.remove_property( collapse_target_ );
				mesh_.remove_property( collapse_target_vh_ );
				mesh_.remove_property( priority_ );
				mesh_.remove_property( heap_position_ );

			}

		public:

			/**
			* @brief Perform a number of collapses on the mesh.
			* @param _n_collapses Desired number of collapses. If zero (default), attempt
			*                     to do as many collapses as possible.
			* @return Number of collapses that were actually performed.
			* @note This operation only marks the removed mesh elements for deletion. In
			*       order to actually remove the decimated elements from the mesh, a
			*       subsequent call to ArrayKernel::garbage_collection() is required.
			*/
			size_t decimate( size_t _n_collapses = 0 ) {

				if( !this->is_initialized() )
					return 0;

				typename Mesh::VertexIter v_it, v_end( mesh_.vertices_end() );
				typename Mesh::VertexHandle vp;
				typename Mesh::VertexHandle vh_target;
				typename Mesh::HalfedgeHandle v0v1;
				typename Mesh::VertexVertexIter vv_it;
				typename Mesh::VertexFaceIter vf_it;
				unsigned int n_collapses( 0 );

				typedef std::set<typename Mesh::VertexHandle> Support;
				typedef typename Support::iterator SupportIterator;

				Support support;
				SupportIterator s_it, s_end;

				// check _n_collapses
				if( !_n_collapses )
					_n_collapses = mesh_.n_vertices();

				// initialize heap
				HeapInterface HI( mesh_, priority_, heap_position_ );

			#if (defined(_MSC_VER) && (_MSC_VER >= 1800)) || __cplusplus > 199711L || defined( __GXX_EXPERIMENTAL_CXX0X__ )
				heap_ = std::unique_ptr<DeciHeap>( new DeciHeap( HI ) );
			#else
				heap_ = std::auto_ptr<DeciHeap>( new DeciHeap( HI ) );
			#endif
				heap_->reserve( mesh_.n_vertices() );

				for( v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it ) {
					heap_->reset_heap_position( *v_it );
					if( !mesh_.status( *v_it ).deleted() )
						heap_vertex( *v_it );
				}

				const bool update_normals = mesh_.has_face_normals();

				// process heap
				while( ( !heap_->empty() ) && ( n_collapses < _n_collapses ) ) {
					// get 1st heap entry
					vp = heap_->front();
					v0v1 = mesh_.property( collapse_target_, vp );
					vh_target = mesh_.property( collapse_target_vh_, vp );
					heap_->pop_front();

					if( !v0v1.is_valid() ) {
						if( mesh_.status( vp ).deleted() || mesh_.status( vh_target ).deleted() ) {
							continue;
						}
						// setup collapse info
						CollapseInfo ci( mesh_, vp, vh_target );

						// store support
						// the two-ring neighborhood is required for updating the element quality correctly
						support.clear();
						for( auto vv : mesh_.vv_range( ci.v0 ) ) {
							support.insert( vv );
							for( auto vvv : mesh_.vv_range( vv ) ) {
								support.insert( vvv );
							}
						}
						for( auto vv : mesh_.vv_range( ci.v1 ) ) {
							support.insert( vv );
							for( auto vvv : mesh_.vv_range( vv ) ) {
								support.insert( vvv );
							}
						}
						support.erase( ci.v0 );

						// pre-processing
						this->preprocess_collapse( ci );

						// perform collapse
						//mesh_.collapse( v0v1 );
						//++n_collapses;
						////////////////////////////////////////// !!!!!!!!!!!!!!!!!!!!!!
						std::vector<std::vector<Mesh::VertexHandle>> faces;
						for( auto vf : mesh_.vf_range( ci.v0 ) ) {
							std::vector<Mesh::VertexHandle> f;
							f.reserve( 3 );
							for( auto fv : vf.vertices() ) {
								f.push_back( fv );
							}
							faces.push_back( f );
						}
						for( auto& f : faces ) {
							std::replace( f.begin(), f.end(), ci.v0, ci.v1 );
						}
						mesh_.delete_vertex( ci.v0, false );
						for( const auto& f : faces ) {
							mesh_.add_face( f );
						}
						++n_collapses;

						if( update_normals ) {
							// update triangle normals
							vf_it = mesh_.vf_iter( ci.v1 );
							for( ; vf_it.is_valid(); ++vf_it )
								if( !mesh_.status( *vf_it ).deleted() )
									mesh_.set_normal( *vf_it, mesh_.calc_face_normal( *vf_it ) );
						}

						// post-process collapse
						this->postprocess_collapse( ci );
						meshPointsNN.buildIndex();

						// update heap (former one ring of decimated vertex)
						for( s_it = support.begin(), s_end = support.end(); s_it != s_end; ++s_it ) {
							assert( !mesh_.status( *s_it ).deleted() );
							heap_vertex( *s_it );
						}

						// notify observer and stop if the observer requests it
						if( !this->notify_observer( n_collapses ) )
							return n_collapses;

						continue;
					}

					// setup collapse info
					CollapseInfo ci( mesh_, v0v1 );

					// check topological correctness AGAIN !
					if( !this->is_collapse_legal( ci ) )
						continue;

					// store support
					// the two-ring neighborhood is required for updating the element quality correctly
					support.clear();
					for( auto vv : mesh_.vv_range( ci.v0 ) ) {
						support.insert( vv );
						for( auto vvv : mesh_.vv_range( vv ) ) {
							support.insert( vvv );
						}
					}
					for( auto vv : mesh_.vv_range( ci.v1 ) ) {
						support.insert( vv );
						for( auto vvv : mesh_.vv_range( vv ) ) {
							support.insert( vvv );
						}
					}
					support.erase( ci.v0 );

					// pre-processing
					this->preprocess_collapse( ci );

					// perform collapse
					mesh_.collapse( v0v1 );
					++n_collapses;

					if( update_normals ) {
						// update triangle normals
						vf_it = mesh_.vf_iter( ci.v1 );
						for( ; vf_it.is_valid(); ++vf_it )
							if( !mesh_.status( *vf_it ).deleted() )
								mesh_.set_normal( *vf_it, mesh_.calc_face_normal( *vf_it ) );
					}

					// post-process collapse
					this->postprocess_collapse( ci );
					meshPointsNN.buildIndex();

					// update heap (former one ring of decimated vertex)
					for( s_it = support.begin(), s_end = support.end(); s_it != s_end; ++s_it ) {
						//assert( !mesh_.status( *s_it ).deleted() );			// TODO change "if" back to assert
						if( !mesh_.status( *s_it ).deleted() )
							heap_vertex( *s_it );
					}

					// notify observer and stop if the observer requests it
					if( !this->notify_observer( n_collapses ) )
						return n_collapses;
				}

				// delete heap
				heap_.reset();

				// duplicate vertices that became non-manifold
				duplicate_non_manifold_vertices();


				// DON'T do garbage collection here! It's up to the application.
				return n_collapses;
			}

			/**
			* @brief Decimate the mesh to a desired target vertex complexity.
			* @param _n_vertices Target complexity, i.e. desired number of remaining
			*                    vertices after decimation.
			* @return Number of collapses that were actually performed.
			* @note This operation only marks the removed mesh elements for deletion. In
			*       order to actually remove the decimated elements from the mesh, a
			*       subsequent call to ArrayKernel::garbage_collection() is required.
			*/
			size_t decimate_to( size_t  _n_vertices ) {
				return ( ( _n_vertices < this->mesh().n_vertices() ) ?
						 decimate( this->mesh().n_vertices() - _n_vertices ) : 0 );
			}

			/**
			* @brief Attempts to decimate the mesh until a desired vertex or face
			*        complexity is achieved.
			* @param _n_vertices Target vertex complexity.
			* @param _n_faces Target face complexity.
			* @return Number of collapses that were actually performed.
			* @note Decimation stops as soon as either one of the two complexity bounds
			*       is satisfied.
			* @note This operation only marks the removed mesh elements for deletion. In
			*       order to actually remove the decimated elements from the mesh, a
			*       subsequent call to ArrayKernel::garbage_collection() is required.
			*/
			size_t decimate_to_faces( size_t  _nv = 0, size_t _nf = 0 ) {

				if( !this->is_initialized() )
					return 0;

				if( _nv >= mesh_.n_vertices() || _nf >= mesh_.n_faces() )
					return 0;

				typename Mesh::VertexIter v_it, v_end( mesh_.vertices_end() );
				typename Mesh::VertexHandle vp;
				typename Mesh::VertexHandle vh_target;
				typename Mesh::HalfedgeHandle v0v1;
				typename Mesh::VertexVertexIter vv_it;
				typename Mesh::VertexFaceIter vf_it;
				size_t nv = mesh_.n_vertices();
				size_t nf = mesh_.n_faces();
				unsigned int n_collapses = 0;

				typedef std::set<typename Mesh::VertexHandle> Support;
				typedef typename Support::iterator SupportIterator;

				Support support;
				SupportIterator s_it, s_end;

				// initialize heap
				HeapInterface HI( mesh_, priority_, heap_position_ );
			#if (defined(_MSC_VER) && (_MSC_VER >= 1800)) || __cplusplus > 199711L || defined( __GXX_EXPERIMENTAL_CXX0X__ )
				heap_ = std::unique_ptr<DeciHeap>( new DeciHeap( HI ) );
			#else
				heap_ = std::auto_ptr<DeciHeap>( new DeciHeap( HI ) );
			#endif
				heap_->reserve( mesh_.n_vertices() );

				for( v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it ) {
					heap_->reset_heap_position( *v_it );
					if( !mesh_.status( *v_it ).deleted() )
						heap_vertex( *v_it );
				}

				const bool update_normals = mesh_.has_face_normals();

				// process heap
				while( ( !heap_->empty() ) && ( _nv < nv ) && ( _nf < nf ) ) {
					// get 1st heap entry
					vp = heap_->front();
					v0v1 = mesh_.property( collapse_target_, vp );
					vh_target = mesh_.property( collapse_target_vh_, vp );
					heap_->pop_front();

					if( !v0v1.is_valid() ) {
						if( mesh_.status( vp ).deleted() || mesh_.status( vh_target ).deleted() ) {
							continue;
						}
						// setup collapse info
						CollapseInfo ci( mesh_, vp, vh_target );

						// store support
						// the two-ring neighborhood is required for updating the element quality correctly
						support.clear();
						for( auto vv : mesh_.vv_range( ci.v0 ) ) {
							support.insert( vv );
							for( auto vvv : mesh_.vv_range( vv ) ) {
								support.insert( vvv );
							}
						}
						for( auto vv : mesh_.vv_range( ci.v1 ) ) {
							support.insert( vv );
							for( auto vvv : mesh_.vv_range( vv ) ) {
								support.insert( vvv );
							}
						}
						support.erase( ci.v0 );

						// pre-processing
						this->preprocess_collapse( ci );

						// perform collapse
						//mesh_.collapse( v0v1 );
						//++n_collapses;
						////////////////////////////////////////// !!!!!!!!!!!!!!!!!!!!!!
						std::vector<std::vector<Mesh::VertexHandle>> faces;
						for( auto vf : mesh_.vf_range( ci.v0 ) ) {
							std::vector<Mesh::VertexHandle> f;
							f.reserve( 3 );
							for( auto fv : vf.vertices() ) {
								f.push_back( fv );
							}
							faces.push_back( f );
						}
						for( auto& f : faces ) {
							std::replace( f.begin(), f.end(), ci.v0, ci.v1 );
						}
						mesh_.delete_vertex( ci.v0, false );
						for( const auto& f : faces ) {
							mesh_.add_face( f );
						}
						++n_collapses;
						--nv;

						if( update_normals ) {
							// update triangle normals
							vf_it = mesh_.vf_iter( ci.v1 );
							for( ; vf_it.is_valid(); ++vf_it )
								if( !mesh_.status( *vf_it ).deleted() )
									mesh_.set_normal( *vf_it, mesh_.calc_face_normal( *vf_it ) );
						}

						// post-process collapse
						this->postprocess_collapse( ci );

						// if internal boundary has only three edges, add face
						if( nf != 1 ) {
							std::set<TriMesh::HalfedgeHandle> bndrs;
							for( auto h : mesh_.voh_range( ci.v0 ) ) {
								if( mesh_.is_boundary( h ) )
									bndrs.insert( h );
								if( mesh_.is_boundary( mesh_.opposite_halfedge_handle( h ) ) )
									bndrs.insert( mesh_.opposite_halfedge_handle( h ) );
							}
							for( auto h : mesh_.voh_range( ci.v1 ) ) {
								if( mesh_.is_boundary( h ) )
									bndrs.insert( h );
								if( mesh_.is_boundary( mesh_.opposite_halfedge_handle( h ) ) )
									bndrs.insert( mesh_.opposite_halfedge_handle( h ) );
							}
							for( auto h : bndrs ) {
								if( !mesh_.is_boundary( h ) )
									continue;

								std::vector<TriMesh::VertexHandle> vhs( 4 );
								vhs[0] = mesh_.from_vertex_handle( mesh_.prev_halfedge_handle( h ) );
								vhs[1] = mesh_.from_vertex_handle( h );
								vhs[2] = mesh_.to_vertex_handle( h );
								vhs[3] = mesh_.to_vertex_handle( mesh_.next_halfedge_handle( h ) );
								if( vhs[0] == vhs[3] ) {
									mesh_.add_face( { vhs[0], vhs[1], vhs[2] } );
									++nf;
								}
							}
						}

						meshPointsNN.buildIndex();

						// update heap (former one ring of decimated vertex)
						for( s_it = support.begin(), s_end = support.end(); s_it != s_end; ++s_it ) {
							assert( !mesh_.status( *s_it ).deleted() );
							heap_vertex( *s_it );
						}

						// notify observer and stop if the observer requests it
						if( !this->notify_observer( n_collapses ) )
							return n_collapses;

						continue;
					}

					// setup collapse info
					CollapseInfo ci( mesh_, v0v1 );

					// check topological correctness AGAIN !
					if( !this->is_collapse_legal( ci ) )
						continue;

					// store support
					// the two-ring neighborhood is required for updating the element quality correctly
					support.clear();
					for( auto vv : mesh_.vv_range( ci.v0 ) ) {
						support.insert( vv );
						for( auto vvv : mesh_.vv_range( vv ) ) {
							support.insert( vvv );
						}
					}
					for( auto vv : mesh_.vv_range( ci.v1 ) ) {
						support.insert( vv );
						for( auto vvv : mesh_.vv_range( vv ) ) {
							support.insert( vvv );
						}
					}
					support.erase( ci.v0 );

					// adjust complexity in advance (need boundary status)
					++n_collapses;
					--nv;
					if( mesh_.is_boundary( ci.v0v1 ) || mesh_.is_boundary( ci.v1v0 ) )
						--nf;
					else
						nf -= 2;

					// pre-processing
					this->preprocess_collapse( ci );

					// perform collapse
					mesh_.collapse( v0v1 );

					// update triangle normals
					if( update_normals ) {
						vf_it = mesh_.vf_iter( ci.v1 );
						for( ; vf_it.is_valid(); ++vf_it )
							if( !mesh_.status( *vf_it ).deleted() )
								mesh_.set_normal( *vf_it, mesh_.calc_face_normal( *vf_it ) );
					}

					// post-process collapse
					this->postprocess_collapse( ci );

					// if internal boundary has only three edges, add face
					if( nf != 1 && (mesh_.is_boundary( ci.v0 ) || mesh_.is_boundary( ci.v0 )) ) {
						std::set<TriMesh::HalfedgeHandle> bndrs;
						for( auto h : mesh_.voh_range( ci.v0 ) ) {
							if( mesh_.is_boundary( h ) )
								bndrs.insert( h );
							if( mesh_.is_boundary( mesh_.opposite_halfedge_handle( h ) ) )
								bndrs.insert( mesh_.opposite_halfedge_handle( h ) );
						}
						for( auto h : mesh_.voh_range( ci.v1 ) ) {
							if( mesh_.is_boundary( h ) )
								bndrs.insert( h );
							if( mesh_.is_boundary( mesh_.opposite_halfedge_handle( h ) ) )
								bndrs.insert( mesh_.opposite_halfedge_handle( h ) );
						}
						for( auto h : bndrs ) {
							if( !mesh_.is_boundary( h ) )
								continue;

							std::vector<TriMesh::VertexHandle> vhs( 4 );
							vhs[0] = mesh_.from_vertex_handle( mesh_.prev_halfedge_handle( h ) );
							vhs[1] = mesh_.from_vertex_handle( h );
							vhs[2] = mesh_.to_vertex_handle( h );
							vhs[3] = mesh_.to_vertex_handle( mesh_.next_halfedge_handle( h ) );
							if( vhs[0] == vhs[3] ) {
								mesh_.add_face( { vhs[0], vhs[1], vhs[2] } );
								++nf;
							}
						}

						// update index if boundary was involved
						meshPointsNN.buildIndex();
					}


					// update heap (former one ring of decimated vertex)
					for( const auto& s : support ) {
						assert( !mesh_.status( s ).deleted() );
						heap_vertex( s );
					}

					// notify observer and stop if the observer requests it
					if( !this->notify_observer( n_collapses ) )
						return n_collapses;
				}

				// delete heap
				heap_.reset();

				// duplicate vertices that became non-manifold
				duplicate_non_manifold_vertices();


				// DON'T do garbage collection here! It's up to the application.
				return n_collapses;
			}

			void duplicate_non_manifold_vertices() {
				for( const auto& v : mesh_.vertices() ) {
					int c = 0;
					for( const auto& e : v.edges() ) {
						if( e.is_boundary() ) ++c;
					}
					if( c == 0 || c == 2 ) {
						continue;
					}

					// solve non-manifoldness
					auto hIt = v.halfedge();
					while( !hIt.is_boundary() ) {
						hIt = hIt.prev().opp();
					}
					hIt = hIt.prev().opp();

					// collect faces on one
					std::vector<OpenMesh::SmartFaceHandle> faceHandles;
					std::vector<std::vector<VertexHandle>> faceVertexHandles;
					while( !hIt.is_boundary() ) {
						faceHandles.push_back( hIt.face() );
						std::vector<VertexHandle> fvhs;
						for( const auto& fv : hIt.face().vertices() ) {
							fvhs.push_back( fv );
						}
						faceVertexHandles.push_back( fvhs );
						hIt = hIt.prev().opp();
					}

					const auto p = mesh_.point( v );
					auto vNew = mesh_.add_vertex( p );
					for( auto& fvhs : faceVertexHandles ) {
						std::replace( fvhs.begin(), fvhs.end(), v, vNew );
					}

					// delete faces and build new ones with new vertex
					for( const auto& f : faceHandles ) {
						mesh_.delete_face( f, false );
					}
					for( const auto& fvhs : faceVertexHandles ) {
						mesh_.add_face( fvhs );
					}
				}
			}
		public:

			typedef typename Mesh::VertexHandle    VertexHandle;
			typedef typename Mesh::HalfedgeHandle  HalfedgeHandle;

			/// Heap interface
			class HeapInterface
			{
			public:

				HeapInterface( Mesh& _mesh,
							   VPropHandleT<float> _prio,
							   VPropHandleT<int>   _pos )
					: mesh_( _mesh ), prio_( _prio ), pos_( _pos ) {}

				inline bool
					less( VertexHandle _vh0, VertexHandle _vh1 ) {
					return mesh_.property( prio_, _vh0 ) < mesh_.property( prio_, _vh1 );
				}

				inline bool
					greater( VertexHandle _vh0, VertexHandle _vh1 ) {
					return mesh_.property( prio_, _vh0 ) > mesh_.property( prio_, _vh1 );
				}

				inline int
					get_heap_position( VertexHandle _vh ) {
					return mesh_.property( pos_, _vh );
				}

				inline void
					set_heap_position( VertexHandle _vh, int _pos ) {
					mesh_.property( pos_, _vh ) = _pos;
				}


			private:
				Mesh& mesh_;
				VPropHandleT<float>  prio_;
				VPropHandleT<int>    pos_;
			};

			typedef Utils::HeapT<VertexHandle, HeapInterface>  DeciHeap;


		private: //---------------------------------------------------- private methods

			/// Insert vertex in heap
			void heap_vertex( VertexHandle _vh ) {
				//   std::clog << "heap_vertex: " << _vh << std::endl;

				float prio, best_prio( FLT_MAX );
				typename Mesh::HalfedgeHandle heh, collapse_target;
				typename Mesh::VertexHandle collapse_target_vh;

				// find best target in one ring
				typename Mesh::VertexOHalfedgeIter voh_it( mesh_, _vh );
				for( ; voh_it.is_valid(); ++voh_it ) {
					heh = *voh_it;
					CollapseInfo ci( mesh_, heh );

					if( this->is_collapse_legal( ci ) ) {
						prio = this->collapse_priority( ci );
						if( prio >= 0.0 && prio < best_prio ) {
							best_prio = prio;
							collapse_target = heh;
						}
					}
				}

				if( mesh_.is_boundary( _vh ) ) {
					// find best target in local neighborhood
					auto vhNN = meshPointsNN.findNearestNeighbor( mesh_.point( _vh ) );
					// only consider if vhNN is not in one-ring
					bool considerNN = true;
					LOG_ASSERT( _vh != vhNN );
					for( const auto& vv : mesh_.vv_range( _vh ) ) {
						if( vv == vhNN ) {
							considerNN = false;
							break;
						}
					}
					
					if( considerNN ) {
						CollapseInfo ciNN( mesh_, _vh, vhNN );
						prio = this->collapse_priority( ciNN );
						if( prio >= 0.0 && prio < best_prio ) {
							best_prio = prio;
							collapse_target_vh = vhNN;
						}
					}
				}


				// target found -> put vertex on heap
				if( collapse_target_vh.is_valid() ) {
					mesh_.property( collapse_target_, _vh ) = OpenMesh::HalfedgeHandle();
					mesh_.property( collapse_target_vh_, _vh ) = collapse_target_vh;
					mesh_.property( priority_, _vh ) = best_prio;
					//LOG_ASSERT( collapse_target_vh != _vh );
				
					if( heap_->is_stored( _vh ) )
						heap_->update( _vh );
					else
						heap_->insert( _vh );
				}
				// target found -> put vertex on heap
				else if( collapse_target.is_valid() ) {
					//     std::clog << "  added|updated" << std::endl;
					mesh_.property( collapse_target_, _vh ) = collapse_target;
					mesh_.property( collapse_target_vh_, _vh ) = collapse_target_vh;
					mesh_.property( priority_, _vh ) = best_prio;

					if( heap_->is_stored( _vh ) )
						heap_->update( _vh );
					else
						heap_->insert( _vh );
				}

				// not valid -> remove from heap
				else {
					//     std::clog << "  n/a|removed" << std::endl;
					if( heap_->is_stored( _vh ) )
						heap_->remove( _vh );

					mesh_.property( collapse_target_, _vh ) = collapse_target;
					mesh_.property( collapse_target_vh_, _vh ) = collapse_target_vh;
					mesh_.property( priority_, _vh ) = -1;
				}
			}

		private: //------------------------------------------------------- private data


			// reference to mesh
			Mesh& mesh_;

			// heap
		#if (defined(_MSC_VER) && (_MSC_VER >= 1800)) || __cplusplus > 199711L || defined( __GXX_EXPERIMENTAL_CXX0X__ )
			std::unique_ptr<DeciHeap> heap_;
		#else
			std::auto_ptr<DeciHeap> heap_;
		#endif

		// vertex properties
			VPropHandleT<HalfedgeHandle>  collapse_target_;
			VPropHandleT<VertexHandle>    collapse_target_vh_;
			VPropHandleT<float>           priority_;
			VPropHandleT<int>             heap_position_;

		};

	} // END_NS_DECIMATER
} // END_NS_OPENMESH
  //=============================================================================
