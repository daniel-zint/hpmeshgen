#pragma once

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <OpenMesh/Core/System/config.h>

struct MyTraits : public OpenMesh::DefaultTraits
{
    VertexTraits
    {
      bool is_feature = false;
    };
};


typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> TriMesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> PolyMesh;