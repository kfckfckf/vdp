#pragma once
#include "MeshDefinition.h"

// Class of mesh subdivision
class Subdivision
{
public:
	Subdivision(const Mesh & m);
	~Subdivision(void);

	// Subdivide the mesh
	void Subdivide(int tarnum);

	// Get the subdivided mesh
	const Mesh & GetMesh(void) const;

	// Get the original index from a subdivided mesh
	int OriginalIndex(int vid) const;

private:
	Mesh mesh; // mesh need to be subdivided
	OpenMesh::VPropHandleT<int> vindex; // record the original indices
	OpenMesh::VPropHandleT<double> curvature; // used to transfer information
};

