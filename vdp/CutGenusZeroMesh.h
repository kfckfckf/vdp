#pragma once
#include "CutClosedMesh.h"

// Class of cutting a genus-zero mesh
class CutGenusZeroMesh :
	public CutClosedMesh
{
public:
	CutGenusZeroMesh(const Mesh& mesh);
	virtual ~CutGenusZeroMesh();

	// main algorithm of cutting a genus-zero mesh
	virtual void Run(void);
};

