#pragma once
#include "MeshDefinition.h"

// Class of mesh cutter
class MeshCut
{
public:
	MeshCut(const Mesh & mesh_);
	virtual ~MeshCut(void);

	// Get edges are boundary or not
	const std::vector<int> & GetCutEdges(void) const { return candidate_seam_edge; }

	// Get vertices are boundary or not
	const std::vector<int> & GetCutVertices(void) const { return candidate_seam_vertex; }

	// Get resulting open mesh
	const Mesh & GetCutedMesh(void) const { return cuted_mesh; }

	// Cutting the mesh
	void mesh_cut(const std::vector<std::vector<int>> & boundary_point);

private:
	// Mark vertices and edges as boundaries
	void MarkVertexEdge(const std::vector<std::vector<int>> & boundarypoints);

	// Compute valence for each cut vertices
	void ComputeValence(void);
	
	const Mesh & mesh; // const reference to the original mesh
	Mesh cuted_mesh;   // resulting open mesh
	std::vector<int> he_to_idx; // to vertex index of a halfedge
	std::vector<int> idx_to_mesh_idx; // vertex indices correspondence
	int boundary_number; // number of boundary edges
	std::vector<int> candidate_seam_vertex; // vertex indices on cut
	std::vector<int> candidate_cut_valence; // vertex valence wrt the cut
	std::vector<int> candidate_seam_edge;   // edge indices on cut
};
