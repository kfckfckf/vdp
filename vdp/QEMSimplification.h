#pragma once
#include "MeshDefinition.h"
#include <Eigen/Eigen>
#include <queue>

// Assistant structure of QEM cost
struct EdgeQEM
{
	Mesh::EdgeHandle idx;
	double qem;
	inline friend bool operator<(const EdgeQEM & lhs, const EdgeQEM & rhs) { return (lhs.qem > rhs.qem); }
	inline friend bool operator>(const EdgeQEM & lhs, const EdgeQEM & rhs) { return (lhs.qem < rhs.qem); }
};

// Class of QEM Simplification
class QEMSimplification
{
public:
	QEMSimplification(const Mesh & m);
	virtual ~QEMSimplification(void);
	
	// Simplify the mesh target number of vertex and a error threshold
	void Simplify(int tar_num_v, double err_threshold);

	// get the original index of a simplified mesh
	int OriginalIndex(int currentid) const;

	// Get the simplified mesh
	const Mesh & GetMesh(void) const;

private:
	// Choose the next edge need to collapsed
	bool ChooseCollapasedEdge(double err_threshold);

	// Initialize edge cost
	void InitialEdgeCost(void);

	// Update edge cost after collapse
	void UpdateEdgeCost(const Mesh::VertexHandle & vh);

	Mesh mesh; // the mesh need to be remeshed
	std::vector<std::vector<int>> index_collapse_v_to_v;
	OpenMesh::VPropHandleT<int> vindex; // initial index of each vertex
	OpenMesh::FPropHandleT<double> area_f;            // area of each face
	OpenMesh::FPropHandleT<Eigen::Matrix4d> QuaMat_f; // matrix of each face
	OpenMesh::VPropHandleT<Eigen::Matrix4d> QuaMat_v; // matrix of each vertex
	OpenMesh::EPropHandleT<double> QuaMat_e;
	OpenMesh::EPropHandleT<OpenMesh::Vec3d> new_point; // new position after collapse
	OpenMesh::EPropHandleT<bool> is_cal_e; // is the edge calculated
	OpenMesh::HPropHandleT<bool> is_cal_he; // is the halfedge calculated
	std::priority_queue<EdgeQEM> qpq; // priority queue according to the edge cost
};

