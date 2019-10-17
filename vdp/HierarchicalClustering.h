#pragma once
#include "MeshDefinition.h"

// Class of hierarchical clustering algorithm
class HierarchicalClustering
{
public:
	// Construct from a mesh
	HierarchicalClustering(const Mesh & m);
	
	// Destructor
	virtual ~HierarchicalClustering(void);
	
	// Compute the distortion of the map to the other mesh
	// Warning: must be called before clustering feature points
	void ComputeDistortion(const Mesh & mesh2);
	
	// Run the algorithm
	void ComputeClusters(const int & size = 10, const bool & usefirst = true);

	// Cluster by persistence
	void ComputeClustersPersistence(const double & tau);

	// Get feature points
	const std::vector<int> & GetFeaturePoints(void) const { return featurepoints; }

private:
	// Pre-compute for computing distortion
	void PrepareComputeDistortion(void);

	// The recursive hierarchical clustering algorithm
	void ClusterRecursively(const std::vector<int> & dataset, const double & threshold);
	
	// Compute the clustering method at one level
	std::vector<std::vector<int>> ClusterOnce(const std::vector<int> & dataset, const double & threshold);
	
	// Compute connected regions for a face set
	std::vector<std::vector<int>> ComputeConnectedRegions(
		const std::vector<int> & dataset, const double & threshold) const;
	
	// Find a feature point in a connected region
	void AddFeaturePointInRegion(const std::vector<int> & region);
	
	// Filter faces by a distortion value
	std::vector<int> FilterDis(const std::vector<int> & dataset,
		const double & filter) const;
	
	// Filter faces by midian distortion
	std::vector<int> FilterMid(const std::vector<int> & dataset) const;
	
	// Filter faces by average distortion
	std::vector<int> FilterAvg(const std::vector<int> & dataset) const;
	
	const Mesh & mesh;// constant reference to the mesh
	std::vector<double> facearea, fpx1, fpx2, fpy2; // for distortion
	std::vector<double> facedistortion; // updated in ComputeDistortion()
	std::vector<int> isfaceselected; // updated in ComputeClusters()
	std::vector<int> featurepoints; // result feature points
	bool isfirstclustering; // updated in ComputeClusters()
};

