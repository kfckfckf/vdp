#pragma once
#include "MeshDefinition.h"

class ShortestPath;

// Class of farthest point sampling
class FarthestPointSampling
{
public:
	FarthestPointSampling(const Mesh& m, const bool& geodesic = false);
	virtual ~FarthestPointSampling(void);

	// Compute n samples
	std::vector<int> ComputeSamples(const int& nsample) const;

private:
	// Compute distance from a point to a set
	double PointSetDistance(const int& vid, const std::vector<std::vector<double>>& distanceset) const;

	// Compute all vertices distance
	std::vector<double> AllVerticesDistance(const int& vid) const;

	const Mesh& mesh; // const reference to the mesh
	bool isgeodesic = false;
	std::unique_ptr<ShortestPath> shortestpath = nullptr;
};

