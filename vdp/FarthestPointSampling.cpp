#include "FarthestPointSampling.h"
#include "ShortestPath.h"

FarthestPointSampling::FarthestPointSampling(const Mesh & m, const bool & geodesic)
	:mesh(m),
	isgeodesic(geodesic)
{
	if (isgeodesic)
	{
		shortestpath = std::make_unique<ShortestPath>(mesh);
	}
}

FarthestPointSampling::~FarthestPointSampling()
{
}

std::vector<int> FarthestPointSampling::ComputeSamples(const int & nsample) const
{
	std::vector<int> selectpts;
	std::vector<std::vector<double>> selectdistance;
	selectpts.push_back(0);
	selectdistance.push_back(AllVerticesDistance(0));
	int nv = static_cast<int>(mesh.n_vertices());
	for (int i = 0; i < nsample; i++)
	{
		int farthestid = 0;
		double dist = 0.0;
		for (int j = 0; j < nv; j++)
		{
			double disttemp = PointSetDistance(j, selectdistance);
			if (disttemp > dist)
			{
				dist = disttemp;
				farthestid = j;
			}
		}
		selectpts.push_back(farthestid);
		selectdistance.push_back(AllVerticesDistance(farthestid));
	}
	return selectpts;
}

double FarthestPointSampling::PointSetDistance(const int & vid, const std::vector<std::vector<double>> & distanceset) const
{
	const auto & p = mesh.point(mesh.vertex_handle(vid));
	double dist = DBL_MAX;
	for (size_t i = 0; i < distanceset.size(); i++)
	{
		double disttemp = distanceset[i][vid];
		dist = disttemp < dist ? disttemp : dist;
	}
	return dist;
}

std::vector<double> FarthestPointSampling::AllVerticesDistance(const int & vid) const
{
	const auto & p = mesh.point(mesh.vertex_handle(vid));
	std::vector<double> dist;
	dist.reserve(mesh.n_vertices());
	if (isgeodesic)
	{
		shortestpath->Compute(vid, std::vector<int>(mesh.n_vertices(), 1), {}, (int)mesh.n_vertices());
		for (const auto & vh : mesh.vertices())
		{
			dist.push_back(shortestpath->GetLength(vh.idx()));
		}
	}
	else
	{
		for (const auto & vh : mesh.vertices())
		{
			dist.push_back((mesh.point(vh) - p).length());
		}
	}
	return dist;
}
