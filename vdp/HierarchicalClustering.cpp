#include "HierarchicalClustering.h"
#include <set>
#include <list>

HierarchicalClustering::HierarchicalClustering(const Mesh& m)
	:mesh(m),
	isfirstclustering(true)
{
	PrepareComputeDistortion();
}

HierarchicalClustering::~HierarchicalClustering(void)
{
}

void HierarchicalClustering::ComputeDistortion(const Mesh& mesh2)
{
	double alpha = 0.5;
	if (mesh2.n_faces() != mesh.n_faces())
	{
		return;
	}
	facedistortion.resize(mesh.n_faces());
	for (const auto& fh : mesh2.faces())
	{
		auto fid = fh.idx();
		auto heh = mesh2.halfedge_handle(fh);
		const auto& p0 = mesh2.point(mesh2.from_vertex_handle(heh));
		const auto& p1 = mesh2.point(mesh2.to_vertex_handle(heh));
		const auto& p2 = mesh2.point(mesh2.to_vertex_handle(mesh2.next_halfedge_handle(heh)));
		auto np = (p1 - p0) % (p2 - p0);
		double area = np.norm() * 0.5;
		np.normalize();
		auto qx1 = (p1 - p0).norm();
		auto pe1 = (p1 - p0).normalized();
		auto pe2 = (np % pe1).normalized();;
		auto qx2 = (p2 - p0) | pe1;
		auto qy2 = (p2 - p0) | pe2;
		double a = qx1 / fpx1[fid];
		double b = qy2 / fpy2[fid];
		double det = a * b;
		if (det > 0.0)
		{
			double c = ((-qx1 * fpx2[fid] + qx2 * fpx1[fid]) / fpx1[fid] / fpy2[fid]);
			facedistortion[fid] = alpha * (a * a + b * b + c * c) / det * 0.5 + (1 - alpha) * (det + 1.0 / det) * 0.5;
			//facedistortion[fid] = (det + 1.0 / det) * 0.5;
			//facedistortion[fid] = 1.0 / det;
		}
		else
		{
			facedistortion[fid] = std::numeric_limits<double>::infinity();
		}
	}
}

void HierarchicalClustering::ComputeClusters(const int& size, const bool& usefirst)
{
	std::vector<int> dataset;
	for (int i = 0; i < facedistortion.size(); i++)
	{
		if (!mesh.is_boundary(mesh.face_handle(i)) && facedistortion[i] > 0)
		{
			dataset.push_back(i);
		}
	}
	isfaceselected.assign(mesh.n_faces(), 0);
	featurepoints.clear();
	isfirstclustering = usefirst;
	ClusterRecursively(dataset, size);
}

void HierarchicalClustering::ComputeClustersPersistence(const double& tau)
{
	// Step 1: compute a point wise function
	std::vector<double> vertdis(mesh.n_vertices());
	std::vector<int> v2v(mesh.n_vertices());
	std::vector<int> v2vi(mesh.n_vertices());
	for (const auto& vh : mesh.vertices())
	{
		double vdis = 0.0;
		double varea = 0.0;
		for (const auto& vfh : mesh.vf_range(vh))
		{
			vdis += facearea[vfh.idx()] * facedistortion[vfh.idx()];
			varea += facearea[vfh.idx()];
		}
		vdis /= varea;
		vertdis[vh.idx()] = vdis;
		v2v[vh.idx()] = vh.idx();
	}

	// Step 2: sort the vertices
	std::sort(v2v.begin(), v2v.end(), [&vertdis](const auto& vv1, const auto& vv2) {
		return vertdis[vv1] > vertdis[vv2];
		});
	for (size_t i = 0; i < v2v.size(); i++)
	{
		v2vi[v2v[i]] = static_cast<int>(i);
	}

	// Step 3: initialize data structure
	std::list<std::pair<size_t, std::set<size_t>>> U;
	std::vector<int> g(mesh.n_vertices(), -1);
	for (size_t i = 0; i < v2v.size(); i++)
	{
		// Step 3.1: compute higher neighbor vertices
		std::vector<int> nei;
		for (const auto& vvh : mesh.vv_range(mesh.vertex_handle(v2v[i])))
		{
			if (v2vi[vvh.idx()] < i)
			{
				nei.push_back(v2vi[vvh.idx()]);
			}
		}
		if (nei.empty())
		{
			// vertex i is a peak of f within mesh
			// create a new entry e in U and attach vertex i to it
			U.push_back(std::make_pair(i, std::set<size_t>({ i })));

		}
		else
		{
			// vertex i is not a peak of f within mesh
			g[i] = nei[0];
			double maxf = vertdis[v2v[g[i]]];
			for (const auto& n : nei)
			{
				if (vertdis[v2v[n]] > maxf)
				{
					maxf = vertdis[v2v[n]];
					g[i] = n;
				}
			}
			auto ei = U.end();
			for (auto eit = U.begin(); eit != U.end(); ++eit)
			{
				if (eit->second.find(g[i]) != eit->second.end())
				{
					ei = eit;
					break;
				}
			}
			ei->second.insert(i);
			for (const auto& j : nei)
			{
				auto e = U.end();
				for (auto eit = U.begin(); eit != U.end(); ++eit)
				{
					if (eit->second.find(j) != eit->second.end())
					{
						e = eit;
						break;
					}
				}
				auto mine = vertdis[v2v[e->first]] < vertdis[v2v[ei->first]] ? e : ei;
				auto maxe = vertdis[v2v[e->first]] < vertdis[v2v[ei->first]] ? ei : e;
				if (e->first != ei->first && vertdis[v2v[mine->first]] < vertdis[v2v[i]] + tau)
				{
					maxe->second.insert(mine->second.begin(), mine->second.end());
					ei = maxe;
					U.erase(mine);
				}
			}
		}
	}
	featurepoints.clear();
	for (const auto& e : U)
	{
		featurepoints.push_back(v2v[e.first]);
	}
}

void HierarchicalClustering::PrepareComputeDistortion(void)
{
	facearea.resize(mesh.n_faces());
	fpx1.resize(mesh.n_faces());
	fpx2.resize(mesh.n_faces());
	fpy2.resize(mesh.n_faces());
	for (const auto& fh : mesh.faces())
	{
		auto heh = mesh.halfedge_handle(fh);
		const auto& p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto& p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto& p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
		auto np = (p1 - p0) % (p2 - p0);
		double area = np.norm() * 0.5;
		np.normalize();
		fpx1[fh.idx()] = (p1 - p0).norm();
		auto pe1 = (p1 - p0).normalized();
		auto pe2 = (np % pe1).normalized();;
		fpx2[fh.idx()] = (p2 - p0) | pe1;
		fpy2[fh.idx()] = (p2 - p0) | pe2;
		facearea[fh.idx()] = area;
	}
}

void HierarchicalClustering::ClusterRecursively(const std::vector<int>& dataset, const double& threshold)
{
	std::vector<std::vector<int>> regions = ComputeConnectedRegions(dataset, threshold);
	for (const auto& r : regions)
	{
		AddFeaturePointInRegion(r);
	}
	for (const auto& r : regions)
	{
		std::vector<int> d_s;
		if (isfirstclustering)
		{
			d_s = FilterDis(r, 2.0);
			isfirstclustering = false;
		}
		else
		{
			d_s = FilterMid(r);
		}
		if (d_s.size() > threshold)
		{
			ClusterRecursively(d_s, threshold);
		}
	}
}

std::vector<std::vector<int>> HierarchicalClustering::ClusterOnce(const std::vector<int>& dataset, const double& threshold)
{
	std::vector<std::vector<int>> regions = ComputeConnectedRegions(dataset, threshold);
	for (const auto& r : regions)
	{
		AddFeaturePointInRegion(r);
	}
	std::vector<std::vector<int>> newdatasets;
	for (const auto& r : regions)
	{
		std::vector<int> d_s;
		if (isfirstclustering)
		{
			d_s = FilterDis(r, 2.0);
			isfirstclustering = false;
		}
		else
		{
			d_s = FilterMid(r);
		}
		if (d_s.size() > threshold)
		{
			newdatasets.push_back(std::move(d_s));
		}
	}
	return newdatasets;
}

std::vector<std::vector<int>> HierarchicalClustering::ComputeConnectedRegions(const std::vector<int>& dataset, const double& threshold) const
{
	std::vector<int> isaccessible(mesh.n_faces(), 0);
	for (const auto& val : dataset)
	{
		isaccessible[val] = 1;
	}
	std::vector<std::vector<int>> regions;
	std::vector<int> isfacevisited(mesh.n_faces(), 0);
	std::vector<int> isvertexvisited(mesh.n_vertices(), 0);
	for (const auto& fid : dataset)
	{
		if (!isfacevisited[fid])
		{
			std::vector<int> r;
			r.push_back(fid);
			isfacevisited[fid] = 1;
			for (size_t j = 0; j < r.size(); j++)
			{
				for (const auto& fvh : mesh.fv_range(mesh.face_handle(r[j])))
				{
					if (!isvertexvisited[fvh.idx()])
					{
						isvertexvisited[fvh.idx()] = 1;
						for (const auto& fvheh : mesh.voh_range(fvh))
						{
							int fid = mesh.face_handle(fvheh).idx();
							if (isaccessible[fid] && !isfacevisited[fid])
							{
								r.push_back(fid);
								isfacevisited[fid] = 1;
							}
						}
					}
				}
			}
			if (r.size() > threshold)
			{
				regions.push_back(std::move(r));
			}
		}
	}
	return regions;
}

void HierarchicalClustering::AddFeaturePointInRegion(const std::vector<int>& region)
{
	double fmaxdis = 0.0;
	int fidmax = 0;
	for (const auto& fid : region)
	{
		if (facedistortion[fid] > fmaxdis)
		{
			fmaxdis = facedistortion[fid];
			fidmax = fid;
		}
	}
	if (!isfaceselected[fidmax])
	{
		isfaceselected[fidmax] = 1;
		auto fh = mesh.face_handle(fidmax);
		int vid = 0;
		double vmaxdis = 0.0;
		for (const auto& fvh : mesh.fv_range(fh))
		{
			double vdis = 0.0;
			int vvalence = 0;
			for (const auto& vfh : mesh.vf_range(fvh))
			{
				vdis += facedistortion[vfh.idx()];
				++vvalence;
			}
			vdis /= vvalence;
			if (vmaxdis < vdis)
			{
				vid = fvh.idx();
				vmaxdis = vdis;
			}
		}
		featurepoints.push_back(vid);
	}
}

std::vector<int> HierarchicalClustering::FilterDis(const std::vector<int>& dataset, const double& filter) const
{
	std::vector<int> newdataset;
	for (const auto& id : dataset)
	{
		if (facedistortion[id] > filter)
		{
			newdataset.push_back(id);
		}
	}
	return newdataset;
}

std::vector<int> HierarchicalClustering::FilterMid(const std::vector<int>& dataset) const
{
	std::vector<double> tempdistortion;
	tempdistortion.reserve(dataset.size());
	for (const auto& id : dataset)
	{
		tempdistortion.push_back(facedistortion[id]);
	}
	std::sort(tempdistortion.begin(), tempdistortion.end());
	return FilterDis(dataset, tempdistortion[dataset.size() / 2]);
}

std::vector<int> HierarchicalClustering::FilterAvg(const std::vector<int>& dataset) const
{
	double avg = 0.0;
	for (const auto& id : dataset)
	{
		avg += facedistortion[id];
	}
	return FilterDis(dataset, avg / dataset.size());
}
