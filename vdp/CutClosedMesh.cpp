#include "CutClosedMesh.h"
#include "QEMSimplification.h"
#include "Subdivision.h"
#include "HierarchicalClustering.h"
#include "MinimalSpanningTree.h"
#include "MeshCut.h"
#include "KPNewton.h"
#include <fstream>
#include <iostream>
#include <set>

CutClosedMesh::CutClosedMesh(const Mesh & mesh)
	:orimesh(mesh)
{
}

CutClosedMesh::~CutClosedMesh(void)
{
}

void CutClosedMesh::Init(const int & s, const bool & issaveintermediate, const bool & isgeodesic, const int & nvsimplified, const int & nrandomcut, const int & nvotes, const int & npostprocess, const std::string & filename)
{
	size = s;
	saveintermediate = issaveintermediate;
	geodesic = isgeodesic;
	nvsimplify = nvsimplified;
	nrandom = nrandomcut;
	nvoting = nvotes;
	npost = npostprocess;
	basefile = filename;
	if (saveintermediate && filename.empty())
	{
		std::cerr << "No base file name specified!" << std::endl;
	}
}

bool CutClosedMesh::Remesh(const int & targetnv)
{
	originalvid.clear();
	if (orimesh.n_vertices() > targetnv)
	{
		QEMSimplification qs(orimesh);
		qs.Simplify(targetnv, 1);
		MeshTools::Reassign(qs.GetMesh(), remesh);
		for (const auto & vh : remesh.vertices())
		{
			originalvid.push_back(qs.OriginalIndex(vh.idx()));
		}
	}
	else
	{
		Subdivision sub(orimesh);
		sub.Subdivide(targetnv);
		MeshTools::Reassign(sub.GetMesh(), remesh);
		for (const auto & vh : remesh.vertices())
		{
			originalvid.push_back(sub.OriginalIndex(vh.idx()));
		}
	}
	hc = std::make_unique<HierarchicalClustering>(remesh);
	InitialCutMesh(remesh);
	return true;
}

bool CutClosedMesh::ClusterLandmarks(const int & size) const
{
	if (hc)
	{
		hc->ComputeDistortion(currentmesh);
		hc->ComputeClusters(size);
		return true;
	}
	return false;
}

bool CutClosedMesh::InitialCutMesh(const Mesh & mesh)
{
	mst = std::make_unique<MinimalSpanningTree>(mesh);
	mc = std::make_unique<MeshCut>(mesh);
	return true;
}

bool CutClosedMesh::CutMesh(const std::vector<int> & landmarks)
{
	if (mc && mst)
	{
		mst->SetNodes(landmarks);
		mst->Prim();
		mc->mesh_cut(mst->GetSpanningTreeEdges());
		currentmesh = mc->GetCutedMesh();
		return true;
	}
	return false;
}

bool CutClosedMesh::Parameterization(void)
{
	KPNewton kpn(currentmesh);
	kpn.Tutte();
	kpn.PrepareDataFree();
	kpn.RunFree(KPNewton::MIPS);
	kpn.RunFree(KPNewton::AMIPS);
	kpn.UpdateMesh();
	double scale = std::sqrt(MeshTools::Area(remesh) / MeshTools::Area(currentmesh));
	for (const auto & vh : currentmesh.vertices())
	{
		currentmesh.point(vh) *= scale;
	}
	return true;
}

std::vector<int> CutClosedMesh::FindOrginalVertexIndex(const std::vector<int>& landmarks) const
{
	std::set<int> orilandmarks;
	for (const auto & id : landmarks)
	{
		orilandmarks.insert(originalvid[id]);
	}
	return std::vector<int>(orilandmarks.begin(), orilandmarks.end());
}

std::vector<std::pair<int, int>> CutClosedMesh::VoteLandmarks(const std::vector<std::vector<int>>& intermediatelandmarks) const
{
	std::vector<std::pair<int, int>> landmarks;
	std::map<int, int> counter;
	for (size_t i = 0; i < intermediatelandmarks.size(); i++)
	{
		for (const auto & k : intermediatelandmarks[i])
		{
			if (counter.find(k) == counter.end())
			{
				counter[k] = 0;
			}
			counter[k]++;
		}
	}
	for (const auto & c : counter)
	{
		if (c.second > nvoting)
		{
			landmarks.push_back({ c.first, c.second });
		}
	}
	return landmarks;
}

std::vector<int> CutClosedMesh::RemoveNearbyLandmarks(const std::vector<std::pair<int, int>>& landmarkwithvote) const
{
	std::map<int, int> landmarkset;
	for (const auto & l : landmarkwithvote)
	{
		bool bfound = false;
		std::set<int> neighborvertices;
		neighborvertices.insert(l.first);
		std::vector<int> onering;
		onering.push_back(l.first);
		for (int j = 0; j < npost; j++)
		{
			std::vector<int> newonering;
			for (const auto & vid : onering)
			{
				for (const auto & vvh : remesh.vv_range(remesh.vertex_handle(vid)))
				{
					if (neighborvertices.find(vvh.idx()) == neighborvertices.end())
					{
						neighborvertices.insert(vvh.idx());
						newonering.push_back(vvh.idx());
						std::map<int, int>::iterator founded;
						if ((founded = landmarkset.find(vvh.idx())) != landmarkset.end())
						{
							if (founded->second < l.second)
							{
								landmarkset.erase(founded);
							}
							else
							{
								bfound = true;
							}
						}
					}
				}
			}
			onering = newonering;
		}
		if (!bfound)
		{
			landmarkset[l.first] = l.second;
		}
	}
	std::vector<int> landmarks;
	for (const auto & l : landmarkset)
	{
		landmarks.push_back(l.first);
	}
	return landmarks;
}

bool CutClosedMesh::SaveCuts(const std::string & filename, const std::vector<int>& iscutvertices, const std::vector<int>& iscutedges)
{
	std::ofstream ofs(filename);
	if (!ofs.is_open())
	{
		std::cerr << "Error: cannot open file " << filename << std::endl;
		return false;
	}
	ofs << "VERTICES" << std::endl;
	for (size_t i = 0; i < iscutvertices.size(); i++)
	{
		if (iscutvertices[i])
		{
			ofs << i << std::endl;
		}
	}
	ofs << "EDGES" << std::endl;
	for (size_t i = 0; i < iscutedges.size(); i++)
	{
		if (iscutedges[i])
		{
			ofs << i << std::endl;
		}
	}
	ofs.close();
	return true;
}

bool CutClosedMesh::SaveLandmarks(const std::string & filename, const std::vector<int>& landmarks)
{
	std::ofstream ofs(filename);
	if (!ofs.is_open())
	{
		std::cerr << "Error: cannot save landmarks to " << filename << std::endl;
		return false;
	}
	for (const auto & l : landmarks)
	{
		ofs << l << std::endl;
	}
	ofs.close();
	return true;
}
