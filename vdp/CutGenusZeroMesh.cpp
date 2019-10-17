#include "CutGenusZeroMesh.h"
#include "FarthestPointSampling.h"
#include "MeshCut.h"
#include "HierarchicalClustering.h"
#include "MinimalSpanningTree.h"
#include <iostream>
#include <fstream>

CutGenusZeroMesh::CutGenusZeroMesh(const Mesh & mesh)
	:CutClosedMesh(mesh)
{
}


CutGenusZeroMesh::~CutGenusZeroMesh()
{
}

void CutGenusZeroMesh::Run(void)
{
	std::cout << "Start remeshing." << std::endl;
	auto t = clock();
	Remesh(nvsimplify);
	std::cout << "Finished. Time: " << (clock() - t) / double(CLOCKS_PER_SEC) << "s." << std::endl;
	if (saveintermediate)
	{
		MeshTools::WriteMesh(remesh, basefile + "_remesh.obj", 16);
	}

	FarthestPointSampling cc(remesh, geodesic);
	auto samplepoints = cc.ComputeSamples(2 * nrandom);

	std::vector<std::vector<int>> intermediatelandmarks;
	for (size_t i = 0; i < samplepoints.size() / 2; i++)
	{
		std::vector<int> landmarks = { samplepoints[2 * i], samplepoints[2 * i + 1] };
		if (saveintermediate)
		{
			if (!SaveLandmarks(basefile + "_landmarks" + std::to_string(i) + ".txt", landmarks))
			{
				std::cout << "Warning: save landmarks failed." << std::endl;
			}
		}
		auto t = clock();
		CutMesh(landmarks);
		std::cout << "Cut the mesh finished. Time: " << (clock() - t) / double(CLOCKS_PER_SEC) << "s." << std::endl;
		if (saveintermediate)
		{
			if (!SaveCuts(basefile + "_cut" + std::to_string(i) + ".txt",
				mc->GetCutVertices(),
				mc->GetCutEdges()))
			{
				std::cout << "Warning: save cut failed." << std::endl;
			}
			MeshTools::WriteMesh(currentmesh, basefile + "_open" + std::to_string(i).c_str() + ".obj", 16);
		}
		t = clock();
		Parameterization();
		std::cout << "Conformal parameterization finished. Time: " << (clock() - t) / double(CLOCKS_PER_SEC) << "s." << std::endl;
		if (saveintermediate)
		{
			MeshTools::WriteMesh(currentmesh, basefile + "_open_conformal" + std::to_string(i).c_str() + ".obj", 16);
		}
		t = clock();
		ClusterLandmarks(size);
		const auto & interlandmarks = hc->GetFeaturePoints();
		std::cout << "Clustering finished. Point size: " << interlandmarks.size()
			<< ". Time: " << (clock() - t) / double(CLOCKS_PER_SEC) << "s." << std::endl;

		intermediatelandmarks.push_back(interlandmarks);
		if (saveintermediate)
		{
			SaveLandmarks(basefile + "_interlandmarks" + std::to_string(i) + ".txt", interlandmarks);
		}
	}

	auto t1 = clock();
	auto landmarks = RemoveNearbyLandmarks(VoteLandmarks(intermediatelandmarks));
	landmarks = FindOrginalVertexIndex(landmarks);
	std::cout << "Remove nearby landmarks finished. Time: " << (clock() - t1) / double(CLOCKS_PER_SEC) << "s." << std::endl;
	SaveLandmarks(basefile + "_landmarks.txt", landmarks);
	InitialCutMesh(orimesh);
	t1 = clock();
	double edgelength = MeshTools::AverageEdgeLength(orimesh);
	std::cout << edgelength << std::endl;
	if (landmarks.size() > 1)
	{
		CutMesh(landmarks);
	}
	else
	{
		std::cout << "Too few landmarks found." << std::endl;
		return;
	}
	std::cout << "Final cut generated. Time: " << (clock() - t1) / double(CLOCKS_PER_SEC) << "s." << std::endl;
	std::cout << "Total time: " << (clock() - t) / double(CLOCKS_PER_SEC) << "s." << std::endl;
	if (!SaveCuts(basefile + "_cut.txt", mc->GetCutVertices(), mc->GetCutEdges()))
	{
		std::cout << "Warning: save cut failed." << std::endl;
	}
	MeshTools::WriteMesh(currentmesh, basefile + "_open.obj", 16);
}
