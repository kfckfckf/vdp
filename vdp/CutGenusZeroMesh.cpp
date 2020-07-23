#include "CutGenusZeroMesh.h"
#include "FarthestPointSampling.h"
#include "MeshCut.h"
#include "HierarchicalClustering.h"
#include "MinimalSpanningTree.h"
#include <iostream>
#include <fstream>

CutGenusZeroMesh::CutGenusZeroMesh(const Mesh& mesh)
	:CutClosedMesh(mesh)
{
}


CutGenusZeroMesh::~CutGenusZeroMesh()
{
}

void CutGenusZeroMesh::Run(void)
{
	if (verbose)
	{
		std::cout << "Start remeshing.\n";
	}
	
	auto tstart = std::chrono::system_clock::now();
	Remesh(nvsimplify);
	if (verbose)
	{
		auto tend = std::chrono::system_clock::now();
		std::cout << "Finished. Time: " << std::chrono::duration<double>(tend - tstart).count() << "s.\n";
	}
	
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
		auto tstart = std::chrono::system_clock::now();
		CutMesh(landmarks);
		if (verbose)
		{
			auto tend = std::chrono::system_clock::now();
			std::cout << "Cut the mesh finished. Time: " << std::chrono::duration<double>(tend - tstart).count() << "s.\n";
		}
		
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
		tstart = std::chrono::system_clock::now();
		Parameterization();
		if (verbose)
		{
			auto tend = std::chrono::system_clock::now();
			std::cout << "Conformal parameterization finished. Time: " << std::chrono::duration<double>(tend - tstart).count() << "s.\n";
		}
		
		if (saveintermediate)
		{
			MeshTools::WriteMesh(currentmesh, basefile + "_open_conformal" + std::to_string(i).c_str() + ".obj", 16);
		}
		tstart = std::chrono::system_clock::now();
		ClusterLandmarks(size);
		const auto& interlandmarks = hc->GetFeaturePoints();
		if (verbose)
		{
			auto tend = std::chrono::system_clock::now();
			std::cout << "Clustering finished. Point size: " << interlandmarks.size() << ". Time: " << std::chrono::duration<double>(tend - tstart).count() << "s.\n";
		}

		intermediatelandmarks.push_back(interlandmarks);
		if (saveintermediate)
		{
			SaveLandmarks(basefile + "_interlandmarks" + std::to_string(i) + ".txt", interlandmarks);
		}
	}

	auto tstart1 = std::chrono::system_clock::now();
	auto landmarks = RemoveNearbyLandmarks(VoteLandmarks(intermediatelandmarks));
	landmarks = FindOrginalVertexIndex(landmarks);
	if (verbose)
	{
		auto tend = std::chrono::system_clock::now();
		std::cout << "Remove nearby landmarks finished. Time: " << std::chrono::duration<double>(tend - tstart1).count() << "s.\n";
	}
	
	SaveLandmarks(basefile + "_landmarks.txt", landmarks);
	InitialCutMesh(orimesh);
	tstart1 = std::chrono::system_clock::now();
	double edgelength = MeshTools::AverageEdgeLength(orimesh);
	if (landmarks.size() > 1)
	{
		CutMesh(landmarks);
	}
	else
	{
		std::cout << "Too few landmarks found." << std::endl;
		return;
	}
	if (verbose)
	{
		auto tend = std::chrono::system_clock::now();
		std::cout << "Final cut generated. Time: " << std::chrono::duration<double>(tend - tstart1).count() << "s.\n";
		std::cout << "Total time: " << std::chrono::duration<double>(tend - tstart).count() << "s.\n";
	}
	if (!SaveCuts(basefile + "_cut.txt", mc->GetCutVertices(), mc->GetCutEdges()))
	{
		std::cout << "Warning: save cut failed." << std::endl;
	}
	MeshTools::WriteMesh(currentmesh, basefile + "_open.obj", 16);
}
