#include "VDPApp.h"
#include "CutGenusZeroMesh.h"
#include "ParameterManager.h"
#include <iostream>
#include <fstream>

VDPApp::VDPApp()
{
	pm = std::make_unique<ParameterManager>();
	pm->AddParameter("inputmesh", "mesh", "", std::string(), false);
	pm->AddParameter("saveintermediate", "intermediate", "a", false);
	pm->AddParameter("usegeodesic", "geodesic", "g", false);
	pm->AddParameter("numvertsimplify", "nsimplify", "p", 13000);
	pm->AddParameter("numrandom", "nrandom", "q", 10);
	pm->AddParameter("numvoting", "nvoting", "r", 2);
	pm->AddParameter("numlocalfeature", "nlocal", "s", 13);
	pm->AddParameter("numpostprocessing", "npost", "t", 5);
	pm->AddParameter("verbose", "verbose", "v", false);
}

VDPApp::~VDPApp()
{
}

VDPApp& VDPApp::GetInstance()
{
	static VDPApp instance;
	return instance;
}

bool VDPApp::ParseCommand(int argc, char* argv[])
{
	if (!pm->ParseCommand(argc, argv))
	{
		std::cerr << "ERROR: wrong argument." << std::endl;
		pm->PrintUsage();
		return false;
	}
	meshfilename = pm->GetParameter<std::string>("inputmesh");
	std::ifstream ifs(meshfilename);
	if (!ifs.is_open())
	{
		std::cerr << "Error: the file " << meshfilename << " does not exist." << std::endl;
		return false;
	}
	ifs.close();

	meshbasename = meshfilename;
	auto point = meshbasename.find_last_of('.');
	if (point != std::string::npos)
	{
		meshbasename = meshbasename.substr(0, point);
	}

	return true;
}

bool VDPApp::Run()
{
	mesh = std::make_unique<Mesh>();
	if (!MeshTools::ReadMesh(*mesh, meshfilename))
	{
		std::cerr << "Error: cannot read mesh from " << meshfilename << std::endl;
		return false;
	}
	if (pm->GetParameter<bool>("verbose"))
	{
		std::cout << "[V, E, F] = [" << mesh->n_vertices() << ", " << mesh->n_edges() << ", " << mesh->n_faces() << "]\n";
	}
	CutMesh();
	return true;
}

void VDPApp::PrintError() const
{
	std::cerr << "ERROR: wrong argument." << std::endl;
	std::cout << "Usage: vdp <mesh> [-a] [-g] [-p <int>] [-q <int>] [-r <int>] [-s <int>] [-t <int>] [-v]\n"
		"Parameters: mesh: file name of a genus zero mesh.\n"
		"            --intermediate, -a: save intermediate results.\n"
		"            --geodesic, -g    : use geodesic distance.\n"
		"            --nsimplify, -p   : vertex number of simplified mesh. default: 13000\n"
		"            --nrandom, -q     : number of random cutting process. default: 10\n"
		"            --nvoting, -r     : voting threshold. default: 2\n"
		"            --nlocal, -s      : local feature size N. default: 13\n"
		"            --npost, -t       : n-ring neighborhood. default: 5\n"
		"            --verbose, -v     : output debug info.\n";
}

void VDPApp::CutMesh()
{
	CutGenusZeroMesh cut(*mesh);
	cut.Init(pm->GetParameter<int>("numlocalfeature"),
		pm->GetParameter<bool>("saveintermediate"),
		pm->GetParameter<bool>("usegeodesic"),
		pm->GetParameter<int>("numvertsimplify"),
		pm->GetParameter<int>("numrandom"),
		pm->GetParameter<int>("numvoting"),
		pm->GetParameter<int>("numpostprocessing"),
		meshbasename,
		pm->GetParameter<bool>("verbose"));
	cut.Run();
}
