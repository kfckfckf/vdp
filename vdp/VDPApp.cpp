#include "VDPApp.h"
#include <iostream>
#include <fstream>
#include "CutGenusZeroMesh.h"

VDPApp::VDPApp(void)
{
}

VDPApp::~VDPApp(void)
{
}

VDPApp & VDPApp::GetInstance(void)
{
	static VDPApp instance;
	return instance;
}

bool VDPApp::ParseCommand(int argc, char * argv[])
{
	std::vector<int> filenameid;
	for (int i = 1; i < argc; i++)
	{
		if (argv[i][0] == '-') // options
		{
			switch (argv[i][1])
			{
			case 'a':
				saveintermediate = true;
				break;
			case 'g': // use geodesic distance (default: use Euclidean)
				geodesic = true;
				break;
			case 'p': // number of simplified mesh (default: 13000)
				i++;
				if (i < argc)
				{
					std::istringstream iss(argv[i]);
					if (iss >> nvsimplify)
					{
						std::cout << nvsimplify << std::endl;
						break;
					}
				}
				PrintError();
				return false;
				break;
			case 'q': // number of random generation process (default: 10)
				i++;
				if (i < argc)
				{
					std::istringstream iss(argv[i]);
					if (iss >> nrandom)
					{
						std::cout << nrandom << std::endl;
						break;
					}
				}
				PrintError();
				return false;
				break;
			case 'r': // number of voting threshold (default: 2)
				i++;
				if (i < argc)
				{
					std::istringstream iss(argv[i]);
					if (iss >> nvoting)
					{
						std::cout << nvoting << std::endl;
						break;
					}
				}
				PrintError();
				return false;
				break;
			case 's': // N
				i++;
				if (i < argc)
				{
					std::istringstream iss(argv[i]);
					if (iss >> size)
					{
						std::cout << size << std::endl;
						break;
					}
				}
				PrintError();
				return false;
				break;
			case 't': // n-ring
				i++;
				if (i < argc)
				{
					std::istringstream iss(argv[i]);
					if (iss >> npost)
					{
						std::cout << npost << std::endl;
						break;
					}
				}
				PrintError();
				return false;
				break;
			default:
				PrintError();
				return false;
				break;
			}
		}
		else
		{
			filenameid.push_back(i);
		}
	}

	if (filenameid.size() == 1)
	{

	}
	else
	{
		PrintError();
		return false;
	}
	meshfilename = argv[filenameid[0]];
	std::ifstream ifs(meshfilename);
	if (!ifs.is_open())
	{
		std::cerr << "Error: the file " << meshfilename << " does not exist." << std::endl;
		return false;
	}
	ifs.close();

	meshpath = ".";
	auto slash = meshfilename.find_last_of('/');
	auto backslash = meshfilename.find_last_of('\\');
	if (slash != std::string::npos && backslash != std::string::npos)
	{
		slash = slash > backslash ? slash : backslash;
		meshpath = meshfilename.substr(0, slash);
	}
	else if (backslash != std::string::npos)
	{
		slash = backslash;
		meshpath = meshfilename.substr(0, slash);
	}
	else if (slash != std::string::npos)
	{
		meshpath = meshfilename.substr(0, slash);
	}

	meshbasename = meshfilename.substr(slash + 1);
	auto point = meshbasename.find_last_of('.');
	if (point != std::string::npos)
	{
		meshbasename = meshbasename.substr(0, point);
	}

	return true;
}

bool VDPApp::Run(void)
{
	if (!MeshTools::ReadMesh(mesh, meshfilename))
	{
		std::cerr << "Error: cannot read mesh from " << meshfilename << std::endl;
		return false;
	}
	std::cout << "[V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]" << std::endl;
	CutMesh();
	return true;
}

void VDPApp::PrintError(void)
{
	std::cerr << "ERROR: wrong argument." << std::endl;
	std::cout << "Usage: ~.exe <mesh1> [-a] [-g] [-p <int>] [-q <int>] [-r <int>][-s <int>]\n"
		"Parameters: mesh1: file name of a genus zero mesh.\n"
		"            -a: save intermediate results.\n"
		"            -g: use geodesic distance.\n"
		"            -p: vertex number of simplified mesh. default: 13000\n"
		"            -q: number of random cutting process. default: 10\n"
		"            -r: voting threshold. default: 2\n"
		"            -s: local feature size N. default: 13\n"
		"            -t: n-ring neighborhood. default: 5" << std::endl;
}

void VDPApp::CutMesh(void)
{
	CutGenusZeroMesh cut(mesh);
	cut.Init(size, saveintermediate, geodesic, nvsimplify, nrandom, nvoting, npost, meshpath + "/" + meshbasename);
	cut.Run();
}
