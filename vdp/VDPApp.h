#pragma once
#include "MeshDefinition.h"

// Class of voting distortion point application
class VDPApp
{
private:
	VDPApp(void);
	~VDPApp(void);
public:
	VDPApp(VDPApp const&) = delete;
	VDPApp(VDPApp &&) = delete;
	VDPApp& operator=(VDPApp const&) = delete;
	VDPApp& operator=(VDPApp &&) = delete;

	// Get application instance
	static VDPApp& GetInstance(void);

	// Parse command line
	bool ParseCommand(int argc, char *argv[]);

	// Run the application
	bool Run(void);

private:
	// Print error information
	void PrintError(void);

	// Cut the mesh
	void CutMesh(void);

	Mesh mesh; // input mesh
	std::string meshfilename; // input mesh file name
	std::string meshpath;     // path to the mesh file
	std::string meshbasename; // base name of the mesh file

	bool saveintermediate = false; // parameter: save intermediate results ?
	int size = 13;                 // parameter: local feature size N
	bool geodesic = false;         // parameter: use geodesic or not
	int nvsimplify = 13000;        // parameter: vertex number of simplified mesh
	int nrandom = 10;              // parameter: number of random cutting process
	int nvoting = 2;               // parameter: voting threshold
	int npost = 5;                 // parameter: n-ring of post-processing
};
