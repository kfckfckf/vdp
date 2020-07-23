#pragma once
#include "MeshDefinition.h"

class ParameterManager;

// Class of voting distortion point application
class VDPApp
{
private:
	VDPApp();
	~VDPApp();

public:
	VDPApp(const VDPApp&) = delete;
	VDPApp(VDPApp&&) = delete;
	VDPApp& operator=(const VDPApp&) = delete;
	VDPApp& operator=(VDPApp&&) = delete;

	// Get application instance
	static VDPApp& GetInstance();

	// Parse command line
	bool ParseCommand(int argc, char* argv[]);

	// Run the application
	bool Run();

private:
	// Print error information
	void PrintError() const;

	// Cut the mesh
	void CutMesh();

	std::unique_ptr<ParameterManager> pm; // parameters

	std::unique_ptr<Mesh> mesh; // input mesh
	std::string meshfilename;   // input mesh file name
	std::string meshbasename;   // base name of the mesh file
};
