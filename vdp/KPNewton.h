#pragma once
#include "MeshDefinition.h"
#include "PardisoSolver.h"
#include <complex>

// Class of KP Newton method
class KPNewton
{
public:
	// Energy type for computing parameterization
	enum EnergyType
	{
		MIPS, AMIPS
	};
private:
	typedef std::vector<int> StdVectori;
	typedef std::vector<StdVectori> StdMatrixi;
	typedef std::vector<double> StdVectord;
	typedef std::vector<StdVectord> StdMatrixd;
	typedef std::vector< std::complex<double> > StdVectorcd;
	typedef std::vector<StdVectorcd> StdMatrixcd;
public:
	// Construct from a mesh
	KPNewton(Mesh & m);
	
	// Destructor
	virtual ~KPNewton(void);

	// Prepare data for fixed boundary parameterization
	void PrepareData(void);

	// Optimize parameterization using a type of energy
	void Run(const EnergyType & etype);

	// Prepare data for free boundary parameterization
	void PrepareDataFree(void);

	// Optimize parameterization using a type of energy
	void RunFree(const EnergyType & etype);

	// Initialize parameterization by Tutte parameterization
	void Tutte(void);

	// Load initial parameterization from a file
	void LoadInitial(const std::string & filename);

	// Update mesh vertices
	void UpdateMesh(void);

private:
	// Get boundary vertices
	std::vector<Mesh::VertexHandle> GetBoundary(void) const;

	// Compute edge length of each face
	StdMatrixd MeshFaceEdgeLen2s(void) const;

	// Compute angles in each triangles
	StdMatrixd MeshAnglesFromFaceEdgeLen2(const StdMatrixd & len2) const;

	// Compute the energy once
	double ComputeEnergy(const StdVectorcd & pos, const StdMatrixcd & D,
		const StdMatrixcd & DC, const StdVectord & area) const;

	// Compute maximum step size in line search
	double ComputeTMax(const StdVectorcd & x, const StdVectorcd & d) const;

	// Compute parameters using MIPS
	static void ComputeMIPS(double & energy, double & alpha1, double & alpha2,
		double & beta1, double & beta2, double & beta3,
		const double & x, const double & y);

	// Compute parameters using AMIPS
	static void ComputeAMIPS(double & energy, double & alpha1, double & alpha2,
		double & beta1, double & beta2, double & beta3,
		const double & x, const double & y);

	// Compute MIPS energy
	static double ComputeEnergyMIPS(const double & x, const double & y);

	// Compute AMIPS energy
	static double ComputeEnergyAMIPS(const double & x, const double & y);
	
private:
	Mesh & mesh; // referece to the mesh
	StdVectori isbv; // is boundary vertices
	StdVectorcd result; // result vertex positions (complex)
	StdMatrixi tri; // indices of 3 vertices in each triangle
	PardisoSolver solver; // Parsido solver
	StdVectori assembleorder; // the order of assemble matrix
	decltype(&ComputeMIPS) computeall; // parameter functional
	decltype(&ComputeEnergyMIPS) computeenergy; // energy functional
};

