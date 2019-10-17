#pragma once
#include <vector>
extern "C"
{
	void pardisoinit(void* pt, int* mtype, int* solver, int* iparm, double* dparm, int*error);
	void pardiso(void* pt, int* maxfct, int* mnum, int* mtype, int* phase, int* n,
		double* a, int* ia, int* ja, int* perm, int* nrhs, int* iparm,
		int* msglvl, double* b, double* x, int* error, double* dparm);
	void pardiso_chkmatrix(int *mtype, int *n, double *a,
		int *ia, int *ja, int *error);
	void pardiso_chkvec(int *n, int *nrhs, double *b, int *error);
	void pardiso_printstats(int *mtype, int *n, double *a, int *ia,
		int *ja, int *nrhs, double *b, int *error);
}

// Class of the Pardiso solver v5.0
class PardisoSolver
{
public:
	PardisoSolver(void);
	~PardisoSolver(void);
	bool Init(int matrixtype);
	bool AnalyzePattern(void);
	bool Factorize(void);
	bool Solve(std::vector<double> &b, std::vector<double> &x);
	void PrintError(void);
public:
	void *pt[64];
	int mtype, solver;
	int iparm[64];
	double dparm[64];
	int error, maxfct, mnum, phase, nrow;
	std::vector<double> a;
	std::vector<int> ia, ja;
	int nrhs, msglvl;
	double ddum; // double dummy
	int idum; // integer dummy
};
