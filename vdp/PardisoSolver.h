#pragma once
#include <vector>

#define USE_MKL_PARDISO

#ifndef USE_MKL_PARDISO
extern "C"
{
	void pardisoinit(void* pt, int* mtype, int* solver, int* iparm, double* dparm, int* error);
	void pardiso(void* pt, int* maxfct, int* mnum, int* mtype, int* phase, int* n,
		double* a, int* ia, int* ja, int* perm, int* nrhs, int* iparm,
		int* msglvl, double* b, double* x, int* error, double* dparm);
	void pardiso_chkmatrix(int* mtype, int* n, double* a,
		int* ia, int* ja, int* error);
	void pardiso_chkvec(int* n, int* nrhs, double* b, int* error);
	void pardiso_printstats(int* mtype, int* n, double* a, int* ia,
		int* ja, int* nrhs, double* b, int* error);
}
#endif

// Class of the Pardiso solver
class PardisoSolver
{
public:
	PardisoSolver(void);
	~PardisoSolver(void);
	bool Init(int matrixtype);
	bool AnalyzePattern(void);
	bool Factorize(void);
	bool Solve(std::vector<double>& b, std::vector<double>& x);
	std::vector<int>& RowIndex(void);
	std::vector<int>& Columns(void);
	std::vector<double>& MatrixValues(void);

private:
	void PrintError(void);

	void* pt[64] = { 0 };  // Solver internal data address pointer
	int mtype = -2;  // Matrix type
	int iparm[64] = { 0 };       // Parameters
#ifndef USE_MKL_PARDISO
	double dparm[64] = { 0.0 };  // Parameters
	int solver = 0;
#endif
	int error = 0;   // Error indicator
	int maxfct = 1;  // Maximal number of factors in memory
	int mnum = 1;    // The number of matrix (from 1 to maxfct) to solve
	int phase = 0;   // Controls the execution of the solver
	int nrow = 0;    // Number of equations in the sparse linear system A*X = B
	std::vector<double> a; // Contains the non-zero elements of the coefficient matrix A
	std::vector<int> ia;   // rowIndex array in CSR3 format
	std::vector<int> ja;   // columns array in CSR3 format
	int nrhs = 1;      // Number of right-hand sides that need to be solved for
	int msglvl = 0;    // Message level information (1 for output)
	double ddum = 0.0; // Double dummy
	int idum = 0;      // Integer dummy
};
