#include "PardisoSolver.h"
#include <iostream>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#else
#if defined(_WIN32)
#pragma comment(lib, "libpardiso600-WIN-X86-64.lib")
#endif
#endif

PardisoSolver::PardisoSolver(void)
{
}

PardisoSolver::~PardisoSolver(void)
{
	if (mtype == 0)
	{
		return;
	}
	phase = -1;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nrow,
		&ddum, ia.data(), ja.data(), &idum, &nrhs, iparm,
		&msglvl, &ddum, &ddum, &error
#ifndef USE_MKL_PARDISO
		, dparm
#endif
	);
}

bool PardisoSolver::Init(int matrixtype)
{
	mtype = matrixtype;
	switch (mtype)
	{
	default:
		std::cerr << "Error: Undefined matrix type." << std::endl;
		return false;
		break;
	case 1:
	case 2:
	case -2:
	case 11:
		break;
	case 3:
	case 4:
	case -4:
	case 6:
	case 13:
		std::cerr << "Error: Complex matrix not supported." << std::endl;
		return false;
		break;
	}

#ifdef USE_MKL_PARDISO
	pardisoinit(pt, &mtype, iparm);
#else
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);
#endif // USE_MKL_PARDISO

	if (error != 0)
	{
		PrintError();
		return false;
	}
#ifdef USE_MKL_PARDISO
	iparm[34] = 0;
#else
	iparm[2] = 8;
#endif
	return true;
}

bool PardisoSolver::AnalyzePattern(void)
{
	if (mtype == 0)
	{
		std::cerr << "Error: Pardiso not initialized." << std::endl;
		return false;
	}
	nrow = static_cast<int>(ia.size()) - 1;
	if (nrow <= 0)
	{
		std::cerr << "Error: Pardiso null matrix." << std::endl;
		return false;
	}
#ifndef USE_MKL_PARDISO
#ifdef _DEBUG
	pardiso_chkmatrix(&mtype, &nrow, a.data(), ia.data(), ja.data(), &error);
	if (error != 0)
	{
		PrintError();
		return false;
	}
#endif
#endif
	phase = 11; // Analysis
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nrow,
		a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm,
		&msglvl, &ddum, &ddum, &error
#ifndef USE_MKL_PARDISO
		, dparm
#endif
	);
	if (error != 0)
	{
		PrintError();
		return false;
	}
	return true;
}

bool PardisoSolver::Factorize(void)
{
	if (mtype == 0)
	{
		std::cerr << "Error: Pardiso not initialized." << std::endl;
		return false;
	}
	phase = 22; // Numerical factorization
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nrow,
		a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm,
		&msglvl, &ddum, &ddum, &error
#ifndef USE_MKL_PARDISO
		, dparm
#endif
	);
	if (error != 0)
	{
		PrintError();
		return false;
	}
	return true;
}

bool PardisoSolver::Solve(std::vector<double>& b, std::vector<double>& x)
{
	if (mtype == 0)
	{
		std::cerr << "Error: Pardiso not initialized." << std::endl;
		return false;
	}
	nrhs = static_cast<int>(b.size()) / nrow;
	if (nrhs <= 0)
	{
		std::cerr << "Error: Pardiso null right-hand-side." << std::endl;
		return false;
	}
#ifndef USE_MKL_PARDISO
#ifdef _DEBUG
	pardiso_chkvec(&nrow, &nrhs, b.data(), &error);
	if (error != 0)
	{
		PrintError();
		return false;
	}
	pardiso_printstats(&mtype, &nrow, a.data(), ia.data(), ja.data(), &nrhs, b.data(), &error);
	if (error != 0)
	{
		PrintError();
		return false;
	}
#endif
#endif // !USE_MKL_PARDISO
	x.resize(nrow);
	phase = 33;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nrow,
		a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm,
		&msglvl, b.data(), x.data(), &error
#ifndef USE_MKL_PARDISO
		, dparm
#endif
	);
	if (error != 0)
	{
		PrintError();
		return false;
	}
	return true;
}

std::vector<int>& PardisoSolver::RowIndex(void)
{
	return ia;
}

std::vector<int>& PardisoSolver::Columns(void)
{
	return ja;
}

std::vector<double>& PardisoSolver::MatrixValues(void)
{
	return a;
}

void PardisoSolver::PrintError(void)
{
	if (error != 0)
	{
		std::cerr << "Error (" << error << "): ";
		switch (error)
		{
		default:
			std::cerr << "Unknown error." << std::endl;
			break;
		case -1:
			std::cerr << "Input inconsistent." << std::endl;
			break;
		case -2:
			std::cerr << "Not enough memory." << std::endl;
			break;
		case -3:
			std::cerr << "Reordering problem." << std::endl;
			break;
		case -4:
			std::cerr << "Zero pivot, numerical fact, or iterative refinement problem." << std::endl;
			break;
		case -5:
			std::cerr << "Unclassified (internal) error." << std::endl;
			break;
		case -6:
			std::cerr << "Preordering failed (matrix type 11, 13 only)." << std::endl;
			break;
		case -7:
			std::cerr << "Diagonal matrix problem." << std::endl;
			break;
		case -8:
			std::cerr << "32-bit integer overflow problem." << std::endl;
			break;
#ifdef USE_MKL_PARDISO
		case -9:
			std::cerr << "Not enough memory for OOC." << std::endl;
			break;
		case -10:
			std::cerr << "Problems with opening OOC temporary files." << std::endl;
			break;
		case -11:
			std::cerr << "Read/write problems with the OOC data file." << std::endl;
			break;
#else
		case -10:
			std::cerr << "No license file pardiso.lic found." << std::endl;
			break;
		case -11:
			std::cerr << "License is expired." << std::endl;
			break;
		case -12:
			std::cerr << "Wrong username or hostname." << std::endl;
			break;
		case -100:
			std::cerr << "Reached maximum number of Krylov-subspace iteration in iterative solver." << std::endl;
			break;
		case -101:
			std::cerr << "No sufficient convergence in Krylov-subspace iteration within 25 iteration." << std::endl;
			break;
		case -102:
			std::cerr << "Error in Krylov-subspace iteration." << std::endl;
			break;
		case -103:
			std::cerr << "Break-Down in Krylov-subspace iteration." << std::endl;
			break;
#endif
		}
	}
}
