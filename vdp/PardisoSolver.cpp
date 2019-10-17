#include "PardisoSolver.h"
#include <iostream>

#if defined(_WIN32)
#pragma comment(lib, "libpardiso600-WIN-X86-64.lib")
#endif

PardisoSolver::PardisoSolver(void)
	:mtype(0),
	solver(0),
	error(0),
	maxfct(1),
	mnum(1),
	nrhs(1),
	msglvl(0)
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
		&msglvl, &ddum, &ddum, &error, dparm);
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
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);
	if (error != 0)
	{
		PrintError();
		return false;
	}
	iparm[2] = 8;
	return true;
}

bool PardisoSolver::AnalyzePattern(void)
{
	if (mtype == 0)
	{
		std::cerr << "Error: Pardiso not initialized." << std::endl;
		return false;
	}
#ifdef _DEBUG
	pardiso_chkmatrix(&mtype, &nrow, a.data(), ia.data(), ja.data(), &error);
	if (error != 0)
	{
		PrintError();
		return false;
	}
#endif
	phase = 11; // Analysis
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nrow,
		a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm,
		&msglvl, &ddum, &ddum, &error, dparm);
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
		&msglvl, &ddum, &ddum, &error, dparm);
	if (error != 0)
	{
		PrintError();
		return false;
	}
	return true;
}

bool PardisoSolver::Solve(std::vector<double> &b, std::vector<double> &x)
{
	if (mtype == 0)
	{
		std::cerr << "Error: Pardiso not initialized." << std::endl;
		return false;
	}
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
	x.resize(nrow);
	phase = 33;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nrow,
		a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm,
		&msglvl, b.data(), x.data(), &error, dparm);
	if (error != 0)
	{
		PrintError();
		return false;
	}
	return true;
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
		}
	}
}
