#include "KPNewton.h"
#include <fstream>
#include <iostream>
#include <unordered_map>

KPNewton::KPNewton(Mesh& m, bool verbose)
	:mesh(m),
	isbv(m.n_vertices(), 0),
	verbose(verbose)
{
	if (!solver.Init(-2))
	{
		std::cerr << "Pardiso initialization failed." << std::endl;
	}
	for (const auto& heh : mesh.halfedges())
	{
		if (mesh.is_boundary(heh))
		{
			auto nextheh = heh;
			do
			{
				isbv[mesh.from_vertex_handle(nextheh).idx()] = 1;
				nextheh = mesh.prev_halfedge_handle(nextheh);
			} while (nextheh != heh);
			break;
		}
	}
}

KPNewton::~KPNewton(void)
{
}

void KPNewton::PrepareData(void)
{
	auto&& ia = solver.RowIndex();
	auto&& ja = solver.Columns();
	int nv = static_cast<int>(mesh.n_vertices());
	ia.clear();
	ja.clear();
	ia.reserve(2 * nv + 1);
	ja.reserve(16 * nv);
	for (const auto& vh : mesh.vertices())
	{
		ia.push_back(static_cast<int>(ja.size()) + 1);
		auto vid = vh.idx();
		if (mesh.is_boundary(vh))
		{
			ja.push_back(vid + 1);
		}
		else
		{
			std::vector<int> rowid;
			rowid.push_back(vid + 1);
			rowid.push_back(vid + nv + 1);
			for (const auto& vvh : mesh.vv_range(vh))
			{
				int vvid = vvh.idx();
				if (!mesh.is_boundary(vvh))
				{
					if (vvid > vid)
					{
						rowid.push_back(vvid + 1);
					}
					rowid.push_back(vvid + nv + 1);
				}
			}
			std::sort(rowid.begin(), rowid.end(), std::less<int>());
			for (const auto& id : rowid)
			{
				ja.push_back(id);
			}
		}
	}
	for (const auto& vh : mesh.vertices())
	{
		ia.push_back(static_cast<int>(ja.size()) + 1);
		auto vid = vh.idx();
		if (mesh.is_boundary(vh))
		{
			ja.push_back(vid + nv + 1);
		}
		else
		{
			std::vector<int> rowid;
			rowid.push_back(vid + nv + 1);
			for (const auto& vvh : mesh.vv_range(vh))
			{
				int vvid = vvh.idx();
				if (!mesh.is_boundary(vvh))
				{
					if (vvid > vid)
					{
						rowid.push_back(vvid + nv + 1);
					}
				}
			}
			std::sort(rowid.begin(), rowid.end(), std::less<int>());
			for (const auto& id : rowid)
			{
				ja.push_back(id);
			}
		}
	}
	ia.push_back(static_cast<int>(ja.size()) + 1);
	solver.MatrixValues().resize(ja.size());
	//solver.nrow = static_cast<int>(2 * mesh.n_vertices());
	if (!solver.AnalyzePattern())
	{
		std::cerr << "Pardiso analyze pattern failed." << std::endl;
		return;
	}
	std::vector<std::unordered_map<int, int>> spid(2 * nv);
	size_t j = 0;
	for (size_t i = 0; i + 1 < ia.size(); i++)
	{
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			spid[i][ja[j - 1] - 1] = j - 1;
		}
	}
	tri.resize(mesh.n_faces());
	for (const auto& fh : mesh.faces())
	{
		auto fid = fh.idx();
		for (const auto& fvh : mesh.fv_range(fh))
		{
			tri[fid].push_back(fvh.idx());
		}
	}
	assembleorder.clear();
	for (size_t fid = 0; fid < tri.size(); fid++)
	{
		const auto& vhs = tri[fid];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (!(isbv[vhs[i]] || isbv[vhs[j]]))
				{
					if (vhs[i] <= vhs[j])
					{
						assembleorder.push_back(spid[vhs[i]][vhs[j]]);
						assembleorder.push_back(spid[vhs[i] + nv][vhs[j] + nv]);
					}
					assembleorder.push_back(spid[vhs[i]][vhs[j] + nv]);
				}
			}
		}
	}
	for (int i = 0; i < nv; i++)
	{
		if (isbv[i])
		{
			assembleorder.push_back(spid[i][i]);
			assembleorder.push_back(spid[i + nv][i + nv]);
		}
	}
}

void KPNewton::Run(const EnergyType& etype)
{
	auto maxiter = 100;
	double Iscale = 1e-8;
	double ls_StepFac = 0.5;
	double ls_GradFac = 0.02;
	double ls_MinStep = 1e-7;
	switch (etype)
	{
	default:
		break;
	case EnergyType::MIPS:
		computeall = ComputeMIPS;
		computeenergy = ComputeEnergyMIPS;
		break;
	case EnergyType::AMIPS:
		computeall = ComputeAMIPS;
		computeenergy = ComputeEnergyAMIPS;
		break;
	}
	if (result.empty())
	{
		std::cerr << "Error: No initial position." << std::endl;
		return;
	}
	auto nV = mesh.n_vertices();
	auto nF = mesh.n_faces();

	StdMatrixd faceElens = MeshFaceEdgeLen2s();
	StdMatrixd faceAngles = MeshAnglesFromFaceEdgeLen2(faceElens);
	StdMatrixcd localframe(nF, StdVectorcd(3));
	StdMatrixcd D(nF, StdVectorcd(3));
	StdMatrixcd DC(nF, StdVectorcd(3));
	StdVectord area(nF);
	double totalarea = 0;
	for (size_t i = 0; i < nF; i++)
	{
		auto&& frame = localframe[i];
		frame[1] = std::sqrt(faceElens[i][2]);
		frame[2] = std::polar(std::sqrt(faceElens[i][1]), faceAngles[i][0]);
		area[i] = frame[1].real() * frame[2].imag() / 2;
		totalarea += area[i];
		DC[i][0] = std::complex<double>(0.0, -0.25) * (frame[1] - frame[2]) / area[i];
		DC[i][1] = std::complex<double>(0.0, -0.25) * (frame[2] - frame[0]) / area[i];
		DC[i][2] = std::complex<double>(0.0, -0.25) * (frame[0] - frame[1]) / area[i];
		D[i][0] = std::conj(DC[i][0]);
		D[i][1] = std::conj(DC[i][1]);
		D[i][2] = std::conj(DC[i][2]);
	}
	for (auto&& a : area)
	{
		a /= totalarea;
	}
	for (int iter = 0; iter < maxiter; iter++)
	{
		auto&& a = solver.MatrixValues();
		a.clear();
		a.resize(solver.Columns().size(), 0);
		double e = 0;
		StdVectord g(2 * nV, 0);
		StdVectord b(2 * nV, 0);
		int nele = 0;
		for (size_t fid = 0; fid < tri.size(); fid++)
		{
			std::complex<double> fz, fzb;
			const auto& vhs = tri[fid];
			for (size_t i = 0; i < 3; i++)
			{
				fz += D[fid][i] * result[vhs[i]];
				fzb += DC[fid][i] * result[vhs[i]];
			}
			double x = std::norm(fz);
			double y = std::norm(fzb);
			double etemp, alpha1, alpha2, beta1, beta2, beta3;
			computeall(etemp, alpha1, alpha2, beta1, beta2, beta3, x, y);
			e += etemp * area[fid];
			for (size_t i = 0; i < 3; i++)
			{
				auto gi = 2.0 * (alpha1 * DC[fid][i] * fz + alpha2 * D[fid][i] * fzb) * area[fid];
				if (isbv[vhs[i]]) continue;
				g[vhs[i]] += gi.real();
				g[vhs[i] + nV] += gi.imag();
			}
			{// fix alpha and beta
				double temp1 = alpha1 + 2 * beta1 * x;
				double temp2 = alpha2 + 2 * beta2 * y;
				double s1 = temp1 + temp2;
				double s2 = temp1 - temp2;
				double lambda3 = s1 + std::sqrt(s2 * s2 + 16 * beta3 * beta3 * x * y);
				double lambda4 = s1 - std::sqrt(s2 * s2 + 16 * beta3 * beta3 * x * y);
				double t1 = (lambda3 - 2 * alpha1 - 4 * beta1 * x) / (4 * beta3 * y);
				double t2 = (lambda4 - 2 * alpha1 - 4 * beta1 * x) / (4 * beta3 * y);
				double e3norm2 = x + y * t1 * t1;
				double e4norm2 = x + y * t2 * t2;
				lambda3 = (lambda3 > 0 ? lambda3 : 0) / e3norm2;
				lambda4 = (lambda4 > 0 ? lambda4 : 0) / e4norm2;
				if (y > 1e-50 && fabs(beta1 * beta2 * beta3) > 1e-50)
				{
					alpha1 = alpha1 > 0 ? alpha1 : 0;
					alpha2 = alpha2 > 0 ? alpha2 : 0;
					beta1 = (lambda3 + lambda4 - alpha1 * 2.0 / x) / 4.0;
					beta2 = (lambda3 * t1 * t1 + lambda4 * t2 * t2 - alpha2 * 2.0 / y) / 4.0;
					beta3 = (lambda3 * t1 + lambda4 * t2) / 4.0;
				}
			}

			double a1 = area[fid] * (2.0 * alpha1 + 4.0 * beta1 * fz.real() * fz.real());
			double a2 = area[fid] * 4.0 * beta1 * fz.real() * fz.imag();
			double a4 = area[fid] * (2.0 * alpha1 + 4.0 * beta1 * fz.imag() * fz.imag());
			double b1 = area[fid] * 4.0 * beta3 * fzb.real() * fz.real();
			double b2 = area[fid] * 4.0 * beta3 * fzb.real() * fz.imag();
			double b3 = area[fid] * 4.0 * beta3 * fzb.imag() * fz.real();
			double b4 = area[fid] * 4.0 * beta3 * fzb.imag() * fz.imag();
			double d1 = area[fid] * (2.0 * alpha2 + 4.0 * beta2 * fzb.real() * fzb.real());
			double d2 = area[fid] * 4.0 * beta2 * fzb.real() * fzb.imag();
			double d4 = area[fid] * (2.0 * alpha2 + 4.0 * beta2 * fzb.imag() * fzb.imag());
			double s1 = a1 + 2.0 * b1 + d1;
			double s2 = a2 + b2 - b3 - d2;
			double s3 = a4 - 2.0 * b4 + d4;
			double s4 = a2 + b2 + b3 + d2;
			double s5 = d1 - a1;
			double s6 = a4 - d4;
			double s7 = -a2 + b2 + b3 - d2;
			double s8 = a4 + 2.0 * b4 + d4;
			double s9 = -a2 + b2 - b3 + d2;
			double s0 = a1 - 2.0 * b1 + d1;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					double RR = D[fid][i].real() * D[fid][j].real();
					double RI = D[fid][i].real() * D[fid][j].imag();
					double IR = D[fid][i].imag() * D[fid][j].real();
					double II = D[fid][i].imag() * D[fid][j].imag();

					if (isbv[vhs[i]] || isbv[vhs[j]])
					{
						if (!isbv[vhs[i]])
						{
							//b[vhs[i]] -= (s1 * RR + s2 * (RI + IR) + s3 * II) * result[vhs[j]].real()
							//	+ (s4 * RR + s5 * RI + s6 * IR + s7 * II) * result[vhs[j]].imag();
							//b[vhs[i] + nV] -= (s4 * RR + s6 * RI + s5 * IR + s7 * II) * result[vhs[j]].real()
							//	+ (s8 * RR + s9 * (RI + IR) + s0 * II) * result[vhs[j]].imag();
						}
					}
					else
					{
						if (vhs[i] <= vhs[j])
						{
							a[assembleorder[nele++]] += s1 * RR + s2 * (RI + IR) + s3 * II;
							a[assembleorder[nele++]] += s8 * RR + s9 * (RI + IR) + s0 * II;
						}
						a[assembleorder[nele++]] += s4 * RR + s5 * RI + s6 * IR + s7 * II;
					}
				}
			}
		}
		for (size_t i = 0; i < nV; i++)
		{
			if (isbv[i])
			{
				a[assembleorder[nele++]] = 1;
				a[assembleorder[nele++]] = 1;
			}
		}
		if (!solver.Factorize())
		{
			std::cerr << "Pardiso factorize failed." << std::endl;
			return;
		}
		double normg = 0.0;
		for (size_t i = 0; i < g.size(); i++)
		{
			b[i] -= g[i];
			normg += g[i] * g[i];
		}
		normg /= std::sqrt(normg);
		StdVectord x;
		if (!solver.Solve(b, x))
		{
			std::cerr << "Pardiso solve x failed." << std::endl;
			return;
		}
		StdVectorcd sDC;
		StdVectorcd newpos;
		double stepDirlen = 0.0;
		double energyGrad = 0.0;
		for (size_t i = 0; i < nV; i++)
		{
			sDC.push_back({ x[i], x[i + nV] });
			newpos.push_back(result[i] + sDC[i]);
			stepDirlen += x[i] * x[i] + x[i + nV] * x[i + nV];
			energyGrad -= b[i] * x[i] + b[i + nV] * x[i + nV];
		}
		energyGrad *= ls_GradFac;
		stepDirlen = std::sqrt(stepDirlen);
		// Compute positive area
		double lsa = ComputeTMax(result, sDC) * 0.9;
		lsa = lsa > 1 ? 1 : lsa;
		//std::cout << lsa << ' ' << stepDirlen << std::endl;
		for (size_t i = 0; i < nV; i++)
		{
			newpos[i] = result[i] + sDC[i] * lsa;
		}
		double newEnergy = ComputeEnergy(newpos, D, DC, area);

		while (stepDirlen * lsa > ls_MinStep && e + lsa * energyGrad < newEnergy)
		{
			lsa = ls_StepFac * lsa;
			for (size_t i = 0; i < nV; i++)
			{
				newpos[i] = result[i] + sDC[i] * lsa;
			}
			newEnergy = ComputeEnergy(newpos, D, DC, area);
		}
		if (verbose)
		{
			std::cout.setf(std::ios::fixed);
			std::cout << iter << ", step: " << lsa << ", \txnorm: " << stepDirlen << ", \tgrad: " << energyGrad << ", \te: " << e << ", \tnewe: " << newEnergy << "\n";
			std::cout.unsetf(std::ios::fixed);
		}
		
		if (newEnergy > e || std::isnan(newEnergy) || std::isinf(newEnergy))
		{
			//std::cout << "Error!" << std::endl;
		}
		else
		{
			for (size_t i = 0; i < nV; i++)
			{
				result[i] = newpos[i];
			}
		}
		if (normg < 1e-4 || (e - newEnergy) / e < 1e-4 || std::isnan(newEnergy) || std::isinf(newEnergy))
		{
			break;
		}
	}
}

void KPNewton::PrepareDataFree(void)
{
	auto&& ia = solver.RowIndex();
	auto&& ja = solver.Columns();
	int nv = static_cast<int>(mesh.n_vertices());
	ia.clear();
	ja.clear();
	ia.reserve(2 * nv + 1);
	ja.reserve(16 * nv);
	for (const auto& vh : mesh.vertices())
	{
		ia.push_back(static_cast<int>(ja.size()) + 1);
		auto vid = vh.idx();
		std::vector<int> rowid;
		rowid.push_back(vid + 1);
		rowid.push_back(vid + nv + 1);
		for (const auto& vvh : mesh.vv_range(vh))
		{
			int vvid = vvh.idx();
			if (vvid > vid)
			{
				rowid.push_back(vvid + 1);
			}
			rowid.push_back(vvid + nv + 1);
		}
		std::sort(rowid.begin(), rowid.end(), std::less<int>());
		for (const auto& id : rowid)
		{
			ja.push_back(id);
		}
	}
	for (const auto& vh : mesh.vertices())
	{
		ia.push_back(static_cast<int>(ja.size()) + 1);
		auto vid = vh.idx();
		std::vector<int> rowid;
		rowid.push_back(vid + nv + 1);
		for (const auto& vvh : mesh.vv_range(vh))
		{
			int vvid = vvh.idx();
			if (vvid > vid)
			{
				rowid.push_back(vvid + nv + 1);
			}
		}
		std::sort(rowid.begin(), rowid.end(), std::less<int>());
		for (const auto& id : rowid)
		{
			ja.push_back(id);
		}
	}
	ia.push_back(static_cast<int>(ja.size()) + 1);
	solver.MatrixValues().resize(ja.size());
	//solver.nrow = static_cast<int>(2 * mesh.n_vertices());
	if (!solver.AnalyzePattern())
	{
		std::cerr << "Pardiso analyze pattern failed." << std::endl;
		return;
	}
	std::vector<std::unordered_map<int, int>> spid(2 * nv);
	size_t j = 0;
	for (size_t i = 0; i + 1 < ia.size(); i++)
	{
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			spid[i][ja[j - 1] - 1] = j - 1;
		}
	}
	tri.resize(mesh.n_faces());
	for (const auto& fh : mesh.faces())
	{
		auto fid = fh.idx();
		for (const auto& fvh : mesh.fv_range(fh))
		{
			tri[fid].push_back(fvh.idx());
		}
	}
	assembleorder.clear();
	for (size_t fid = 0; fid < tri.size(); fid++)
	{
		const auto& vhs = tri[fid];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (vhs[i] <= vhs[j])
				{
					assembleorder.push_back(spid[vhs[i]][vhs[j]]);
					assembleorder.push_back(spid[vhs[i] + nv][vhs[j] + nv]);
				}
				assembleorder.push_back(spid[vhs[i]][vhs[j] + nv]);
			}
		}
	}
}

void KPNewton::RunFree(const EnergyType& etype)
{
	auto maxiter = 100;
	double Iscale = 1e-8;
	double ls_StepFac = 0.5;
	double ls_GradFac = 0.02;
	double ls_MinStep = 1e-7;
	switch (etype)
	{
	default:
		break;
	case EnergyType::MIPS:
		computeall = ComputeMIPS;
		computeenergy = ComputeEnergyMIPS;
		break;
	case EnergyType::AMIPS:
		computeall = ComputeAMIPS;
		computeenergy = ComputeEnergyAMIPS;
		break;
	}
	if (result.empty())
	{
		std::cerr << "Error: No initial position." << std::endl;
		return;
	}
	// prepare inputs
	auto nV = mesh.n_vertices();
	auto nF = mesh.n_faces();

	StdMatrixd faceElens = MeshFaceEdgeLen2s();
	StdMatrixd faceAngles = MeshAnglesFromFaceEdgeLen2(faceElens);
	StdMatrixcd localframe(nF, StdVectorcd(3));
	StdMatrixcd D(nF, StdVectorcd(3));
	StdMatrixcd DC(nF, StdVectorcd(3));
	StdVectord area(nF);
	double totalarea = 0;
	for (size_t i = 0; i < nF; i++)
	{
		auto&& frame = localframe[i];
		frame[1] = std::sqrt(faceElens[i][2]);
		frame[2] = std::polar(std::sqrt(faceElens[i][1]), faceAngles[i][0]);
		area[i] = frame[1].real() * frame[2].imag() / 2;
		totalarea += area[i];
		DC[i][0] = std::complex<double>(0.0, -0.25) * (frame[1] - frame[2]) / area[i];
		DC[i][1] = std::complex<double>(0.0, -0.25) * (frame[2] - frame[0]) / area[i];
		DC[i][2] = std::complex<double>(0.0, -0.25) * (frame[0] - frame[1]) / area[i];
		D[i][0] = std::conj(DC[i][0]);
		D[i][1] = std::conj(DC[i][1]);
		D[i][2] = std::conj(DC[i][2]);
	}
	for (auto&& a : area)
	{
		a /= totalarea;
	}
	for (int iter = 0; iter < maxiter; iter++)
	{
		auto&& a = solver.MatrixValues();
		a.clear();
		a.resize(solver.Columns().size(), 0);
		double e = 0;
		StdVectord g(2 * nV, 0);
		StdVectord b(2 * nV, 0);
		int nele = 0;
		for (size_t fid = 0; fid < tri.size(); fid++)
		{
			std::complex<double> fz, fzb;
			const auto& vhs = tri[fid];
			for (size_t i = 0; i < 3; i++)
			{
				fz += D[fid][i] * result[vhs[i]];
				fzb += DC[fid][i] * result[vhs[i]];
			}
			double x = std::norm(fz);
			double y = std::norm(fzb);
			double etemp, alpha1, alpha2, beta1, beta2, beta3;
			computeall(etemp, alpha1, alpha2, beta1, beta2, beta3, x, y);
			e += etemp * area[fid];
			for (size_t i = 0; i < 3; i++)
			{
				auto gi = 2.0 * (alpha1 * DC[fid][i] * fz + alpha2 * D[fid][i] * fzb) * area[fid];
				g[vhs[i]] += gi.real();
				g[vhs[i] + nV] += gi.imag();
			}
			{// fix alpha and beta
				double temp1 = alpha1 + 2 * beta1 * x;
				double temp2 = alpha2 + 2 * beta2 * y;
				double s1 = temp1 + temp2;
				double s2 = temp1 - temp2;
				double lambda3 = s1 + std::sqrt(s2 * s2 + 16 * beta3 * beta3 * x * y);
				double lambda4 = s1 - std::sqrt(s2 * s2 + 16 * beta3 * beta3 * x * y);
				double t1 = (lambda3 - 2 * alpha1 - 4 * beta1 * x) / (4 * beta3 * y);
				double t2 = (lambda4 - 2 * alpha1 - 4 * beta1 * x) / (4 * beta3 * y);
				double e3norm2 = x + y * t1 * t1;
				double e4norm2 = x + y * t2 * t2;
				lambda3 = (lambda3 > 0 ? lambda3 : 0) / e3norm2;
				lambda4 = (lambda4 > 0 ? lambda4 : 0) / e4norm2;
				if (y > 1e-50 && fabs(beta1 * beta2 * beta3) > 1e-50)
				{
					alpha1 = alpha1 > 0 ? alpha1 : 0;
					alpha2 = alpha2 > 0 ? alpha2 : 0;
					beta1 = (lambda3 + lambda4 - alpha1 * 2.0 / x) / 4.0;
					beta2 = (lambda3 * t1 * t1 + lambda4 * t2 * t2 - alpha2 * 2.0 / y) / 4.0;
					beta3 = (lambda3 * t1 + lambda4 * t2) / 4.0;
				}
			}

			double a1 = area[fid] * (2.0 * alpha1 + 4.0 * beta1 * fz.real() * fz.real());
			double a2 = area[fid] * 4.0 * beta1 * fz.real() * fz.imag();
			double a4 = area[fid] * (2.0 * alpha1 + 4.0 * beta1 * fz.imag() * fz.imag());
			double b1 = area[fid] * 4.0 * beta3 * fzb.real() * fz.real();
			double b2 = area[fid] * 4.0 * beta3 * fzb.real() * fz.imag();
			double b3 = area[fid] * 4.0 * beta3 * fzb.imag() * fz.real();
			double b4 = area[fid] * 4.0 * beta3 * fzb.imag() * fz.imag();
			double d1 = area[fid] * (2.0 * alpha2 + 4.0 * beta2 * fzb.real() * fzb.real());
			double d2 = area[fid] * 4.0 * beta2 * fzb.real() * fzb.imag();
			double d4 = area[fid] * (2.0 * alpha2 + 4.0 * beta2 * fzb.imag() * fzb.imag());
			double s1 = a1 + 2.0 * b1 + d1;
			double s2 = a2 + b2 - b3 - d2;
			double s3 = a4 - 2.0 * b4 + d4;
			double s4 = a2 + b2 + b3 + d2;
			double s5 = d1 - a1;
			double s6 = a4 - d4;
			double s7 = -a2 + b2 + b3 - d2;
			double s8 = a4 + 2.0 * b4 + d4;
			double s9 = -a2 + b2 - b3 + d2;
			double s0 = a1 - 2.0 * b1 + d1;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					double RR = D[fid][i].real() * D[fid][j].real();
					double RI = D[fid][i].real() * D[fid][j].imag();
					double IR = D[fid][i].imag() * D[fid][j].real();
					double II = D[fid][i].imag() * D[fid][j].imag();

					if (vhs[i] <= vhs[j])
					{
						a[assembleorder[nele++]] += s1 * RR + s2 * (RI + IR) + s3 * II;
						a[assembleorder[nele++]] += s8 * RR + s9 * (RI + IR) + s0 * II;
					}
					a[assembleorder[nele++]] += s4 * RR + s5 * RI + s6 * IR + s7 * II;
				}
			}
		}
		if (!solver.Factorize())
		{
			std::cerr << "Pardiso factorize failed." << std::endl;
			return;
		}
		double normg = 0.0;
		for (size_t i = 0; i < g.size(); i++)
		{
			b[i] -= g[i];
			normg += g[i] * g[i];
		}
		normg /= std::sqrt(normg);
		StdVectord x;
		if (!solver.Solve(b, x))
		{
			std::cerr << "Pardiso solve x failed." << std::endl;
			return;
		}
		StdVectorcd sDC;
		StdVectorcd newpos;
		double stepDirlen = 0.0;
		double energyGrad = 0.0;
		for (size_t i = 0; i < nV; i++)
		{
			sDC.push_back({ x[i], x[i + nV] });
			newpos.push_back(result[i] + sDC[i]);
			stepDirlen += x[i] * x[i] + x[i + nV] * x[i + nV];
			energyGrad -= b[i] * x[i] + b[i + nV] * x[i + nV];
		}
		energyGrad *= ls_GradFac;
		stepDirlen = std::sqrt(stepDirlen);
		// Compute positive area
		double lsa = ComputeTMax(result, sDC) * 0.9;
		lsa = lsa > 1 ? 1 : lsa;
		//std::cout << lsa << ' ' << stepDirlen << std::endl;
		for (size_t i = 0; i < nV; i++)
		{
			newpos[i] = result[i] + sDC[i] * lsa;
		}
		double newEnergy = ComputeEnergy(newpos, D, DC, area);

		while (stepDirlen * lsa > ls_MinStep && e + lsa * energyGrad < newEnergy)
		{
			lsa = ls_StepFac * lsa;
			for (size_t i = 0; i < nV; i++)
			{
				newpos[i] = result[i] + sDC[i] * lsa;
			}
			newEnergy = ComputeEnergy(newpos, D, DC, area);
		}
		if (verbose)
		{
			std::cout.setf(std::ios::fixed);
			std::cout << iter << ", step: " << lsa << ", \txnorm: " << stepDirlen << ", \tgrad: " << energyGrad << ", \te: " << e << ", \tnewe: " << newEnergy << "\n";
			std::cout.unsetf(std::ios::fixed);
		}
		if (newEnergy > e || std::isnan(newEnergy) || std::isinf(newEnergy))
		{
			//std::cout << "Error!" << std::endl;
		}
		else
		{
			for (size_t i = 0; i < nV; i++)
			{
				result[i] = newpos[i];
			}
		}
		if (normg < 1e-4 || (e - newEnergy) / e < 1e-4 || std::isnan(newEnergy) || std::isinf(newEnergy))
		{
			break;
		}
	}
}

void KPNewton::Tutte(void)
{
	if (verbose)
	{
		std::cout << "Tutte start ... ";
	}
	auto boundaryvhs = GetBoundary();

	double delta_angle = 2 * M_PI / boundaryvhs.size();
	double area_1_factor = 0.5 * M_2_SQRTPI;
	auto nv = mesh.n_vertices();
	std::vector<double> posx(nv), posy(nv);
	for (size_t i = 0; i < boundaryvhs.size(); ++i)
	{
		posx[boundaryvhs[i].idx()] = area_1_factor * cos(i * delta_angle);
		posy[boundaryvhs[i].idx()] = area_1_factor * sin(i * delta_angle);
	}

	auto&& pardiso_it = solver.RowIndex();
	auto&& pardiso_jt = solver.Columns();
	auto&& pardiso_t = solver.MatrixValues();
	std::vector<double> pardiso_tu;
	std::vector<double> pardiso_tv;

	pardiso_it.reserve(nv + 1);
	pardiso_jt.reserve(6 * nv);
	pardiso_t.reserve(6 * nv);
	pardiso_tu.resize(nv, 0.0);
	pardiso_tv.resize(nv, 0.0);
	for (const auto& vh : mesh.vertices())
	{
		pardiso_it.push_back(static_cast<int>(pardiso_jt.size()) + 1);
		auto vid = vh.idx();
		if (mesh.is_boundary(vh))
		{
			pardiso_jt.push_back(vid + 1);
			pardiso_t.push_back(1);
			pardiso_tu[vh.idx()] = posx[vid];
			pardiso_tv[vh.idx()] = posy[vid];
		}
		else
		{
			pardiso_jt.push_back(vid + 1);
			pardiso_t.push_back(mesh.valence(vh));
			std::vector<int> row_id;
			row_id.reserve(mesh.valence(vh));
			for (const auto& vvh : mesh.vv_range(vh))
			{
				int vvid = vvh.idx();
				if (mesh.is_boundary(vvh))
				{
					pardiso_tu[vid] += posx[vvid];
					pardiso_tv[vid] += posy[vvid];
				}
				else
				{
					if (vvid > vid)
					{
						row_id.push_back(vvid);
					}
				}
			}
			std::sort(row_id.begin(), row_id.end(), std::less<int>());
			for (size_t j = 0; j < row_id.size(); j++)
			{
				pardiso_jt.push_back(row_id[j] + 1);
				pardiso_t.push_back(-1);
			}
		}
	}
	pardiso_it.push_back(static_cast<int>(pardiso_jt.size()) + 1);
	//solver.nrow = static_cast<int>(nv);
	if (!solver.AnalyzePattern())
	{
		std::cerr << "Pardiso analyze pattern ... failed." << std::endl;
		return;
	}
	if (!solver.Factorize())
	{
		std::cerr << "Pardiso factorize ... failed." << std::endl;
		return;
	}
	if (!solver.Solve(pardiso_tu, posx))
	{
		std::cerr << "Pardiso solve posx ... failed." << std::endl;
		return;
	}
	if (!solver.Solve(pardiso_tv, posy))
	{
		std::cerr << "Pardiso solve posy ... failed." << std::endl;
		return;
	}
	result.clear();
	result.reserve(nv);
	for (size_t i = 0; i < nv; i++)
	{
		result.push_back({ posx[i], posy[i] });
	}
	//UpdateMesh();
	if (verbose)
	{
		std::cout << "finished!" << std::endl;
	}
}

void KPNewton::LoadInitial(const std::string& filename)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open())
	{
		std::cerr << "Error: cannot open file " << filename << std::endl;
		return;
	}
	auto vhs = GetBoundary();
	std::string line;
	std::string keyWrd;
	Mesh::Point::value_type x, y, z;
	std::stringstream stream;
	result.clear();
	result.reserve(mesh.n_vertices());
	while (ifs && !ifs.eof())
	{
		std::getline(ifs, line);
		if (ifs.bad())
		{
			std::cerr << "  Warning! Could not read file properly!\n";
			return;
		}
		auto start = line.find_first_not_of(" \t\r\n");
		auto end = line.find_last_not_of(" \t\r\n");
		if ((std::string::npos == start) || (std::string::npos == end))
			line = "";
		else
			line = line.substr(start, end - start + 1);

		if (line.size() == 0 || line[0] == '#' || isspace(line[0]))
		{
			continue;
		}
		stream.str(line);
		stream.clear();
		stream >> keyWrd;
		if (keyWrd == "v")
		{
			stream >> x >> y >> z;
			if (!stream.fail())
			{
				result.push_back({ x, y });
			}
		}
	}
	ifs.close();
}

std::vector<Mesh::VertexHandle> KPNewton::GetBoundary(void) const
{
	std::vector<Mesh::VertexHandle> vhs;
	for (const auto& heh : mesh.halfedges())
	{
		if (mesh.is_boundary(heh))
		{
			auto nextheh = heh;
			do
			{
				vhs.push_back(mesh.from_vertex_handle(nextheh));
				nextheh = mesh.prev_halfedge_handle(nextheh);
			} while (nextheh != heh);
			break;
		}
	}
	return vhs;
}

KPNewton::StdMatrixd KPNewton::MeshFaceEdgeLen2s(void) const
{
	StdMatrixd len;
	len.reserve(mesh.n_faces());
	for (const auto& fh : mesh.faces())
	{
		auto fvit = mesh.fv_iter(fh);
		const auto& p0 = mesh.point(*fvit);
		fvit++;
		const auto& p1 = mesh.point(*fvit);
		fvit++;
		const auto& p2 = mesh.point(*fvit);
		len.push_back({ (p2 - p1).sqrnorm(), (p0 - p2).sqrnorm(), (p1 - p0).sqrnorm() });
	}
	return len;
}

KPNewton::StdMatrixd KPNewton::MeshAnglesFromFaceEdgeLen2(const StdMatrixd& len2) const
{
	StdMatrixd ang;
	ang.reserve(len2.size());
	for (const auto& fel : len2) // cos(a) = (b^2 + c^2 - a^2) / sqrt(bc) / 2
	{
		ang.push_back({
			std::acos((fel[1] + fel[2] - fel[0]) / std::sqrt(fel[1] * fel[2]) / 2.0),
			std::acos((fel[2] + fel[0] - fel[1]) / std::sqrt(fel[2] * fel[0]) / 2.0),
			std::acos((fel[0] + fel[1] - fel[2]) / std::sqrt(fel[0] * fel[1]) / 2.0)
			});
	}
	return ang;
}

double KPNewton::ComputeEnergy(const StdVectorcd& pos, const StdMatrixcd& D, const StdMatrixcd& DC, const StdVectord& area) const
{
	double energy = 0;
	for (size_t fid = 0; fid < tri.size(); fid++)
	{
		std::complex<double> fz = 0;
		std::complex<double> fzb = 0;
		for (size_t i = 0; i < 3; i++)
		{
			fz += D[fid][i] * pos[tri[fid][i]];
			fzb += DC[fid][i] * pos[tri[fid][i]];
		}
		double x = std::norm(fz);
		double y = std::norm(fzb);
		if (x < y) return std::numeric_limits<double>::infinity();
		energy += computeenergy(x, y) * area[fid];
	}
	return energy;
}

double KPNewton::ComputeTMax(const StdVectorcd& x, const StdVectorcd& d) const
{
	double temp_t = std::numeric_limits<double>::infinity();
	double a, b, c, b1, b2, tt, tt1, tt2;
	int V_N = (int)mesh.n_vertices();
	for (const auto& f : tri)
	{
		const auto& f0 = f[0];
		const auto& f1 = f[1];
		const auto& f2 = f[2];
		const auto& x0 = x[f0].real();
		const auto& x1 = x[f1].real();
		const auto& x2 = x[f2].real();
		const auto& x3 = x[f0].imag();
		const auto& x4 = x[f1].imag();
		const auto& x5 = x[f2].imag();
		const auto& d0 = d[f0].real();
		const auto& d1 = d[f1].real();
		const auto& d2 = d[f2].real();
		const auto& d3 = d[f0].imag();
		const auto& d4 = d[f1].imag();
		const auto& d5 = d[f2].imag();

		a = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1 = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2 = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b = b1 - b2;
		c = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt = std::numeric_limits<double>::infinity();
		//tt = 10000;
		if (b * b - 4 * a * c >= 0)
		{
			tt1 = 1 / (2 * a) * (-b + sqrt(b * b - 4 * a * c));
			tt2 = 1 / (2 * a) * (-b - sqrt(b * b - 4 * a * c));
			if (tt1 > 0 && tt2 > 0)
			{
				tt = std::min(tt1, tt2);
			}
			if (tt1 > 0 && tt2 < 0)
			{
				tt = tt1;
			}
			if (tt1 < 0 && tt2 > 0)
			{
				tt = tt2;
			}
		}
		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}
	return temp_t;
}

void KPNewton::UpdateMesh(void)
{
	for (const auto& vh : mesh.vertices())
	{
		mesh.set_point(vh, { result[vh.idx()].real(), result[vh.idx()].imag(), 0.0 });
	}
}

void KPNewton::ComputeMIPS(double& energy, double& alpha1, double& alpha2, double& beta1, double& beta2, double& beta3, const double& x, const double& y)
{
	energy = (x + y) / (x - y);
	alpha1 = -2 * y / (x - y) / (x - y);
	alpha2 = 2 * x / (x - y) / (x - y);
	beta1 = -2 * alpha1 / (x - y);
	beta2 = 2 * alpha2 / (x - y);
	beta3 = (alpha1 - alpha2) / (x - y);
}

void KPNewton::ComputeAMIPS(double& energy, double& alpha1, double& alpha2, double& beta1, double& beta2, double& beta3, const double& x, const double& y)
{
	energy = std::exp((x + y) / (x - y));
	alpha1 = -2 * y / (x - y) / (x - y) * energy;
	alpha2 = 2 * x / (x - y) / (x - y) * energy;
	beta1 = -2 * x * alpha1 / (x - y) / (x - y);
	beta2 = (4 * x - 2 * y) * alpha2 / (x - y) / (x - y);
	beta3 = ((2 * x - y) * alpha1 - x * alpha2) / (x - y) / (x - y);
}

double KPNewton::ComputeEnergyMIPS(const double& x, const double& y)
{
	return (x + y) / (x - y);
}

double KPNewton::ComputeEnergyAMIPS(const double& x, const double& y)
{
	return std::exp((x + y) / (x - y));
}
