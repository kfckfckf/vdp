#include "QEMSimplification.h"
#include <iostream>

QEMSimplification::QEMSimplification(const Mesh& m)
	:mesh(m)
{
	mesh.add_property(vindex);
	for (const auto& vh : mesh.vertices())
	{
		mesh.property(vindex, vh) = vh.idx();
	}
	mesh.add_property(QuaMat_f);
	mesh.add_property(QuaMat_v);
	mesh.add_property(QuaMat_e);
	mesh.add_property(valence);
	mesh.add_property(area_f);
	mesh.add_property(new_point);
	mesh.add_property(is_cal_e);
	for (const auto& eh : mesh.edges())
	{
		mesh.property(is_cal_e, eh) = false;
	}
}

QEMSimplification::~QEMSimplification(void)
{
}

void QEMSimplification::Simplify(int tar_num_v, double threshold)
{
	InitialEdgeCost();
	index_collapse_v_to_v.clear();
	index_collapse_v_to_v.reserve(mesh.n_vertices());
	size_t counter = 0;
	auto nv = mesh.n_vertices();
	while (counter + tar_num_v < nv && ChooseCollapasedEdge(threshold))
	{
		++counter;
	}
	mesh.garbage_collection(true);
}

const std::vector<int>& QEMSimplification::OriginalIndex(void) const
{
	return mesh.property(vindex).data_vector();
}

const Mesh& QEMSimplification::GetMesh(void) const
{
	return mesh;
}

bool QEMSimplification::ChooseCollapasedEdge(double err_threshhold)
{
	while (!qpq.empty() && ((mesh.status(mesh.edge_handle(mesh.halfedge_handle(qpq.top().idx, 0))).deleted())
		|| (mesh.property(QuaMat_e, qpq.top().idx) != qpq.top().qem)
		|| (mesh.property(valence, mesh.from_vertex_handle(mesh.halfedge_handle(qpq.top().idx, 0))) + mesh.property(valence, mesh.to_vertex_handle(mesh.halfedge_handle(qpq.top().idx, 0))) > 12)))
	{
		qpq.pop();
	}
	if (qpq.empty())
	{
		return false;
	}
	if (qpq.top().qem < err_threshhold)
	{
		auto idx = mesh.halfedge_handle(qpq.top().idx, 0);
		index_collapse_v_to_v.push_back({
			mesh.property(vindex, mesh.from_vertex_handle(idx)),
			mesh.property(vindex, mesh.to_vertex_handle(idx)),
			mesh.property(vindex, mesh.to_vertex_handle(mesh.next_halfedge_handle(idx))),
			mesh.property(vindex, mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(idx))))
			});
		mesh.set_point(mesh.to_vertex_handle(idx), mesh.property(new_point, qpq.top().idx));
		mesh.collapse(idx);
		qpq.pop();
		UpdateEdgeCost(mesh.to_vertex_handle(idx));
		return true;
	}
	else
	{
		std::cout << "can not collapsed anymore !   " << "Current minimize qem is : " << qpq.top().qem << std::endl;
		return false;
	}
}

void QEMSimplification::InitialEdgeCost(void)
{
	for (const auto& fh : mesh.faces())
	{
		auto heh = mesh.halfedge_handle(fh);
		const auto& p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto& p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto& p2 = mesh.point(mesh.opposite_vh(heh));
		Mesh::Point n = (p1 - p0) % (p2 - p0);
		double area = n.norm();
		mesh.property(area_f, fh) = area;
		n = (area == 0.0) ? Mesh::Point(0) : (n / area);
		Eigen::Vector4d vec(n[0], n[1], n[2], -(n | p0));
		mesh.property(QuaMat_f, fh) = vec * vec.transpose();
	}

	for (const auto& vh : mesh.vertices())
	{
		mesh.property(QuaMat_v, vh) = Eigen::Matrix4d::Zero();
		mesh.property(valence, vh) = 0;
		for (const auto& voheh : mesh.voh_range(vh))
		{
			auto fh = mesh.face_handle(voheh);
			if (fh.is_valid())
			{
				mesh.property(QuaMat_v, vh) += mesh.property(QuaMat_f, fh);
			}
			++mesh.property(valence, vh);
		}
	}
	std::vector<EdgeQEM> initialcost;
	initialcost.reserve(mesh.n_edges());
	for (const auto& eh : mesh.edges())
	{
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh0 = mesh.to_vertex_handle(heh);
		auto vh1 = mesh.from_vertex_handle(heh);
		auto mat = mesh.property(QuaMat_v, vh0) + mesh.property(QuaMat_v, vh1);
		if (std::fabs(mat.topLeftCorner(3, 3).determinant()) < 1e-4)
		{
			Mesh::Point newp = (mesh.point(vh0) + mesh.point(vh1)) * 0.5;
			mesh.property(new_point, eh) = newp;
			if (mesh.is_collapse_ok(heh))
			{
				double err = mat(3, 3) + Eigen::RowVector3d(newp[0], newp[1], newp[2]) * (mat.topLeftCorner(3, 3) * Eigen::Vector3d(newp[0], newp[1], newp[2]) + 2.0 * mat.topRightCorner(3, 1));
				mesh.property(QuaMat_e, eh) = err * mesh.property(valence, vh0) * mesh.property(valence, vh1);
			}
			else
			{
				mesh.property(QuaMat_e, eh) = DBL_MAX;
			}
		}
		else
		{
			Eigen::Vector3d newpe = mat.topLeftCorner(3, 3).selfadjointView<Eigen::Upper>().ldlt().solve(mat.topRightCorner(3, 1));
			mesh.property(new_point, eh) = { -newpe[0], -newpe[1], -newpe[2] };
			if (mesh.is_collapse_ok(heh))
			{
				double err = mat(3, 3) - Eigen::Vector3d(mat.topRightCorner(3, 1)).dot(newpe);
				mesh.property(QuaMat_e, eh) = err * mesh.property(valence, vh0) * mesh.property(valence, vh1);
			}
			else
			{
				mesh.property(QuaMat_e, eh) = DBL_MAX;
			}
		}
		initialcost.push_back({ eh, mesh.property(QuaMat_e, eh) });
	}
	qpq = std::priority_queue<EdgeQEM>(initialcost.begin(), initialcost.end());
}

void QEMSimplification::UpdateEdgeCost(const Mesh::VertexHandle& vh)
{
	mesh.property(QuaMat_v, vh) = Eigen::Matrix4d::Zero();
	mesh.property(valence, vh) = 0;
	const auto& p = mesh.point(vh);
	for (const auto& voheh : mesh.voh_range(vh))
	{
		auto vfh = mesh.face_handle(voheh);
		if (vfh.is_valid())
		{
			const auto& p1 = mesh.point(mesh.to_vertex_handle(voheh));
			const auto& p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(voheh)));
			Mesh::Point n = (p1 - p) % (p2 - p);
			double area = n.norm();
			mesh.property(area_f, vfh) = area;
			n /= area;
			Eigen::Vector4d vec(n[0], n[1], n[2], -(n | p));
			mesh.property(QuaMat_f, vfh) = vec * vec.transpose();
			mesh.property(QuaMat_v, vh) += mesh.property(QuaMat_f, vfh);
		}
		++mesh.property(valence, vh);
	}

	for (const auto& vvh : mesh.vv_range(vh))
	{
		mesh.property(QuaMat_v, vvh) = Eigen::Matrix4d::Zero();
		mesh.property(valence, vvh) = 0;
		for (const auto& voheh : mesh.voh_range(vvh))
		{
			auto fh = mesh.face_handle(voheh);
			if (fh.is_valid())
			{
				mesh.property(QuaMat_v, vvh) += mesh.property(QuaMat_f, fh);
			}
			++mesh.property(valence, vvh);
		}
	}

	for (const auto& vvh : mesh.vv_range(vh))
	{
		for (const auto& veh : mesh.ve_range(vvh))
		{
			if (!mesh.property(is_cal_e, veh))
			{
				auto heh = mesh.halfedge_handle(veh, 0);
				auto vh0 = mesh.to_vertex_handle(heh);
				auto vh1 = mesh.from_vertex_handle(heh);
				auto mat = mesh.property(QuaMat_v, vh0) + mesh.property(QuaMat_v, vh1);
				if (std::fabs(mat.topLeftCorner(3, 3).determinant()) < 1e-4)
				{
					Mesh::Point newp = (mesh.point(vh0) + mesh.point(vh1)) * 0.5;
					mesh.property(new_point, veh) = newp;
					if (mesh.is_collapse_ok(heh))
					{
						double err = mat(3, 3) + Eigen::RowVector3d(newp[0], newp[1], newp[2]) * (mat.topLeftCorner(3, 3) * Eigen::Vector3d(newp[0], newp[1], newp[2]) + 2.0 * mat.topRightCorner(3, 1));
						mesh.property(QuaMat_e, veh) = err * mesh.property(valence, vh0) * mesh.property(valence, vh1);
						qpq.push({ veh, mesh.property(QuaMat_e, veh) });
					}
					else
					{
						mesh.property(QuaMat_e, veh) = DBL_MAX;
					}
				}
				else
				{
					Eigen::Vector3d newpe = mat.topLeftCorner(3, 3).selfadjointView<Eigen::Upper>().ldlt().solve(mat.topRightCorner(3, 1));
					mesh.property(new_point, veh) = { -newpe[0], -newpe[1], -newpe[2] };
					if (mesh.is_collapse_ok(heh))
					{
						double err = mat(3, 3) - Eigen::Vector3d(mat.topRightCorner(3, 1)).dot(newpe);
						mesh.property(QuaMat_e, veh) = err * mesh.property(valence, vh0) * mesh.property(valence, vh1);
						qpq.push({ veh, mesh.property(QuaMat_e, veh) });
					}
					else
					{
						mesh.property(QuaMat_e, veh) = DBL_MAX;
					}
				}
				mesh.property(is_cal_e, veh) = true;
			}
		}
	}
	for (const auto& vvh : mesh.vv_range(vh))
	{
		for (const auto& veh : mesh.ve_range(vvh))
		{
			mesh.property(is_cal_e, veh) = false;
		}
	}
}
