#include "QEMSimplification.h"
#include <iostream>

QEMSimplification::QEMSimplification(const Mesh & m)
	:mesh(m)
{
	mesh.add_property(vindex);
	for (const auto & vh : mesh.vertices())
	{
		mesh.property(vindex, vh) = vh.idx();
	}
	mesh.add_property(QuaMat_f);
	mesh.add_property(QuaMat_v);
	mesh.add_property(QuaMat_e);
	mesh.add_property(area_f);
	mesh.add_property(new_point);
	mesh.add_property(is_cal_e);
	mesh.add_property(is_cal_he);
	for (const auto & eh : mesh.edges())
	{
		mesh.property(is_cal_e, eh) = false;
		mesh.property(is_cal_he, mesh.halfedge_handle(eh, 0)) = false;
		mesh.property(is_cal_he, mesh.halfedge_handle(eh, 1)) = false;
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
	auto counter = 0;
	auto nv = mesh.n_vertices();
	while (ChooseCollapasedEdge(threshold))
	{
		counter++;
		if (nv - counter <= tar_num_v)
		{
			mesh.garbage_collection(true);
			break;
		}
	}
	std::cout << "[V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]" << std::endl;
	std::cout << "Simplify is over!  the iter number is:  " << counter << std::endl;
}

int QEMSimplification::OriginalIndex(int currentid) const
{
	return mesh.property(vindex, mesh.vertex_handle(currentid));
}

const Mesh & QEMSimplification::GetMesh(void) const
{
	return mesh;
}

bool QEMSimplification::ChooseCollapasedEdge(double err_threshhold)
{
	EdgeQEM currentEdge;
	int valence = 0;
	do
	{
		if (qpq.size() == 0)
		{
			mesh.garbage_collection(true);
			return false;
		}
		else
		{
			currentEdge = qpq.top();
			qpq.pop();
			valence = mesh.valence(mesh.from_vertex_handle(mesh.halfedge_handle(currentEdge.idx, 0)));
			valence += mesh.valence(mesh.to_vertex_handle(mesh.halfedge_handle(currentEdge.idx, 0)));
		}
	} while ((mesh.status(mesh.edge_handle(mesh.halfedge_handle(currentEdge.idx, 0))).deleted()) || (mesh.property(QuaMat_e, currentEdge.idx) != currentEdge.qem) || (valence > 12));
	if (currentEdge.qem < err_threshhold)
	{
		auto idx = mesh.halfedge_handle(currentEdge.idx, 0);
		std::vector<int> re_in = {
			mesh.property(vindex, mesh.from_vertex_handle(idx)),
			mesh.property(vindex, mesh.to_vertex_handle(idx)),
			mesh.property(vindex, mesh.opposite_vh(idx)),
			mesh.property(vindex, mesh.opposite_he_opposite_vh(idx))
		};
		index_collapse_v_to_v.push_back(re_in);
		mesh.set_point(mesh.to_vertex_handle(idx), mesh.property(new_point, currentEdge.idx));
		mesh.collapse(idx);
		UpdateEdgeCost(mesh.to_vertex_handle(idx));
		return true;
	}
	else
	{
		mesh.garbage_collection(true);
		std::cout << "can not collapsed anymore !   " << "Current minimize qem is : " << currentEdge.qem << std::endl;
		return false;
	}
}

void QEMSimplification::InitialEdgeCost(void)
{
	for (const auto & fh : mesh.faces())
	{
		auto heh = mesh.halfedge_handle(fh);
		const auto & p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto & p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto & p2 = mesh.point(mesh.opposite_vh(heh));
		auto n = (p1 - p0) % (p2 - p0);
		double area = n.norm();
		mesh.property(area_f, fh) = area;
		n = (area == 0.0) ? Mesh::Point(0) : (n / area);

		double a = n[0];
		double b = n[1];
		double c = n[2];
		double d = -(n | p0);
		mesh.property(QuaMat_f, fh) <<
			a * a, a * b, a * c, a * d,
			b * a, b * b, b * c, b * d,
			c * a, c * b, c * c, c * d,
			d * a, d * b, d * c, d * d;
	}

	for (const auto & vh : mesh.vertices())
	{
		mesh.property(QuaMat_v, vh) = Eigen::Matrix4d::Zero();
		for (const auto & vfh : mesh.vf_range(vh))
		{
			mesh.property(QuaMat_v, vh) += mesh.property(QuaMat_f, vfh);
		}
	}
	std::vector<EdgeQEM> initialcost;
	initialcost.reserve(mesh.n_edges());
	for (const auto & eh : mesh.edges())
	{
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh0 = mesh.to_vertex_handle(heh);
		auto vh1 = mesh.from_vertex_handle(heh);
		auto mat = mesh.property(QuaMat_v, vh0) + mesh.property(QuaMat_v, vh1);

		Eigen::Matrix4d mid_m = mat;
		mid_m(3, 0) = mid_m(3, 1) = mid_m(3, 2) = 0.0;
		mid_m(3, 3) = 1.0;
		Eigen::Vector4d b(0.0, 0.0, 0.0, 1.0);
		Eigen::Vector4d v;
		OpenMesh::Vec3d new_p;
		if (fabs(mid_m.determinant()) < 1e-4)
		{
			new_p = (mesh.point(vh0) + mesh.point(vh1)) / 2;
			v << new_p[0], new_p[1], new_p[2], 1.0;
		}
		else
		{
			v = mid_m.colPivHouseholderQr().solve(b);
			new_p = { v[0], v[1], v[2] };
		}
		mesh.property(new_point, eh) = new_p;
		double err = v.transpose() * mat * v;
		err *= (mesh.valence(vh0) * mesh.valence(vh1));
		err = mesh.is_collapse_ok(heh) ? err : DBL_MAX;
		mesh.property(QuaMat_e, eh) = err;
		initialcost.push_back({ eh, err });
	}
	qpq = std::priority_queue<EdgeQEM>(initialcost.begin(), initialcost.end());
}

void QEMSimplification::UpdateEdgeCost(const Mesh::VertexHandle & vh)
{
	mesh.property(QuaMat_v, vh) = Eigen::Matrix4d::Zero();
	for (const auto & vfh : mesh.vf_range(vh))
	{
		auto heh = mesh.halfedge_handle(vfh);
		const auto & p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto & p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto & p2 = mesh.point(mesh.opposite_vh(heh));
		auto n = (p1 - p0) % (p2 - p0);
		double area = n.norm();
		mesh.property(area_f, vfh) = area;
		n /= area;
		double a = n[0];
		double b = n[1];
		double c = n[2];
		double d = -(n | p0);
		Eigen::Matrix4d mat;
		mat << a * a, a*b, a*c, a*d,
			b*a, b*b, b*c, b*d,
			c*a, c*b, c*c, c*d,
			d*a, d*b, d*c, d*d;
		mesh.property(QuaMat_f, vfh) = mat;
		mesh.property(QuaMat_v, vh) += mat;
	}

	for (const auto & vvh : mesh.vv_range(vh))
	{
		mesh.property(QuaMat_v, vvh) = Eigen::Matrix4d::Zero();
		for (const auto & vvfh : mesh.vf_range(vvh))
		{
			mesh.property(QuaMat_v, vvh) += mesh.property(QuaMat_f, vvfh);
		}
	}


	for (const auto & vvh : mesh.vv_range(vh))
	{
		for (const auto & veh : mesh.ve_range(vvh))
		{
			if (!mesh.property(is_cal_e, veh))
			{
				auto heh = mesh.halfedge_handle(veh, 0);
				auto vh0 = mesh.to_vertex_handle(heh);
				auto vh1 = mesh.from_vertex_handle(heh);
				auto mat = mesh.property(QuaMat_v, vh0) + mesh.property(QuaMat_v, vh1);
				Eigen::Matrix4d mid_m = mat;
				mid_m(3, 0) = mid_m(3, 1) = mid_m(3, 2) = 0.0;
				mid_m(3, 3) = 1.0;
				Eigen::Vector4d b(0.0, 0.0, 0.0, 1.0);
				Eigen::Vector4d v;
				OpenMesh::Vec3d new_p;
				if (fabs(mid_m.determinant()) < 1e-4)
				{
					new_p = (mesh.point(vh0) + mesh.point(vh1)) / 2;
					v << new_p[0], new_p[1], new_p[2], 1.0;
				}
				else
				{
					v = mid_m.colPivHouseholderQr().solve(b);
					new_p = { v[0], v[1], v[2] };
				}
				mesh.property(new_point, veh) = new_p;
				double err = v.transpose() * mat * v;
				err *= (mesh.valence(vh0) * mesh.valence(vh1));
				if (mesh.is_collapse_ok(heh) && mesh.is_collapse_ok(mesh.opposite_halfedge_handle(heh)))
				{

					mesh.property(QuaMat_e, veh) = err;
					qpq.push({ veh, err });
				}
				else
				{
					mesh.property(QuaMat_e, veh) = DBL_MAX;
				}
				mesh.property(is_cal_e, veh) = true;
			}
		}
	}
	for (const auto & vvh : mesh.vv_range(vh))
	{
		for (const auto & veh : mesh.ve_range(vvh))
		{
			mesh.property(is_cal_e, veh) = false;
		}
	}
}
