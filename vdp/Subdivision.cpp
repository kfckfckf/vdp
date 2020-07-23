#include "Subdivision.h"
#include <queue>

Subdivision::Subdivision(const Mesh& m)
	:mesh(m)
{
	mesh.add_property(vindex);
	mesh.add_property(curvature);
	for (const auto& vh : mesh.vertices())
	{
		mesh.property(vindex, vh) = vh.idx();
		double c = 0;
		for (const auto& viheh : mesh.vih_range(vh))
		{
			c += mesh.calc_sector_angle(viheh);
		}
		mesh.property(curvature, vh) = 1.0 - c / 2.0 / M_PI;
	}
}

Subdivision::~Subdivision(void)
{
}

void Subdivision::Subdivide(int tarnum)
{
	typedef std::pair<Mesh::EdgeHandle, double> EdgeInfo;
	auto cmp = [](const EdgeInfo& left, const EdgeInfo& right) { return left.second < right.second; };
	std::priority_queue<EdgeInfo, std::vector<EdgeInfo>, decltype(cmp)> pq(cmp);
	for (const auto& eh : mesh.edges())
	{
		pq.push({ eh, mesh.calc_edge_sqr_length(eh) });
	}
	while (mesh.n_vertices() < tarnum)
	{
		auto e0 = pq.top();
		pq.pop();
		auto heh = mesh.halfedge_handle(e0.first, 0);
		const auto& p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto& p1 = mesh.point(mesh.to_vertex_handle(heh));
		auto newvh = mesh.split(e0.first, (p0 + p1) / 2);
		mesh.property(curvature, newvh) = -DBL_MAX;
		mesh.property(vindex, newvh) = -1;
		for (const auto& voheh : mesh.voh_range(newvh))
		{
			pq.push({ mesh.edge_handle(voheh),mesh.calc_edge_sqr_length(voheh) });
			if (mesh.property(curvature, mesh.to_vertex_handle(voheh)) >
				mesh.property(curvature, newvh))
			{
				mesh.property(curvature, newvh) = mesh.property(curvature, mesh.to_vertex_handle(voheh));
				mesh.property(vindex, newvh) = mesh.property(vindex, mesh.to_vertex_handle(voheh));
			}
		}
		if (mesh.property(vindex, newvh) == -1)
		{
			throw std::exception("Error in splitting vertex.");
		}
	}
	mesh.garbage_collection();
}

const Mesh& Subdivision::GetMesh(void) const
{
	return mesh;
}

const std::vector<int>& Subdivision::OriginalIndex(void) const
{
	return mesh.property(vindex).data_vector();
}
