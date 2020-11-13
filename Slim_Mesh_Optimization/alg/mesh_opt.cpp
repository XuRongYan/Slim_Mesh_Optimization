//
// Created by 徐溶延 on 2020/11/10.
//

#include "mesh_opt.h"
#include "grid_meshing/grid_hex_meshing.h"


mesh_opt::mesh_opt(const arguments &args) {
	MD = Mesh_Domain();
	MD.mesh_entire.type = Mesh_type::Hex;
	MD.mesh_subA.type = Mesh_type::Hex;
	MD.mesh_subB.type = Mesh_type::Hex;
	STOP_EXTENT_MIN = STOP_EXTENT_MAX = 15;
	if (args.edge_length_ratio != 0)
		STOP_EXTENT_MIN = STOP_EXTENT_MAX = args.edge_length_ratio;
	weight_opt = args.weight_opt;
}

void mesh_opt::hex_mesh_opt(const Mesh &hex_mesh_in, Mesh &hex_mesh_out) {

	MD.mesh_entire = hex_mesh_in;
	deformation(MD);
	hex_mesh_out = MD.mesh_entire;
}

//deformation
bool mesh_opt::deformation(Mesh_Domain &md) {

	auto &m = md.mesh_entire;
	// get boundary vertices
	std::vector<Hybrid_V> boundary_Vs;
	for (const auto &v : m.Vs) {
		if (v.boundary) {
			boundary_Vs.emplace_back(v);
		}
	}

	scaled_jacobian(m, mq);
	std::cout << "before deformation: minimum scaled J: " << mq.min_Jacobian << " average scaled J: " << mq.ave_Jacobian
			  << endl;

	vector<bool> Huntangle_flag(m.Hs.size(), false), H_flag(m.Hs.size(), false), H_inout_tag(m.Hs.size(), true);
	vector<uint32_t> Hids;
	for (uint32_t i = 0; i < m.Hs.size(); i++)
		Hids.push_back(i);

	double MESHRATIO = 1;
	ts = Tetralize_Set();
	ts.V = m.V.transpose();
	Vector4i t;
	if (m.type == Hex) {
		ts.T.resize(m.Hs.size() * 8, 4);
		for (auto &h : m.Hs) {
			for (uint32_t i = 0; i < 8; i++) {
				for (uint32_t j = 0; j < 4; j++) t[j] = h.vs[hex_tetra_table[i][j]];
				ts.T.row(h.id * 8 + i) = t;
			}
		}
	} else if (m.type == Tet) {
		ts.T.resize(m.Hs.size(), 4);
		for (size_t i = 0; i < m.Hs.size(); i++) {
			for (size_t j = 0; j < 4; j++) {
				t[j] = m.Hs[i].vs[j];
			}
			ts.T.row(i) = t;
		}
	}

	fc.b.resize(boundary_Vs.size());
	fc.bc.resize(boundary_Vs.size(), 3);
	for (size_t i = 0; i < boundary_Vs.size(); i++) {
		fc.b[i] = boundary_Vs[i].id;
		const auto &pos = boundary_Vs[i].v;
		fc.bc.row(i) << pos[0], pos[1], pos[2];
	}

	ts.energy_type = SYMMETRIC_DIRICHLET;
	ts.UV = ts.V;
	ts.fc = fc;
	ts.projection = false;
	ts.global = true;
	ts.glue = false;
	ts.lamda_glue = 0;
	ts.lamda_region = 0;
	ts.record_Sequence = false;
	ts.b = fc.b;
	ts.bc = fc.bc;
	ts.mesh_type = m.type;
	ts.lamda_b = 0;

	optimization opt;
	opt.weight_opt = weight_opt;

	if (m.type == Hex) {
		compute_referenceMesh(ts.V, m.Hs, H_inout_tag, Hids, ts.RT, true);
	}

	opt.slim_m_opt(ts, 30, -1);

	m.V = ts.UV.transpose();
	ts.V = ts.UV;

	for (auto &v : m.Vs) {
		v.v[0] = m.V(0, v.id);
		v.v[1] = m.V(1, v.id);
		v.v[2] = m.V(2, v.id);
	}
	for (auto &v : md.mesh_subA.Vs) {
		v.v = md.mesh_entire.Vs[md.V_map_reverse[v.id]].v;
		md.mesh_subA.V.col(v.id) = md.mesh_entire.V.col(md.V_map_reverse[v.id]);
	}
	return true;
}

