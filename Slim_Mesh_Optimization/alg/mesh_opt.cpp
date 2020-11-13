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

bool mesh_opt::feature(Mesh &mesh) {
	mf.angle_threshold = 0;
	if (!mf.read_from_file) {
		std::cout << "no feature file, detect based on angle" << std::endl;
		mf.angle_threshold = 0;
		mf.orphan_curve = true;
		mf.orphan_curve_single = true;
	}

	mf.tri = mesh;
	triangle_mesh_feature(mf);
	build_feature_graph(mf, fg);
	return false;
}

void mesh_opt::build_aabb_tree(Mesh &tmi, Treestr &a_tree, bool is_tri) {
	a_tree.TriV = tmi.V.transpose();
	a_tree.TriF.resize(tmi.Fs.size(), 3);

	for (uint32_t i = 0; i < tmi.Fs.size(); i++) {
		a_tree.TriF(i, 0) = tmi.Fs[i].vs[0];
		a_tree.TriF(i, 1) = tmi.Fs[i].vs[1];
		a_tree.TriF(i, 2) = tmi.Fs[i].vs[2];
	}
	if (!is_tri) {
		// Precompute signed distance AABB tree
		a_tree.tree.init(a_tree.TriV, a_tree.TriF);
		// Precompute vertex,edge and face normals
		igl::per_face_normals(a_tree.TriV, a_tree.TriF, a_tree.TriFN);
		igl::per_vertex_normals(a_tree.TriV, a_tree.TriF, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, a_tree.TriFN,
								a_tree.TriVN);
		igl::per_edge_normals(a_tree.TriV, a_tree.TriF, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, a_tree.TriFN,
							  a_tree.TriEN, a_tree.TriE, a_tree.TriEMAP);
	}
}

bool mesh_opt::hausdorff_ratio_check(Mesh &m0, Mesh &m1, double &hausdorff_dis_threshold) {
	std::function<void(Mesh &, Mesh &, vector<bool> &, int &)> re_indexing = [&](Mesh &M, Mesh &m, vector<bool> &V_flag,
																				 int &N) -> void {
		m.V.resize(3, N);
		N = 0;
		vector<int> v_map(M.Vs.size(), 0);
		for (uint32_t i = 0; i < V_flag.size(); i++)
			if (V_flag[i]) {
				m.V.col(N) = M.V.col(i);
				v_map[i] = N++;
			}
		for (auto f : M.Fs) {
			if (!f.boundary) continue;
			vector<uint32_t> bvs;
			for (auto vid : f.vs)if (V_flag[vid]) bvs.push_back(v_map[vid]);
			if (bvs.size() == 3) {
				Hybrid_F hf;
				hf.vs = bvs;
				m.Fs.push_back(hf);
			} else if (bvs.size() == 4) {
				Hybrid_F hf;
				hf.vs.push_back(bvs[0]);
				hf.vs.push_back(bvs[1]);
				hf.vs.push_back(bvs[2]);
				m.Fs.push_back(hf);
				hf.vs.clear();
				hf.vs.push_back(bvs[2]);
				hf.vs.push_back(bvs[3]);
				hf.vs.push_back(bvs[0]);
				m.Fs.push_back(hf);
			}
		}
	};
	vector<bool> V_flag;
	int bvN = 0, bbvN = 0;

	Mesh Mglobal;
	V_flag.resize(m1.Vs.size());
	std::fill(V_flag.begin(), V_flag.end(), false);
	bvN = 0;
	for (uint32_t i = 0; i < m1.Vs.size(); i++)
		if (m1.Vs[i].boundary) {
			V_flag[i] = true;
			bvN++;
		}
	re_indexing(m1, Mglobal, V_flag, bvN);

	if (!compute(m0, Mglobal, hausdorff_ratio, hausdorff_ratio_threshould, hausdorff_dis_threshold)) {
		cout << "too large hausdorff distance" << endl;
		return false;
	}
	return true;
}

//boundary projection & feature capture
bool mesh_opt::surface_mapping(Mesh &tmi, Mesh_Domain &md) {
	Mesh_Feature mf_temp;
	vector<bool> Corner_tag;
	vector<vector<uint32_t>> circle2curve_map;
	break_circles(mf_temp, Corner_tag, circle2curve_map);

	//build surface feture graph
	Feature_Graph Tfg;

	md.clear_Qfg();
	md.clear_mf_quad();
	Mesh htri;
	htri.type = Mesh_type::Tri;
	auto &hqua = md.quad_tree.mesh;
	hqua.type = Mesh_type::Qua;

	extract_surface_conforming_mesh(md.mesh_subA, htri, md.TV_map, md.TV_map_reverse, md.TF_map, md.TF_map_reverse);
	extract_surface_conforming_mesh(md.mesh_subA, hqua, md.QV_map, md.QV_map_reverse, md.QF_map, md.QF_map_reverse);
	//aabb stree for surface of the hex-mesh
	build_aabb_tree(htri, md.quad_tree);

	cout << "node mapping" << endl;
	if (!node_mapping(mf_temp, Tfg, md))
		return false;
	cout << "curve mapping" << endl;
	curve_mapping(Tfg, md);
	if (!reunion_circles(mf_temp, Corner_tag, circle2curve_map, md))
		return false;
	cout << "patch mapping" << endl;
	if (!patch_mapping(md))
		return false;
	patch_trees(md);

	Q_final_fg = md.Qfg;
	GRAPH_MATCHES = md.graph_matches;
	return true;
}

void
mesh_opt::break_circles(Mesh_Feature &mf_temp, vector<bool> &Corner_tag, vector<vector<uint32_t>> &circle2curve_map) {
	mf_temp.corners = mf.corners;
	mf_temp.corner_curves = mf.corner_curves;
	mf_temp.curve_vs = mf.curve_vs;

	Corner_tag.resize(mf.corners.size(), true);
	circle2curve_map.resize(mf.circles.size());
	for (uint32_t i = 0; i < mf.circles.size(); i++) {
		if (mf.circles[i]) {
			int lenv = mf.curve_vs[i].size();
			int vid0 = mf.curve_vs[i][0], vid1 = mf.curve_vs[i][lenv / 2];
			mf_temp.corners.push_back(vid0);
			Corner_tag.push_back(false);
			mf_temp.corners.push_back(vid1);
			Corner_tag.push_back(false);
			vector<uint32_t> curve_vs;
			for (auto vid : mf.curve_vs[i])
				if (vid == vid1) {
					curve_vs.push_back(vid);
					mf_temp.curve_vs[i] = curve_vs;
					curve_vs.clear();
					curve_vs.push_back(vid);
				} else curve_vs.push_back(vid);
			curve_vs.push_back(vid0);
			mf_temp.curve_vs.push_back(curve_vs);
			vector<uint32_t> curve_ids;
			curve_ids.push_back(i);
			curve_ids.push_back(mf_temp.curve_vs.size() - 1);
			mf_temp.corner_curves.push_back(curve_ids);
			mf_temp.corner_curves.push_back(curve_ids);
			circle2curve_map[i] = curve_ids;
		} else circle2curve_map[i].push_back(i);
	}
}

bool mesh_opt::node_mapping(Mesh_Feature &mf_temp, Feature_Graph &Tfg, Mesh_Domain &md) {
	auto &quad_tree = md.quad_tree;
	auto &TF_map_reverse = md.TF_map_reverse;
	auto &QV_map_reverse = md.QV_map_reverse;
	auto &QF_map = md.QF_map;
	auto &Qfg = md.Qfg;
	//find corresponding vs of corners on the triangle mesh
	MatrixXd Ps(mf_temp.corners.size(), 3);
	int num_corner = 0;
	for (auto vid : mf_temp.corners) { Ps.row(num_corner++) = mf.tri.V.col(vid).transpose(); }

	VectorXd signed_dis;
	VectorXi ids;
	MatrixXd VT, NT;
	signed_distance_pseudonormal(Ps, quad_tree.TriV, quad_tree.TriF, quad_tree.tree, quad_tree.TriFN, quad_tree.TriVN,
								 quad_tree.TriEN, quad_tree.TriEMAP, signed_dis, ids, VT, NT);
	//build nodes of the graph for both triangle mesh and quad mesh
	vector<bool> TE_tag(mf.tri.Es.size(), false), QE_tag(quad_tree.mesh.Es.size(), false);
	vector<bool> V_flag(quad_tree.mesh.Vs.size(), true);

	for (uint32_t i = 0; i < mf_temp.corners.size(); i++) {
		Vector3d v = mf.tri.V.col(mf_temp.corners[i]);
		int qid = QF_map[TF_map_reverse[ids[i]]];
		int vid = -1, vid_min = -1;
		double min_dis = std::numeric_limits<double>::max();
		for (uint32_t j = 0; j < quad_tree.mesh.Fs[qid].vs.size(); j++) {
			int vid_ = quad_tree.mesh.Fs[qid].vs[j];
			if (!V_flag[vid_])continue;//already taken by other corners
			double dis = (v - quad_tree.mesh.V.col(vid_)).norm();
			if (j == 0 || min_dis > dis) {
				vid_min = vid_;
				min_dis = dis;
				if (quad_tree.mesh.Vs[vid_].neighbor_es.size() >= mf_temp.corner_curves[i].size())
					vid = vid_;
			}
		}
		Feature_Corner c;
		c.id = Qfg.Cs.size();
		c.original = v;
		c.projected = VT.row(i);

		if (vid == -1) {
			if (vid_min == -1) {
				vid_min = quad_tree.mesh.Fs[qid].vs[0];
			}
			for (auto nfid : quad_tree.mesh.Vs[vid_min].neighbor_fs) {
				for (uint32_t j = 0; j < quad_tree.mesh.Fs[nfid].vs.size(); j++) {
					int vid_ = quad_tree.mesh.Fs[nfid].vs[j];
					if (!V_flag[vid_])continue;

					double dis = (v - quad_tree.mesh.V.col(vid_)).norm();
					if (quad_tree.mesh.Vs[vid_].neighbor_es.size() >= mf_temp.corner_curves[i].size() &&
						(j == 0 || min_dis > dis)) {
						vid = vid_;
						min_dis = dis;
					}
				}
			}
			if (vid == -1) {
				auto g_id = md.V_map_reverse[QV_map_reverse[vid_min]];
				vector<uint32_t> nvs = md.mesh_entire.Vs[g_id].neighbor_vs;
				tb_subdivided_cells.insert(tb_subdivided_cells.end(), nvs.begin(), nvs.end());
			} else c.vs.push_back(vid);
		} else c.vs.push_back(vid);

		for (auto vid : c.vs) {
			V_flag[vid] = false;
		}
		Qfg.Cs.push_back(c);
	}

	if (tb_subdivided_cells.size()) return false;
//Tfg
	for (uint32_t i = 0; i < mf_temp.corners.size(); i++) {
		Feature_Corner tc;
		tc.id = Tfg.Cs.size();
		tc.vs.push_back(mf_temp.corners[i]);
		//neighbor vs
		vector<uint32_t> nts;
		for (auto vid : tc.vs)
			nts.insert(nts.end(), mf.tri.Vs[vid].neighbor_fs.begin(), mf.tri.Vs[vid].neighbor_fs.end());
		std::sort(nts.begin(), nts.end());
		nts.erase(unique(nts.begin(), nts.end()), nts.end());
		face_soup_info(mf.tri, nts, TE_tag, tc.ring_vs);

		Tfg.Cs.push_back(tc);
	}
	Tfg.Ls.resize(mf_temp.curve_vs.size());
	for (uint32_t i = 0; i < Tfg.Cs.size(); i++) {
		Feature_Corner &tc = Tfg.Cs[i];
		tc.neighbor_ls = mf_temp.corner_curves[i];
		for (auto cid : mf_temp.corner_curves[i]) Tfg.Ls[cid].cs.push_back(tc.id);
	}
	for (int i = 0; i < Tfg.Ls.size(); i++) {
		Tfg.Ls[i].id = i;
		Tfg.Ls[i].vs = mf_temp.curve_vs[i];
	}
	//corner mismatch check
	for (uint32_t i = 0; i < mf_temp.corners.size(); i++) {
		Feature_Corner &c = Qfg.Cs[i];
		//neighbor vs
		vector<uint32_t> nqs;
		for (auto vid : c.vs)
			nqs.insert(nqs.end(), quad_tree.mesh.Vs[vid].neighbor_fs.begin(), quad_tree.mesh.Vs[vid].neighbor_fs.end());
		std::sort(nqs.begin(), nqs.end());
		nqs.erase(unique(nqs.begin(), nqs.end()), nqs.end());
		vector<uint32_t> ring_vs;
		face_soup_info(quad_tree.mesh, nqs, QE_tag, ring_vs);

		vector<int> wrong_corners;
		for (auto vid : ring_vs) {
			bool is_neighbor = false;
			for (auto vid_ : c.vs) {
				if (find(quad_tree.mesh.Vs[vid_].neighbor_vs.begin(), quad_tree.mesh.Vs[vid_].neighbor_vs.end(), vid) !=
					quad_tree.mesh.Vs[vid_].neighbor_vs.end()) {
					is_neighbor = true;
					break;
				}
			}
			if (!is_neighbor) continue;
			//conflicting corners, ring_vs with other corners, ring_vs of other corners
			bool wrong_connection = false;
			if (!V_flag[vid]) {
				for (const auto &c_ : Qfg.Cs)
					if (c_.vs[0] == vid) {
						auto es0 = Tfg.Cs[c.id].neighbor_ls, es1 = Tfg.Cs[c_.id].neighbor_ls;
						vector<uint32_t> sharedes;
						std::sort(es0.begin(), es0.end());
						std::sort(es1.begin(), es1.end());
						set_intersection(es0.begin(), es0.end(), es1.begin(), es1.end(), back_inserter(sharedes));
						if (!sharedes.size()) {
							wrong_corners.push_back(c_.vs[0]);
							wrong_connection = true;// break;
						}
					} else {
						for (auto ovid : c_.ring_vs)
							if (vid == ovid) {
								auto es0 = Tfg.Cs[c.id].neighbor_ls, es1 = Tfg.Cs[c_.id].neighbor_ls;
								vector<uint32_t> sharedes;
								std::sort(es0.begin(), es0.end());
								std::sort(es1.begin(), es1.end());
								set_intersection(es0.begin(), es0.end(), es1.begin(), es1.end(),
												 back_inserter(sharedes));
								if (!sharedes.size()) {
									wrong_corners.push_back(c_.vs[0]);
									wrong_connection = true; //break;
								}
							}
					}
			}
			if (wrong_connection) continue;
			c.ring_vs.push_back(vid);
		}
		if (c.ring_vs.size() < mf_temp.corner_curves[i].size()) {
			auto g_id = md.V_map_reverse[QV_map_reverse[c.vs[0]]];
			tb_subdivided_cells.push_back(g_id);
			for (auto wc : wrong_corners) {
				auto g_id = md.V_map_reverse[QV_map_reverse[wc]];
				tb_subdivided_cells.push_back(g_id);
			}
		}
		for (auto vid : c.ring_vs)V_flag[vid] = false;
	}
	if (tb_subdivided_cells.size()) return false;

	//corner direction match
	Qfg.Ls.resize(mf_temp.curve_vs.size());
	for (uint32_t i = 0; i < Tfg.Cs.size(); i++) {
		Feature_Corner &tc = Tfg.Cs[i];
		tc.ring_vs_tag.resize(tc.ring_vs.size(), -1);
		tc.neighbor_ls = mf_temp.corner_curves[i];

		vector<uint32_t> vs;
		for (auto cid : mf_temp.corner_curves[i]) {
			vector<uint32_t> &curve = mf_temp.curve_vs[cid];
			if (curve[0] == mf_temp.corners[i]) vs.push_back(curve[1]);
			else if (curve[curve.size() - 1] == mf_temp.corners[i]) vs.push_back(curve[curve.size() - 2]);
		}
		vector<Vector3d> Tdirs;
		vector<int> lids;
		for (uint32_t j = 0; j < tc.ring_vs.size(); j++) {
			if (find(vs.begin(), vs.end(), tc.ring_vs[j]) != vs.end()) {
				tc.ring_vs_tag[j] = mf_temp.corner_curves[i][find(vs.begin(), vs.end(), tc.ring_vs[j]) - vs.begin()];
				lids.push_back(tc.ring_vs_tag[j]);
				Tdirs.push_back((mf.tri.V.col(tc.ring_vs[j]) - mf.tri.V.col(tc.vs[0])).normalized());
			}
		}
		Feature_Corner &qc = Qfg.Cs[i];
		vector<Vector3d> Qdirs;
		for (uint32_t j = 0; j < qc.ring_vs.size(); j++)
			Qdirs.push_back((quad_tree.mesh.V.col(qc.ring_vs[j]) - quad_tree.mesh.V.col(qc.vs[0])).normalized());
		//compare Tdirs and Qdirs
		vector<bool> Qdirs_tag(Qdirs.size(), false);
		std::fill(Qdirs_tag.end() - Tdirs.size(), Qdirs_tag.end(), true);
		vector<int> best_set;
		double max_alignment = -std::numeric_limits<double>::infinity();
		do {
			vector<int> cur_set;
			for (uint32_t j = 0; j < Qdirs_tag.size(); j++) if (Qdirs_tag[j]) cur_set.push_back(j);

			int start_id = -1;
			for (uint32_t j = 0; j < cur_set.size(); j++) {
				double cost = 0;
				for (uint32_t k = 0; k < cur_set.size(); k++)
					cost += Tdirs[k].dot(Qdirs[cur_set[(k + j) % cur_set.size()]]);
				if (cost > max_alignment) {
					start_id = j;
					max_alignment = cost;
				}
			}
			if (start_id != -1) {
				best_set.clear();
				for (uint32_t k = 0; k < cur_set.size(); k++)
					best_set.push_back(cur_set[(k + start_id) % cur_set.size()]);
			}
		} while (std::next_permutation(Qdirs_tag.begin(), Qdirs_tag.end()));

		qc.ring_vs_tag.resize(qc.ring_vs.size(), -1);
		for (uint32_t j = 0; j < best_set.size(); j++)qc.ring_vs_tag[best_set[j]] = lids[j];
	}
	//corner deep mismatch check
	fill(V_flag.begin(), V_flag.end(), true);
	for (const auto &c : Qfg.Cs) {
		for (const auto &vid : c.vs) V_flag[vid] = false;
		for (const auto &vid : c.ring_vs) V_flag[vid] = false;
	}
	for (const auto &c : Qfg.Cs) {
		vector<int> wrong_corners;
		for (int i = 0; i < c.ring_vs_tag.size(); i++) {
			if (c.ring_vs_tag[i] == -1)continue;
			int lid = c.ring_vs_tag[i];
			const auto vid = c.ring_vs[i];
			bool wrong_connection = false;
			if (!V_flag[vid]) {
				for (const auto &c_ : Qfg.Cs)
					if (c_.vs[0] == vid) {
						if (std::find(c_.ring_vs_tag.begin(), c_.ring_vs_tag.end(), lid) != c_.ring_vs_tag.end()) {
							auto id = std::find(c_.ring_vs_tag.begin(), c_.ring_vs_tag.end(), lid) -
									  c_.ring_vs_tag.begin();
							if (c_.ring_vs[id] != c.vs[0]) {
								wrong_corners.push_back(c_.vs[0]);
								wrong_connection = true;// break;
							}
						} else {
							wrong_corners.push_back(c_.vs[0]);
							wrong_connection = true;// break;
						}
					} else {
						for (auto ovid : c_.ring_vs)
							if (vid == ovid) {
								if (std::find(c_.ring_vs_tag.begin(), c_.ring_vs_tag.end(), lid) !=
									c_.ring_vs_tag.end()) {
									auto id = std::find(c_.ring_vs_tag.begin(), c_.ring_vs_tag.end(), lid) -
											  c_.ring_vs_tag.begin();
									if (c_.ring_vs[id] != c.ring_vs[i]) {
										wrong_corners.push_back(c_.vs[0]);
										wrong_connection = true;// break;
									}
								} else {
									wrong_corners.push_back(c_.vs[0]);
									wrong_connection = true;// break;
								}
							}
					}
			}
		}
		if (wrong_corners.size()) {
			auto g_id = md.V_map_reverse[QV_map_reverse[c.vs[0]]];
			tb_subdivided_cells.push_back(g_id);
			for (auto wc : wrong_corners) {
				auto g_id = md.V_map_reverse[QV_map_reverse[wc]];
				tb_subdivided_cells.push_back(g_id);
			}
		}
	}
	if (tb_subdivided_cells.size()) return false;
	return true;
}

bool mesh_opt::curve_mapping(Feature_Graph &Tfg, Mesh_Domain &md) {
	auto &quad_tree = md.quad_tree;
	auto &TF_map_reverse = md.TF_map_reverse;
	auto &QF_map = md.QF_map;
	auto &Qfg = md.Qfg;
//setup before line trace
	vector<bool> V_flag(quad_tree.mesh.Vs.size(), true);
	for (auto &c : Qfg.Cs) {
		for (auto vid : c.vs) V_flag[vid] = false;
		for (uint32_t i = 0; i < c.ring_vs_tag.size(); i++) if (c.ring_vs_tag[i] != -1) V_flag[c.ring_vs[i]] = false;
	}

	adjacency_list_t adjacency_list(quad_tree.mesh.Vs.size()), adjacency_list_w;
	std::vector<weight_t> min_distance;
	std::vector<vertex_t> previous;
	for (auto &v : quad_tree.mesh.Vs) {
		for (auto &nv : v.neighbor_vs) adjacency_list[v.id].push_back(neighbor(nv, 1));
	}
	adjacency_list_w = adjacency_list;
	for (auto &v : quad_tree.mesh.Vs) {
		if (V_flag[v.id]) continue;
		for (auto &nv : v.neighbor_vs) adjacency_list[v.id].clear();
	}

	//trace lines
	for (uint32_t i = 0; i < Tfg.Ls.size(); i++) {
		Qfg.Ls[i].id = i;
		Qfg.Ls[i].cs = Tfg.Ls[i].cs;
		//start to trace a curve
		vector<uint32_t> &curve = Tfg.Ls[i].vs;
		int start_cid = Tfg.Ls[i].cs[0], start2_vid = -1, end2_vid = -1, end_cid = Tfg.Ls[i].cs[1];

		if (curve[0] == Tfg.Cs[end_cid].vs[0]) swap(start_cid, end_cid);
		int which_id = find(Qfg.Cs[start_cid].ring_vs_tag.begin(), Qfg.Cs[start_cid].ring_vs_tag.end(), i) -
					   Qfg.Cs[start_cid].ring_vs_tag.begin();
		start2_vid = Qfg.Cs[start_cid].ring_vs[which_id];
		which_id = find(Qfg.Cs[end_cid].ring_vs_tag.begin(), Qfg.Cs[end_cid].ring_vs_tag.end(), i) -
				   Qfg.Cs[end_cid].ring_vs_tag.begin();
		end2_vid = Qfg.Cs[end_cid].ring_vs[which_id];

		if (start2_vid == Qfg.Cs[end_cid].vs[0] || Qfg.Cs[start_cid].vs[0] == end2_vid) {
			Qfg.Ls[i].vs.push_back(Qfg.Cs[start_cid].vs[0]);
			Qfg.Ls[i].vs.push_back(Qfg.Cs[end_cid].vs[0]);

			for (auto vid : Qfg.Ls[i].vs)adjacency_list[vid].clear();
			continue;
		} else if (start2_vid == end2_vid) {
			Qfg.Ls[i].vs.push_back(Qfg.Cs[start_cid].vs[0]);
			Qfg.Ls[i].vs.push_back(start2_vid);
			Qfg.Ls[i].vs.push_back(Qfg.Cs[end_cid].vs[0]);

			for (auto vid : Qfg.Ls[i].vs)adjacency_list[vid].clear();
			continue;
		}
		//sample line, compute weight of the mesh, update adjacency_list, run dijstra
		for (auto nv : quad_tree.mesh.Vs[start2_vid].neighbor_vs) adjacency_list[start2_vid].push_back(neighbor(nv, 1));
		for (auto nv : quad_tree.mesh.Vs[end2_vid].neighbor_vs) adjacency_list[end2_vid].push_back(neighbor(nv, 1));

		double sample_density = voxel_size * pow(2, STOP_EXTENT_MAX) * 0.5;
		if (!args.octree)
			sample_density = voxel_size * 0.5;

		vector<Vector3d> samples;
		int l_len = Tfg.Ls[i].vs.size() - 1;
		double len_cur = 0;
		int num_segment = 0;
		for (uint32_t j = 0; j < l_len; j++) {
			Vector3d v0 = mf.tri.V.col(Tfg.Ls[i].vs[j]), v1 = mf.tri.V.col(Tfg.Ls[i].vs[j + 1]);
			double dis = (v0 - v1).norm(), t = 0;
			if (dis < 1.0e-15)continue;

			while (len_cur <= num_segment * sample_density && len_cur + dis > num_segment * sample_density) {
				t = (len_cur + dis - num_segment * sample_density) / dis;
				samples.push_back(v0 + (1 - t) * (v1 - v0));
				num_segment++;
			}
			len_cur += dis;
		}
		Qfg.Ls[i].guiding_v = samples;

		MatrixXd Ps(samples.size(), 3);
		VectorXd signed_dis;
		VectorXi ids;
		MatrixXd VT, NT;
		for (uint32_t j = 0; j < samples.size(); j++) Ps.row(j) = samples[j].transpose();
		signed_distance_pseudonormal(Ps, quad_tree.TriV, quad_tree.TriF, quad_tree.tree, quad_tree.TriFN,
									 quad_tree.TriVN, quad_tree.TriEN, quad_tree.TriEMAP, signed_dis, ids, VT, NT);
		vector<int> l_source;
		for (uint32_t j = 0; j < samples.size(); j++) {
			Vector3d &v = samples[j];
			int qid = QF_map[TF_map_reverse[ids[j]]];
			int vid = -1;
			double min_dis;
			for (uint32_t k = 0; k < quad_tree.mesh.Fs[qid].vs.size(); k++) {
				int vid_ = quad_tree.mesh.Fs[qid].vs[k];
				double dis = (v - quad_tree.mesh.V.col(vid_)).norm();
				if (k == 0 || min_dis > dis) {
					vid = vid_;
					min_dis = dis;
				}
			}
			l_source.push_back(vid);
		}
		Qfg.Ls[i].guiding_vs = l_source;
		DijkstraComputePaths(l_source, adjacency_list_w, min_distance, previous);
		for (uint32_t j = 0; j < adjacency_list.size(); j++)
			for (uint32_t k = 0; k < adjacency_list[j].size(); k++) {
				adjacency_list[j][k].weight = 2 * max(min_distance[j], min_distance[adjacency_list[j][k].target]);
			}

		DijkstraComputePaths(start2_vid, adjacency_list, min_distance, previous);
		vector<int> path = DijkstraGetShortestPathTo(end2_vid, previous);

		Qfg.Ls[i].vs.push_back(Qfg.Cs[start_cid].vs[0]);
		Qfg.Ls[i].vs.insert(Qfg.Ls[i].vs.end(), path.begin(), path.end());

		Qfg.Ls[i].vs.push_back(Qfg.Cs[end_cid].vs[0]);

		for (auto vid : Qfg.Ls[i].vs)adjacency_list[vid].clear();
	}
	return true;
}

bool mesh_opt::reunion_circles(
		Mesh_Feature &mf_temp, vector<bool> &Corner_tag, vector<vector<uint32_t>> &circle2curve_map, Mesh_Domain &md) {
	auto &quad_tree = md.quad_tree;
	auto &TF_map_reverse = md.TF_map_reverse;
	auto &QV_map_reverse = md.QV_map_reverse;
	auto &QF_map = md.QF_map;
	auto &Qfg = md.Qfg;
	auto &Qfg_sur = md.Qfg_sur;
	auto &mf_quad = md.mf_quad;
	//compute Qfg
	for (uint32_t i = 0; i < Qfg.Cs.size(); i++) if (Corner_tag[i]) mf_quad.corners.push_back(Qfg.Cs[i].vs[0]);
	mf_quad.corner_curves = mf.corner_curves;
	mf_quad.curve_vs.resize(mf.curve_vs.size());
	mf_quad.curve_es.resize(mf.curve_es.size());
	mf_quad.broken_curves.resize(mf.curve_es.size());
	std::fill(mf_quad.broken_curves.begin(), mf_quad.broken_curves.end(), false);
	mf_quad.circles = mf.circles;

	mf_quad.ave_length = average_edge_length(quad_tree.mesh);
	mf_quad.tri = quad_tree.mesh;

	function<bool(Mesh &, vector<uint32_t> &, vector<uint32_t> &, const bool &)> VSES = [&](Mesh &mesh,
																							vector<uint32_t> &vs,
																							vector<uint32_t> &es,
																							const bool &circle) -> bool {
		int vlen = vs.size(), vlen_ = vlen;
		if (!circle) vlen--;
		es.clear();
		bool broken_curve = false;
		for (uint32_t i = 0; i < vlen; i++) {
			int v0 = vs[i], v1 = vs[(i + 1) % vlen_];
			vector<uint32_t> sharedes, es0 = mesh.Vs[v0].neighbor_es, es1 = mesh.Vs[v1].neighbor_es;
			std::sort(es0.begin(), es0.end());
			std::sort(es1.begin(), es1.end());
			set_intersection(es0.begin(), es0.end(), es1.begin(), es1.end(), back_inserter(sharedes));
			if (sharedes.size()) es.push_back(sharedes[0]);
			else {
				broken_curve = true;
			}
		}
		return broken_curve;
	};

	for (uint32_t i = 0; i < circle2curve_map.size(); i++) {
		if (circle2curve_map[i].size() == 2) {
			//cout << "circular curve " << i << endl;
			vector<uint32_t> vs0, es0, vs1, es1;
			vs0 = Qfg.Ls[circle2curve_map[i][0]].vs;
			vs1 = Qfg.Ls[circle2curve_map[i][1]].vs;
			if (VSES(quad_tree.mesh, vs0, es0, false)) mf_quad.broken_curves[i] = true;
			if (VSES(quad_tree.mesh, vs1, es1, false)) mf_quad.broken_curves[i] = true;
			mf_quad.curve_vs[i] = vs0;
			mf_quad.curve_es[i] = es0;
			if (vs0[0] == vs1[0]) {
				for (uint32_t j = vs1.size() - 2; j > 1; j--) mf_quad.curve_vs[i].push_back(vs1[j]);
				reverse(es0.begin(), es0.end());
				mf_quad.curve_es[i].insert(mf_quad.curve_es[i].end(), es1.begin(), es1.end());
			} else if (vs0[0] == vs1[vs1.size() - 1]) {
				for (uint32_t j = 1; j < vs1.size() - 1; j++) mf_quad.curve_vs[i].push_back(vs1[j]);
				mf_quad.curve_es[i].insert(mf_quad.curve_es[i].end(), es1.begin(), es1.end());
			} else cout << "Bug!!!" << endl;
		} else {
			//cout << "non-circular curve " << i << endl;
			mf_quad.curve_vs[i] = Qfg.Ls[i].vs;
			if (VSES(quad_tree.mesh, mf_quad.curve_vs[i], mf_quad.curve_es[i], false))
				mf_quad.broken_curves[i] = true;
		}
	}

	Qfg_sur = Qfg;
	Qfg.Cs.clear();
	Qfg.Ls.clear();
	Qfg.Rs.clear();

	build_feature_graph(mf_quad, Qfg);

	//handle broken curves
	vector<bool> V_flag(md.mesh_subA.Vs.size(), false);
	for (int i = 0; i < mf_quad.broken_curves.size(); i++) {
		if (!mf_quad.broken_curves[i]) continue;

		for (auto vid : mf_quad.curve_vs[i])V_flag[QV_map_reverse[vid]] = true;

		vector<int32_t> guiding_vs;
		for (auto lid: circle2curve_map[i])
			guiding_vs.insert(guiding_vs.end(), Qfg_sur.Ls[lid].guiding_vs.begin(), Qfg_sur.Ls[lid].guiding_vs.end());

		for (auto vid : guiding_vs) {
			if (!V_flag[QV_map_reverse[vid]]) {
				auto g_id = md.V_map_reverse[QV_map_reverse[vid]];
				tb_subdivided_cells.push_back(g_id);
			}
		}
	}
	if (tb_subdivided_cells.size()) {
		cout << "broen curves!!!  " << tb_subdivided_cells.size() << endl;
		return false;
	}

	if (tb_subdivided_cells.size())return false;

	//update mf_quad
	for (auto &vid: mf_quad.corners) vid = md.V_map_reverse[QV_map_reverse[vid]];
	for (auto &vs:mf_quad.curve_vs) for (auto &vid : vs) vid = md.V_map_reverse[QV_map_reverse[vid]];
	for (uint32_t i = 0; i < mf_quad.curve_es.size(); i++)
		VSES(md.mesh_entire, mf_quad.curve_vs[i], mf_quad.curve_es[i], mf_quad.circles[i]);

	return true;
}

bool mesh_opt::patch_mapping(Mesh_Domain &md) {
	auto &Qfg = md.Qfg;
	auto &Qfg_sur = md.Qfg_sur;
	auto &graph_matches = md.graph_matches;
	//matching graph fg & Qfg!
	graph_matches.clear();
	graph_matches.resize(Qfg.Rs.size());
	//1. compare edges and corners
	for (auto r : Qfg.Rs) {
		//find its corresponding patches in fg
		vector<uint32_t> ls = r.ls;

		if (!ls.size()) {
			graph_matches[r.id].push_back(0);
			continue;
		}
		std::sort(ls.begin(), ls.end());
		bool found = false;
		for (auto l : ls) {
			for (auto trid : fg.Ls[l].neighbor_rs) {
				vector<uint32_t> tls = fg.Rs[trid].ls;
				if (tls.size() != ls.size()) continue;
				std::sort(tls.begin(), tls.end());
				if (std::equal(ls.begin(), ls.end(), tls.begin())) {
					graph_matches[r.id].push_back(trid);
					found = true;
				}
			}
		}
		sort(graph_matches[r.id].begin(), graph_matches[r.id].end());
		graph_matches[r.id].erase(unique(graph_matches[r.id].begin(), graph_matches[r.id].end()),
								  graph_matches[r.id].end());
	}

	cout << "start to find bad lines!" << endl;
	//kill miss match because of shared edges!
	vector<uint32_t> lines, regions;
	for (auto r : Qfg.Rs) {
		if (graph_matches[r.id].size())continue;
		lines.insert(lines.end(), r.ls.begin(), r.ls.end());
		regions.push_back(r.id);
		std::cout << "missed regions: " << r.id << std::endl;
	}
	Feature_Graph tfg_temp = fg, tfg_temp2 = fg;
	tfg_temp.Ls.clear();
	tfg_temp.Rs.clear();
	tfg_temp2.Ls.clear();
	tfg_temp2.Rs.clear();
	vector<bool> r_flag(fg.Rs.size(), false);
	for (auto m : graph_matches) for (auto rid : m)r_flag[rid] = true;

	vector<uint32_t> tlines, tregions;
	for (auto r : fg.Rs)
		if (!r_flag[r.id]) {
			tlines.insert(tlines.end(), r.ls.begin(), r.ls.end());
			tregions.push_back(r.id);
		}
	vector<bool> tl_flag(fg.Ls.size(), false);
	for (auto lid : tlines)
		if (tl_flag[lid])tfg_temp2.Ls.push_back(fg.Ls[lid]);
		else {
			tfg_temp.Ls.push_back(fg.Ls[lid]);
			tl_flag[lid] = true;
		}

	for (auto l : tfg_temp2.Ls)
		for (auto vid : Qfg_sur.Ls[l.id].guiding_vs) {
			auto g_id = md.V_map_reverse[md.QV_map_reverse[vid]];
			tb_subdivided_cells.push_back(g_id);
		}
	std::cout << "line->tb cells: " << tb_subdivided_cells.size() << std::endl;
	if (tb_subdivided_cells.size()) return false;
	//kill miss match because of connecting patches
	vector<bool> V_flag(Qfg.mf.tri.Vs.size(), false);
	for (auto c : Qfg.Cs)V_flag[c.id] = true;
	for (auto l : Qfg.Ls)for (auto lvid : l.vs)V_flag[lvid] = true;
	for (auto rid : regions) {
		for (auto fid : Qfg.Rs[rid].tris)
			for (auto vid : Qfg.mf.tri.Fs[fid].vs)
				if (!V_flag[vid]) {
					auto g_id = md.V_map_reverse[md.QV_map_reverse[vid]];
					tb_subdivided_cells.push_back(g_id);
					V_flag[vid] = true;
				}
	}
	std::cout << "patch->tb cells: " << tb_subdivided_cells.size() << std::endl;
	if (tb_subdivided_cells.size())return false;

	for (auto r : Qfg.Rs) {
		//find its corresponding patches in fg
		if (graph_matches[r.id].size() > 1) {
			cout << r.id << " ambiguous matching!" << endl;
			continue;
		}
	}
	cout << "patch matched!" << endl;
	return true;
}

void mesh_opt::patch_trees(Mesh_Domain &md) {

	auto &graph_matches = md.graph_matches;
	auto &tri_Trees = md.tri_Trees;
	// build aabbtrees
	tri_Trees.clear();
	tri_Trees.resize(graph_matches.size());

	for (uint32_t j = 0; j < graph_matches.size(); j++) {
		auto &R = graph_matches[j];
		vector<bool> V_tag(mf.tri.Vs.size(), false), T_tag(mf.tri.Fs.size(), false);
		for (auto rid : R)
			for (auto tid : fg.Rs[rid].tris) {
				T_tag[tid] = true;
				for (auto vid : mf.tri.Fs[tid].vs)V_tag[vid] = true;
			}
		int vn = 0;
		vector<int> v_map(mf.tri.Vs.size(), -1);
		for (uint32_t i = 0; i < V_tag.size(); i++)if (V_tag[i])v_map[i] = vn++;
		Mesh tri_mesh;
		tri_mesh.V.resize(3, vn);
		vn = 0;
		for (uint32_t i = 0; i < V_tag.size(); i++)if (V_tag[i]) tri_mesh.V.col(vn++) = mf.tri.V.col(i);
		for (uint32_t i = 0; i < T_tag.size(); i++)
			if (T_tag[i]) {
				tri_mesh.Fs.push_back(mf.tri.Fs[i]);
				for (auto &vid : tri_mesh.Fs[tri_mesh.Fs.size() - 1].vs) vid = v_map[vid];
			}

		tri_Trees[j].v_map = v_map;
		build_aabb_tree(tri_mesh, tri_Trees[j]);
	}
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
	ts.T.resize(m.Hs.size() * 8, 4);
	Vector4i t;
	for (auto &h : m.Hs) {
		for (uint32_t i = 0; i < 8; i++) {
			for (uint32_t j = 0; j < 4; j++) t[j] = h.vs[hex_tetra_table[i][j]];
			ts.T.row(h.id * 8 + i) = t;
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

void mesh_opt::projection_smooth(const Mesh &hmi, Feature_Constraints &fc) {
}

bool mesh_opt::smoothing(Mesh_Domain &md) {

	auto &m = md.mesh_entire;
	scaled_jacobian(m, mq);
	std::cout << "before deformation: minimum scaled J: " << mq.min_Jacobian << " average scaled J: " << mq.ave_Jacobian
			  << endl;

	vector<bool> Huntangle_flag(m.Hs.size(), false), H_flag(m.Hs.size(), false), H_inout_tag(m.Hs.size(), true);
	vector<uint32_t> Hids;
	for (uint32_t i = 0; i < m.Hs.size(); i++)Hids.push_back(i);

	ts = Tetralize_Set();
	ts.V = m.V.transpose();
	ts.T.resize(m.Hs.size() * 8, 4);
	Vector4i t;
	for (auto &h : m.Hs) {
		for (uint32_t i = 0; i < 8; i++) {
			for (uint32_t j = 0; j < 4; j++) t[j] = h.vs[hex_tetra_table[i][j]];
			ts.T.row(h.id * 8 + i) = t;
		}
	}
	//smooth feature only
	fc = Feature_Constraints();
//	fc.V_types.resize(m.Vs.size());
//	fill(fc.V_types.begin(), fc.V_types.end(), Feature_V_Type::INTERIOR);
//	fc.V_ids.resize(m.Vs.size());
//	fc.RV_type.resize(m.Vs.size());
//	fill(fc.RV_type.begin(), fc.RV_type.end(), true);
//	fc.lamda_C = fc.lamda_L = fc.lamda_T = 0;

	ts.energy_type = SYMMETRIC_DIRICHLET;
	ts.UV = ts.V;
	ts.fc = fc;
	ts.projection = false;
	ts.global = true;
	ts.glue = false;
	ts.lamda_glue = 0;
	ts.lamda_region = 0;
	ts.record_Sequence = false;
	if (args.scaffold_type == 2 || args.scaffold_type == 3) {
		ts.known_value_post = true;
		//md.post_index = ts.V.rows();
		ts.post_index = md.post_index;
		ts.post_Variables.resize(ts.V.rows() - md.post_index, 3);
		for (int i = md.post_index; i < m.Vs.size(); i++) {
			auto &v = m.Vs[i];
			if (v.boundary)
				ts.post_Variables.row(i - md.post_index) = ts.V.row(i);
		}
	}
	ts.b.resize(0);
	ts.bc.resize(0, 3);
	ts.bc.setZero();
	ts.lamda_b = 0;

	optimization opt;
	opt.weight_opt = weight_opt;
	improve_Quality_after_Untangle_Iter_MAX = 3;
	for (uint32_t i = 0; i < improve_Quality_after_Untangle_Iter_MAX; i++) {

		compute_referenceMesh(ts.V, m.Hs, H_inout_tag, Hids, ts.RT, true);
		projection_smooth(m, fc);
		ts.fc = fc;

		opt.slim_opt_igl(ts, 5);
		m.V = ts.V.transpose();
	}

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

bool mesh_opt::feature_alignment(Mesh_Domain &md) {

	auto &m = md.mesh_entire;
	auto &m_ = md.mesh_subA;

	vector<bool> H_flag(m.Hs.size(), true);
	//hex2tet
	ts = Tetralize_Set();

	ts.V = m.V.transpose();
	ts.T.resize(m.Hs.size() * 8, 4);
	Vector4i t;
	for (auto &h : m.Hs) {
		for (uint32_t i = 0; i < 8; i++) {
			for (uint32_t j = 0; j < 4; j++) t[j] = h.vs[hex_tetra_table[i][j]];
			ts.T.row(h.id * 8 + i) = t;
		}
	}
	vector<uint32_t> Hids;
	for (uint32_t i = 0; i < m.Hs.size(); i++)Hids.push_back(i);

	double MESHRATIO = 1;
	double LAMDA_FEATURE_PROJECTION = MESHRATIO * args.feature_weight;//0.05
	double LAMDA_FEATURE_PROJECTION_BOUND = MESHRATIO * 1e+16;

	ts.b.resize(0);
	ts.bc.resize(0, 3);
	ts.bc.setZero();
	ts.lamda_b = MESHRATIO * 1e+3;
	ts.regionb.resize(0);
	ts.regionbc.resize(0, 3);
	ts.regionbc.setZero();
	ts.lamda_region = 0;
//	fc.lamda_C = fc.lamda_L = fc.lamda_T = LAMDA_FEATURE_PROJECTION;

	ts.energy_type = SYMMETRIC_DIRICHLET;
	ts.UV = ts.V;
	ts.fc = fc;
	ts.projection = false;
	ts.global = true;
	ts.glue = false;
	ts.lamda_glue = 0;
	ts.record_Sequence = false;

	if (args.scaffold_type == 2 || args.scaffold_type == 3) {
		ts.known_value_post = true;
		//md.post_index = ts.V.rows();
		ts.post_index = md.post_index;
		ts.post_Variables.resize(ts.V.rows() - md.post_index, 3);
		for (int i = md.post_index; i < m.Vs.size(); i++) {
			auto &v = m.Vs[i];
			if (v.boundary)
				ts.post_Variables.row(i - md.post_index) = ts.V.row(i);
		}
	}

	optimization opt;
	opt.weight_opt = weight_opt * 0.001;

	Mesh htri;
	htri.type = Mesh_type::Tri;
	extract_surface_conforming_mesh(m_, htri, md.TV_map, md.TV_map_reverse, md.TF_map, md.TF_map_reverse);
	ave_Hausdorff_dises.clear();
	ratio_ave_Hausdorff.clear();
	ratio_max_Hausdorff.clear();

	//dirty feature
	dirty_region_identification(md, Dirty_Vertices);
	//Dirty_Vertices.clear();
	dirty_local_feature_update(md, Dirty_Vertices);
	//loop
	bool Max_dis_satisified = false, Ave_dis_satisfied = false;
	improve_Quality_after_Untangle_Iter_MAX = 30;
	double energy_pre = -1;
	for (uint32_t i = 0; i < improve_Quality_after_Untangle_Iter_MAX; i++) {

		dirty_graph_projection(md, fc, Dirty_Vertices);

		ts.fc = fc;

		compute_referenceMesh(ts.UV, m.Hs, md.H_flag, Hids, ts.RT, true);
		opt.slim_m_opt(ts, 1, -1);
		m.V = ts.UV.transpose();
		ts.V = ts.UV;

		if (!Max_dis_satisified) {
//			std::cout << "energy_quality " << opt.engery_quality << "; soft: " << opt.engery_soft << " lamda_c "
//					  << fc.lamda_C << std::endl;
//			fc.lamda_T = fc.lamda_L = fc.lamda_C = std::min(
//					LAMDA_FEATURE_PROJECTION * std::max(opt.engery_quality / opt.engery_soft * fc.lamda_C, 1.0),
//					LAMDA_FEATURE_PROJECTION_BOUND);
		}

		if (i <= 1)
			energy_pre = opt.energy;
		else if (i > 1 && energy_pre != opt.energy)
			energy_pre = opt.energy;
		else if (i > 1 && energy_pre == opt.energy)
			break;

		if (stop_criterior_satisfied(md, i, mf.tri, htri, Max_dis_satisified, Ave_dis_satisfied)) {

			break;
		}
	}
	for (auto &v : m.Vs) {
		v.v[0] = m.V(0, v.id);
		v.v[1] = m.V(1, v.id);
		v.v[2] = m.V(2, v.id);
	}

	for (auto &v : md.mesh_subA.Vs) {
		v.v = md.mesh_entire.Vs[md.V_map_reverse[v.id]].v;
		md.mesh_subA.V.col(v.id) = md.mesh_entire.V.col(md.V_map_reverse[v.id]);
	}

	scaled_jacobian(m, mq);
	std::cout << "after: V, H, minimum scaled J: " << m.Vs.size() << " " << m.Hs.size() << " " << mq.min_Jacobian
			  << " average scaled J: " << mq.ave_Jacobian << endl;
	scaled_jacobian(m_, mq);
	std::cout << "after: V, H, minimum scaled J: " << m_.Vs.size() << " " << m_.Hs.size() << " " << mq.min_Jacobian
			  << " average scaled J: " << mq.ave_Jacobian << endl;
	mq.ave_hausdorff = ratio_ave_Hausdorff;
	mq.max_hausdorff = ratio_max_Hausdorff;

	double J_bound = 0.0;
	if (mq.min_Jacobian < J_bound) {
		for (int i = 0; i < m.Hs.size(); i++) {
			for (int j = 0; j < 8; j++)
				if (mq.V_Js[i * 8 + j] < J_bound) {
					tb_subdivided_cells.push_back(m.Hs[i].vs[j]);
				}
		}
		return false;
	}
	return true;
}

void mesh_opt::dirty_local_feature_update(Mesh_Domain &md, std::vector<int> &Dirty_Vs) {
	auto &hmi = md.mesh_entire;
	auto &mf_quad = md.mf_quad;

//	fc.V_types.resize(hmi.Vs.size());
//	fill(fc.V_types.begin(), fc.V_types.end(), Feature_V_Type::INTERIOR);
//	fc.V_ids.resize(hmi.Vs.size());
//	fc.RV_type.resize(hmi.Vs.size());
//	fill(fc.RV_type.begin(), fc.RV_type.end(), true);
//	//matrices
//
//	for (uint32_t i = 0; i < mf_quad.corners.size(); i++) {
//		fc.V_types[mf_quad.corners[i]] = Feature_V_Type::CORNER;
//		fc.V_ids[mf_quad.corners[i]] = mf.corners[i];
//	}
//	for (uint32_t i = 0; i < mf_quad.curve_vs.size(); i++)
//		for (auto vid : mf_quad.curve_vs[i]) {
//			if (fc.V_types[vid] == Feature_V_Type::CORNER)continue;
//			fc.V_types[vid] = Feature_V_Type::LINE;
//			fc.V_ids[vid] = i;
//		}
//	for (auto v : hmi.Vs) {
//		auto i = v.id;
//		if (v.on_medial_surface && fc.V_types[i] == Feature_V_Type::INTERIOR) {
//			fc.V_types[i] = Feature_V_Type::REGULAR;
//		}
//	}
//
//	for (const auto &v : Dirty_Vs) {
//		fc.V_types[v] = Feature_V_Type::REGULAR;
//	}
//
//	uint32_t num_corners = 0, num_lines = 0, num_regulars = 0;
//	for (const auto &type : fc.V_types)
//		if (type == Feature_V_Type::CORNER)num_corners++;
//		else if (type == Feature_V_Type::LINE)num_lines++;
//		else if (type == Feature_V_Type::REGULAR)num_regulars++;
//
//	fc.ids_C.resize(num_corners);
//	fc.C.resize(num_corners, 3);
//	fc.num_a = num_lines;
//	fc.ids_L.resize(num_lines);
//	fc.Axa_L.resize(num_lines, 3);
//	fc.origin_L.resize(num_lines, 3);
//	fc.ids_T.resize(num_regulars);
//	fc.normal_T.resize(num_regulars, 3);
//	fc.dis_T.resize(num_regulars);
//	fc.V_T.resize(num_regulars, 3);
//	num_corners = num_lines = num_regulars = 0;
}

void mesh_opt::dirty_region_identification(Mesh_Domain &md, std::vector<int> &Dirty_Vs) {

	auto &Qfg = md.Qfg;
	auto &mf_quad = md.mf_quad;
	auto &QV_map_reverse = md.QV_map_reverse;
	std::vector<int> C_flag(md.mesh_entire.Vs.size(), -1);

	std::vector<bool> V_flag(md.mesh_entire.Vs.size(), false);
	for (const auto &tc : fg.Cs) {

		std::fill(C_flag.begin(), C_flag.end(), -1);
		for (const auto &l : Qfg.Ls)for (const auto &vid : l.vs)C_flag[md.V_map_reverse[QV_map_reverse[vid]]] = l.id;

		std::vector<double> angles;
		int pre = -1, cur = -1;
		bool Has_Tiny = false;
		for (int i = 0; i < tc.ring_vs_tag.size(); i++) {
			if (tc.ring_vs_tag[i] >= 0) {
				if (pre == -1 && cur == -1) {
					pre = tc.ring_vs[i];
				} else if (pre != -1 && cur == -1) {
					cur = tc.ring_vs[i];
				} else if (pre != -1 && cur != -1) {
					pre = cur;
					cur = tc.ring_vs[i];
				}

				if (pre != -1 && cur != -1) {
					Vector3d l0 = (mf.tri.V.col(pre) - mf.tri.V.col(tc.vs[0])).normalized();
					Vector3d l1 = (mf.tri.V.col(cur) - mf.tri.V.col(tc.vs[0])).normalized();

					angles.push_back(std::abs(std::acos(l0.dot(l1))));
					if (angles[angles.size() - 1] < Tiny_angle_threshold) {
						Has_Tiny = true;
						break;
					}
				}
			}
		}
		if (Has_Tiny) {

			for (int i = 0; i < tc.ring_vs_tag.size(); i++) {
				if (tc.ring_vs_tag[i] >= 0)
					for (const auto &vid : Qfg.Ls[tc.ring_vs_tag[i]].vs)C_flag[md.V_map_reverse[QV_map_reverse[vid]]] = -1;
			}

			int qvid = Qfg.Cs[tc.id].vs[0];
			vector<int> lregion_vs;
			lregion_vs.push_back(qvid);

			int N_ring = 3;
			for (int j = 0; j < N_ring; j++) {
				vector<int> lregion_vs_;
				for (int k = 0; k < lregion_vs.size(); k++) {
					const auto &nfs = mf_quad.tri.Vs[lregion_vs[k]].neighbor_fs;
					for (const auto &fid : nfs) {
						const auto &nvs = mf_quad.tri.Fs[fid].vs;
						for (const auto &vid : nvs)
							if (C_flag[md.V_map_reverse[QV_map_reverse[vid]]] != -1)continue;
							else
								V_flag[md.V_map_reverse[QV_map_reverse[vid]]] = true;
						lregion_vs_.insert(lregion_vs_.end(), nvs.begin(), nvs.end());
					}
				}
				lregion_vs = lregion_vs_;
			}
		}
	}
	Dirty_Vs.clear();
	for (int i = 0; i < V_flag.size(); i++)if (V_flag[i])Dirty_Vs.push_back(i);
}

bool mesh_opt::stop_criterior_satisfied(Mesh_Domain &md, const int iter_after_untangle, const Mesh &tm0, Mesh &tm1,
										bool &max_dis_satisified, bool &ave_dis_satisfied) {

	double bbox_diagonal, max_hausdorff_dis, ave_hausdorff_dis;
	auto &hmi = md.mesh_entire;
	for (int i = 0; i < md.TV_map_reverse.size(); i++)tm1.V.col(i) = hmi.V.col(md.V_map_reverse[md.TV_map_reverse[i]]);

	compute(tm0, tm1, bbox_diagonal, max_hausdorff_dis, ave_hausdorff_dis);

	cout << "ave_hausdorff_dis: " << ave_hausdorff_dis << endl;

	ave_Hausdorff_dises.push_back(ave_hausdorff_dis);

	if (iter_after_untangle >= start_ave_hausdorff_count_Iter) {
		double ratio_max = max_hausdorff_dis / bbox_diagonal;
		int cur_ind = (int) ave_Hausdorff_dises.size() - 1;
		if (cur_ind < 1) return false;

		double ratio_ave =
				(ave_Hausdorff_dises[cur_ind - 1] - ave_Hausdorff_dises[cur_ind]) / ave_Hausdorff_dises[cur_ind];

		cout << "ratio_max: " << ratio_max << endl;
		cout << "ratio_ave: " << ratio_ave << endl;

		max_HR = ratio_max;

		ratio_ave_Hausdorff.push_back(ratio_ave);
		ratio_max_Hausdorff.push_back(ratio_max);

		if (ratio_max < 0.8 * HR) max_dis_satisified = true;
		else max_dis_satisified = false;
		if (ratio_ave < 0.8 * STOP_AVE_HAUSDORFF_THRESHOLD) ave_dis_satisfied = true;
		else ave_dis_satisfied = false;
		if (ratio_ave < STOP_AVE_HAUSDORFF_THRESHOLD && ratio_max < HR) return true;
		if (ave_dis_satisfied);// return true;

	}
	return false;
}

void mesh_opt::dirty_graph_projection(Mesh_Domain &md, Feature_Constraints &fc, std::vector<int> &Dirty_Vs) {

	auto &Qfg = md.Qfg;
	auto &hmi = md.mesh_entire;
	auto &tri_Trees = md.tri_Trees;

	uint32_t num_corners = 0, num_lines = 0, num_regulars = 0;
//	for (uint32_t i = 0; i < fc.V_types.size(); i++) {
//		if (fc.V_types[i] == Feature_V_Type::CORNER) {
//			fc.ids_C[num_corners] = i;
//			fc.C.row(num_corners) = mf.tri.V.col(fc.V_ids[i]);
//			num_corners++;
//		} else if (fc.V_types[i] == Feature_V_Type::LINE) {
//			Vector3d v = hmi.V.col(i).transpose();
//			uint32_t curve_id = fc.V_ids[i];
//			vector<uint32_t> &curve = mf.curve_vs[curve_id];
//			Vector3d tangent(1, 0, 0);
//			uint32_t curve_len = curve.size();
//
//			if (!mf.circles[curve_id]) curve_len--;
//
//			vector<Vector3d> pvs, tangents;
//			vector<pair<double, uint32_t>> dis_ids;
//			Vector3d pv;
//			for (uint32_t j = 0; j < curve_len; j++) {
//				uint32_t pos_0 = curve[j], pos_1 = curve[(j + 1) % curve.size()];
//				double t, precision_here = 1.0e1;
//				point_line_projection(mf.tri.V.col(pos_0), mf.tri.V.col(pos_1), v, pv, t);
//				tangent = (mf.tri.V.col(pos_1) - mf.tri.V.col(pos_0)).normalized();
//				dis_ids.push_back(make_pair((v - pv).norm(), pvs.size()));
//				pvs.push_back(pv);
//				tangents.push_back(tangent);
//			}
//			sort(dis_ids.begin(), dis_ids.end());
//
//			if (dis_ids.size()) {
//				uint32_t cloestid = dis_ids[0].second;
//				pv = pvs[cloestid];
//				tangent = tangents[cloestid];
//			} else {
//				//brute-force search
//				for (uint32_t j = 0; j < curve.size(); j++) {
//					double dis = (mf.tri.V.col(curve[j]) - v).norm();
//					dis_ids.push_back(make_pair(dis, j));
//				}
//				sort(dis_ids.begin(), dis_ids.end());
//
//				int pos = dis_ids[0].second;
//				pv = mf.tri.V.col(curve[pos]);
//
//				curve_len = curve.size();
//				if (mf.circles[curve_id] || (!mf.circles[curve_id] && pos != 0 && pos != curve_len - 1)) {
//					uint32_t pos_0 = (pos - 1 + curve_len) % curve_len, pos_1 = (pos + 1) % curve_len;
//					tangent += (mf.tri.V.col(curve[pos]) - mf.tri.V.col(curve[pos_0])).normalized();
//					tangent += (mf.tri.V.col(curve[pos_1]) - mf.tri.V.col(curve[pos])).normalized();
//				} else if (!mf.circles[curve_id] && pos == 0) {
//					uint32_t pos_1 = (pos + 1) % curve_len;
//					tangent += (mf.tri.V.col(curve[pos_1]) - mf.tri.V.col(curve[pos])).normalized();
//				} else if (!mf.circles[curve_id] && pos == curve_len - 1) {
//					uint32_t pos_0 = (pos - 1 + curve_len) % curve_len;
//					tangent += (mf.tri.V.col(curve[pos]) - mf.tri.V.col(curve[pos_0])).normalized();
//				}
//				tangent.normalize();
//			}
//			fc.ids_L[num_lines] = i;
//			fc.origin_L.row(num_lines) = pv;
//			fc.Axa_L.row(num_lines) = tangent;
//			num_lines++;
//		}
//	}
	//cout << "graph_projection patch projection" << endl;
	vector<bool> V_flag(hmi.Vs.size(), false);
	num_regulars = 0;

	for (const auto &vid : Dirty_Vs) V_flag[vid] = true;
	std::cout << "dirty_vs size: " << Dirty_Vs.size() << std::endl;
	if (Dirty_Vs.size()) {
		MatrixXd Ps(Dirty_Vs.size(), 3);
		int num_regulars_ = 0;
		for (const auto &vid : Dirty_Vs) {
//			fc.ids_T(num_regulars + num_regulars_) = vid;
			Ps.row(num_regulars_++) = hmi.V.col(vid).transpose();
		}

		VectorXd signed_dis;
		VectorXi ids;
		MatrixXd V_T(Dirty_Vs.size(), 3), normal_T(Dirty_Vs.size(), 3);
		signed_distance_pseudonormal(Ps, tri_tree.TriV, tri_tree.TriF, tri_tree.tree, tri_tree.TriFN, tri_tree.TriVN,
									 tri_tree.TriEN,
									 tri_tree.TriEMAP, signed_dis, ids, V_T, normal_T);

		for (uint32_t j = 0; j < Dirty_Vs.size(); j++) {
//			fc.normal_T.row(num_regulars + j) = normal_T.row(j);
//			fc.V_T.row(num_regulars + j) = V_T.row(j);
//			fc.dis_T[num_regulars + j] = normal_T.row(j).dot(V_T.row(j));
//			fc.V_ids[Dirty_Vs[j]] = -1;
		}

		num_regulars = Dirty_Vs.size();
	}

	for (uint32_t i = 0; i < Qfg.Rs.size(); i++) {
		auto &R = Qfg.Rs[i];
		vector<int> vs;
		for (const auto &vid : R.vs)if (!V_flag[vid]) vs.push_back(vid);

		if (!vs.size()) {
			continue;
		}

		MatrixXd Ps(vs.size(), 3);
		int num_regulars_ = 0;
		for (auto vid : vs) {
//			fc.ids_T(num_regulars + num_regulars_) = vid;
			Ps.row(num_regulars_++) = hmi.V.col(vid).transpose();
		}

		VectorXd signed_dis;
		VectorXi ids;
		MatrixXd V_T(vs.size(), 3), normal_T(vs.size(), 3);
		signed_distance_pseudonormal(Ps, tri_Trees[i].TriV, tri_Trees[i].TriF, tri_Trees[i].tree, tri_Trees[i].TriFN,
									 tri_Trees[i].TriVN, tri_Trees[i].TriEN,
									 tri_Trees[i].TriEMAP, signed_dis, ids, V_T, normal_T);

		for (uint32_t j = 0; j < vs.size(); j++) {
//			fc.normal_T.row(num_regulars + j) = normal_T.row(j);
//			fc.V_T.row(num_regulars + j) = V_T.row(j);
//			fc.dis_T[num_regulars + j] = normal_T.row(j).dot(V_T.row(j));
//			fc.V_ids[vs[j]] = -1;
		}
		num_regulars += vs.size();
	}
}

