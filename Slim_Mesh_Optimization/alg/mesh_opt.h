//
// Created by 徐溶延 on 2020/11/10.
//

#ifndef SLIM_MESH_OPTIMIZATION_MESH_OPT_H
#define SLIM_MESH_OPTIMIZATION_MESH_OPT_H

#include "optimization.h"
#include "../../common/global_functions.h"
#include "../../io/io.h"
#include "grid_meshing/octree.h"
#include "metro_hausdorff.h"
using namespace common;

class mesh_opt {
	enum NV_type {
		IV = 0,//inner vertex
		EV,//edge vertex
		DEV,//edge vertex shared by two base cells
		NRV,//new regular edge vertex
		NV//not new vertex
	};

public:
	explicit mesh_opt(const arguments &args);

	void hex_mesh_opt(const Mesh &hex_mesh_in, Mesh &hex_mesh_out);

private:
	//deformation
	bool deformation(Mesh_Domain &md);

public:
#define	untangling_Iter_DEFAULT 30;

	int num_voxels = 1048576;
	bool simplification_bool = false;
	double HR = 0.005;
	double weight_opt = 1;
	int STOP_EXTENT_MIN;
	int STOP_EXTENT_MAX;
	vector<int> tb_subdivided_cells;
	vector<int> hex2Octree_map;
	int local_refinement_ringN = 2;
	int region_Size_untangle = 2;
	double Ratio_grow = 0.1;
	int start_ave_hausdorff_count_Iter = 5;
	int improve_Quality_after_Untangle_Iter_MAX = 10;
	std::vector<double>ave_Hausdorff_dises;
	std::vector<double>ratio_max_Hausdorff;
	std::vector<double>ratio_ave_Hausdorff;
	double STOP_AVE_HAUSDORFF_THRESHOLD = 0.01;
	bool re_Octree_Meshing = false;

	double hausdorff_ratio_threshould = 0.01;
	double hausdorff_ratio = 0;
	OctreeGrid octree;
	h_io io;
	Mesh_Quality mq;
	Mesh_Domain MD;
	double voxel_size = 1;
	bool graded = false;
	bool paired = true;
	//cleaning
	vector<NV_type> nv_TAG;
	//build AABB tree for a triangle mesh
	Treestr tri_tree;
	//optimization
	Feature_Constraints fc;
	Tetralize_Set ts;
	double Tiny_angle_threshold = 30.0 / 180.0 * PAI;
	std::vector<int> Dirty_Vertices;
	std::vector<int> Dirty_Quads;
	//local padding
	vector<vector<uint32_t>> Q_regions;
	vector<vector<uint32_t>> V_layers, H_Tlayers, H_Nlayers;
	vector<int32_t> V2V_Map_Reverse_final;//
	//pocket
	vector<int> Pocket_vs;
	//target volume of reference shape
	double VOLUME_input;
	vector<double> Target_Vols;
	//for debugging
	int file_id = 0;
	double max_HR = 0;
};


#endif //SLIM_MESH_OPTIMIZATION_MESH_OPT_H
