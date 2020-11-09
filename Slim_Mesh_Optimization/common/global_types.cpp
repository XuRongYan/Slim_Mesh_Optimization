#include "global_types.h"


namespace common {

	arguments args;
	//meshes
	//Mesh hyrbid_mesh;
	Mesh_Feature mf;
	GEO::Mesh M_i;
	Feature_Graph fg;
	Feature_Graph Q_final_fg;
	vector<vector<uint32_t>> GRAPH_MATCHES;


	char path_out_[300];
	int32_t GRAIN_SIZE;

	double diagonal_len;
}

