//
// Created by 徐溶延 on 2020/11/8.
//

#include <iostream>
#include <string>

#include "parse_config_file.h"
#include "../timer.h"
#include "../io/io.h"
#include "../alg/meshing.h"
#include "../alg/mesh_opt.h"




using namespace std;

h_io io;

int grid_method(common::arguments &args);

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "can't execute without configure file." << endl;
	}

	// parse configure file
	common::arguments args;
	parse_config_file(DATA_SHARED_PATH "/config/" + string(argv[1]), args);
	grid_method(args);
	return 1;
}

int grid_method(common::arguments &args) {
	Mesh meshIn, meshOut;
	size_t suffix_idx = args.input.rfind('.');
	std::string path_without_suffix = args.input.substr(0, suffix_idx);
	std::string input_format = args.input.substr(suffix_idx);

	GEO::initialize();

	if (args.choice == "HEX") {
		meshIn.type = Mesh_type::Hex;
	} else if (args.choice == "TET") {
		meshIn.type = Mesh_type::Tet;
	} else if (args.choice == "TRI") {
		meshIn.type = Mesh_type::Tri;
	} else {
		cout << "not support this mesh type" << endl;
	}

	if (input_format == ".obj" || input_format == ".OBJ") {
		io.read_hybrid_mesh_OBJ(meshIn, args.input);
	} else if (input_format == ".vtk" || input_format == ".VTK") {
		io.read_hybrid_mesh_VTK(meshIn, args.input);
	} else if (input_format == ".mesh" || input_format == ".MESH") {
		io.read_hybrid_mesh_MESH(meshIn, args.input);
	} else {
		cout << "not support such input format!" << endl;
		return 0;
	}

	build_connectivity(meshIn);

	for (size_t i = 0; i < meshIn.Vs.size(); i++) {
		meshIn.Vs[i].v.resize(3);
		meshIn.Vs[i].v[0] = meshIn.V(0, i);
		meshIn.Vs[i].v[1] = meshIn.V(1, i);
		meshIn.Vs[i].v[2] = meshIn.V(2, i);
	}

	int nprocess = -1;
	tbb::task_scheduler_init init(nprocess == -1 ? tbb::task_scheduler_init::automatic : nprocess);
	Eigen::initParallel();
	int n = 1;
	Eigen::setNbThreads(n);
	meshing m;

	Timer<> timer;
	timer.beginStage("START MESH OPT");
	cout << endl;

	mesh_opt meshOpt(args);
	meshOpt.hex_mesh_opt(meshIn, meshOut);

	timer.endStage("END MESH OPT");
	std::cout << "TIMING: " << timer.value() << "ms" << endl;

	if (args.output.empty()) {
		string path_out;
		path_out = path_without_suffix + ".mesh";
		cout << "writing mesh to path_out: " << path_out << endl;
		io.write_hybrid_mesh_MESH(meshOut, path_out);
		path_out = path_without_suffix + ".vtk";
		cout << "writing mesh to path_out: " << path_out << endl;
		io.write_hybrid_mesh_VTK(meshOut, path_out);

		path_out = path_without_suffix + "_subB.mesh";
		cout << "writing mesh to path_out: " << path_out << endl;
		io.write_hybrid_mesh_MESH(meshOut, path_out);
		path_out = path_without_suffix + "_subB.vtk";
		cout << "writing mesh to path_out: " << path_out << endl;
		io.write_hybrid_mesh_VTK(meshOut, path_out);
	} else {
		string out_format;
		suffix_idx = args.output.rfind('.');
		out_format = args.output.substr(suffix_idx);
		if (out_format == ".vtk" || out_format == ".VTK"){
			io.write_hybrid_mesh_VTK(meshOut, args.output);
			cout << "writing mesh to path_out: " << args.output << endl;
		}
		else if (out_format == ".mesh" || out_format == ".MESH") {
			path_without_suffix = args.output.substr(0, suffix_idx);
			string path_out;
			path_out = path_without_suffix + ".mesh";
			cout << "writing mesh to path_out: " << path_out << endl;
			io.write_hybrid_mesh_MESH(meshOut, path_out);
			path_out = path_without_suffix + ".vtk";
			cout << "writing mesh to path_out: " << path_out << endl;
			io.write_hybrid_mesh_VTK(meshOut, path_out);

			path_out = path_without_suffix + "_subB.mesh";
			cout << "writing mesh to path_out: " << path_out << endl;
			io.write_hybrid_mesh_MESH(meshOut, path_out);
			path_out = path_without_suffix + "_subB.vtk";
			cout << "writing mesh to path_out: " << path_out << endl;
			io.write_hybrid_mesh_VTK(meshOut, path_out);
		} else {
			return false;
		}
	}
	return 1;
}

