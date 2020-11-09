//
// Created by 徐溶延 on 2020/11/8.
//

#include <iostream>
#include <string>

#include "parse_config_file.h"
#include "../timer.h"
#include "../io/io.h"




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

	if (args.choice == "GRID") {
		grid_method(args);
	}
	return 1;
}

int grid_method(common::arguments &args) {
	size_t suffix_idx = args.input.rfind('.');
	std::string path_without_suffix = args.input.substr(0, suffix_idx);
	std::string input_format = args.input.substr(suffix_idx);
	if (!(input_format == ".obj" || input_format == ".OBJ")) {
		cout << "not support such input format!" << endl;
		return 0;
	}
	// 最后输出的mesh
	Mesh mesh_out;

	Timer<> timer;
	timer.beginStage("START MESH OPT");
	cout << endl;

	// TODO: Mesh Opt

	timer.endStage("END MESH OPT");
	std::cout << "TIMING: " << timer.value() << "ms" << endl;

	if (args.output.empty()) {
		string path_out;
		path_out = path_without_suffix + ".mesh";
		cout << "writing mesh to path_out: " << path_out << endl;
		io.write_hybrid_mesh_MESH(mesh_out, path_out);
		path_out = path_without_suffix + ".vtk";
		cout << "writing mesh to path_out: " << path_out << endl;
		io.write_hybrid_mesh_VTK(mesh_out, path_out);

		path_out = path_without_suffix + "_subB.mesh";
		cout << "writing mesh to path_out: " << path_out << endl;
		io.write_hybrid_mesh_MESH(mesh_out, path_out);
		path_out = path_without_suffix + "_subB.vtk";
		cout << "writing mesh to path_out: " << path_out << endl;
		io.write_hybrid_mesh_VTK(mesh_out, path_out);
	} else {
		string out_format;
		suffix_idx = args.output.rfind('.');
		out_format = args.output.substr(suffix_idx);
		if (out_format == ".vtk" || out_format == ".VTK")
			io.write_hybrid_mesh_VTK(mesh_out, args.output);
		else if (out_format == ".mesh" || out_format == ".MESH") {
			path_without_suffix = args.output.substr(0, suffix_idx);
			string path_out;
			path_out = path_without_suffix + ".mesh";
			cout << "writing mesh to path_out: " << path_out << endl;
			io.write_hybrid_mesh_MESH(mesh_out, path_out);
			path_out = path_without_suffix + ".vtk";
			cout << "writing mesh to path_out: " << path_out << endl;
			io.write_hybrid_mesh_VTK(mesh_out, path_out);

			path_out = path_without_suffix + "_subB.mesh";
			cout << "writing mesh to path_out: " << path_out << endl;
			io.write_hybrid_mesh_MESH(mesh_out, path_out);
			path_out = path_without_suffix + "_subB.vtk";
			cout << "writing mesh to path_out: " << path_out << endl;
			io.write_hybrid_mesh_VTK(mesh_out, path_out);
		} else {
			return false;
		}
	}
	return 1;
}

