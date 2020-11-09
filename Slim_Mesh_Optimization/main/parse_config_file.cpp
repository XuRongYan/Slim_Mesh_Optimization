//
// Created by 徐溶延 on 2020/11/8.
//

#include "parse_config_file.h"
#include <fstream>

#include <boost/program_options.hpp>

using namespace boost;


int parse_config_file(const std::string &filename, common::arguments &args) {
	std::ifstream ifs(filename);
	if (!ifs.is_open()) {
		return __LINE__;
	}
	program_options::options_description config("command parameters");
	config.add_options()
			("help, h", "help message")
			("choice, s", program_options::value<std::string>(), "functionality choice")
			("in, s", program_options::value<std::string>(), "Input mesh")
			("out, s", program_options::value<std::string>(), "Output mesh")
			("octree, i", program_options::value<int>()->default_value(1), "whether using octree meshing")
			("num_cell, i", program_options::value<int>()->default_value(30), "num_cells for voxel meshing.")
			("hausdorff, d", program_options::value<double>()->default_value(0.005), "Hausdorff_ratio_t")
			("edge_length_ratio, d", program_options::value<double>()->default_value(15), "edge_length_ratio")
			("weight_opt, d", program_options::value<double>()->default_value(1), "weight_opt")
			("feature_weight, d", program_options::value<double>()->default_value(0.05), "feature_weight_opt")
			("pca_oobb, i", program_options::value<int>()->default_value(1), "pca_oobb")
			("scaffold_type, i", program_options::value<int>()->default_value(1), "scaffold_type")
			("hard_feature, i", program_options::value<int>()->default_value(1), "Hard_Feature")
			("iteration_base, i", program_options::value<int>()->default_value(3), "optimization Iteration_Base");
	program_options::variables_map vm;
	program_options::store(program_options::parse_config_file(ifs, config), vm);
	ifs.close();
	program_options::notify(vm);
	if (vm.count("help")) {
		std::cout << config << std::endl;
		return __LINE__;
	}
	args.choice = vm["choice"].as<std::string>();
	args.input = DATA_SHARED_PATH "/objs/" + vm["in"].as<std::string>();
	args.output = DATA_SHARED_PATH "/output/" + vm["out"].as<std::string>();
	args.octree = vm["octree"].as<int>();
	args.num_cell = vm["num_cell"].as<int>();
	args.Hausdorff_ratio_t = vm["hausdorff"].as<double>();
	args.edge_length_ratio = vm["edge_length_ratio"].as<double>();
	args.weight_opt = vm["weight_opt"].as<double>();
	args.feature_weight = vm["feature_weight"].as<double>();
	args.pca_oobb = vm["pca_oobb"].as<int>();
	args.scaffold_type = vm["scaffold_type"].as<int>();
	args.Hard_Feature = vm["hard_feature"].as<int>();
	args.Iteration_Base = vm["iteration_base"].as<int>();

#ifndef NDEBUG
	std::cout << args << std::endl;
#endif

	return 0;
}
