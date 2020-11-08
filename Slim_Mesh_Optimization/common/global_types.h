//
// Created by 徐溶延 on 2020/11/8.
//

#ifndef SLIM_MESH_OPTIMIZATION_GLOBAL_TYOES_H
#define SLIM_MESH_OPTIMIZATION_GLOBAL_TYOES_H

#endif //SLIM_MESH_OPTIMIZATION_GLOBAL_TYOES_H

#include <string>
#include <ostream>

namespace common {
	using namespace std;

	struct arguments {
		string choice;
		double Hausdorff_ratio_t = 0.005;
		double edge_length_ratio = 15;
		double edge_length = 0;
		int Iteration_Base = 3;
		bool Hard_Feature = true;
		string input = "", output = "";
		bool octree = true;
		bool whole_domain = false;
		int num_cell = 30;
		int scaffold_layer = 3;
//
		bool pca_oobb = true;
		int scaffold_type = 1;//1 box scaffold free boundary; 2 layered scaffold free boundary; 3 box scaffold fixed boundary; 4 layered scaffold fixed boundary.
		double weight_opt = 1;
		double feature_weight = 0.05;

		friend ostream &operator<<(ostream &os, const arguments &arguments) {
			os << "choice: " << "\t" << arguments.choice << " Hausdorff_ratio_t: " << arguments.Hausdorff_ratio_t
			   << std::endl
			   << " edge_length_ratio: " << "\t" << arguments.edge_length_ratio << " edge_length: "
			   << arguments.edge_length << std::endl
			   << " Iteration_Base: " << "\t" << arguments.Iteration_Base << " Hard_Feature: " << arguments.Hard_Feature
			   << std::endl
			   << " input: " << "\t" << arguments.input << " output: " << arguments.output << " octree: "
			   << arguments.octree << std::endl
			   << " whole_domain: " << "\t" << arguments.whole_domain << " num_cell: " << arguments.num_cell
			   << std::endl
			   << " scaffold_layer: " << "\t" << arguments.scaffold_layer << " pca_oobb: " << arguments.pca_oobb
			   << std::endl
			   << " scaffold_type: " << "\t" << arguments.scaffold_type << " weight_opt: " << arguments.weight_opt
			   << std::endl
			   << " feature_weight: " << "\t" << arguments.feature_weight << std::endl;
			return os;
		}
	};
}