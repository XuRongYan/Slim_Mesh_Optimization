//
// Created by 徐溶延 on 2020/11/8.
//

#ifndef SLIM_MESH_OPTIMIZATION_PARSE_CONFIG_FILE_H
#define SLIM_MESH_OPTIMIZATION_PARSE_CONFIG_FILE_H
#include <iostream>
#include "../common/global_types.h"

int parse_config_file(const std::string &filename,
					  common::arguments &args);


#endif //SLIM_MESH_OPTIMIZATION_PARSE_CONFIG_FILE_H
