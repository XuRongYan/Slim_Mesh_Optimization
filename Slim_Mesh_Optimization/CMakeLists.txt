# CMakeList.txt: Slim_Mesh_Optimization 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("Slim_Mesh_Optimization")

# 将源代码添加到此项目的可执行文件。
add_executable (Slim_Mesh_Optimization "Slim_Mesh_Optimization.cpp" "Slim_Mesh_Optimization.h")

# TODO: 如有需要，请添加测试并安装目标。
