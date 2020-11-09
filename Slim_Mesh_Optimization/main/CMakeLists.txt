file(GLOB header *.h)
file(GLOB source *.cpp)

set(DATA_SHARED_PATH ${PROJECT_SOURCE_DIR}/data CACHE PATH "location of shared data resources")

add_executable(slim_mesh_opt_exe ${source} ${header} ${SRCFILES})
target_link_libraries(slim_mesh_opt_exe
        boost_program_options
        ${Boost_LIBRARIES}
        io
        )

target_compile_definitions(slim_mesh_opt_exe PUBLIC "-DDATA_SHARED_PATH=\"${DATA_SHARED_PATH}\"")
target_include_directories(slim_mesh_opt_exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

