cmake_minimum_required (VERSION 3.8)

project (Slim_Mesh_Optimization)

file(GLOB SRCFILES *.cpp)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(MSVC)
    # Enable parallel compilation
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj /openmp")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
    # Libigl requires a modern C++ compiler that supports c++11
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
endif()

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)

include(geogram.cmake)
set_property(TARGET geogram PROPERTY FOLDER "dependencies")

# Compilation flags: adapt to your needs
if(MSVC)
    # Enable parallel compilation
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj /openmp")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
    # Libigl requires a modern C++ compiler that supports c++11
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
endif()

if(WIN32)
    LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}") # Add the directory where FindSuiteSparse.cmake module can be
    set(SuiteSparse_USE_LAPACK_BLAS ON)
    find_package(SuiteSparse QUIET NO_MODULE)  # 1st: Try to locate the *config.cmake file.
    if(NOT SuiteSparse_FOUND)
        find_package(SuiteSparse REQUIRED) # 2nd: Use FindSuiteSparse.cmake module
        include_directories(${SuiteSparse_INCLUDE_DIRS})
    else()
        message(STATUS "Find SuiteSparse : include(${USE_SuiteSparse})")
        include(${USE_SuiteSparse})
    endif()
    MESSAGE(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")
    # ------------------------------------------------------------------
    #   End of SuiteSparse detection
    # ------------------------------------------------------------------
else()
    # ------------------------------------------------------------------
    # Detect SuiteSparse libraries:
    # If not found automatically, set SuiteSparse_DIR in CMake to the
    # directory where SuiteSparse was built.
    # ------------------------------------------------------------------
    LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/..") # Add the directory where FindSuiteSparse.cmake module can be found.

    find_package(SuiteSparse REQUIRED)  # 1st: Try to locate the *config.cmake file.
    if(CSPARSE_FOUND AND CHOLMOD_FOUND)
        include_directories(${CHOLMOD_INCLUDE_DIR} ${CSPARSE_INCLUDE_DIR})
    else()
        #message(FATAL_ERROR "UNABLE to Find SuiteSpars, pass SUITE_SPARSE_ROOT to cmake")
        message("UNABLE to Find SuiteSpars, pass SUITE_SPARSE_ROOT to cmake")
    endif()
    MESSAGE(STATUS "SuiteSparse_LIBS: ${CHOLMOD_LIBRARY} ${CHOLMOD_INCLUDE_DIR}")
    # ------------------------------------------------------------------
    #   End of SuiteSparse detection
    # ------------------------------------------------------------------
endif()

find_package(Boost COMPONENTS system thread)
if (Boost_FOUND)
    message("find boost by find_package(), " PATH=${Boost_INCLUDE_DIRS})
    include_directories(${Boost_INCLUDE_DIRS})
else (Boost_FOUND)
    if (UNIX)
        message("find boost by UNIX, " PATH=$ENV{HOME}/boost/)
        link_directories($ENV{HOME}/boost/)
    elseif (WIN32)
        message("find boost by WIN32")
        include_directories($ENV{HOME}/boost/boost_1_58_0)
        link_directories($ENV{HOME}/boost/boost_1_58_0/lib64-msvc-12.0)
        include_directories($ENV{HOME}/boost_1_64_0)
        link_directories($ENV{HOME}/boost_1_64_0/lib64-msvc-14.1)
    elseif (APPLE)
        message("find boost by APPLE, " PATH=$ENV{HOME})
    endif ()
endif ()

find_package(tbb)
if (tbb_FOUND)
    message("find tbb by find_package(), " PATH=${tbb_INCLUDE_DIRS})
    include_directories(${tbb_INCLUDE_DIRS})
    link_libraries(/usr/local/lib/libtbb.a)
else(tbb_FOUND)
    message("find tbb at /usr/local/include")
    include_directories(/usr/local/include)
    link_libraries(/usr/local/lib/libtbb.a)
endif()

# libigl options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)

# add a customizable menu bar
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   ON)

# libigl options: choose your dependencies (by default everything is OFF except opengl)
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_GLFW        "Use GLFW"           ON)
option(LIBIGL_WITH_BBW         "Use BBW"            OFF)
option(LIBIGL_WITH_EMBREE      "Use Embree"         ON)
option(LIBIGL_WITH_PNG         "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       ON)
option(LIBIGL_WITH_XML         "Use XML"            OFF)
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF) # This option is not supported yet
option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF) # This option is not supported yet
option(LIBIGL_WITH_CGAL        "Use CGAL"           OFF)
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, the cgal build system fails without it
    find_package(CGAL REQUIRED)
    set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
    include(${CGAL_USE_FILE})
endif()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/libigl/shared/cmake" "libigl")

# Build VCG
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/vcg external_build/vcg)
set_property(TARGET vcg PROPERTY FOLDER "dependencies")

# Prepare the build environment
include_directories(
        ${LIBIGL_INCLUDE_DIRS}
        # Intel Thread Building Blocks
        ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb/include
        # VCG
        ${CMAKE_CURRENT_SOURCE_DIR}/external/vcg
        ${CMAKE_CURRENT_SOURCE_DIR}/external/geogram/src/lib
)
add_definitions(${LIBIGL_DEFINITIONS})
add_definitions(-DENABLE_SERIALIZATION)

add_definitions(-DNOMINMAX)

# libigl information
message("libigl includes: ${LIBIGL_INCLUDE_DIRS}")
message("libigl libraries: ${LIBIGL_LIBRARIES}")
message("libigl extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libigl extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libigl definitions: ${LIBIGL_DEFINITIONS}")

# Build geogram
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
include(geogram.cmake)
set_property(TARGET geogram PROPERTY FOLDER "dependencies")

if (APPLE)
    INCLUDE_DIRECTORIES(
            /usr/local/include
    )

else ()
    INCLUDE_DIRECTORIES(
            $ENV{HOME}/usr/include/
    )
endif ()


if (APPLE)
    link_directories(
            /usr/local/lib
    )

else ()
    link_directories(
            $ENV{HOME}/usr/lib
            $ENV{HOME}/usr/Linux/64/gcc4/lib/
    )
endif ()


# TODO: 如有需要，请添加测试并安装目标。
add_subdirectory(common)
add_subdirectory(io)
add_subdirectory(main)
