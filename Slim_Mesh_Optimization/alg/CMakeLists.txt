file(GLOB HEADERS *.h *.hpp grid_meshing/*.h grid_meshing/*.hpp)
file(GLOB SOURCES *.cpp grid_meshing/*.cpp)

if (UNIX)
    add_library(alg SHARED ${SOURCES} ${HEADERS})
elseif (WIN32)
    add_library(alg STATIC ${SOURCES} ${HEADERS})
elseif (APPLE)
    add_library(alg SHARED ${SOURCES} ${HEADERS})
endif ()

target_link_libraries(alg PUBLIC
        geogram
        vcg
        ${LIBIGL_LIBRARIES}
        ${LIBIGL_EXTRA_LIBRARIES}
        common
        io
        )