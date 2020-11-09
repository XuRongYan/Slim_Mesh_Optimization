file(GLOB HEADERS *.h *.hpp)
file(GLOB SOURCES *.cpp)

if (UNIX)
    add_library(common SHARED ${SOURCES} ${HEADERS})
elseif (WIN32)
    add_library(common STATIC ${SOURCES} ${HEADERS})
elseif (APPLE)
    add_library(common SHARED ${SOURCES} ${HEADERS})
endif ()

target_link_libraries(common PUBLIC
        geogram
        vcg
        ${LIBIGL_LIBRARIES}
        ${LIBIGL_EXTRA_LIBRARIES}
        )