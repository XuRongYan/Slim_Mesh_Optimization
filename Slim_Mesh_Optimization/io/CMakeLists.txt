file(GLOB HEADERS *.h *.hpp)
file(GLOB SOURCES *.cpp)

if (UNIX)
    add_library(io SHARED ${SOURCES} ${HEADERS})
elseif (WIN32)
    add_library(io STATIC ${SOURCES} ${HEADERS})
elseif (APPLE)
    add_library(io SHARED ${SOURCES} ${HEADERS})
endif ()

target_link_libraries(io
        boost_program_options
        ${Boost_LIBRARIES}
        ${LIBIGL_LIBRARIES}
        ${LIBIGL_EXTRA_LIBRARIES}
        )