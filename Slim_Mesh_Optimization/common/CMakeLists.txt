file(GLOB HEADERS *.h *.hpp)
file(GLOB SOURCES *.cpp)

# Build geogram
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
include(../geogram.cmake)
set_property(TARGET geogram PROPERTY FOLDER "dependencies")

if (UNIX)
    add_library(common SHARED ${SOURCES} ${HEADERS})
elseif (WIN32)
    add_library(common STATIC ${SOURCES} ${HEADERS})
elseif (APPLE)
    add_library(common SHARED ${SOURCES} ${HEADERS})
endif ()

target_link_libraries(common PUBLIC
        boost_program_options
        ${Boost_LIBRARIES}
        geogram
        vcg
        ${LIBIGL_LIBRARIES}
        ${LIBIGL_EXTRA_LIBRARIES}
        )