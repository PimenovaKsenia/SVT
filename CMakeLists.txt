cmake_minimum_required(VERSION 3.10)
project(VectorFieldOnMesh)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(INMOST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../2/INMOST-install")

find_package(INMOST REQUIRED)

include_directories(${INMOST_INCLUDE_DIRS})
link_directories(${INMOST_LIBRARY_DIRS})
add_definitions(${INMOST_DEFINITIONS})

add_executable(compute_tangent mesh.cpp)
target_link_libraries(compute_tangent ${INMOST_LIBRARIES})
