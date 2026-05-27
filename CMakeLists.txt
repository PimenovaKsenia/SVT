cmake_minimum_required(VERSION 3.10)
project(Poisson2D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Путь к установленной библиотеке INMOST (измените при необходимости)
set(INMOST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/INMOST-install")

find_package(INMOST REQUIRED)

include_directories(${INMOST_INCLUDE_DIRS})
add_executable(solver Poisson_solver.cpp)
target_link_libraries(solver ${INMOST_LIBRARIES})
