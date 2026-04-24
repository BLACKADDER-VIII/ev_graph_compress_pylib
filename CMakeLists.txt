cmake_minimum_required(VERSION 3.15)
project(nd_compress CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)

# Request shared libs only to avoid LTO-compiled static igraph causing
# R_X86_64_TPOFF32 relocation errors on its hidden TLS symbols.
set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
pkg_check_modules(IGRAPH REQUIRED igraph)
pkg_check_modules(LIBXML2 REQUIRED libxml-2.0)

pybind11_add_module(nd_compress compress_bindings.cpp)

target_include_directories(nd_compress PRIVATE
    ${IGRAPH_INCLUDE_DIRS}
    ${LIBXML2_INCLUDE_DIRS}
)

target_compile_options(nd_compress PRIVATE
    ${IGRAPH_CFLAGS_OTHER}
    ${LIBXML2_CFLAGS_OTHER}
    -O2
)

target_link_directories(nd_compress PRIVATE
    ${IGRAPH_LIBRARY_DIRS}
    ${LIBXML2_LIBRARY_DIRS}
)

target_link_libraries(nd_compress PRIVATE
    OpenMP::OpenMP_CXX
    ${IGRAPH_LIBRARIES}
    ${LIBXML2_LIBRARIES}
)

install(TARGETS nd_compress DESTINATION ev_graph_compress)
