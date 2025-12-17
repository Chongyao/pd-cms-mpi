project(pd-cms-mpi)
cmake_minimum_required(VERSION 3.5.0)

set(CPM_USE_NAMED_CACHE_DIRECTORIES ON)

set(CPM_SOURCE_CACHE ${CMAKE_SOURCE_DIR}/../3rd/CPMCache)
include(cmake/CPM.cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")


#O3 optimization
message("cmake build type: ${CMAKE_BUILD_TYPE}")
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")
endif()
#enable all instructions

set(DISABLED_WARNING "-Wno-unused-but-set-variable  -Wno-unused-parameter -Wno-c++11-narrowing -Wno-unused-variable  -Wno-array-bounds -Wno-tautological-compare -Wno-deprecated-literal-operator")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native ${DISABLED_WARNING}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wno-everything")



# PIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)



# List of all includes
set(INCLUDES
  blas
  openmp
  suitesparse
  eigen
  spectra
  fmt
  spdlog
  mpi
)

foreach(INCLUDE IN LISTS INCLUDES)
  include(cmake/include_${INCLUDE}.cmake)
  message("Load ${INCLUDE} done.")
endforeach()













include_directories(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(src)
add_subdirectory(apps)

