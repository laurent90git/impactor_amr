cmake_minimum_required(VERSION 3.16)
project(eulerian_hypervelocity_impact CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(samurai CONFIG REQUIRED)

# one executable per mesh dimension
foreach(_d 1 2 3)
  add_executable(hvi${_d}d src/app/main.cpp)
  target_link_libraries(hvi${_d}d samurai::samurai)
  target_compile_definitions(hvi${_d}d PRIVATE HVI_DIM=${_d})
endforeach()

# tests of the physics kernel (no samurai needed)
add_executable(test_physics tests/test_physics.cpp)
target_include_directories(test_physics PRIVATE src)
enable_testing()
add_test(NAME physics COMMAND test_physics)

# OpenMP build: cmake -S . -B build-omp -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DSAMURAI_WITH_OPENMP=ON
