cmake_minimum_required(VERSION 3.24)

project(
  CFD_LBM
  VERSION 0.1.0
  DESCRIPTION "General-purpose high-performance Lattice Boltzmann framework"
  LANGUAGES CXX
)

option(LBM_BUILD_EXAMPLES "Build LBM example executables" ON)
option(LBM_BUILD_TESTS "Build LBM tests" ON)

add_library(lbm_core INTERFACE)
add_library(lbm::core ALIAS lbm_core)

target_compile_features(lbm_core INTERFACE cxx_std_20)
target_include_directories(
  lbm_core
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

if(MSVC)
  target_compile_options(lbm_core INTERFACE /W4 /permissive-)
else()
  target_compile_options(lbm_core INTERFACE -Wall -Wextra -Wpedantic)
endif()

if(LBM_BUILD_EXAMPLES)
  add_executable(lbm_periodic_smoke examples/periodic_smoke.cpp)
  target_link_libraries(lbm_periodic_smoke PRIVATE lbm::core)

  add_executable(lbm_poiseuille examples/poiseuille_flow.cpp)
  target_link_libraries(lbm_poiseuille PRIVATE lbm::core)

  add_executable(lbm_couette examples/couette_flow.cpp)
  target_link_libraries(lbm_couette PRIVATE lbm::core)
endif()

include(CTest)

if(LBM_BUILD_TESTS AND BUILD_TESTING)
  add_executable(lbm_core_tests tests/lbm_core_tests.cpp)
  target_link_libraries(lbm_core_tests PRIVATE lbm::core)
  add_test(NAME lbm_core_tests COMMAND lbm_core_tests)
endif()
