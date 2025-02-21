cmake_minimum_required(VERSION 3.21)

project(
  atmosphericphysics
  VERSION 0.1.0
  LANGUAGES Fortran C CXX
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_MODULE_PATH}/SetDefaults.cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ---------------------------------------------------------------------------------------------
# NOTE: If 'CCPP_ENABLE_MUSICA_TESTS' is enabled, this is no longer a stand-alone CMake project.
#       The MUSICA CCPP wrapper requires both the MUSICA library and ccpp-framework/src.
#       To enable 'CCPP_ENABLE_MUSICA_TESTS', you can either build a CMake project using
#       'docker/Dockerfile.musica' or follow the build instructions outlined in that file.
#       The following '$ENV' variables are configured by the Dockerfile.
# ---------------------------------------------------------------------------------------------
option(CCPP_ENABLE_MUSICA_TESTS "Build the MUSICA tests" OFF)
option(CCPP_ENABLE_MEMCHECK "Enable memory checks in tests" OFF)

set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

if (CCPP_ENABLE_MUSICA_TESTS)
  set(MUSICA_SRC_PATH ${CMAKE_SOURCE_DIR}/../schemes/musica)
  set(TO_BE_CCPPIZED_SRC_PATH ${CMAKE_SOURCE_DIR}/../to_be_ccppized)
  set(CCPP_SRC_PATH ${CMAKE_SOURCE_DIR}/$ENV{CCPP_SRC_PATH})
  set(CCPP_TEST_SRC_PATH ${CMAKE_SOURCE_DIR}/include)

  include(TestUtils)
  include(CTest)
  enable_testing()

  add_subdirectory(musica)
endif()
