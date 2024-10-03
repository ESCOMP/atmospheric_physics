cmake_minimum_required(VERSION 3.21)

project(
  atmosphericphysics
  VERSION 0.0.0
  LANGUAGES Fortran C CXX
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_MODULE_PATH}/SetDefaults.cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# --------------------------------------------------------------------------------
# NOTE: If 'CCPP_ENABLE_MUSICA_TESTS' on, this is not a stand-alone cmake project anymore.
#       MUSICA CCPP wrapper needs MUSICA library and ccpp-framework/src.
#       To 'CCPP_ENABLE_MUSICA_TESTS', you either build a cmake project through
#       'docker/Dockerfile.musica' or follow the build instructions in the file.
#       The following '$ENV' variables are set by the docker file.
# --------------------------------------------------------------------------------
option(CCPP_ENABLE_MUSICA_TESTS "Build the MUSICA tests" OFF)
option(CCPP_ENABLE_MEMCHECK "Enable memory checks in tests" OFF)

set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

if (CCPP_ENABLE_MUSICA_TESTS)
  set(MUSICA_VERSION $ENV{MUSICA_VERSION})
  set(MUSICA_SRC_PATH ${CMAKE_SOURCE_DIR}/../schemes/musica)
  set(CCPP_SRC_PATH ${CMAKE_SOURCE_DIR}/$ENV{CCPP_SRC_PATH})
  set(CCPP_TEST_SRC_PATH ${CMAKE_SOURCE_DIR}/include)

  add_subdirectory(musica)

  enable_testing()
endif()
