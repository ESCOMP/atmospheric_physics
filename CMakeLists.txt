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
# NOTE: If 'ENABLE_MICM_TEST' on, this is not a stand-alone cmake project anymore. 
#       MICM CCPP wrapper needs MUSICA library and ccpp-framework/src.
#       To 'ENABLE_MICM_TEST', you either build a cmake project through
#       'docker/Dockerfile.musica' or follow the build instructions in the file.
#       The following '$ENV' variables are set by the docker file.
# --------------------------------------------------------------------------------
option(ENABLE_MICM_TEST "Build the micm test" OFF)

if (ENABLE_MICM_TEST)
  set(MUSICA_VERSION $ENV{MUSICA_VERSION})
  set(MICM_SRC_PATH ${CMAKE_SOURCE_DIR}/musica/micm)
  set(CCPP_SRC_PATH ${CMAKE_SOURCE_DIR}/$ENV{CCPP_SRC_PATH})
  set(CCPP_TEST_SRC_PATH ${CMAKE_SOURCE_DIR}/test/include)

  add_subdirectory(test/musica/micm)

  enable_testing()
endif()