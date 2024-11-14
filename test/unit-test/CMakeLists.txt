cmake_minimum_required(VERSION 3.17)

project(atmospheric_physics VERSION 0.0.1 LANGUAGES Fortran)

find_package(PFUNIT REQUIRED)

if(NOT ATMOSPHERIC_PHYSICS_IS_TOP_LEVEL)
  message(WARNING "atmospheric-physics is not integrated into the CMake build of any top level "
                  "project yet and this CMake is for testing purposes only.  "
                  "Making a change to this project's CMake will not impact the build of "
                  "a parent project at this time.")
endif()

option(ATMOSPHERIC_PHYSICS_ENABLE_TESTS "Run pFUnit unit tests" OFF)
option(ATMOSPHERIC_PHYSICS_ENABLE_CODE_COVERAGE "Run code coverage tool" OFF)

if(ATMOSPHERIC_PHYSICS_ENABLE_CODE_COVERAGE)
  add_compile_options(-O0 --coverage)
  add_link_options(--coverage)
endif()

set(CMAKE_BUILD_TYPE Debug)

set(UTILITIES_SRC
    ../../schemes/utilities/state_converters.F90
    ../../schemes/utilities/static_energy.F90
    ../../schemes/utilities/physics_tendency_updaters.F90
    include/ccpp_kinds.F90
)

add_library(utilities ${UTILITIES_SRC})
target_compile_options(utilities PRIVATE -ffree-line-length-none)
target_include_directories(utilities PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if(ATMOSPHERIC_PHYSICS_ENABLE_TESTS OR ATMOSPHERIC_PHYSICS_ENABLE_CODE_COVERAGE)
  enable_testing()
  add_subdirectory(tests)
endif()

