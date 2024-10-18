include(FetchContent)
include(TestUtils)

FetchContent_Declare(musica
  GIT_REPOSITORY https://github.com/NCAR/musica.git
  GIT_TAG        $ENV{MUSICA_GIT_TAG}  # Set by docker
)
message(STATUS "Using MUSICA commit: $ENV{MUSICA_GIT_TAG}")

set(MUSICA_BUILD_C_CXX_INTERFACE ON)
set(MUSICA_BUILD_FORTRAN_INTERFACE ON)
set(MUSICA_ENABLE_TESTS OFF)
set(MUSICA_ENABLE_INSTALL OFF)

FetchContent_MakeAvailable(musica)

# ---------------------------------------------------------
# Create a test for MUSICA CCPP wrapper
# ---------------------------------------------------------

add_executable(test_musica_api test_musica_api.F90 musica_ccpp_namelist.F90)

file(GLOB MUSICA_CCPP_SOURCES
      ${MUSICA_SRC_PATH}/*.F90
      ${MUSICA_SRC_PATH}/micm/*.F90
      ${MUSICA_SRC_PATH}/tuvx/*.F90
)

target_sources(test_musica_api
  PUBLIC
    ${MUSICA_CCPP_SOURCES}
    ${CCPP_SRC_PATH}/ccpp_constituent_prop_mod.F90
    ${CCPP_SRC_PATH}/ccpp_hash_table.F90
    ${CCPP_SRC_PATH}/ccpp_hashable.F90
    ${CCPP_SRC_PATH}/ccpp_types.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_musica_api
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_musica_api
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

include(CTest)

add_test(
  NAME test_musica_api
  COMMAND $<TARGET_FILE:test_musica_api>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_musica_api $<TARGET_FILE:test_musica_api> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_custom_target(
  copy_micm_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/micm/configs/chapman ${CMAKE_BINARY_DIR}/chapman
)

add_custom_target(
  copy_tuvx_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/tuvx/configs ${CMAKE_BINARY_DIR}/configs
)

add_subdirectory(micm)
add_subdirectory(tuvx)

# Test metdadata
find_package(Python REQUIRED)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/metadata_test)
add_custom_target(
  copy_metadata_test_files ALL ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/../../schemes/musica/musica_ccpp.meta
  ${CMAKE_CURRENT_SOURCE_DIR}/../../schemes/musica/musica_ccpp.F90
  ${CMAKE_BINARY_DIR}/metadata_test
)

add_test(
  NAME test_metadata
  COMMAND ${Python_EXECUTABLE} ${CMAKE_BINARY_DIR}/../$ENV{CCPP_FORTRAN_TOOLS_PATH}/offline_check_fortran_vs_metadata.py
  --directory ${CMAKE_BINARY_DIR}/metadata_test
)