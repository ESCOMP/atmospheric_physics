include(FetchContent)
include(TestUtils)

set(MUSICA_GIT_TAG 92b3365c1b567fbdfd135d8643b7c4cf81933bf4 CACHE STRING "Git tag for the musica repository")

FetchContent_Declare(musica
  GIT_REPOSITORY https://github.com/NCAR/musica.git
  GIT_TAG        ${MUSICA_GIT_TAG}
)

set(MUSICA_BUILD_C_CXX_INTERFACE OFF)
set(MUSICA_BUILD_FORTRAN_INTERFACE ON)
set(MUSICA_ENABLE_TESTS OFF)
set(MUSICA_ENABLE_INSTALL OFF)

FetchContent_MakeAvailable(musica)

# ---------------------------------------------------------
# Create a test for MUSICA CCPP wrapper
# ---------------------------------------------------------

add_executable(test_musica_api test_musica_api.F90 musica_ccpp_namelist.F90)

target_sources(test_musica_api
  PUBLIC
    ${MUSICA_SRC_PATH}/micm/musica_ccpp_micm.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx.F90
    ${MUSICA_SRC_PATH}/musica_ccpp.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
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