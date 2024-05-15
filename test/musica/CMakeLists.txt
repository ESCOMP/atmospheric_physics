include(FetchContent)
include(TestUtils)

FetchContent_Declare(musica
  GIT_REPOSITORY https://github.com/NCAR/musica.git
#  GIT_TAG        v0.6.1
  GIT_TAG        2a5eeaac982a3eb80b96d1e2087b91b301d1e748
)

set(MUSICA_BUILD_C_CXX_INTERFACE OFF)
set(MUSICA_BUILD_FORTRAN_INTERFACE ON)
set(MUSICA_ENABLE_TESTS OFF)
set(MUSICA_ENABLE_INSTALL OFF)

FetchContent_MakeAvailable(musica)

# ---------------------------------------------------------
# Create a test for MUSICA CCPP wrapper
# ---------------------------------------------------------

add_executable(test_musica_api test_musica_api.F90 musica_namelist.F90)

target_sources(test_musica_api
  PUBLIC
    ${MUSICA_SRC_PATH}/micm/micm.F90
    ${MUSICA_SRC_PATH}/musica.F90
    ${MUSICA_SRC_PATH}/util.F90
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