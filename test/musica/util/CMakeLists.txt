# MUSICA grid test
add_executable(test_musica_grid test_musica_grid.F90)

target_sources(test_musica_grid
  PUBLIC
    ${MUSICA_CCPP_SOURCES}
    ${TO_BE_CCPPIZED_SRC_PATH}/ccpp_const_utils.F90
    ${TO_BE_CCPPIZED_SRC_PATH}/ccpp_tuvx_utils.F90
    ${CCPP_SRC_PATH}/ccpp_constituent_prop_mod.F90
    ${CCPP_SRC_PATH}/ccpp_hash_table.F90
    ${CCPP_SRC_PATH}/ccpp_hashable.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
    ../musica_ccpp_namelist.F90
)

target_link_libraries(test_musica_grid
  PRIVATE
    musica-fortran musica yaml-cpp stdc++ netcdff
)

set_target_properties(test_musica_grid
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_musica_grid
  COMMAND $<TARGET_FILE:test_musica_grid>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_musica_grid $<TARGET_FILE:test_musica_grid> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
