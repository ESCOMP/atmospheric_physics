# MUSICA stub aerosol model tests
add_executable(test_stub_aerosol_model test_stub_aerosol_model.F90)

target_sources(test_stub_aerosol_model
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

target_link_libraries(test_stub_aerosol_model
  PRIVATE
    musica-fortran musica yaml-cpp stdc++ netcdff
)

set_target_properties(test_stub_aerosol_model
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_stub_aerosol_model
  COMMAND $<TARGET_FILE:test_stub_aerosol_model>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_stub_aerosol_model $<TARGET_FILE:test_stub_aerosol_model> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})