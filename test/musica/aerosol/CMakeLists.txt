# MUSICA integration test for the abstract aerosol model
add_executable(test_musica_aerosol_model_mock_host mock_host.F90)

target_sources(test_musica_aerosol_model_mock_host
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

target_link_libraries(test_musica_aerosol_model_mock_host
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_musica_aerosol_model_mock_host
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_musica_aerosol_model_mock_host
  COMMAND $<TARGET_FILE:test_musica_aerosol_model_mock_host>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_musica_aerosol_model_mock_host $<TARGET_FILE:test_musica_aerosol_model_mock_host> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})