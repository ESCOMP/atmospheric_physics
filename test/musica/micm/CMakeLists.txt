
add_executable(test_micm_util test_micm_util.F90)

target_sources(test_micm_util
  PUBLIC
    ${MUSICA_SRC_PATH}/micm/micm_util.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_micm_util
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_micm_util
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_micm_util
  COMMAND $<TARGET_FILE:test_micm_util>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_micm_util $<TARGET_FILE:test_micm_util> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})