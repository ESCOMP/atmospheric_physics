add_executable(test_micm_util test_micm_util.F90 test_micm_mock_mods.F90)

target_sources(test_micm_util
  PUBLIC
    ${MUSICA_SRC_PATH}/micm/musica_ccpp_micm_util.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_micm_util
  PRIVATE
    musica-fortran musica mechanism_configuration yaml-cpp stdc++ netcdff
)

add_test(
  NAME test_micm_util
  COMMAND $<TARGET_FILE:test_micm_util>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_micm_util $<TARGET_FILE:test_micm_util> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})