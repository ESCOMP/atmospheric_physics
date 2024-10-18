# Height grid
add_executable(test_tuvx_height_grid test_tuvx_height_grid.F90)

target_sources(test_tuvx_height_grid
  PUBLIC
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_height_grid.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_tuvx_height_grid
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_tuvx_height_grid
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_height_grid
  COMMAND $<TARGET_FILE:test_tuvx_height_grid>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_height_grid $<TARGET_FILE:test_tuvx_height_grid> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Temperature
add_executable(test_tuvx_temperature test_tuvx_temperature.F90)

target_sources(test_tuvx_temperature
  PUBLIC
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_height_grid.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_temperature.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_tuvx_temperature
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_tuvx_temperature
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_temperature
  COMMAND $<TARGET_FILE:test_tuvx_temperature>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_temperature $<TARGET_FILE:test_tuvx_temperature> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Surface albedo
add_executable(test_tuvx_surface_albedo test_tuvx_surface_albedo.F90)

target_sources(test_tuvx_surface_albedo
  PUBLIC
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_wavelength_grid.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_surface_albedo.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_tuvx_surface_albedo
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_tuvx_surface_albedo
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_surface_albedo
  COMMAND $<TARGET_FILE:test_tuvx_surface_albedo>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_surface_albedo $<TARGET_FILE:test_tuvx_surface_albedo> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})