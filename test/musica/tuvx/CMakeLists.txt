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

# Wavelength grid
add_executable(test_tuvx_wavelength_grid test_tuvx_wavelength_grid.F90)

target_sources(test_tuvx_wavelength_grid
  PUBLIC
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_wavelength_grid.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_tuvx_wavelength_grid
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_tuvx_wavelength_grid
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_wavelength_grid
  COMMAND $<TARGET_FILE:test_tuvx_wavelength_grid>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_wavelength_grid $<TARGET_FILE:test_tuvx_wavelength_grid> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

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

# Extraterrestrial flux
add_executable(test_tuvx_extraterrestrial_flux test_tuvx_extraterrestrial_flux.F90)

target_sources(test_tuvx_extraterrestrial_flux
  PUBLIC
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_wavelength_grid.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_extraterrestrial_flux.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
    ${TO_BE_CCPPIZED_SRC_PATH}/ccpp_tuvx_utils.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_tuvx_extraterrestrial_flux
  PRIVATE
    musica-fortran
)

set_target_properties(test_tuvx_extraterrestrial_flux
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_extraterrestrial_flux
  COMMAND $<TARGET_FILE:test_tuvx_extraterrestrial_flux>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_extraterrestrial_flux $<TARGET_FILE:test_tuvx_extraterrestrial_flux> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Cloud optics
add_executable(test_tuvx_cloud_optics test_tuvx_cloud_optics.F90)

target_sources(test_tuvx_cloud_optics
  PUBLIC
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_height_grid.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_wavelength_grid.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_cloud_optics.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
    ${TO_BE_CCPPIZED_SRC_PATH}/ccpp_tuvx_utils.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_tuvx_cloud_optics
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_tuvx_cloud_optics
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_cloud_optics
  COMMAND $<TARGET_FILE:test_tuvx_cloud_optics>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_cloud_optics $<TARGET_FILE:test_tuvx_cloud_optics> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Aerosol optics
add_executable(test_tuvx_aerosol_optics test_tuvx_aerosol_optics.F90)

target_sources(test_tuvx_aerosol_optics
  PUBLIC
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_height_grid.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_wavelength_grid.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_aerosol_optics.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_tuvx_aerosol_optics
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_tuvx_aerosol_optics
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_aerosol_optics
  COMMAND $<TARGET_FILE:test_tuvx_aerosol_optics>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_aerosol_optics $<TARGET_FILE:test_tuvx_aerosol_optics> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# TUV-x gas species profiles
add_executable(test_tuvx_gas_species test_tuvx_gas_species.F90)

target_sources(test_tuvx_gas_species
  PUBLIC
    ${MUSICA_CCPP_SOURCES}
    ${TO_BE_CCPPIZED_SRC_PATH}/ccpp_const_utils.F90
    ${TO_BE_CCPPIZED_SRC_PATH}/ccpp_tuvx_utils.F90
    ${CCPP_SRC_PATH}/ccpp_constituent_prop_mod.F90
    ${CCPP_SRC_PATH}/ccpp_hash_table.F90
    ${CCPP_SRC_PATH}/ccpp_hashable.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/../musica_ccpp_namelist.F90
)

target_link_libraries(test_tuvx_gas_species
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_tuvx_gas_species
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_gas_species
  COMMAND $<TARGET_FILE:test_tuvx_gas_species>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_gas_species $<TARGET_FILE:test_tuvx_gas_species> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Configure TUV-x species
add_executable(test_tuvx_load_species test_tuvx_load_species.F90)

target_sources(test_tuvx_load_species
  PUBLIC
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_load_species.F90
    ${MUSICA_SRC_PATH}/tuvx/musica_ccpp_tuvx_gas_species.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_species.F90
    ${MUSICA_SRC_PATH}/musica_ccpp_util.F90
    ${TO_BE_CCPPIZED_SRC_PATH}/ccpp_const_utils.F90
    ${CCPP_SRC_PATH}/ccpp_constituent_prop_mod.F90
    ${CCPP_SRC_PATH}/ccpp_hash_table.F90
    ${CCPP_SRC_PATH}/ccpp_hashable.F90
    ${CCPP_TEST_SRC_PATH}/ccpp_kinds.F90
)

target_link_libraries(test_tuvx_load_species
  PRIVATE
    musica::musica-fortran
)

set_target_properties(test_tuvx_load_species
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

add_test(
  NAME test_tuvx_load_species
  COMMAND $<TARGET_FILE:test_tuvx_load_species>
  WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
)

add_memory_check_test(test_tuvx_load_species $<TARGET_FILE:test_tuvx_load_species> "" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})