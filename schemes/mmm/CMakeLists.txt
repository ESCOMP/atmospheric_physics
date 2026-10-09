cmake_minimum_required(VERSION 3.21)

# `mmm_physics_compat` has not been integrated into the CMake build of any top level projects yet,
# and this CMakeLists.txt file is currently for unit testing purposes only.
# Therefore, making a change to this CMakeLists.txt file will not impact the build of a parent project at this time.
project(mmm_physics_compat
    VERSION
        0.1.0
    DESCRIPTION
        "MMM physics compatibility layer for CCPP"
    LANGUAGES
        Fortran
)

# Build the actual MMM physics code.
add_library(mmm_physics)
target_sources(mmm_physics
    PRIVATE
        ccpp_kind_types.F90
        mmm_physics/bl_gwdo.F90
        mmm_physics/bl_mynn.F90
        mmm_physics/bl_mynn_subroutines.F90
        mmm_physics/cu_ntiedtke.F90
        mmm_physics/mynn_shared.F90
        mmm_physics/sf_mynn.F90
)
target_include_directories(mmm_physics
    INTERFACE
        ${CMAKE_CURRENT_BINARY_DIR}
)
target_link_libraries(mmm_physics
    PUBLIC
        ccpp_framework_stub
)

add_library(mmm_physics_compat)
target_sources(mmm_physics_compat
    PRIVATE
        bl_gwdo_compat.F90
        bl_mynn_compat.F90
        cldfra_xu_randall.F90
        cu_ntiedtke_compat.F90
        mmm_physics_compat.F90
        sf_mynn_compat.F90
)
target_include_directories(mmm_physics_compat
    INTERFACE
        ${CMAKE_CURRENT_BINARY_DIR}
)
target_compile_options(mmm_physics_compat
    PRIVATE
        $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:GNU>>:-fbacktrace -fcheck=all -ffpe-trap=invalid,overflow,zero -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -std=f2018 -Wall -Wextra -Wpedantic>
        $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:Intel>>:-check all -fp-model=precise -stand f18 -traceback -warn all>
        $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:IntelLLVM>>:-check all -fp-model=precise -stand f18 -traceback -warn all>
)
target_link_libraries(mmm_physics_compat
    PRIVATE
        mmm_physics
    PUBLIC
        ccpp_framework_stub
)
