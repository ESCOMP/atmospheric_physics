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

add_library(mmm_physics_compat)
target_sources(mmm_physics_compat
    PRIVATE
        ../../test/unit-test/include/ccpp_kinds.F90
        ccpp_kind_types.F90
        mmm_physics_compat.F90
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
