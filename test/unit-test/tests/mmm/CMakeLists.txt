add_pfunit_ctest(mmm_physics_compat_tests
    TEST_SOURCES
        mmm_physics_compat_tests.pf
    OTHER_SOURCES
        ../../include/ccpp_kinds.F90
    LINK_LIBRARIES
        mmm_physics_compat
)
target_compile_options(mmm_physics_compat_tests
    PRIVATE
        $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:GNU>>:-fbacktrace -fcheck=all -ffpe-trap=invalid,overflow,zero -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -std=f2018 -Wall -Wextra -Wpedantic>
        $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:Intel>>:-check all -fp-model=precise -stand f18 -traceback -warn all>
        $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:IntelLLVM>>:-check all -fp-model=precise -stand f18 -traceback -warn all>
)
