add_pfunit_ctest(mmm_physics_compat_tests
    TEST_SOURCES
        mmm_physics_compat_tests.pf
    OTHER_SOURCES
        ../../include/ccpp_kinds.F90
    LINK_LIBRARIES
        mmm_physics_compat
)
