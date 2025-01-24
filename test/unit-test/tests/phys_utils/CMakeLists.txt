add_pfunit_ctest(phys_utils_tests
  TEST_SOURCES test_atmos_pbl_utils.pf
  LINK_LIBRARIES phys_utils
)
# target_compile_options(phys_utils_tests PRIVATE -std=f2008 -fmax-identifier-length=64)