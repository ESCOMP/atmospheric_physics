add_pfunit_ctest(rrtmgp_tests
  TEST_SOURCES test_radiation_utils.pf
  TEST_SOURCES test_radiation_tools.pf
  TEST_SOURCES test_rrtmgp_sw_solar_var_setup.pf
  TEST_SOURCES test_rrtmgp_dry_static_energy_tendency.pf
  TEST_SOURCES test_calculate_net_heating.pf
  LINK_LIBRARIES rrtmgp
)
