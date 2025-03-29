! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_test_data

  implicit none

  private
  public :: get_wavelength_edges

contains

  !> Returns the wavelength edges used in the CAM-Chem photolysis rate constant lookup table
  !! These are the values that will be used in CAM-SIMA and correspond to the wavelength
  !! bins used in the CAM-Chem photolysis rate constant lookup table.
  !!
  !! We're using the actual values here because several of the TS1/TSMLT photolysis
  !! rate constant configurations are sensitive to the wavelength grid.
  subroutine get_wavelength_edges(edges)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), dimension(:), intent(out) :: edges ! m

    edges = (/ &
      120.0e-9_kind_phys, &
      121.4e-9_kind_phys, &
      121.9e-9_kind_phys, &
      123.5e-9_kind_phys, &
      124.3e-9_kind_phys, &
      125.5e-9_kind_phys, &
      126.3e-9_kind_phys, &
      127.1e-9_kind_phys, &
      130.1e-9_kind_phys, &
      131.1e-9_kind_phys, &
      135.0e-9_kind_phys, &
      140.0e-9_kind_phys, &
      145.0e-9_kind_phys, &
      150.0e-9_kind_phys, &
      155.0e-9_kind_phys, &
      160.0e-9_kind_phys, &
      165.0e-9_kind_phys, &
      168.0e-9_kind_phys, &
      171.0e-9_kind_phys, &
      173.0e-9_kind_phys, &
      174.4e-9_kind_phys, &
      175.4e-9_kind_phys, &
      177.0e-9_kind_phys, &
      178.6e-9_kind_phys, &
      180.2e-9_kind_phys, &
      181.8e-9_kind_phys, &
      183.5e-9_kind_phys, &
      185.2e-9_kind_phys, &
      186.9e-9_kind_phys, &
      188.7e-9_kind_phys, &
      190.5e-9_kind_phys, &
      192.3e-9_kind_phys, &
      194.2e-9_kind_phys, &
      196.1e-9_kind_phys, &
      198.0e-9_kind_phys, &
      200.0e-9_kind_phys, &
      202.0e-9_kind_phys, &
      204.1e-9_kind_phys, &
      206.2e-9_kind_phys, &
      208.0e-9_kind_phys, &
      211.0e-9_kind_phys, &
      214.0e-9_kind_phys, &
      217.0e-9_kind_phys, &
      220.0e-9_kind_phys, &
      223.0e-9_kind_phys, &
      226.0e-9_kind_phys, &
      229.0e-9_kind_phys, &
      232.0e-9_kind_phys, &
      235.0e-9_kind_phys, &
      238.0e-9_kind_phys, &
      241.0e-9_kind_phys, &
      244.0e-9_kind_phys, &
      247.0e-9_kind_phys, &
      250.0e-9_kind_phys, &
      253.0e-9_kind_phys, &
      256.0e-9_kind_phys, &
      259.0e-9_kind_phys, &
      263.0e-9_kind_phys, &
      267.0e-9_kind_phys, &
      271.0e-9_kind_phys, &
      275.0e-9_kind_phys, &
      279.0e-9_kind_phys, &
      283.0e-9_kind_phys, &
      287.0e-9_kind_phys, &
      291.0e-9_kind_phys, &
      295.0e-9_kind_phys, &
      298.5e-9_kind_phys, &
      302.5e-9_kind_phys, &
      305.5e-9_kind_phys, &
      308.5e-9_kind_phys, &
      311.5e-9_kind_phys, &
      314.5e-9_kind_phys, &
      317.5e-9_kind_phys, &
      322.5e-9_kind_phys, &
      327.5e-9_kind_phys, &
      332.5e-9_kind_phys, &
      337.5e-9_kind_phys, &
      342.5e-9_kind_phys, &
      347.5e-9_kind_phys, &
      350.0e-9_kind_phys, &
      355.0e-9_kind_phys, &
      360.0e-9_kind_phys, &
      365.0e-9_kind_phys, &
      370.0e-9_kind_phys, &
      375.0e-9_kind_phys, &
      380.0e-9_kind_phys, &
      385.0e-9_kind_phys, &
      390.0e-9_kind_phys, &
      395.0e-9_kind_phys, &
      400.0e-9_kind_phys, &
      405.0e-9_kind_phys, &
      410.0e-9_kind_phys, &
      415.0e-9_kind_phys, &
      420.0e-9_kind_phys, &
      430.0e-9_kind_phys, &
      440.0e-9_kind_phys, &
      450.0e-9_kind_phys, &
      500.0e-9_kind_phys, &
      550.0e-9_kind_phys, &
      600.0e-9_kind_phys, &
      650.0e-9_kind_phys, &
      700.0e-9_kind_phys, &
      750.0e-9_kind_phys &
    /)
 
  end subroutine get_wavelength_edges

end module musica_test_data