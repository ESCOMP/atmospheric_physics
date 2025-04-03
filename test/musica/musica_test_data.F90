! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_test_data

  implicit none

  private
  public :: get_wavelength_edges, get_extrarterrestrial_fluxes

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

  !> file NECKEL.FLX  extraterrestrial solar irradiance 330-1250nm
  !! from: H.Neckel and D.Labs, "The Solar Radiation Between 3300 and 12500 A",
  !!       Solar Physics v.90, pp.205-258 (1984).
  !! Resolution:  1 nm between 330.5 and 529.5 nm
  !!              2 nm between 631.0 and 709.0 nm
  !!              5 nm between 872.5 and 1247.4 nm
  !! Units: converted from 1.0E-6 Watts cm-2 A-1 to photons cm-2 s-1 nm-1
  !! Values are for the whole disk irradiance, at mean Earth-Sun distance
  !! wavelengths are for the center of the intervals.
  !!
  !! The flux values are extracted from the specified file for the first 102 entries.
  subroutine get_extrarterrestrial_fluxes(fluxes)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), dimension(:), intent(out) :: fluxes ! photons cm-2 s-1 nm-1

    fluxes = (/ &
      1.67555e+14_kind_phys, &
      1.61720e+14_kind_phys, &
      1.54341e+14_kind_phys, &
      1.52118e+14_kind_phys, &
      1.58469e+14_kind_phys, &
      1.66036e+14_kind_phys, &
      1.29769e+14_kind_phys, &
      1.47316e+14_kind_phys, &
      1.56274e+14_kind_phys, &
      1.60325e+14_kind_phys, &
      1.70225e+14_kind_phys, &
      1.61097e+14_kind_phys, &
      1.71742e+14_kind_phys, &
      1.70514e+14_kind_phys, &
      1.24876e+14_kind_phys, &
      1.68376e+14_kind_phys, &
      1.60490e+14_kind_phys, &
      1.57979e+14_kind_phys, &
      1.66505e+14_kind_phys, &
      1.52378e+14_kind_phys, &
      1.97811e+14_kind_phys, &
      1.75901e+14_kind_phys, &
      1.54751e+14_kind_phys, &
      1.98793e+14_kind_phys, &
      2.02567e+14_kind_phys, &
      1.89537e+14_kind_phys, &
      1.68532e+14_kind_phys, &
      1.60546e+14_kind_phys, &
      1.13346e+14_kind_phys, &
      2.05967e+14_kind_phys, &
      1.77864e+14_kind_phys, &
      1.62888e+14_kind_phys, &
      2.14804e+14_kind_phys, &
      1.75501e+14_kind_phys, &
      1.86444e+14_kind_phys, &
      2.32774e+14_kind_phys, &
      2.30828e+14_kind_phys, &
      2.24982e+14_kind_phys, &
      2.02219e+14_kind_phys, &
      2.47972e+14_kind_phys, &
      2.00891e+14_kind_phys, &
      2.44825e+14_kind_phys, &
      1.99913e+14_kind_phys, &
      1.57765e+14_kind_phys, &
      1.65729e+14_kind_phys, &
      2.16080e+14_kind_phys, &
      2.09073e+14_kind_phys, &
      2.45738e+14_kind_phys, &
      2.55917e+14_kind_phys, &
      1.91251e+14_kind_phys, &
      2.47308e+14_kind_phys, &
      2.10889e+14_kind_phys, &
      1.41346e+14_kind_phys, &
      1.32255e+14_kind_phys, &
      1.98997e+14_kind_phys, &
      1.85347e+14_kind_phys, &
      2.08594e+14_kind_phys, &
      1.88650e+14_kind_phys, &
      1.78574e+14_kind_phys, &
      2.41000e+14_kind_phys, &
      2.40832e+14_kind_phys, &
      2.75942e+14_kind_phys, &
      1.88910e+14_kind_phys, &
      9.70730e+13_kind_phys, &
      2.19069e+14_kind_phys, &
      2.74779e+14_kind_phys, &
      1.29952e+14_kind_phys, &
      2.08327e+14_kind_phys, &
      3.08963e+14_kind_phys, &
      3.33270e+14_kind_phys, &
      3.32895e+14_kind_phys, &
      3.63642e+14_kind_phys, &
      3.65966e+14_kind_phys, &
      3.37217e+14_kind_phys, &
      3.26648e+14_kind_phys, &
      3.41746e+14_kind_phys, &
      3.32766e+14_kind_phys, &
      3.17377e+14_kind_phys, &
      3.75741e+14_kind_phys, &
      3.52127e+14_kind_phys, &
      3.10827e+14_kind_phys, &
      3.77464e+14_kind_phys, &
      3.72567e+14_kind_phys, &
      3.66392e+14_kind_phys, &
      3.63313e+14_kind_phys, &
      3.63562e+14_kind_phys, &
      3.87293e+14_kind_phys, &
      3.50809e+14_kind_phys, &
      3.55652e+14_kind_phys, &
      3.60092e+14_kind_phys, &
      3.73018e+14_kind_phys, &
      3.82393e+14_kind_phys, &
      3.37355e+14_kind_phys, &
      3.65658e+14_kind_phys, &
      3.78703e+14_kind_phys, &
      3.63957e+14_kind_phys, &
      3.65457e+14_kind_phys, &
      3.38550e+14_kind_phys, &
      3.43225e+14_kind_phys, &
      3.19808e+14_kind_phys, &
      2.46645e+14_kind_phys, &
      3.67134e+14_kind_phys &
    /)
  end subroutine get_extrarterrestrial_fluxes

end module musica_test_data