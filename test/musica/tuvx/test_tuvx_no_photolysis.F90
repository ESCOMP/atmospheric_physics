program test_tuvx_surface_albedo

   use musica_ccpp_tuvx_no_photolysis_rate

   implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

   call test_calculate_NO_photolysis_rate()

contains

   subroutine test_calculate_NO_photolysis_rate()
      use ccpp_kinds, only: kind_phys
      implicit none

      real(kind_phys) :: solar_zenith_angle
      real(kind_phys), dimension(8) :: extraterrestrial_flux 
      real(kind_phys), dimension(2, 4) :: constituents ! (layers, constituents)
      real(kind_phys), dimension(2) :: height_at_interfaces ! (layers)
      real(kind_phys), dimension(2) :: dry_air_density
      integer :: N2_index, O2_index, O3_index, NO_index
      real(kind_phys) :: molar_mass_N2, molar_mass_O2, molar_mass_O3, molar_mass_NO
      real(kind_phys), dimension(4) :: jNO

      ! Initialize test data
      solar_zenith_angle = 60.0_kind_phys
      extraterrestrial_flux = (/ 1.5e13_kind_phys, 1.4e13_kind_phys, 1.3e13_kind_phys, 1.2e13_kind_phys, &
                                1.1e13_kind_phys, 1.0e13_kind_phys, 9.0e12_kind_phys, 8.0e12_kind_phys /)
      constituents = reshape([0.21_kind_phys, 0.79_kind_phys, 1.0e-4_kind_phys, 1.0e-9_kind_phys, &
                              0.21_kind_phys, 0.79_kind_phys, 2.0e-4_kind_phys, 2.0e-9_kind_phys], shape(constituents))
      height_at_interfaces = (/ 4.0_kind_phys, 3.0_kind_phys /)
      dry_air_density = (/ 1.2_kind_phys, 1.1_kind_phys /)
      N2_index = 1
      O2_index = 2
      O3_index = 3
      NO_index = 4
      molar_mass_N2 = 0.0280134_kind_phys
      molar_mass_O2 = 0.0319988_kind_phys
      molar_mass_O3 = 0.0479982_kind_phys
      molar_mass_NO = 0.0300061_kind_phys

      ! Call the function to test
      jNO = calculate_NO_photolysis_rate(solar_zenith_angle, extraterrestrial_flux, constituents, height_at_interfaces, &
                                         dry_air_density, N2_index, O2_index, O3_index, NO_index, molar_mass_N2, molar_mass_O2, molar_mass_O3, molar_mass_NO)

      ! Validate the results
      ASSERT(size(jNO) == 4)
      print *, jNO
      ! ASSERT_NEAR(jNO(1), 1.0e-5_kind_phys, 1.0e-6_kind_phys)
      ! ASSERT_NEAR(jNO(2), 1.0e-5_kind_phys, 1.0e-6_kind_phys)
      ! ASSERT_NEAR(jNO(3), 1.0e-5_kind_phys, 1.0e-6_kind_phys)
      ! ASSERT_NEAR(jNO(4), 1.0e-5_kind_phys, 1.0e-6_kind_phys)
   end subroutine test_calculate_NO_photolysis_rate

end program test_tuvx_surface_albedo
