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
      real(kind_phys), dimension(4) :: extraterrestrial_flux 
      real(kind_phys), dimension(1,2,4) :: constituents ! (column, layers, constituents)  kg kg-1
      real(kind_phys), dimension(3) :: height_at_interfaces ! (layers + 1) km
      real(kind_phys), dimension(1,2) :: dry_air_density ! (column, layers) kg m-3
      integer :: N2_index, O2_index, O3_index, NO_index
      real(kind_phys) :: molar_mass_N2, molar_mass_O2, molar_mass_O3, molar_mass_NO ! kg mol-1
      real(kind_phys), dimension(4) :: jNO

      ! Some of this initialization data corresponds to values listed in the Minschwaner and Siskind (1993) paper
      ! see the musica_ccpp_tuvx_no_photolysis_rate.F90 file for the full citation


      ! Initialize test data
      solar_zenith_angle = 60.0_kind_phys
      extraterrestrial_flux = (/ 1.5e13_kind_phys, 1.4e13_kind_phys, 1.3e13_kind_phys, 1.2e13_kind_phys /)
      height_at_interfaces = reshape([40.0_kind_phys, 30.0_kind_phys, 20.0_kind_phys], shape(height_at_interfaces))
      dry_air_density = reshape([1.2_kind_phys, 1.1_kind_phys], shape(dry_air_density))
      N2_index = 1
      O2_index = 2
      O3_index = 3
      NO_index = 4
      molar_mass_N2 = 0.0280134_kind_phys
      molar_mass_O2 = 0.0319988_kind_phys
      molar_mass_O3 = 0.0479982_kind_phys
      molar_mass_NO = 0.3000610_kind_phys

      constituents(1,:,N2_index) = 0.78_kind_phys
      constituents(1,:,O2_index) = 0.21_kind_phys
      constituents(1,:,O3_index) = 20.0e-9_kind_phys
      constituents(1,:,NO_index) = 0.3e-9_kind_phys

      ! print out the consitutents by index to ensure they are set correctly
      print *, "N2: ", constituents(1,1,N2_index)
      print *, "O2: ", constituents(1,1,O2_index)
      print *, "O3: ", constituents(1,1,O3_index)
      print *, "NO: ", constituents(1,1,NO_index)

      ! Call the function to test
      jNO = calculate_NO_photolysis_rate(size(constituents, dim=2), solar_zenith_angle, extraterrestrial_flux, constituents(1,:,:), height_at_interfaces, &
                                         dry_air_density(1,:), N2_index, O2_index, O3_index, NO_index, molar_mass_N2, molar_mass_O2, molar_mass_O3, molar_mass_NO)

      ! Validate the results
      ASSERT(size(jNO) == 4)
      print *, jNO
      ! ASSERT_NEAR(jNO(1), 1.0e-5_kind_phys, 1.0e-6_kind_phys)
      ! ASSERT_NEAR(jNO(2), 1.0e-5_kind_phys, 1.0e-6_kind_phys)
      ! ASSERT_NEAR(jNO(3), 1.0e-5_kind_phys, 1.0e-6_kind_phys)
      ! ASSERT_NEAR(jNO(4), 1.0e-5_kind_phys, 1.0e-6_kind_phys)
   end subroutine test_calculate_NO_photolysis_rate

end program test_tuvx_surface_albedo
