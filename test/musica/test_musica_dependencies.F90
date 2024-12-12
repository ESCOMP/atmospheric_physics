! Copyright (C) 2024 National Science Foundation-National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_musica_ccpp_dependencies

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  implicit none

  call test_dependencies()

contains

  subroutine test_dependencies()

    use ccpp_kinds, only: kind_phys
    use musica_ccpp_dependencies, only: musica_ccpp_dependencies_init

    integer :: photolysis_wavelength_grid_section_dimension
    integer :: photolysis_wavelength_grid_interface_dimension
    real(kind_phys) :: surface_albedo
    real(kind_phys), allocatable :: photolysis_wavelength_grid_interfaces(:)
    real(kind_phys), allocatable :: extraterrestrial_flux(:)

    call musica_ccpp_dependencies_init(photolysis_wavelength_grid_section_dimension, &
                                       photolysis_wavelength_grid_interface_dimension, &
                                       surface_albedo, photolysis_wavelength_grid_interfaces, &
                                       extraterrestrial_flux)    

    ASSERT(size(photolysis_wavelength_grid_interfaces) == photolysis_wavelength_grid_interface_dimension)
    ASSERT(size(extraterrestrial_flux) == photolysis_wavelength_grid_section_dimension)
    ASSERT(photolysis_wavelength_grid_interface_dimension == photolysis_wavelength_grid_section_dimension + 1)
    ASSERT(all(photolysis_wavelength_grid_interfaces > 1.0e-10))
    ASSERT(all(photolysis_wavelength_grid_interfaces < 1.0e-4))
    ASSERT(surface_albedo >= 0.0)

  end subroutine test_dependencies

end program test_musica_ccpp_dependencies