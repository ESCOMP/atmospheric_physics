! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
program test_micm_surface

  use ccpp_kinds,               only: kind_phys
  use musica_ccpp_micm_surface, only: surface_reaction_parameters_from_mode, &
                                      find_surface_reactions, mode_from_label
  use musica_ccpp_util,         only: PI

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) if( (abs(a - b) >= abs_error) .and. (abs(a - b) /= 0.0) ) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a, b"; stop 1; endif

  call test_known_answer()
  call test_identity_with_cam_rate()
  call test_tropopause_mask()
  call test_no_surface_or_diameter()
  call test_find_surface_reactions()
  call test_find_surface_reactions_missing_partner()
  call test_mode_from_label()

contains

  ! Hand-computed values for one cell:
  ! sad = 1e-6 cm2 cm-3, d = 2e-5 cm (0.2 um)
  ! r = 0.5 * 2e-5 cm * 1e-2 m cm-1 = 1e-7 m
  ! N = 1e-6 cm2 cm-3 * 1e2 / (4 pi (1e-7 m)^2) = 1e-4 / (4 pi 1e-14) = 7.957747154594767e8 m-3
  subroutine test_known_answer()
    real(kind_phys) :: sad(1,1), diameter(1,1), radius(1,1), number(1,1)
    integer         :: tropopause_layer(1)

    sad(1,1) = 1.0e-6_kind_phys
    diameter(1,1) = 2.0e-5_kind_phys
    tropopause_layer(1) = 0

    call surface_reaction_parameters_from_mode(sad, diameter, tropopause_layer, radius, number)

    ASSERT_NEAR(radius(1,1), 1.0e-7_kind_phys, 1.0e-20_kind_phys)
    ASSERT_NEAR(number(1,1) / 7.957747154594767e8_kind_phys, 1.0_kind_phys, 1.0e-12_kind_phys)
  end subroutine test_known_answer

  ! For a spread of realistic per-mode values, MICM's 4 pi r^2 N and r must
  ! recover CAM's surface area density and d/2 to 1e-10 relative (the
  ! mam_surfarea_diagnostics identity check), and MICM's rate expression must
  ! then equal CAM's hetrxtrate term for the same gas parameters.
  subroutine test_identity_with_cam_rate()
    integer, parameter :: NUM_COLUMNS = 3
    integer, parameter :: NUM_LAYERS = 2
    real(kind_phys), parameter :: RTOL = 1.0e-10_kind_phys
    real(kind_phys), parameter :: DIFFUSION_CM2_S = 0.1_kind_phys     ! CAM dg
    real(kind_phys), parameter :: SPEED_CM_S = 2.4e4_kind_phys        ! mean molecular speed
    real(kind_phys), parameter :: GAMMA = 0.02_kind_phys              ! reaction probability
    real(kind_phys) :: sad(NUM_COLUMNS,NUM_LAYERS)       ! cm2 cm-3
    real(kind_phys) :: diameter(NUM_COLUMNS,NUM_LAYERS)  ! cm
    real(kind_phys) :: radius(NUM_COLUMNS,NUM_LAYERS)    ! m
    real(kind_phys) :: number(NUM_COLUMNS,NUM_LAYERS)    ! m-3
    real(kind_phys) :: recovered_sad, recovered_radius, expected_radius, cam_rate, micm_rate
    integer         :: tropopause_layer(NUM_COLUMNS), i_column, i_layer

    ! Aitken-, accumulation- and coarse-like magnitudes
    sad(:,1) = (/ 3.2e-8_kind_phys, 4.7e-7_kind_phys, 1.1e-6_kind_phys /)
    sad(:,2) = (/ 8.9e-9_kind_phys, 2.5e-6_kind_phys, 6.0e-7_kind_phys /)
    diameter(:,1) = (/ 3.1e-6_kind_phys, 1.8e-5_kind_phys, 2.2e-4_kind_phys /)
    diameter(:,2) = (/ 2.4e-6_kind_phys, 3.5e-5_kind_phys, 4.9e-4_kind_phys /)
    tropopause_layer(:) = 0

    call surface_reaction_parameters_from_mode(sad, diameter, tropopause_layer, radius, number)

    do i_layer = 1, NUM_LAYERS
      do i_column = 1, NUM_COLUMNS
        ! back to CAM units: m2 m-3 -> cm2 cm-3 is x1e-2, m -> cm is x1e2
        recovered_sad = 4.0_kind_phys * PI * radius(i_column, i_layer)**2 &
                        * number(i_column, i_layer) * 1.0e-2_kind_phys
        recovered_radius = radius(i_column, i_layer) * 1.0e2_kind_phys
        expected_radius = 0.5_kind_phys * diameter(i_column, i_layer)
        ASSERT(abs(recovered_sad - sad(i_column, i_layer)) <= RTOL * sad(i_column, i_layer))
        ASSERT(abs(recovered_radius - expected_radius) <= RTOL * expected_radius)

        ! both rates in CAM's cgs units, s-1
        cam_rate = sad(i_column, i_layer) &
                   / (0.5_kind_phys * diameter(i_column, i_layer) / DIFFUSION_CM2_S &
                      + 4.0_kind_phys / (SPEED_CM_S * GAMMA))
        micm_rate = recovered_sad &
                    / (recovered_radius / DIFFUSION_CM2_S &
                       + 4.0_kind_phys / (SPEED_CM_S * GAMMA))
        ASSERT(abs(micm_rate - cam_rate) <= RTOL * cam_rate)
      end do
    end do
  end subroutine test_identity_with_cam_rate

  ! Layers at and above the tropopause layer index must be zero for both
  ! parameters, matching CAM's tropospheric-only surface area density, even
  ! though the wet diameter is nonzero there.
  subroutine test_tropopause_mask()
    integer, parameter :: NUM_COLUMNS = 3
    integer, parameter :: NUM_LAYERS = 4
    real(kind_phys) :: sad(NUM_COLUMNS,NUM_LAYERS), diameter(NUM_COLUMNS,NUM_LAYERS)
    real(kind_phys) :: radius(NUM_COLUMNS,NUM_LAYERS), number(NUM_COLUMNS,NUM_LAYERS)
    integer         :: tropopause_layer(NUM_COLUMNS), i_column, i_layer

    sad(:,:) = 1.0e-6_kind_phys
    diameter(:,:) = 2.0e-5_kind_phys
    tropopause_layer(:) = (/ 0, 2, 4 /)

    call surface_reaction_parameters_from_mode(sad, diameter, tropopause_layer, radius, number)

    do i_layer = 1, NUM_LAYERS
      do i_column = 1, NUM_COLUMNS
        if (i_layer <= tropopause_layer(i_column)) then
          ASSERT(radius(i_column, i_layer) == 0.0_kind_phys)
          ASSERT(number(i_column, i_layer) == 0.0_kind_phys)
        else
          ASSERT(radius(i_column, i_layer) > 0.0_kind_phys)
          ASSERT(number(i_column, i_layer) > 0.0_kind_phys)
        end if
      end do
    end do
  end subroutine test_tropopause_mask

  ! A mode with no surface (e.g. primary carbon, whose surface surf_area_dens
  ! leaves at zero) or no diameter must give zeros, never 0/0 or division by zero.
  subroutine test_no_surface_or_diameter()
    real(kind_phys) :: sad(3,1), diameter(3,1), radius(3,1), number(3,1)
    integer         :: tropopause_layer(3)

    sad(:,1) = (/ 0.0_kind_phys, 1.0e-6_kind_phys, 0.0_kind_phys /)
    diameter(:,1) = (/ 2.0e-5_kind_phys, 0.0_kind_phys, 0.0_kind_phys /)
    tropopause_layer(:) = 0

    call surface_reaction_parameters_from_mode(sad, diameter, tropopause_layer, radius, number)

    ASSERT(all(radius(:,1) == 0.0_kind_phys))
    ASSERT(all(number(:,1) == 0.0_kind_phys))
  end subroutine test_no_surface_or_diameter

  ! Surface reactions are found from their effective radius parameter and
  ! paired with their particle number concentration, regardless of order and
  ! interleaved with other rate parameters.
  subroutine test_find_surface_reactions()
    character(len=64) :: names(5)
    integer           :: indices(5)
    character(len=:), allocatable :: labels(:)
    integer,          allocatable :: radius_index(:), number_index(:)
    character(len=512) :: errmsg
    integer            :: errcode

    names(1) = 'PHOTO.O3_1'
    names(2) = 'SURF.N2O5_a2.particle number concentration [# m-3]'
    names(3) = 'SURF.N2O5_a1.effective radius [m]'
    names(4) = 'SURF.N2O5_a1.particle number concentration [# m-3]'
    names(5) = 'SURF.N2O5_a2.effective radius [m]'
    indices(:) = (/ 1, 2, 3, 4, 5 /)

    call find_surface_reactions(names, indices, labels, radius_index, number_index, &
                                errmsg, errcode)

    ASSERT(errcode == 0)
    ASSERT(size(labels) == 2)
    ASSERT(trim(labels(1)) == 'SURF.N2O5_a1')
    ASSERT(radius_index(1) == 3)
    ASSERT(number_index(1) == 4)
    ASSERT(trim(labels(2)) == 'SURF.N2O5_a2')
    ASSERT(radius_index(2) == 5)
    ASSERT(number_index(2) == 2)

    ! no surface reactions is not an error
    call find_surface_reactions(names(1:1), indices(1:1), labels, radius_index, &
                                number_index, errmsg, errcode)
    ASSERT(errcode == 0)
    ASSERT(size(labels) == 0)
  end subroutine test_find_surface_reactions

  subroutine test_find_surface_reactions_missing_partner()
    character(len=64) :: names(1)
    integer           :: indices(1)
    character(len=:), allocatable :: labels(:)
    integer,          allocatable :: radius_index(:), number_index(:)
    character(len=512) :: errmsg
    integer            :: errcode

    names(1) = 'SURF.N2O5_a1.effective radius [m]'
    indices(1) = 1

    call find_surface_reactions(names, indices, labels, radius_index, number_index, &
                                errmsg, errcode)

    ASSERT(errcode /= 0)
    ASSERT(index(errmsg, 'SURF.N2O5_a1') > 0)
  end subroutine test_find_surface_reactions_missing_partner

  ! The '_a<mm>' suffix convention: CAM's per-mode naming
  subroutine test_mode_from_label()
    ASSERT(mode_from_label('SURF.N2O5_hydrolysis_a1') == 1)
    ASSERT(mode_from_label('SURF.N2O5_hydrolysis_a04') == 4)
    ASSERT(mode_from_label('SURF.N2O5_hydrolysis_a12   ') == 12)
    ASSERT(mode_from_label('SURF.N2O5_hydrolysis') == 0)
    ASSERT(mode_from_label('SURF.N2O5_hydrolysis_a') == 0)
    ASSERT(mode_from_label('SURF.N2O5_hydrolysis_b1') == 0)
    ASSERT(mode_from_label('SURF.N2O5_hydrolysis1') == 0)
    ASSERT(mode_from_label('SURF.N2O5_a1234') == 0)
    ASSERT(mode_from_label('1') == 0)
    ASSERT(mode_from_label('') == 0)
  end subroutine test_mode_from_label

end program test_micm_surface
