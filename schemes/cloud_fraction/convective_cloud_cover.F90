! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Computes deep and shallow convective cloud fractions
! CCPP-ized: Haipeng Lin, February 2025
module convective_cloud_cover

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: convective_cloud_cover_init
  public :: convective_cloud_cover_run

  ! Module variables for tuning parameters
  real(kind_phys) :: sh1 ! Shallow convection tuning parameter 1 [fraction]
  real(kind_phys) :: sh2 ! Shallow convection tuning parameter 2 [m2 s kg-1]
  real(kind_phys) :: dp1 ! Deep convection    tuning parameter 1 [fraction]
  real(kind_phys) :: dp2 ! Deep convection    tuning parameter 2 [m2 s kg-1]

contains

!> \section arg_table_convective_cloud_cover_init Argument Table
!! \htmlinclude convective_cloud_cover_init.html
  subroutine convective_cloud_cover_init( &
    amIRoot, iulog, &
    sh1_in, sh2_in, dp1_in, dp2_in, errmsg, errflg)

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog      ! log output unit
    real(kind_phys),    intent(in)  :: sh1_in     ! Shallow convection parameter 1
    real(kind_phys),    intent(in)  :: sh2_in     ! Shallow convection parameter 2
    real(kind_phys),    intent(in)  :: dp1_in     ! Deep convection parameter 1
    real(kind_phys),    intent(in)  :: dp2_in     ! Deep convection parameter 2
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Set module variables from namelist
    sh1 = sh1_in
    sh2 = sh2_in
    dp1 = dp1_in
    dp2 = dp2_in

    if(amIRoot) then
      write(iulog,*) 'tuning parameters convective_cloud_cover: dp1 ',dp1,' dp2 ',dp2,' sh1 ',sh1,' sh2 ',sh2
    endif

  end subroutine convective_cloud_cover_init

  ! Compute convective cloud cover (deep and shallow)
  ! Should produce typical numbers of 20%
  ! Shallow and deep convective cloudiness are evaluated separately and summed
  ! because the processes are evaluated separately.
!> \section arg_table_convective_cloud_cover_run Argument Table
!! \htmlinclude convective_cloud_cover_run.html
  subroutine convective_cloud_cover_run( &
    ncol, pver, &
    top_lev_cloudphys, &
    use_shfrc, shfrc, &
    cmfmc_total, cmfmc_sh, &
    shallowcu, deepcu, concld, &
    errmsg, errflg)

    ! Input arguments
    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: top_lev_cloudphys ! Top vertical level for cloud physics [index]

    logical,          intent(in)  :: use_shfrc         ! Use cloud area fraction provided by shallow convection? [flag]
    real(kind_phys),  intent(in)  :: shfrc(:, :)       ! Input shallow cloud fraction [fraction]

    real(kind_phys),  intent(in)  :: cmfmc_total(:, :) ! atmosphere_convective_mass_flux_due_to_all_convection [kg m-2 s-1]
    real(kind_phys),  intent(in)  :: cmfmc_sh(:, :)    ! atmosphere_convective_mass_flux_due_to_shallow_convection [kg m-2 s-1]

    ! Output arguments
    real(kind_phys),  intent(out) :: shallowcu(:, :)   ! Shallow convective cloud fraction [fraction]
    real(kind_phys),  intent(out) :: deepcu(:, :)      ! Deep convective cloud fraction [fraction]
    real(kind_phys),  intent(out) :: concld(:, :)      ! Convective cloud cover [fraction]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    concld(:ncol, :) = 0.0_kind_phys
    do k = top_lev_cloudphys, pver
      do i = 1, ncol
        if (.not. use_shfrc) then
          ! Compute the shallow convection cloud cover using the shallow convective mass flux.
          shallowcu(i, k) = max(0.0_kind_phys, &
                                min(sh1*log(1.0_kind_phys + sh2*cmfmc_sh(i, k + 1)), 0.30_kind_phys))
        else
          ! If use shallow convective calculated clouds, then just assign
          shallowcu(i, k) = shfrc(i, k)
        end if

        ! REMOVECAM: This could be changed to use deep convective mass flux
        ! since it is independently available in CCPP, once CAM is retired
        deepcu(i, k) = max(0.0_kind_phys, &
                           min(dp1 * &
                               log(1.0_kind_phys + &
                                   dp2 * (cmfmc_total(i, k + 1) - cmfmc_sh(i, k + 1)) &
                                ), &
                               0.60_kind_phys))

        ! Estimate of local convective cloud cover based on convective mass flux
        ! Modify local large-scale relative humidity to account for presence of
        ! convective cloud when evaluating relative humidity based layered cloud amount
        concld(i, k) = min(shallowcu(i, k) + deepcu(i, k), 0.80_kind_phys)
      end do
    end do

  end subroutine convective_cloud_cover_run

end module convective_cloud_cover
