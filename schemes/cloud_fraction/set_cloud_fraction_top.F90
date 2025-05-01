! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Stub scheme to set top of cloud physics to below top cloud level.
! Used for all macrophysical schemes except RK.
module set_cloud_fraction_top
  implicit none
  private

  public :: set_cloud_fraction_top_init

contains

!> \section arg_table_set_cloud_fraction_top_init Argument Table
!! \htmlinclude set_cloud_fraction_top_init.html
  subroutine set_cloud_fraction_top_init(trop_cloud_top_lev, top_lev, errmsg, errflg)
    integer,            intent(in)  :: trop_cloud_top_lev ! Troposphere cloud physics top level [index]
    integer,            intent(out) :: top_lev            ! Cloud fraction top level [index]
    character(len=512), intent(out) :: errmsg             ! Error message
    integer,            intent(out) :: errflg             ! Error flag

    ! Initialize error handling
    errmsg = ''
    errflg = 0

    ! Limit CAM5+ cloud physics to below top cloud level.
    top_lev = trop_cloud_top_lev

 end subroutine set_cloud_fraction_top_init

end module set_cloud_fraction_top
