module check_energy_scaling
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: check_energy_scaling_run

contains

  ! CCPP routine to get scaling factor for conversion of temperature increment.
  ! This is extracted to a separate subroutine so that scaling_dycore can be passed
  ! directly to the CCPP-ized check_energy_chng_run subroutine from CAM with subcolumns.
  !
  ! When subcolumns are removed from CAM, this dummy scheme can be removed, and
  ! scaling_dycore can just be calculated in check_energy_chng. (hplin, 9/5/24)
!> \section arg_table_check_energy_scaling_run Argument Table
!! \htmlinclude arg_table_check_energy_scaling_run.html
  subroutine check_energy_scaling_run( &
       ncol,                           &
       cp_or_cv_dycore, cpairv,        &
       scaling_dycore,                 &
       errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol                  ! number of atmospheric columns
    real(kind_phys),    intent(in)    :: cp_or_cv_dycore(:,:)  ! cp or cv from dycore [J kg-1 K-1]
    real(kind_phys),    intent(in)    :: cpairv(:,:)           ! specific heat of dry air at constant pressure [J kg-1 K-1]

    ! Output arguments
    real(kind_phys),    intent(out)   :: scaling_dycore(:,:)   ! scaling for conversion of temperature increment [1]
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    scaling_dycore(:ncol,:) = cpairv(:ncol,:) / cp_or_cv_dycore(:ncol,:)

  end subroutine check_energy_scaling_run

end module check_energy_scaling
