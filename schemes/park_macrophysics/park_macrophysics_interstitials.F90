! Host-side interstitials for the Park (CAM5) cloud macrophysics.
! These reproduce the parts of CAM's macrop_driver_tend that are not part of the
! portable park_macrophysics scheme itself (the remaining "CAM interface").
module park_macrophysics_interstitials

  implicit none
  private

  public :: park_macrophysics_save_qtlcwat_run

contains

  ! Save the post-macrophysics equilibrium state for use as the "previous timestep"
  ! reference state at the next step (CAM macrop_driver.F90).
  ! Only levels top_lev:pver are updated, matching CAM; levels above top_lev retain
  ! their incoming values (hence intent(inout) on the saved fields).
!> \section arg_table_park_macrophysics_save_qtlcwat_run Argument Table
!! \htmlinclude park_macrophysics_save_qtlcwat_run.html
  subroutine park_macrophysics_save_qtlcwat_run( &
    ncol, pver, top_lev, &
    t, q_wv, cldice, cldliq, numliq, numice, &
    tcwat, qcwat, lcwat, iccwat, nlwat, niwat, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: top_lev
    real(kind_phys),    intent(in)    :: t(:, :)        ! Temperature [K]
    real(kind_phys),    intent(in)    :: q_wv(:, :)     ! Water vapor mixing ratio [kg kg-1]
    real(kind_phys),    intent(in)    :: cldice(:, :)   ! Cloud ice mixing ratio [kg kg-1]
    real(kind_phys),    intent(in)    :: cldliq(:, :)   ! Cloud liquid water mixing ratio [kg kg-1]
    real(kind_phys),    intent(in)    :: numliq(:, :)   ! Cloud liquid droplet number concentration [kg-1]
    real(kind_phys),    intent(in)    :: numice(:, :)   ! Cloud ice crystal number concentration [kg-1]

    ! Output arguments (saved equilibrium reference state for the next timestep)
    real(kind_phys),    intent(inout) :: tcwat(:, :)    ! Equilibrium reference temperature [K]
    real(kind_phys),    intent(inout) :: qcwat(:, :)    ! Equilibrium reference water vapor [kg kg-1]
    real(kind_phys),    intent(inout) :: lcwat(:, :)    ! Equilibrium reference total cloud water (liquid + ice) [kg kg-1]
    real(kind_phys),    intent(inout) :: iccwat(:, :)   ! Equilibrium reference cloud ice [kg kg-1]
    real(kind_phys),    intent(inout) :: nlwat(:, :)    ! Equilibrium reference liquid droplet number [kg-1]
    real(kind_phys),    intent(inout) :: niwat(:, :)    ! Equilibrium reference ice crystal number [kg-1]
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables
    integer :: k

    errmsg = ''
    errflg = 0

    do k = top_lev, pver
      tcwat(:ncol, k)  = t(:ncol, k)
      qcwat(:ncol, k)  = q_wv(:ncol, k)
      lcwat(:ncol, k)  = cldliq(:ncol, k) + cldice(:ncol, k)
      iccwat(:ncol, k) = cldice(:ncol, k)
      nlwat(:ncol, k)  = numliq(:ncol, k)
      niwat(:ncol, k)  = numice(:ncol, k)
    end do

  end subroutine park_macrophysics_save_qtlcwat_run

end module park_macrophysics_interstitials
