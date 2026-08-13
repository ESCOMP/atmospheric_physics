module micro_pumas_ccpp_optics_limiter

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: micro_pumas_ccpp_optics_limiter_run

  !Convective size distribution parameters, used in place of the stratiform
  !parameters where the stratiform cloud fraction is negligible.
  !Mirrors CAM micro_pumas_cam.F90.
  real(kind_phys), parameter :: dcon   = 25.e-6_kind_phys !Convective size distribution effective radius (m)
  real(kind_phys), parameter :: mucon  = 5.3_kind_phys    !Convective size distribution shape parameter (1)
  real(kind_phys), parameter :: deicon = 50._kind_phys    !Convective ice effective diameter (um)

  !Stratiform cloud fraction below which the stratiform size distribution is
  !discarded. PUMAS divides cloud water by the cloud fraction to obtain the
  !in-cloud values these parameters are derived from, so as the cloud fraction
  !vanishes the slope and diameters it produces grow without bound.
  real(kind_phys), parameter :: cldfrc_min = 1.e-4_kind_phys

contains
  !> \section arg_table_micro_pumas_ccpp_optics_limiter_run Argument Table
  !! \htmlinclude micro_pumas_ccpp_optics_limiter_run.html
  subroutine micro_pumas_ccpp_optics_limiter_run(ncol, nlev, trop_cloud_top_lev,   &
                   strat_cldfrc, pgamrad, lamcrad, deffi, errmsg, errcode)

    integer,            intent(in)    :: ncol               !Number of horizontal columns (count)
    integer,            intent(in)    :: nlev               !Number of vertical layers (count)
    integer,            intent(in)    :: trop_cloud_top_lev !Index of the top model level for which
                                                            !tropospheric cloud physics is run (index)
    real(kind_phys),    intent(in)    :: strat_cldfrc(:, :) !Stratiform cloud area fraction (fraction)
    real(kind_phys),    intent(inout) :: pgamrad(:, :)      !Size distribution shape parameter (1)
    real(kind_phys),    intent(inout) :: lamcrad(:, :)      !Slope of droplet distribution (m-1)
    real(kind_phys),    intent(inout) :: deffi(:, :)        !Effective diameter of cloud ice particle (um)
    character(len=*),   intent(out)   :: errmsg             !CCPP error message (none)
    integer,            intent(out)   :: errcode            !CCPP error code (1)

    integer :: i, k

    errmsg  = ''
    errcode = 0

    !lamcrad is an inverse length; it is declared dimensionless in the metadata
    !to match the RRTMGP and host declarations, so dcon must stay in meters for
    !the value handed to radiation to be consistent with CAM.
    do k = trop_cloud_top_lev, nlev
       do i = 1, ncol
          if (strat_cldfrc(i, k) < cldfrc_min) then
             pgamrad(i, k) = mucon
             lamcrad(i, k) = (mucon + 1._kind_phys) / dcon
             deffi(i, k)   = deicon
          end if
       end do
    end do

  end subroutine micro_pumas_ccpp_optics_limiter_run

end module micro_pumas_ccpp_optics_limiter
