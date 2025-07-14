!> \file rrtmgp_sw_main.F90
!!

!> This module contains the call to the RRTMGP-sw radiation routine
module rrtmgp_sw_main
  implicit none
  private

  public rrtmgp_sw_main_run
contains

!> \section arg_table_rrtmgp_sw_main_run Argument Table
!! \htmlinclude rrtmgp_sw_main_run.html
!!
   subroutine rrtmgp_sw_main_run(doswrad, doswclrsky, doswallsky, nday, iter_num, rrtmgp_phys_blksz, sw_optical_props, &
                                 sw_optical_props_clouds, top_at_1, aersw, coszen, toa_src_sw,      &
                                 sfc_alb_dir, sfc_alb_dif, flux_clrsky, flux_allsky, errmsg, errflg)
    use machine,                  only: kind_phys
    use mo_rte_sw,                only: rte_sw
    use ccpp_optical_props,       only: ty_optical_props_2str_ccpp
    use ccpp_fluxes_byband,       only: ty_fluxes_byband_ccpp
    use ccpp_fluxes,              only: ty_fluxes_broadband_ccpp
    use radiation_tools,          only: check_error_msg

    ! Inputs
    logical, intent(in) :: doswrad                                              !< Flag to perform shortwave calculation
    logical, intent(in) :: doswclrsky                                           !< Flag to compute clear-sky fluxes
    logical, intent(in) :: doswallsky                                           !< Flag to compute all-sky fluxes
    logical, intent(in) :: top_at_1                                             !< Flag for vertical ordering convention

    integer, intent(in) :: nday                                                 !< Number of horizontal daylight points
    integer, intent(in) :: iter_num                                             !< Radiation subcycle iteration number
    integer, intent(in) :: rrtmgp_phys_blksz                                    !< Number of horizontal points to process at once

    real(kind_phys), dimension(:,:),   intent(in) :: toa_src_sw                 !< Top-of-atmosphere flux on g-points [W m-2]
    real(kind_phys), dimension(:,:),   intent(in) :: sfc_alb_dir                !< Albedo direct at surface [fraction]
    real(kind_phys), dimension(:,:),   intent(in) :: sfc_alb_dif                !< Albedo diffuse at surface [fraction]
    real(kind_phys), dimension(:),     intent(in) :: coszen                     !< Cosine of solar zenith angle for daytime points

    ! Outputs
    class(ty_fluxes_byband_ccpp),      intent(inout) :: flux_allsky             !< All-sky flux [W m-2]
    class(ty_fluxes_broadband_ccpp),   intent(inout) :: flux_clrsky             !< Clear-sky flux [W m-2]
    class(ty_optical_props_2str_ccpp), intent(inout) :: aersw                   !< Aerosol optical properties object
    class(ty_optical_props_2str_ccpp), intent(inout) :: sw_optical_props        !< Clear-sky optical properties object
    class(ty_optical_props_2str_ccpp), intent(inout) :: sw_optical_props_clouds !< Cloud optical properties object

    character(len=*), intent(out) :: errmsg                                     !< CCPP error message
    integer,          intent(out) :: errflg                                     !< CCPP error flag

    ! Local variables
    integer :: iCol, iCol2 

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doswrad) return

    iCol = ((iter_num - 1) * rrtmgp_phys_blksz) + 1
    iCol2 = min(iCol + rrtmgp_phys_blksz - 1, nday)

    ! ###################################################################################
    !
    ! Compute clear-sky fluxes (gaseous+aerosol)
    !
    ! ###################################################################################
    ! Increment optics (always)
    errmsg = aersw%optical_props%increment(sw_optical_props%optical_props)
    call check_error_msg('rrtmgp_sw_main_increment_aerosol_to_clrsky', errmsg)
    if (len_trim(errmsg) /= 0) then
        errflg = 1
        return
    end if

    ! Optionally compute clear-sky fluxes
    if (doswclrsky) then
       call check_error_msg('rrtmgp_sw_main_rte_sw_clrsky',rte_sw(     &
                  sw_optical_props%optical_props,    & ! IN  - optical-properties
                  top_at_1,                          & ! IN  - veritcal ordering flag
                  coszen(iCol:iCol2),                      & ! IN  - Cosine of solar zenith angle
                  toa_src_sw,                              & ! IN  - incident solar flux at TOA
                  sfc_alb_dir,                             & ! IN  - Shortwave surface albedo (direct)
                  sfc_alb_dif,                             & ! IN  - Shortwave surface albedo (diffuse)
                  flux_clrsky%fluxes))                       ! OUT - Fluxes, clear-sky, 3D (1,nLay,nBand)
       if (len_trim(errmsg) /= 0) then
           errflg = 1
           return
       end if
    end if

    ! ###################################################################################
    !
    ! All-sky fluxes (clear-sky + clouds + precipitation)
    ! 
    ! ###################################################################################

    if (doswallsky) then
       ! Delta scale
       !call check_error_msg('rrtmgp_sw_main_delta_scale',sw_optical_props_clouds%delta_scale())

       ! Increment
       call check_error_msg('rrtmgp_sw_main_increment_clouds_to_clrsky', &
            sw_optical_props_clouds%optical_props%increment(sw_optical_props%optical_props))

       ! Compute fluxes
       call check_error_msg('rrtmgp_sw_main_rte_sw_allsky',rte_sw(     &
            sw_optical_props%optical_props,  & ! IN  - optical-properties
            top_at_1,                        & ! IN  - veritcal ordering flag
            coszen(iCol:iCol2),              & ! IN  - Cosine of solar zenith angle
            toa_src_sw,                      & ! IN  - incident solar flux at TOA
            sfc_alb_dir,                     & ! IN  - Shortwave surface albedo (direct)
            sfc_alb_dif,                     & ! IN  - Shortwave surface albedo (diffuse)
            flux_allsky%fluxes))               ! OUT - Fluxes, clear-sky, 3D (1,nLay,nBand)

       if (len_trim(errmsg) /= 0) then
          errflg = 1
       end if
    end if

  end subroutine rrtmgp_sw_main_run
end module rrtmgp_sw_main
