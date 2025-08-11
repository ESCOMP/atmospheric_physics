!-------------------------------------------------------------------------------
! This module uses the solar irradiance data 
! to provide a spectral scaling factor
! to approximate the spectral distribution of irradiance
! when the radiation scheme might use a different solar source function
!-------------------------------------------------------------------------------
! peverwhee - dependencies = radiation_utils, mo_util
module rrtmgp_sw_solar_var

  use ccpp_kinds,        only : kind_phys

  implicit none
  save

  private
  public :: rrtmgp_sw_solar_var_init
  public :: rrtmgp_sw_solar_var_run

  real(kind_phys), allocatable :: irrad(:)           ! solar irradiance at model timestep in each band

  real(kind_phys), allocatable :: radbinmax(:)
  real(kind_phys), allocatable :: radbinmin(:)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> \section arg_table_rrtmgp_sw_solar_var_init Argument Table
!! \htmlinclude rrtmgp_sw_solar_var_init.html
!!
  subroutine rrtmgp_sw_solar_var_init(nswbands, do_spctrl_scaling, has_spectrum, errmsg, errflg)
    use radiation_utils,   only : get_sw_spectral_boundaries_ccpp
    integer, intent(in) :: nswbands            ! number of shortwave bands
    logical, intent(in) :: do_spctrl_scaling   ! flag to do spectral scaling
    logical, intent(in) :: has_spectrum        ! flag for whether solar input file has irradiance spectrum
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: radmax_loc
    character(len=256) :: alloc_errmsg

    if ( do_spctrl_scaling ) then

       if ( .not.has_spectrum ) then
          write(errmsg, *) 'rrtmgp_sw_solar_var_init: solar input fil must have irradiance spectrum'
          errflg = 1
          return
       endif

       allocate (radbinmax(nswbands),stat=errflg,errmsg=alloc_errmsg)
       if (errflg /= 0) then
          write(errmsg,*) 'rrtmgp_sw_solar_var_init: Error allocating space for radbinmax - message: ', alloc_errmsg
          return
       end if

       allocate (radbinmin(nswbands),stat=errflg,errmsg=alloc_errmsg)
       if (errflg /= 0) then
          write(errmsg,*) 'rrtmgp_sw_solar_var_init: Error allocating space for radbinmin - message: ', alloc_errmsg
          return
       end if

       allocate (irrad(nswbands), stat=errflg, errmsg=alloc_errmsg)
       if (errflg /= 0) then
          write(errmsg,*) 'rrtmgp_sw_solar_var_init: Error allocating space for irrad - message: ', alloc_errmsg
          return
       end if

       call get_sw_spectral_boundaries_ccpp(radbinmin, radbinmax, 'nm', errmsg, errflg)
       if (errflg /= 0) then
          return
       end if

       ! Make sure that the far-IR is included, even if radiation grid does not
       ! extend that far down. 10^5 nm corresponds to a wavenumber of
       ! 100 cm^-1.
       radmax_loc = maxloc(radbinmax,1)
       radbinmax(radmax_loc) = max(100000._kind_phys,radbinmax(radmax_loc))

    endif

  end subroutine rrtmgp_sw_solar_var_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!> \section arg_table_rrtmgp_sw_solar_var_run Argument Table
!! \htmlinclude rrtmgp_sw_solar_var_run.html
!!
  subroutine rrtmgp_sw_solar_var_run(toa_flux, band2gpt_sw, nswbands, sol_irrad, we, nbins, sol_tsi, do_spctrl_scaling, &
                                     sfac, eccf, errmsg, errflg)

     ! Arguments 
     real(kind_phys),    intent(inout) :: toa_flux(:,:)         ! top-of-atmosphere flux to be scaled (columns,gpts)
     real(kind_phys),    intent(in)    :: sol_tsi               ! total solar irradiance
     real(kind_phys),    intent(in)    :: sol_irrad(:)          ! solar irradiance
     real(kind_phys),    intent(in)    :: we(:)                 ! wavelength endpoints
     integer,            intent(in)    :: nbins                 ! number of bins
     integer,            intent(in)    :: band2gpt_sw(:,:)      ! array for converting shortwave band limits to g-points
     integer,            intent(in)    :: nswbands              ! number of shortwave bands
     logical,            intent(in)    :: do_spctrl_scaling     ! flag to do spectral scaling
     real(kind_phys),    intent(in)    :: eccf                  ! eccentricity factor
     real(kind_phys),    intent(out)   :: sfac(:,:)             ! scaling factors (columns,gpts)
     character(len=512), intent(out)   :: errmsg
     integer,            intent(out)   :: errflg

     ! Local variables 
     integer :: i, j, gpt_start, gpt_end, ncols
     real(kind_phys), allocatable :: scale(:)
     character(len=256)          :: alloc_errmsg
     character(len=*), parameter :: sub = 'rrtmgp_sw_solar_var_run'

     ! Initialize error variables
     errflg = 0
     errmsg = ''
    
     if (do_spctrl_scaling) then 

        ! Determine target irradiance for each band
        call integrate_spectrum(nbins, nswbands, we, radbinmin, radbinmax, sol_irrad, irrad)

        ncols = size(toa_flux, 1)
        allocate(scale(ncols), stat=errflg, errmsg=alloc_errmsg)
        if (errflg /= 0) then
           write(errmsg,*) sub, ': Error allocating "scale", message - ', alloc_errmsg
           errflg = 1
           return
        end if

        do i = 1, nswbands 
           gpt_start = band2gpt_sw(1,i) 
           gpt_end   = band2gpt_sw(2,i) 
           scale = spread(irrad(i), 1, ncols) / sum(toa_flux(:, gpt_start:gpt_end), dim=2)
           do j = gpt_start, gpt_end
              sfac(:,j) = scale
           end do
        end do

     else 
        sfac(:,:) = sol_tsi / spread(sum(toa_flux, 2), 2, size(toa_flux, 2))
     end if

     toa_flux = toa_flux * sfac * eccf

  end subroutine rrtmgp_sw_solar_var_run


!-------------------------------------------------------------------------------
! private method.........
!-------------------------------------------------------------------------------

  subroutine integrate_spectrum( nsrc, ntrg, src_x, min_trg, max_trg, src, trg )

    use mo_util, only : rebin

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer,  intent(in)  :: nsrc                  ! dimension source array
    integer,  intent(in)  :: ntrg                  ! dimension target array
    real(kind_phys), intent(in)  :: src_x(nsrc+1)         ! source coordinates
    real(kind_phys), intent(in)  :: max_trg(ntrg)         ! target coordinates
    real(kind_phys), intent(in)  :: min_trg(ntrg)         ! target coordinates
    real(kind_phys), intent(in)  :: src(nsrc)             ! source array
    real(kind_phys), intent(out) :: trg(ntrg)             ! target array
 
    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    real(kind_phys) :: trg_x(2), targ(1)         ! target coordinates
    integer  :: i

    do i = 1, ntrg

       trg_x(1) = min_trg(i)
       trg_x(2) = max_trg(i)

       call rebin( nsrc, 1, src_x, trg_x, src(1:nsrc), targ(:) )
       ! W/m2/nm --> W/m2
       trg( i ) = targ(1)*(trg_x(2)-trg_x(1))

    enddo


  end subroutine integrate_spectrum

end module rrtmgp_sw_solar_var
