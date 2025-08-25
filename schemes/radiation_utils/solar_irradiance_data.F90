!-------------------------------------------------------------------------------
! This module uses the solar irradiance data 
! to provide a spectral scaling factor
! to approximate the spectral distribution of irradiance
! when the radiation scheme might use a different solar source function
!-------------------------------------------------------------------------------
! peverwhee - dependencies = radiation_utils, mo_util
module solar_irradiance_data

  use ccpp_kinds,        only : kind_phys

  implicit none
  save

  private
  public :: solar_irradiance_data_init
  public :: solar_irradiance_data_run

  real(kind_phys), allocatable :: irrad(:)           ! solar irradiance at model timestep in each band

  real(kind_phys), allocatable :: radbinmax(:)
  real(kind_phys), allocatable :: radbinmin(:)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

  subroutine solar_irradiance_data_init(irrad_file_path, nswbands, do_spctrl_scaling, has_spectrum, errmsg, errflg)
    use radiation_utils,   only : get_sw_spectral_boundaries_ccpp
    integer, intent(in) :: nswbands            ! number of shortwave bands
    logical, intent(in) :: do_spctrl_scaling   ! flag to do spectral scaling
    logical, intent(in) :: has_spectrum        ! flag for whether solar input file has irradiance spectrum
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: radmax_loc
    character(len=256) :: alloc_errmsg



  end subroutine solar_irradiance_data_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  subroutine solar_irradiance_data_run(toa_flux, band2gpt_sw, nswbands, sol_irrad, we, nbins, sol_tsi, do_spctrl_scaling, &
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
     character(len=*), parameter :: sub = 'solar_irradiance_data_run'

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

  end subroutine solar_irradiance_data_run

end module solar_irradiance_data
