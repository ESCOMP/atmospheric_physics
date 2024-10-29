module musica_ccpp_tuvx_extraterrestrial_flux
    implicit none
  
    private :: rebin
    public :: create_extraterrestrial_flux_profile, set_extraterrestrial_flux_values, &
              extraterrestrial_flux_label, extraterrestrial_flux_unit 
  
    !> Label for extraterrestrial_flux in TUV-x
    character(len=*), parameter :: extraterrestrial_flux_label = "extraterrestrial flux"
    !> Unit for extraterrestrial_flux in TUV-x
    character(len=*), parameter :: extraterrestrial_flux_unit = "photon cm-2 s-1"
    !> Default value of number of wavelength bins
    integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
    !> Number of wavelength bins
    integer, protected :: num_wavelength_bins_ = DEFAULT_NUM_WAVELENGTH_BINS

contains

  !> Regrids normalized flux to TUV-x wavelength grid
  subroutine rebin( source_dimension, target_dimension, source_coordinates, &
                    target_coordinates, source, target )
    use ccpp_kinds, only: kind_phys

    integer,         intent(in)  :: source_dimension
    integer,         intent(in)  :: target_dimension
    real(kind_phys), intent(in)  :: source_coordinates(source_dimension+1)
    real(kind_phys), intent(in)  :: target_coordinates(target_dimension+1)
    real(kind_phys), intent(in)  :: source(source_dimension)
    real(kind_phys), intent(out) :: target(target_dimension)

    ! local variables
    integer         :: i, si, si1, sil, siu
    real(kind_phys) :: y, sl, su, tl, tu

    do i = 1, target_dimension
      tl = target_coordinates(i)
      if( tl < source_coordinates( source_dimension + 1) ) then
        do sil = 1, source_dimension+1
          if( tl <= source_coordinates( sil ) ) then
            exit
          end if
        end do
        tu = target_coordinates( i + 1 )
        do siu = 1, source_dimension + 1
          if( tu <= source_coordinates( siu ) ) then
            exit
          end if
        end do
        y   = 0._kind_phys
        sil = max( sil,2 )
        siu = min( siu, source_dimension + 1 )
        do si = sil, siu
          si1 = si - 1
          sl  = max( tl, source_coordinates( si1 ) )
          su  = min( tu,source_coordinates( si ) )
          y   = y + ( su - sl ) * source( si1 )
        end do
        target(i) = y / (target_coordinates( i + 1 ) - target_coordinates( i ) )
      else
        target(i) = 0._kind_phys
      end if
    end do

  end subroutine rebin

  !> Creates a TUV-x extraterrestrial flux profile from the host-model wavelength grid
  function create_extraterrestrial_flux_profile( wavelength_grid, errmsg, errcode ) & 
      result( profile )
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(grid_t),     intent(in)  :: wavelength_grid
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode
    type(profile_t),  pointer     :: profile

    ! local variables
    type(error_t) :: error

    num_wavelength_bins_ = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    profile => profile_t( extraterrestrial_flux_label, extraterrestrial_flux_unit, &
                          wavelength_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end function create_extraterrestrial_flux_profile

  !> Sets TUV-x extraterrestrial flux midpoints
  !
  ! Extraterrestrial flux is read from data files and interpolated to the
  ! TUV-x wavelength grid. CAM ET Flux values are multiplied by the
  ! width of the wavelength bins to get the TUV-x units of photon cm-2 s-1
  !
  ! TUV-x only uses mid-point values for extraterrestrial flux
  subroutine set_extraterrestrial_flux( profile, wavelength_grid, errmsg, errcode )
    use ccpp_kinds,          only: kind_phys
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t
    use mo_util,             only: rebin
    use solar_irrad_data,    only: nbins,  & ! number of wavelength bins
                                   we,     & ! wavelength bin edges
                                   sol_etf   ! extraterrestrial flux, photon cm-2 nm-1 s-1

    type(profile_t),  intent(inout) :: profile
    type(grid_t),     intent(in)    :: wavelength_grid
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: wavelegnth_interfaces(num_wavelength_bins_ + 1) 
    real(kind_phys) :: extraterrestrial_flux(nbins)
    real(kind_phys) :: midpoints(num_wavelength_bins_) ! photon cm-2 nm-1 s-1
    integer         :: i_bin

    if (num_wavelength_bins_ <= DEFAULT_NUM_WAVELENGTH_BINS) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x wavelength bins."
      errcode = 1
      return
    end if
    
    call wavelength_grid%get_edge_values(wavelegnth_interfaces, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    ! regrid normalized flux to TUV-x wavelength grid
    extraterrestrial_flux(:) = sol_etf(:)
    call rebin( nbins, num_wavelength_bins_, we, wavelegnth_interfaces, &
                extraterrestrial_flux, midpoints )

    ! convert normalized flux to flux on TUV-x wavelength grid
    midpoints(:) = midpoints(:) * ( wavelegnth_interfaces(2 : num_wavelength_bins + 1) &
                 - wavelegnth_interfaces(1 :num_wavelength_bins_) )

    call profile%set_midpoints( midpoints, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) return

  end subroutine set_extraterrestrial_flux_values

end module musica_ccpp_tuvx_extraterrestrial_flux