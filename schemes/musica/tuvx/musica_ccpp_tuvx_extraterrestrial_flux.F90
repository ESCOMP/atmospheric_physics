module musica_ccpp_tuvx_extraterrestrial_flux
    use ccpp_kinds, only: kind_phys

    implicit none
  
    private :: rebin
    public :: create_extraterrestrial_flux_profile, set_extraterrestrial_flux_values
  
    !> Label for extraterrestrial_flux in TUV-x
    character(len=*), parameter, public :: extraterrestrial_flux_label = "extraterrestrial flux"
    !> Unit for extraterrestrial_flux in TUV-x
    character(len=*), parameter, public :: extraterrestrial_flux_unit = "photon cm-2 s-1"
    !> Wavelength grid interface values
    real(kind_phys), protected, allocatable :: wavelength_grid_interfaces_(:) ! nm
    !> Default value of number of wavelength grid bins
    integer, parameter :: DEFAULT_NUM_WAVELENGTH_BINS = 0
    !> Number of wavelength grid bins
    integer, protected :: num_wavelength_bins_ = DEFAULT_NUM_WAVELENGTH_BINS


contains

  !> Regrids normalized flux data to match a specified wavelength grid
  !! This function is copied from CAM/src/chemistry/utils/mo_util.F90
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
        do sil = 1, source_dimension + 1
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
        sil = max( sil, 2 )
        siu = min( siu, source_dimension + 1 )
        do si = sil, siu
          si1 = si - 1
          sl  = max( tl, source_coordinates( si1 ) )
          su  = min( tu, source_coordinates( si ) )
          y   = y + ( su - sl ) * source( si1 )
        end do
        target(i) = y / (target_coordinates( i + 1 ) - target_coordinates( i ) )
      else
        target(i) = 0._kind_phys
      end if
    end do

  end subroutine rebin

  !> Creates a TUV-x extraterrestrial flux profile from the host-model wavelength grid
  function create_extraterrestrial_flux_profile(wavelength_grid, errmsg, errcode) &
      result( profile )
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_grid,    only: grid_t
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t

    type(grid_t),     intent(inout) :: wavelength_grid
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode
    type(profile_t),  pointer       :: profile

    ! local variables
    type(error_t) :: error

    num_wavelength_bins_ = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    allocate(wavelength_grid_interfaces_(num_wavelength_bins_ + 1))

    call wavelength_grid%get_edges(wavelength_grid_interfaces_, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) then
      deallocate(wavelength_grid_interfaces_)
      return
    end if

    wavelength_grid_interfaces_ = wavelength_grid_interfaces_ ! nm

    profile => profile_t( extraterrestrial_flux_label, extraterrestrial_flux_unit, &
                          wavelength_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) then
      deallocate(wavelength_grid_interfaces_)
      return
    end if

  end function create_extraterrestrial_flux_profile

  !> Sets TUV-x extraterrestrial flux midpoints
  !
  ! Extraterrestrial flux is read from data files and interpolated to the
  ! TUV-x wavelength grid. CAM ET Flux values are multiplied by the
  ! width of the wavelength bins to get the TUV-x units of photon cm-2 s-1
  !
  ! TUV-x only uses mid-point values for extraterrestrial flux
  subroutine set_extraterrestrial_flux_values(profile, from_data_num_wavelength_grid_bins, &
      from_data_wavelength_grid_interfaces, from_data_extraterrestrial_flux, errmsg, errcode)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t
    use ccpp_kinds,          only: kind_phys

    type(profile_t),  intent(inout) :: profile
    integer,          intent(in)    :: from_data_num_wavelength_grid_bins      ! (count)
    real(kind_phys),  intent(in)    :: from_data_wavelength_grid_interfaces(:) ! nm
    real(kind_phys),  intent(in)    :: from_data_extraterrestrial_flux(:)      ! photons cm-2 s-1 nm-1
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: midpoints(num_wavelength_bins_)

    if (num_wavelength_bins_ <= DEFAULT_NUM_WAVELENGTH_BINS) then
      errmsg = "[MUSICA Error] Invalid size of TUV-x wavelength bins."
      errcode = 1
      return
    end if

    if (.not. allocated(wavelength_grid_interfaces_)) then
      errmsg = "[MUSICA Error] Failed to allocate the TUV-x wavelength grid interface array"
      errcode = 1
      return
    end if

    ! Regrid normalized flux to TUV-x wavelength grid
    call rebin( from_data_num_wavelength_grid_bins, num_wavelength_bins_,     &
                from_data_wavelength_grid_interfaces, wavelength_grid_interfaces_, &
                from_data_extraterrestrial_flux, midpoints )

    ! Convert normalized flux to flux on TUV-x wavelength grid
    midpoints = midpoints * ( wavelength_grid_interfaces_(2 : num_wavelength_bins_ + 1) &
                 - wavelength_grid_interfaces_(1 :num_wavelength_bins_) )

    call profile%set_midpoint_values( midpoints, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    deallocate(wavelength_grid_interfaces_)

  end subroutine set_extraterrestrial_flux_values

end module musica_ccpp_tuvx_extraterrestrial_flux