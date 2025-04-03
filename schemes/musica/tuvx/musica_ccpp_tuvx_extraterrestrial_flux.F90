! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_tuvx_extraterrestrial_flux
  use ccpp_kinds, only: kind_phys

  implicit none

  private
  save

  public :: create_extraterrestrial_flux_profile, set_extraterrestrial_flux_values, &
            cleanup_photolysis_wavelength_grid_interfaces

  type, private :: wavelength_grid_interfaces_t
    real(kind_phys), allocatable :: interfaces(:) ! nm
    integer                      :: size = 0
  contains
    procedure :: deallocate_data => wavelength_grid_interfaces_t_deallocate_data
  end type wavelength_grid_interfaces_t

  interface wavelength_grid_interfaces_t
    procedure wavelength_grid_interfaces_t_constructor
  end interface wavelength_grid_interfaces_t

  !> Label for extraterrestrial_flux in TUV-x
  character(len=*), parameter, public :: extraterrestrial_flux_label = "extraterrestrial flux"
  !> Unit for extraterrestrial_flux in TUV-x
  character(len=*), parameter, public :: extraterrestrial_flux_unit = "photon cm-2 s-1" ! photon cm-2 s-1 nm-1
  !> Wavelength grid interface values
  type(wavelength_grid_interfaces_t), protected, allocatable, public :: host_wavelength_grid_interfaces_ ! nm
  type(wavelength_grid_interfaces_t), protected, allocatable, public :: tuvx_wavelength_grid_interfaces_ ! nm

contains

  !> Constructor for wavelength grid interface object
  function wavelength_grid_interfaces_t_constructor(wavelength_grid_interfaces, size) &
      result( this )

    real(kind_phys), intent(in)        :: wavelength_grid_interfaces(:) ! nm
    integer,         intent(in)        :: size
    type(wavelength_grid_interfaces_t) :: this

    allocate( this%interfaces( size ) )
    this%interfaces(:) = wavelength_grid_interfaces(:)
    this%size = size

  end function wavelength_grid_interfaces_t_constructor

  !> Deallocates memory for interface array of wavelength grid interface object
  subroutine wavelength_grid_interfaces_t_deallocate_data(this)
    class(wavelength_grid_interfaces_t), intent(inout) :: this

    if (allocated(this%interfaces)) then
      deallocate(this%interfaces)
    end if

  end subroutine wavelength_grid_interfaces_t_deallocate_data

  !> Deallocates memory associated with photolysis wavelength grid interfaces
  subroutine cleanup_photolysis_wavelength_grid_interfaces()

    if (allocated( host_wavelength_grid_interfaces_ )) then
      call host_wavelength_grid_interfaces_%deallocate_data()
      deallocate( host_wavelength_grid_interfaces_ )
    end if

    if (allocated( tuvx_wavelength_grid_interfaces_ )) then
      call tuvx_wavelength_grid_interfaces_%deallocate_data()
      deallocate( tuvx_wavelength_grid_interfaces_ )
    end if

  end subroutine cleanup_photolysis_wavelength_grid_interfaces

  !> Creates a TUV-x extraterrestrial flux profile based on the TUV-x wavelength grid
  !  and initializes photolysis wavelength grid interfaces for each CAM-SIMA and TUV-x
  function create_extraterrestrial_flux_profile(wavelength_grid, &
      photolysis_wavelength_grid_interfaces, errmsg, errcode) result( profile )
    use musica_util,                      only: error_t
    use musica_ccpp_util,                 only: has_error_occurred
    use musica_tuvx_grid,                 only: grid_t
    use musica_tuvx_profile,              only: profile_t

    type(grid_t),     intent(inout) :: wavelength_grid
    real(kind_phys),  intent(in)    :: photolysis_wavelength_grid_interfaces(:) ! nm
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode
    type(profile_t),  pointer       :: profile

    ! local variables
    real(kind_phys), allocatable :: interfaces(:) ! nm
    integer                      :: num_wavelength_grid_sections
    type(error_t)                :: error

    profile => profile_t( extraterrestrial_flux_label, extraterrestrial_flux_unit, &
                          wavelength_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    num_wavelength_grid_sections = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    allocate( interfaces( num_wavelength_grid_sections + 1) )

    call wavelength_grid%get_edges( interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    allocate(tuvx_wavelength_grid_interfaces_)
    tuvx_wavelength_grid_interfaces_ = wavelength_grid_interfaces_t_constructor( &
      interfaces, num_wavelength_grid_sections + 1 )

    allocate(host_wavelength_grid_interfaces_)
    host_wavelength_grid_interfaces_ = wavelength_grid_interfaces_t_constructor( &
      photolysis_wavelength_grid_interfaces, size( photolysis_wavelength_grid_interfaces ) )

    deallocate( interfaces )

  end function create_extraterrestrial_flux_profile

  !> Sets TUV-x extraterrestrial flux midpoints
  !
  ! Extraterrestrial flux is read from data files and interpolated to the
  ! TUV-x wavelength grid. CAM extraterrestrial flux values are multiplied by the
  ! width of the wavelength bins to get the TUV-x units of photon cm-2 s-1
  !
  ! TUV-x only uses mid-point values for extraterrestrial flux
  subroutine set_extraterrestrial_flux_values(profile, extraterrestrial_flux, errmsg, errcode)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t
    use ccpp_kinds,          only: kind_phys
    use ccpp_tuvx_utils,     only: rebin

    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: extraterrestrial_flux(:) ! photons cm-2 s-1 nm-1
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    real(kind_phys), allocatable :: midpoints(:)
    type(error_t)                :: error

    if (.not. allocated(host_wavelength_grid_interfaces_)) then
      errmsg = "[MUSICA Error] Failed to allocate the host model wavelength grid interfaces"
      errcode = 1
      call cleanup_photolysis_wavelength_grid_interfaces()
      return
    end if

    if (.not. allocated(tuvx_wavelength_grid_interfaces_)) then
      errmsg = "[MUSICA Error] Failed to allocate the TUV-x wavelength grid interfaces"
      errcode = 1
      call cleanup_photolysis_wavelength_grid_interfaces()
      return
    end if

    allocate( midpoints( tuvx_wavelength_grid_interfaces_%size - 1 ))

    ! Regrid normalized flux to TUV-x wavelength grid
    ! This function fills the TUV-x midpoints.
    call rebin( host_wavelength_grid_interfaces_%size - 1,   &
                tuvx_wavelength_grid_interfaces_%size - 1,   &
                host_wavelength_grid_interfaces_%interfaces, &
                tuvx_wavelength_grid_interfaces_%interfaces, &
                extraterrestrial_flux, midpoints )

    ! Convert normalized flux to flux on TUV-x wavelength grid
    midpoints(:) = midpoints(:) * &
      ( tuvx_wavelength_grid_interfaces_%interfaces(2 : tuvx_wavelength_grid_interfaces_%size) - &
        tuvx_wavelength_grid_interfaces_%interfaces(1 : tuvx_wavelength_grid_interfaces_%size - 1) )

    call profile%set_midpoint_values( midpoints, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) then
      call cleanup_photolysis_wavelength_grid_interfaces()
      return
    end if

    deallocate( midpoints )

  end subroutine set_extraterrestrial_flux_values

end module musica_ccpp_tuvx_extraterrestrial_flux