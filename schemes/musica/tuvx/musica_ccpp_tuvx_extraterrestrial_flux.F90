! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module musica_ccpp_tuvx_extraterrestrial_flux
  use ccpp_kinds, only: kind_phys

  implicit none

  private
  save

  public :: create_extraterrestrial_flux_profile, set_extraterrestrial_flux_values, &
            deallocate_photolysis_wavelength_grid_interfaces

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

  !> TODO(Jiwon)
  type, private :: wavelength_grid_interfaces_t
    real(kind_phys), allocatable :: wavelength_grid_interfaces(:) ! nm
    integer                      :: size = 0
  contains
    ! Deallocates memory associated with this musica species object
    procedure :: deallocate_data => wavelength_grid_interfaces_t_deallocate
  end type wavelength_grid_interfaces_t

  interface wavelength_grid_interfaces_t
    procedure wavelength_grid_interfaces_constructor
  end interface wavelength_grid_interfaces_t

  type(wavelength_grid_interfaces_t), protected, allocatable :: host_wavelength_grid_interfaces ! nm
  type(wavelength_grid_interfaces_t), protected, allocatable :: tuvx_wavelength_grid_interfaces ! nm


contains

  function wavelength_grid_interfaces_constructor(wavelength_grid_interfaces, size) &
      result( this )

    real(kind_phys), intent(in) :: wavelength_grid_interfaces(:) ! nm
    integer,         intent(in) :: size
    type(wavelength_grid_interfaces_t)  :: this

    allocate( this%wavelength_grid_interfaces( size ) )
    this%wavelength_grid_interfaces(:) = wavelength_grid_interfaces(:)
    this%size = size

  end function wavelength_grid_interfaces_constructor

  !> TODO(jiwon) - Deallocates memory associated with this musica species object
  subroutine wavelength_grid_interfaces_t_deallocate(this)
    class(wavelength_grid_interfaces_t), intent(inout) :: this

      if (allocated(this%wavelength_grid_interfaces)) then
       deallocate(this%wavelength_grid_interfaces)
      end if

  end subroutine wavelength_grid_interfaces_t_deallocate

  !> Deallocates photolysis wavelength grid interfaces
  subroutine deallocate_photolysis_wavelength_grid_interfaces()

    if (allocated( host_wavelength_grid_interfaces )) then
      call host_wavelength_grid_interfaces%deallocate_data()
      deallocate( host_wavelength_grid_interfaces )
    end if

    if (allocated( tuvx_wavelength_grid_interfaces )) then
      call tuvx_wavelength_grid_interfaces%deallocate_data()
      deallocate( tuvx_wavelength_grid_interfaces )
    end if

  end subroutine deallocate_photolysis_wavelength_grid_interfaces

  !> Creates a TUV-x extraterrestrial flux profile based on the TUV-x wavelength grid
  !  and initializes photolysis wavelength grid interfaces, which are configured differently
  !  for each CAM-SIMA and TUV-x.
  function create_extraterrestrial_flux_profile(wavelength_grid, &
      photolysis_wavelength_grid_interfaces, errmsg, errcode) result( profile )
    use musica_util,                      only: error_t
    use musica_ccpp_util,                 only: has_error_occurred
    use musica_ccpp_tuvx_wavelength_grid, only: m_to_nm
    use musica_tuvx_grid,                 only: grid_t
    use musica_tuvx_profile,              only: profile_t

    type(grid_t),     intent(inout) :: wavelength_grid
    real(kind_phys),  intent(in)    :: photolysis_wavelength_grid_interfaces(:) ! nm
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode
    type(profile_t),  pointer       :: profile

    ! local variables
    real(kind_phys), allocatable :: interfaces(:) ! nm
    type(error_t)                :: error

    profile => profile_t( extraterrestrial_flux_label, extraterrestrial_flux_unit, &
                          wavelength_grid, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    num_wavelength_grid_sections = wavelength_grid%number_of_sections( error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    allocate( interfaces( num_wavelength_grid_sections ) )

    call wavelength_grid%get_edges( interfaces, error )
    if ( has_error_occurred( error, errmsg, errcode ) ) return

    allocate(host_wavelength_grid_interfaces)
    host_wavelength_grid_interfaces = wavelength_grid_interfaces_constructor( &
      photolysis_wavelength_grid_interfaces, size( photolysis_wavelength_grid_interfaces ) )

    allocate(tuvx_wavelength_grid_interfaces)
    tuvx_wavelength_grid_interfaces = wavelength_grid_interfaces_constructor( &
      interfaces, num_wavelength_grid_sections )

    deallocate( interfaces )

  end function create_extraterrestrial_flux_profile

  !> Sets TUV-x extraterrestrial flux midpoints
  !
  ! Extraterrestrial flux is read from data files and interpolated to the
  ! TUV-x wavelength grid. CAM extraterrestrial flux values are multiplied by the
  ! width of the wavelength bins to get the TUV-x units of photon cm-2 s-1
  !
  ! TUV-x only uses mid-point values for extraterrestrial flux
  subroutine set_extraterrestrial_flux_values(profile, photolysis_wavelength_grid_interfaces, &
      extraterrestrial_flux, errmsg, errcode)
    use musica_ccpp_util,    only: has_error_occurred
    use musica_tuvx_profile, only: profile_t
    use musica_util,         only: error_t
    use ccpp_kinds,          only: kind_phys
    use ccpp_tuvx_utils,     only: rebin

    type(profile_t),  intent(inout) :: profile
    real(kind_phys),  intent(in)    :: photolysis_wavelength_grid_interfaces(:) ! nm
    real(kind_phys),  intent(in)    :: extraterrestrial_flux(:)                 ! photons cm-2 s-1 nm-1
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errcode

    ! local variables
    type(error_t)   :: error
    real(kind_phys) :: midpoints(num_wavelength_bins_)

    if (.not. allocated(host_wavelength_grid_interfaces)) then
      errmsg = "[MUSICA Error] Failed to allocate the host model wavelength grid interfaces"
      errcode = 1
      return
    end if

    if (.not. allocated(tuvx_wavelength_grid_interfaces)) then
      errmsg = "[MUSICA Error] Failed to allocate the TUV-x wavelength grid interfaces"
      errcode = 1
      return
    end if

    ! Regrid normalized flux to TUV-x wavelength grid
    ! This function fills the TUV-x midpoints.
    call rebin( host_wavelength_grid_interfaces%size - 1,                   & ! todo(jiwon) don't know why -1
                tuvx_wavelength_grid_interfaces%size,                       &
                host_wavelength_grid_interfaces%wavelength_grid_interfaces, &
                tuvx_wavelength_grid_interfaces%wavelength_grid_interfaces, &
                extraterrestrial_flux, midpoints )

    ! Convert normalized flux to flux on TUV-x wavelength grid
    midpoints = midpoints * &
      ( tuvx_wavelength_grid_interfaces(2 : tuvx_wavelength_grid_interfaces%size + 1) &
      - tuvx_wavelength_grid_interfaces(1 : tuvx_wavelength_grid_interfaces%size) )

    call profile%set_midpoint_values( midpoints, error)
    if ( has_error_occurred( error, errmsg, errcode ) ) then
      call deallocate_photolysis_wavelength_grid_interfaces()
      return
    end if

  end subroutine set_extraterrestrial_flux_values

end module musica_ccpp_tuvx_extraterrestrial_flux