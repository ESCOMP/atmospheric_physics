! Clamp number concentration constituent tendencies
! so that after tendency application, values remain within [qmin, qmax].
!
! this subroutine updates const_q directly because const_q >> qmin
! using tendencies in this subroutine would result in catastrophic cancellation
!
! This scheme must be placed immediately after apply_constituent_tendencies
! in the SDF to replicate the same behavior as in physics_types.F90 in CAM.
module clamp_number_concentrations

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: clamp_number_concentrations_init
  public :: clamp_number_concentrations_run

  ! Number of number-concentration species to clamp
  integer, parameter :: num_species = 4

  ! Constituent indices (-1 = not present; looked up in _init)
  integer :: ix_numliq  = -1
  integer :: ix_numrai  = -1
  integer :: ix_numice  = -1
  integer :: ix_numsno  = -1

  ! Clamp bounds
  real(kind_phys), parameter :: qmin = 1.0e-12_kind_phys ! minimum number concentration [kg-1]
  real(kind_phys), parameter :: qmax = 1.0e10_kind_phys  ! maximum number concentration [kg-1]

  ! Standard names for each species (order matches ix_ variables above)
  character(len=*), parameter :: std_names(num_species) = &
    [character(len=128) :: 'mass_number_concentration_of_cloud_liquid_water_droplets_in_moist_air_and_condensed_water', &
                           'mass_number_concentration_of_rain_drops_in_moist_air_and_condensed_water', &
                           'mass_number_concentration_of_cloud_ice_water_crystals_in_moist_air_and_condensed_water', &
                           'mass_number_concentration_of_snow_crystals_in_moist_air_and_condensed_water']

contains

!> \section arg_table_clamp_number_concentrations_init Argument Table
!! \htmlinclude clamp_number_concentrations_init.html
  subroutine clamp_number_concentrations_init(errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    integer :: ix_species(num_species)
    integer :: n

    errmsg = ''
    errflg = 0

    ix_species(:) = -1
    ! Look up constituent indices. A missing species is skipped
    do n = 1, num_species
      call ccpp_constituent_index(trim(std_names(n)), ix_species(n), errmsg=errmsg, errcode=errflg)
      if(errflg /= 0) return
      if (ix_species(n) < 0) then
        ! Constituent not found — mark as inactive
        ix_species(n) = -1
      end if
    end do

    ! Store in module variables
    ix_numliq = ix_species(1)
    ix_numrai = ix_species(2)
    ix_numice = ix_species(3)
    ix_numsno = ix_species(4)

  end subroutine clamp_number_concentrations_init

!> \section arg_table_clamp_number_concentrations_run Argument Table
!! \htmlinclude clamp_number_concentrations_run.html
  subroutine clamp_number_concentrations_run( &
    ncol, pver, &
    const_q, const_tend, &
    errmsg, errflg)

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(inout) :: const_q(:, :, :)    ! constituent mixing ratios
    real(kind_phys),    intent(inout) :: const_tend(:, :, :) ! constituent tendencies
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables
    integer :: ix_species(num_species)
    integer :: n, i, k, ix
    real(kind_phys) :: q_projected  ! projected value after tendency application [kg-1]

    errmsg = ''
    errflg = 0

    ix_species = [ix_numliq, ix_numrai, ix_numice, ix_numsno]

    do n = 1, num_species
      ix = ix_species(n)
      if (ix <= 0) cycle

      do k = 1, pver
        do i = 1, ncol
          if (const_q(i, k, ix) < qmin) then
            const_q(i, k, ix)    = qmin
            const_tend(i, k, ix) = 0.0_kind_phys
          else if (const_q(i, k, ix) > qmax) then
            const_q(i, k, ix)    = qmax
            const_tend(i, k, ix) = 0.0_kind_phys
          end if
        end do
      end do
    end do
  end subroutine clamp_number_concentrations_run

end module clamp_number_concentrations
