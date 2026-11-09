!-------------------------------------------------------------------------------
! This module uses namelist variables to set prescribed concentrations
! for radiatively-active gases.
!
! This scheme "owns" (provides) the gases listed in the prescribed_ghg_list
! namelist variable, and thus is responsible for registering non-advected
! constituents for each.
!
! at init:
!  values for these GHGs are filled at init with a spatially uniform value
!  from the prescribed_*_vmr namelist variable.
! at timestep init:
!  if these GHGs have LBCs specified in flbc_list then they are
!  updated with time-interpolated uniform global means read by the LBC scheme.
!-------------------------------------------------------------------------------
module prescribe_radiative_gas_concentrations
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: prescribe_radiative_gas_concentrations_register
  public :: prescribe_radiative_gas_concentrations_init
  public :: prescribe_radiative_gas_concentrations_timestep_init

  ! Gases this scheme can prescribe:
  character(len=*), parameter :: provided_gases(6) = &
       [character(len=8) :: 'O2', 'CO2', 'N2O', 'CH4', 'CFC11', 'CFC12']

  ! Binding of prescribed GHGs to a LBC entry if available.
  type :: lbc_binding
    integer         :: handle = 0                 ! species handle in prescribe_lower_boundary_conditions
    integer         :: const_idx = -1             ! constituent the gas fills
    real(kind_phys) :: mmr_factor = -1._kind_phys ! conversion factor: vmr -> mass mixing ratio w.r.t. dry air
  end type lbc_binding

  integer                        :: n_bindings = 0
  type(lbc_binding), allocatable :: bindings(:)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> \section arg_table_prescribe_radiative_gas_concentrations_register Argument Table
!! \htmlinclude prescribe_radiative_gas_concentrations_register.html
!!
  subroutine prescribe_radiative_gas_concentrations_register(prescribed_ghg_list, dyn_consts, &
                                               errmsg, errcode)

    ! Use statements
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use ccpp_kinds,                only: kind_phys

    ! Input arguments
    character(len=*), intent(in) :: prescribed_ghg_list(:) ! (namelist) gases this scheme prescribes

    ! Output arguments
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: dyn_consts(:) ! Runtime constituent properties
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! Local variables
    character(len=256) :: alloc_errmsg
    integer            :: const_idx, num_consts, ierr

    errmsg = ''
    errcode = 0

    ! Count and validate the requested gases.
    num_consts = 0
    count_loop: do const_idx = 1, size(prescribed_ghg_list)

       if ( len_trim(prescribed_ghg_list(const_idx)) == 0 ) then
          exit count_loop
       end if

       if (.not. any(provided_gases == prescribed_ghg_list(const_idx))) then
          write(errmsg, *) 'prescribe_radiative_gas_concentrations_register: prescribed_ghg_list entry "', &
               trim(prescribed_ghg_list(const_idx)), '" is not a gas this scheme can prescribe (allowed: ', &
               provided_gases, ')'
          errcode = 1
          return
       end if

       num_consts = num_consts + 1

    end do count_loop

    ! Allocate the dynamic constituents array
    allocate(dyn_consts(num_consts), stat=ierr, errmsg=alloc_errmsg)
    if (ierr /= 0) then
       write(errmsg, *) 'prescribe_radiative_gas_concentrations_register: Unable to allocate dyn_consts - message: ', &
            trim(alloc_errmsg)
       errcode = 1
       return
    end if

    ! Register a non-advected constituent for each prescribed gas.
    do const_idx = 1, num_consts
       call dyn_consts(const_idx)%instantiate(               &
          std_name = trim(prescribed_ghg_list(const_idx)),   &
          long_name = trim(prescribed_ghg_list(const_idx)),  &
          units = 'kg kg-1',                                 &
          vertical_dim = 'vertical_layer_dimension',         &
          min_value = 0.0_kind_phys,                         &
          advected = .false.,                                &
          diag_name = trim(prescribed_ghg_list(const_idx)),  &
          water_species = .false.,                           &
          mixing_ratio_type = 'dry',                         &
          errcode = errcode,                                 &
          errmsg = errmsg)
       if (errcode /= 0) then
          return
       end if
    end do

  end subroutine prescribe_radiative_gas_concentrations_register

!> \section arg_table_prescribe_radiative_gas_concentrations_init Argument Table
!! \htmlinclude prescribe_radiative_gas_concentrations_init.html
  subroutine prescribe_radiative_gas_concentrations_init(prescribed_ghg_list, &
                                               ch4_vmr,   co2_vmr, cfc11_vmr,   &
                                               cfc12_vmr, n2o_vmr, o2_vmr, const_array, &
                                               errmsg, errcode)

    ! Use statements
    use ccpp_kinds, only: kind_phys

    ! Input arguments
    character(len=*), intent(in) :: prescribed_ghg_list(:) ! (namelist) gases this scheme prescribes
    real(kind_phys),  intent(in) :: ch4_vmr
    real(kind_phys),  intent(in) :: co2_vmr
    real(kind_phys),  intent(in) :: cfc11_vmr
    real(kind_phys),  intent(in) :: cfc12_vmr
    real(kind_phys),  intent(in) :: n2o_vmr
    real(kind_phys),  intent(in) :: o2_vmr

    ! Input/output arguments
    real(kind_phys), intent(inout) :: const_array(:,:,:) ! Constituents array

    ! Output arguments
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    ! Gases not in prescribed_ghg_list are left untouched: they belong to
    ! another provider (e.g. prognostic chemistry).
    call fill_prescribed_gas('CH4',   ch4_vmr,   prescribed_ghg_list, const_array, errmsg, errcode)
    if (errcode /= 0) then
      return
    end if

    call fill_prescribed_gas('CO2',   co2_vmr,   prescribed_ghg_list, const_array, errmsg, errcode)
    if (errcode /= 0) then
      return
    end if

    call fill_prescribed_gas('CFC11', cfc11_vmr, prescribed_ghg_list, const_array, errmsg, errcode)
    if (errcode /= 0) then
      return
    end if

    call fill_prescribed_gas('CFC12', cfc12_vmr, prescribed_ghg_list, const_array, errmsg, errcode)
    if (errcode /= 0) then
      return
    end if

    call fill_prescribed_gas('N2O',   n2o_vmr,   prescribed_ghg_list, const_array, errmsg, errcode)
    if (errcode /= 0) then
      return
    end if

    call fill_prescribed_gas('O2',    o2_vmr,    prescribed_ghg_list, const_array, errmsg, errcode)
    if (errcode /= 0) then
      return
    end if

    ! Bind to LBCs if any are provided:
    call bind_lbc_gases(prescribed_ghg_list, errmsg, errcode)
    if (errcode /= 0) then
      return
    end if

  end subroutine prescribe_radiative_gas_concentrations_init

!> \section arg_table_prescribe_radiative_gas_concentrations_timestep_init Argument Table
!! \htmlinclude arg_table_prescribe_radiative_gas_concentrations_timestep_init.html
  subroutine prescribe_radiative_gas_concentrations_timestep_init(ncol, const_array, &
                                               errmsg, errcode)

    use prescribe_lower_boundary_conditions, only: lbc_global_mean_vmr

    integer,            intent(in)    :: ncol
    real(kind_phys),    intent(inout) :: const_array(:,:,:)
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errcode

    ! Local variables
    integer         :: g
    real(kind_phys) :: vmr_gmean

    errmsg = ''
    errcode = 0

    ! Update prescribed GHGs with global mean VMRs from the LBC module if available:
    do g = 1, n_bindings
      call lbc_global_mean_vmr(bindings(g)%handle, ncol, vmr_gmean, errmsg, errcode)
      if (errcode /= 0) then
        return
      end if
      const_array(:, :, bindings(g)%const_idx) = vmr_gmean*bindings(g)%mmr_factor
    end do

  end subroutine prescribe_radiative_gas_concentrations_timestep_init

  ! Bind each prescribed GHG to an entry in the LBC list if available.
  ! note on CFC11:
  !  it accepts the CMIP equivalent entry CFC11eq (all minor halogens expressed as CFC11)
  !  if CFC11 itself is not listed (simple chemistry configurations use CFC11eq instead of CFC11)
  subroutine bind_lbc_gases(prescribed_ghg_list, errmsg, errcode)

    use ccpp_constituent_prop_mod, only: int_unassigned
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use radiation_utils,           only: get_molar_mass_ratio
    use prescribe_lower_boundary_conditions, only: lbc_activate_species

    ! Input arguments
    character(len=*), intent(in)  :: prescribed_ghg_list(:)

    ! Output arguments
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errcode

    ! Local variables
    character(len=256) :: alloc_errmsg
    integer            :: g, handle, const_idx
    real(kind_phys)    :: dry_air_to_const_molar_mass_ratio

    errmsg = ''
    errcode = 0
    n_bindings = 0

    if (allocated(bindings)) then
      deallocate(bindings)
    end if
    allocate (bindings(size(prescribed_ghg_list)), stat=errcode, errmsg=alloc_errmsg)
    if (errcode /= 0) then
      write(errmsg, *) 'prescribe_radiative_gas_concentrations (bind_lbc_gases): ', &
           'Unable to allocate bindings - message: ', trim(alloc_errmsg)
      return
    end if

    gas_loop: do g = 1, size(prescribed_ghg_list)

      if (len_trim(prescribed_ghg_list(g)) == 0) then
        exit gas_loop
      end if

      call lbc_activate_species(trim(prescribed_ghg_list(g)), handle)
      if (handle == 0 .and. trim(prescribed_ghg_list(g)) == 'CFC11') then
        ! allow fallback to CFC11eq for simple chemistry configurations to match CAM.
        call lbc_activate_species('CFC11eq', handle)
      end if
      if (handle == 0) then
        cycle gas_loop
      end if

      ! Find the constituent registered by this scheme's register phase:
      call ccpp_constituent_index(trim(prescribed_ghg_list(g)), const_idx, errcode, errmsg)
      if (errcode /= 0) then
        return
      end if
      if (const_idx == int_unassigned) then
        write(errmsg, *) 'prescribe_radiative_gas_concentrations (bind_lbc_gases): ', &
             'no constituent registered for prescribed gas "', trim(prescribed_ghg_list(g)), '"'
        errcode = 1
        return
      end if

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio(trim(prescribed_ghg_list(g)), dry_air_to_const_molar_mass_ratio, &
           errmsg, errcode)
      if (errcode /= 0) then
        return
      end if

      n_bindings = n_bindings + 1
      bindings(n_bindings)%handle = handle
      bindings(n_bindings)%const_idx = const_idx
      bindings(n_bindings)%mmr_factor = 1._kind_phys/dry_air_to_const_molar_mass_ratio

    end do gas_loop

  end subroutine bind_lbc_gases

  ! If gas_name is in prescribed_ghg_list, convert its namelist-provided
  ! number/mole fraction into a mass mixing ratio w.r.t. dry air and fill
  ! its constituent (uniform over all columns and levels).
  subroutine fill_prescribed_gas(gas_name, gas_vmr, prescribed_ghg_list, const_array, &
                                 errmsg, errcode)

    ! Use statements
    use ccpp_constituent_prop_mod, only: int_unassigned
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_kinds,                only: kind_phys

    ! Use statement from RRTMGP,
    ! which should hopefully be replaced once
    ! molar masses for these species are included
    ! in the constituents properties themselves:
    use radiation_utils, only: get_molar_mass_ratio

    ! Input arguments
    character(len=*), intent(in) :: gas_name
    real(kind_phys),  intent(in) :: gas_vmr
    character(len=*), intent(in) :: prescribed_ghg_list(:)

    ! Input/output arguments
    real(kind_phys), intent(inout) :: const_array(:,:,:) ! Constituents array

    ! Output arguments
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! Local variables
    integer         :: const_idx                         ! Constituents object index
    real(kind_phys) :: dry_air_to_const_molar_mass_ratio ! Ratio of dry air molar mass to constituent molar mass

    errmsg = ''
    errcode = 0

    if (.not. any(prescribed_ghg_list == gas_name)) then
      return
    end if

    ! Find the constituent registered by this scheme's register phase:
    call ccpp_constituent_index(gas_name, const_idx, errcode, errmsg)
    if (errcode /= 0) then
      return
    end if
    if (const_idx == int_unassigned) then
      write(errmsg, *) 'prescribe_radiative_gas_concentrations_init: no constituent registered for prescribed gas "', &
           trim(gas_name), '"'
      errcode = 1
      return
    end if

    !+++++++++++++++++++++++++++++++++++
    ! Convert number/mole fraction into
    ! mass mixing ratio w.r.t dry air
    !+++++++++++++++++++++++++++++++++++

    ! Get ratio of molar mass of dry air / constituent molar mass
    call get_molar_mass_ratio(gas_name, dry_air_to_const_molar_mass_ratio, errmsg, errcode)
    if (errcode /= 0) then
      return
    end if

    ! Convert namelist-provided number/mole fraction to
    ! mass mixing ratio w.r.t. dry air, and set constituents
    ! array to new converted value:
    const_array(:,:,const_idx) = gas_vmr/dry_air_to_const_molar_mass_ratio

  end subroutine fill_prescribed_gas

end module prescribe_radiative_gas_concentrations
