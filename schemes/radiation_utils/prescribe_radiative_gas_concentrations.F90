!-------------------------------------------------------------------------------
! This module uses namelist variables to set prescribed concentrations
! for radiatively-active gases.

! The scheme owns the gases listed in the prescribed_ghg_list namelist
! variable: it registers a non-advected constituent for each and fills it at
! initialization with a spatially uniform value from the corresponding
! prescribed_*_vmr namelist variable.  A gas provided by another scheme
! (e.g. a prognostic chemistry constituent) should not be listed.

! Eventually this module should be replaced with a more comprehensive atmospheric
! composition/chemistry system, but is fine to use for now when running low-top,
! non-exoplanet CAM-SIMA configurations with minimal chemistry.
!-------------------------------------------------------------------------------
module prescribe_radiative_gas_concentrations

  implicit none

  private
  public :: prescribe_radiative_gas_concentrations_register
  public :: prescribe_radiative_gas_concentrations_init

  ! Gases this scheme can prescribe. Ozone is not included because its
  ! constituent is registered by the scheme providing it (e.g. prescribed_ozone).
  character(len=8), parameter :: provided_gases(6) = &
       [character(len=8) :: 'O2', 'CO2', 'N2O', 'CH4', 'CFC11', 'CFC12']

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
    character(len=*), intent(out) :: errmsg
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
       call dyn_consts(const_idx)%instantiate(     &
          std_name = trim(prescribed_ghg_list(const_idx)),   &
          long_name = trim(prescribed_ghg_list(const_idx)),  &
          units = 'kg kg-1',                            &
          vertical_dim = 'vertical_layer_dimension', &
          min_value = 0.0_kind_phys,                 &
          advected = .false.,                         &
          diag_name = trim(prescribed_ghg_list(const_idx)), &
          water_species = .false.,                    &
          mixing_ratio_type = 'dry',                 &
          errcode = errcode,                         &
          errmsg = errmsg)
       if (errcode /= 0) then
          return
       end if
    end do

  end subroutine prescribe_radiative_gas_concentrations_register

!> \section arg_table_prescribe_radiative_gas_concentrations_init Argument Table
!! \htmlinclude prescribe_radiative_gas_concentrations_init.html
!!
  subroutine prescribe_radiative_gas_concentrations_init(prescribed_ghg_list, &
                                               ch4_vmr,   co2_vmr, cfc11_vmr,   &
                                               cfc12_vmr, n2o_vmr, o2_vmr, const_array, &
                                               errmsg, errcode)

    ! Use statements
    use ccpp_kinds, only: kind_phys

    ! Input arguments
    character(len=*), intent(in) :: prescribed_ghg_list(:) ! (namelist) gases this scheme prescribes
    real(kind_phys), intent(in) :: ch4_vmr
    real(kind_phys), intent(in) :: co2_vmr
    real(kind_phys), intent(in) :: cfc11_vmr
    real(kind_phys), intent(in) :: cfc12_vmr
    real(kind_phys), intent(in) :: n2o_vmr
    real(kind_phys), intent(in) :: o2_vmr

    ! Input/output arguments
    real(kind_phys), intent(inout) :: const_array(:,:,:) ! Constituents array

    ! Output arguments
    character(len=*), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errmsg = ''
    errcode = 0

    !+++++++++++++++++++++++++++++++++++
    ! Convert number/mole fraction into
    ! mass mixing ratio w.r.t dry air
    !+++++++++++++++++++++++++++++++++++

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

  end subroutine prescribe_radiative_gas_concentrations_init

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
    character(len=*), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! Local variables
    integer         :: const_idx                         !Constituents object index
    real(kind_phys) :: dry_air_to_const_molar_mass_ratio !Ratio of dry air molar mass to constituent molar mass

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
