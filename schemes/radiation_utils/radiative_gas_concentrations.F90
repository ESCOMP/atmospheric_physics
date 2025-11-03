!-------------------------------------------------------------------------------
! This module uses namelist variable to set prescribed concentrations
! for radiatively-active gases.

! Eventually this module should be replaced with a more comprehensive atmospheric
! composition/chemistry system, but is fine to use for now when running low-top,
! non-exoplanet CAM-SIMA configurations with minimal chemistry.
!-------------------------------------------------------------------------------
module radiative_gas_concentrations

  implicit none

  private
  public :: radiative_gas_concentrations_init


!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> \section arg_table_radiative_gas_concentrations_init Argument Table
!! \htmlinclude radiative_gas_concentrations_init.html
!!
  subroutine radiative_gas_concentrations_init(ch4_vmr, co2_vmr, cfc11_vmr, &
                                               cfc12_vmr, n2o_vmr,          &
                                               const_array, errmsg, errcode)

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
    real(kind_phys), intent(in) :: ch4_vmr
    real(kind_phys), intent(in) :: co2_vmr
    real(kind_phys), intent(in) :: cfc11_vmr
    real(kind_phys), intent(in) :: cfc12_vmr
    real(kind_phys), intent(in) :: n2o_vmr

    ! Input/output arguments
    real(kind_phys), intent(inout) :: const_array(:,:,:) ! Constituents array

    ! Output arguments:
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! Local variables
    integer :: cnst_idx !Constituents object index
    real(kind_phys) :: dry_air_to_const_molar_mass_ratio !Ratio of dry air molar mass to constituent molar mass
    real(kind_phys) :: gas_mmr                           !Mass mixing ratio wrt dry air for specified gas

    ! Convert number/mole fraction into mass mixing ratio wrt dry air:

    !----
    ! CH4:
    !----

    ! Check if CH4 is present in constituents object:
    call ccpp_constituent_index('CH4', const_idx, errflg, errmsg)
    if (errflg /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('CH4', dry_air_to_const_molar_mass_ratio, errmsg, errflg)
      if (errflg /= 0) then
        return
      end if

      ! If so, then convert namelist-provided number/mole fraction to
      !mass mixing ratio w.r.t. dry air:
      gas_mmr = chv_vmr*dry_air_to_const_molar_mass_ratio

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new coonverted value:
      const_array(:,: cnst_idx) = chv_vmr*dry_air_to_const_molar_mass_ratio

    end if

    !----
    ! CO2:
    !----

    ! Check if CO2 is present in constituents object:
    call ccpp_constituent_index('CO2', const_idx, errflg, errmsg)
    if (errflg /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('CO2', dry_air_to_const_molar_mass_ratio, errmsg, errflg)
      if (errflg /= 0) then
        return
      end if

      ! If so, then convert namelist-provided number/mole fraction to
      !mass mixing ratio w.r.t. dry air:
      gas_mmr = chv_vmr*dry_air_to_const_molar_mass_ratio

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new coonverted value:
      const_array(:,: cnst_idx) = chv_vmr*dry_air_to_const_molar_mass_ratio

    end if

    !------
    ! CFC11:
    !------

    ! Check if CFC-11 is present in constituents object:
    call ccpp_constituent_index('CFC11', const_idx, errflg, errmsg)
    if (errflg /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('CFC11', dry_air_to_const_molar_mass_ratio, errmsg, errflg)
      if (errflg /= 0) then
        return
      end if

      ! If so, then convert namelist-provided number/mole fraction to
      !mass mixing ratio w.r.t. dry air:
      gas_mmr = chv_vmr*dry_air_to_const_molar_mass_ratio

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new coonverted value:
      const_array(:,: cnst_idx) = chv_vmr*dry_air_to_const_molar_mass_ratio

    end if

    !------
    ! CFC12:
    !------

    ! Check if CFC-12 is present in constituents object:
    call ccpp_constituent_index('CFC12', const_idx, errflg, errmsg)
    if (errflg /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('CFC12', dry_air_to_const_molar_mass_ratio, errmsg, errflg)
      if (errflg /= 0) then
        return
      end if

      ! If so, then convert namelist-provided number/mole fraction to
      !mass mixing ratio w.r.t. dry air:
      gas_mmr = chv_vmr*dry_air_to_const_molar_mass_ratio

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new coonverted value:
      const_array(:,: cnst_idx) = chv_vmr*dry_air_to_const_molar_mass_ratio

    end if

    !----
    ! N2O:
    !----

    ! Check if N2O is present in constituents object:
    call ccpp_constituent_index('N2O', const_idx, errflg, errmsg)
    if (errflg /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('N2O', dry_air_to_const_molar_mass_ratio, errmsg, errflg)
      if (errflg /= 0) then
        return
      end if

      ! If so, then convert namelist-provided number/mole fraction to
      !mass mixing ratio w.r.t. dry air:
      gas_mmr = chv_vmr*dry_air_to_const_molar_mass_ratio

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new coonverted value:
      const_array(:,: cnst_idx) = chv_vmr*dry_air_to_const_molar_mass_ratio

    end if

    ! Set error variables
    errmsg = ''
    errflg = 0

  end subroutine radiative_gas_concentrations_init

end module radiative_gas_concentrations
