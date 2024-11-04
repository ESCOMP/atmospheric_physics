module convect_shallow_zm_conv_evap_prep
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: convect_shallow_zm_conv_evap_prep_run

contains

!> \section arg_table_convect_shallow_zm_conv_evap_prep_run Argument Table
!! \htmlinclude convect_shallow_zm_conv_evap_prep_run.html
  subroutine convect_shallow_zm_conv_evap_prep_run( &
      ncol, pver, &
      cldice_sh, cldliq_sh, &
      prec_sh, prec_dp_in)     ! renaming prec_sh -> prec_deep for inout in zm_conv_evap

    ! Input arguments
    integer,         intent(in)    :: ncol               ! number of atmospheric columns
    integer,         intent(in)    :: pver               ! number of vertical levels
    real(kind_phys), intent(in)    :: prec_sh(:,:)       ! rename cmfdqr -> prec_dp_in [kg kg-1 s-1]

    ! Output arguments
    real(kind_phys), intent(out)   :: cldice_sh(:,:)     ! shallow convection ice cloud mixing ratio [kg kg-1]
    real(kind_phys), intent(out)   :: cldliq_sh(:,:)     ! shallow convection liquid cloud mixing ratio [kg kg-1]
    real(kind_phys), intent(out)   :: prec_dp_in(:,:)    ! input to deep zm_conv_evap: -> prdprec [kg kg-1 s-1]


    ! eventually set cldice_sh and cldliq_sh to zero
    ! as " clouds have no water... :)" from convect_shallow.F90
    cldice_sh(:ncol,:) = 0._kind_phys
    cldliq_sh(:ncol,:) = 0._kind_phys

    ! eventually rename and prepare physical quantities for zm_conv_evap
    ! cmfdqr/rprdsh (shallow convective rainout) -> input to deep zm_conv_evap
    !
    ! devnote hplin: I think this may not be actually necessary if zm_conv_evap can rename
    ! tendency_of_precipitation_wrt_moist_air_and_condensed_water_due_to_deep_convection_excluding_subcloud_evaporation
    ! to use tendency_of_precipitation_wrt_moist_air_and_condensed_water_due_to_shallow_convection_excluding_subcloud_evaporation instead. to check where other places this is used. if used in deep conv as well then
    ! maybe a more generic name can be adopted so we do not have to do this interstitial for renaming
    prec_dp_in(:ncol,:) = prec_sh(:ncol,:)

  end subroutine convect_shallow_zm_conv_evap_prep_run


end module convect_shallow_zm_conv_evap_prep
