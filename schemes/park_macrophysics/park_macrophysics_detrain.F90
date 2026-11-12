! Park cloud macrophysics: detrainment of convective condensate
! Extracted from cldwat2m_macro.F90 and macrop_driver.F90.
! Park, S., Bretherton, C. S., and Rasch, P. J.: Integrating Cloud Processes
! in the Commmunity Atmosphere Model, Version 5
! https://doi.org/10.1175/JCLI-D-14-00087.1
! Original author: Sungsu Park
module park_macrophysics_detrain
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: park_macrophysics_detrain_run

  ! 'cu_det_st' : If .true. (.false.), detrain cumulus liquid condensate into the pre-existing liquid stratus
  !               (environment) without (with) macrophysical evaporation. If there is no pre-esisting stratus,
  !               evaporate cumulus liquid condensate. This option only influences the treatment of cumulus
  !               liquid condensate, not cumulus ice condensate.
  logical,         parameter :: cu_det_st = .false.

contains

  ! Detrainment of convective condensate into the environment or stratiform cloud
!> \section arg_table_park_macrophysics_detrain_run Argument Table
!! \htmlinclude park_macrophysics_detrain_run.html
  subroutine park_macrophysics_detrain_run( &
    ncol, pver, top_lev, &
    latice, cpair, gravit, &
    do_detrain, &
    t, pdel, dlf, dlf2, &
    tend_cldliq, tend_cldice, tend_numliq, tend_numice, tend_s, &
    det_s, det_ice, &
    dpdlfliq, dpdlfice, shdlfliq, shdlfice, dpdlft, shdlft, &
    dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni, &
    scheme_name, errmsg, errflg)

    ! Input arguments
    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    integer,         intent(in)  :: top_lev
    real(kind_phys), intent(in)  :: latice                ! Latent heat of fusion [J kg-1]
    real(kind_phys), intent(in)  :: cpair                 ! Specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys), intent(in)  :: gravit                ! Gravitational acceleration [m s-2]
    logical,         intent(in)  :: do_detrain            ! If .true., compute detrainment of convective condensate [flag]
    real(kind_phys), intent(in)  :: t(:, :)               ! Temperature [K]
    real(kind_phys), intent(in)  :: pdel(:, :)            ! Layer pressure thickness [Pa]
    real(kind_phys), intent(in)  :: dlf(:, :)             ! Detraining cloud water from deep+shallow convection [kg kg-1 s-1]
    real(kind_phys), intent(in)  :: dlf2(:, :)            ! Detraining cloud water from shallow convection only [kg kg-1 s-1]

    ! Output arguments
    real(kind_phys), intent(out) :: tend_cldliq(:, :)     ! Tendency of grid-mean cloud liquid from detrainment [kg kg-1 s-1]
    real(kind_phys), intent(out) :: tend_cldice(:, :)     ! Tendency of grid-mean cloud ice from detrainment [kg kg-1 s-1]
    real(kind_phys), intent(out) :: tend_numliq(:, :)     ! Tendency of cloud liquid droplet number from detrainment [# kg-1 s-1]
    real(kind_phys), intent(out) :: tend_numice(:, :)     ! Tendency of cloud ice crystal number from detrainment [# kg-1 s-1]
    real(kind_phys), intent(out) :: tend_s(:, :)          ! Energy tendency from ice detrainment (latent heat of fusion) [J kg-1 s-1]
    real(kind_phys), intent(out) :: det_s(:)              ! Column-integrated energy from ice detrainment [W m-2]
    real(kind_phys), intent(out) :: det_ice(:)            ! Column-integrated detrained ice, liquid water equivalent [m s-1]
    real(kind_phys), intent(out) :: dpdlfliq(:, :)        ! Deep convective detraining cloud liquid [kg kg-1 s-1]
    real(kind_phys), intent(out) :: dpdlfice(:, :)        ! Deep convective detraining cloud ice [kg kg-1 s-1]
    real(kind_phys), intent(out) :: shdlfliq(:, :)        ! Shallow convective detraining cloud liquid [kg kg-1 s-1]
    real(kind_phys), intent(out) :: shdlfice(:, :)        ! Shallow convective detraining cloud ice [kg kg-1 s-1]
    real(kind_phys), intent(out) :: dpdlft(:, :)          ! Deep convective detraining temperature tendency [K s-1]
    real(kind_phys), intent(out) :: shdlft(:, :)          ! Shallow convective detraining temperature tendency [K s-1]
    real(kind_phys), intent(out) :: dlf_T(:, :)           ! Targeted stratus detrainment: temperature forcing [K s-1]
    real(kind_phys), intent(out) :: dlf_qv(:, :)          ! Targeted stratus detrainment: water vapor forcing [kg kg-1 s-1]
    real(kind_phys), intent(out) :: dlf_ql(:, :)          ! Targeted stratus detrainment: cloud liquid forcing [kg kg-1 s-1]
    real(kind_phys), intent(out) :: dlf_qi(:, :)          ! Targeted stratus detrainment: cloud ice forcing [kg kg-1 s-1]
    real(kind_phys), intent(out) :: dlf_nl(:, :)          ! Targeted stratus detrainment: liquid number forcing [# kg-1 s-1]
    real(kind_phys), intent(out) :: dlf_ni(:, :)          ! Targeted stratus detrainment: ice number forcing [# kg-1 s-1]
    character(len=64),  intent(out) :: scheme_name
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer         :: i, k
    real(kind_phys) :: dum1

    errmsg = ''
    errflg = 0
    scheme_name = 'park_macrophysics_detrain'

    ! Initialize all outputs to zero
    tend_cldliq(:, :) = 0._kind_phys
    tend_cldice(:, :) = 0._kind_phys
    tend_numliq(:, :) = 0._kind_phys
    tend_numice(:, :) = 0._kind_phys
    tend_s(:, :) = 0._kind_phys

    det_s(:) = 0._kind_phys
    det_ice(:) = 0._kind_phys

    dpdlfliq = 0._kind_phys
    dpdlfice = 0._kind_phys
    shdlfliq = 0._kind_phys
    shdlfice = 0._kind_phys
    dpdlft = 0._kind_phys
    shdlft = 0._kind_phys

    dlf_T(:, :) = 0._kind_phys
    dlf_qv(:, :) = 0._kind_phys
    dlf_ql(:, :) = 0._kind_phys
    dlf_qi(:, :) = 0._kind_phys
    dlf_nl(:, :) = 0._kind_phys
    dlf_ni(:, :) = 0._kind_phys

    do k = top_lev, pver
    do i = 1, ncol
      if (t(i, k) > 268.15_kind_phys) then
        dum1 = 0.0_kind_phys
      elseif (t(i, k) < 238.15_kind_phys) then
        dum1 = 1.0_kind_phys
      else
        dum1 = (268.15_kind_phys - t(i, k))/30._kind_phys
      end if

      ! If detrainment was done elsewhere, still update the variables used for output
      ! assuming that the temperature split between liquid and ice is the same as assumed
      ! here.
      if (do_detrain) then
        tend_cldliq(i, k) = dlf(i, k)*(1._kind_phys - dum1)
        tend_cldice(i, k) = dlf(i, k)*dum1
        ! dum2              = dlf(i,k) * ( 1._kind_phys - dum1 )
        tend_numliq(i, k) = 3._kind_phys*(max(0._kind_phys, (dlf(i, k) - dlf2(i, k)))*(1._kind_phys - dum1))/ &
                            (4._kind_phys*3.14_kind_phys*8.e-6_kind_phys**3*997._kind_phys) + & ! Deep    Convection
                            3._kind_phys*(dlf2(i, k)*(1._kind_phys - dum1))/ &
                            (4._kind_phys*3.14_kind_phys*10.e-6_kind_phys**3*997._kind_phys)     ! Shallow Convection
        ! dum2              = dlf(i,k) * dum1
        tend_numice(i, k) = 3._kind_phys*(max(0._kind_phys, (dlf(i, k) - dlf2(i, k)))*dum1)/ &
                            (4._kind_phys*3.14_kind_phys*25.e-6_kind_phys**3*500._kind_phys) + & ! Deep    Convection
                            3._kind_phys*(dlf2(i, k)*dum1)/ &
                            (4._kind_phys*3.14_kind_phys*50.e-6_kind_phys**3*500._kind_phys)     ! Shallow Convection
        tend_s(i, k) = dlf(i, k)*dum1*latice
      else
        tend_cldliq(i, k) = 0._kind_phys
        tend_cldice(i, k) = 0._kind_phys
        tend_numliq(i, k) = 0._kind_phys
        tend_numice(i, k) = 0._kind_phys
        tend_s(i, k) = 0._kind_phys
      end if

      ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
      !   track of the integrals of ice and static energy that is effected from conversion to ice
      !   so that the energy checker doesn't complain.
      det_s(i) = det_s(i) + tend_s(i, k)*pdel(i, k)/gravit
      det_ice(i) = det_ice(i) - tend_cldice(i, k)*pdel(i, k)/gravit

      ! Targetted detrainment of convective liquid water either directly into the
      ! existing liquid stratus or into the environment.
      if (cu_det_st) then
        dlf_T(i, k) = tend_s(i, k)/cpair
        dlf_qv(i, k) = 0._kind_phys
        dlf_ql(i, k) = tend_cldliq(i, k)
        dlf_qi(i, k) = tend_cldice(i, k)
        dlf_nl(i, k) = tend_numliq(i, k)
        dlf_ni(i, k) = tend_numice(i, k)
        tend_cldliq(i, k) = 0._kind_phys
        tend_cldice(i, k) = 0._kind_phys
        tend_numliq(i, k) = 0._kind_phys
        tend_numice(i, k) = 0._kind_phys
        tend_s(i, k) = 0._kind_phys
        dpdlfliq(i, k) = 0._kind_phys
        dpdlfice(i, k) = 0._kind_phys
        shdlfliq(i, k) = 0._kind_phys
        shdlfice(i, k) = 0._kind_phys
        dpdlft(i, k) = 0._kind_phys
        shdlft(i, k) = 0._kind_phys
      else
        dpdlfliq(i, k) = (dlf(i, k) - dlf2(i, k))*(1._kind_phys - dum1)
        dpdlfice(i, k) = (dlf(i, k) - dlf2(i, k))*(dum1)
        dpdlft(i, k) = (dlf(i, k) - dlf2(i, k))*dum1*latice/cpair

        shdlfliq(i, k) = dlf2(i, k)*(1._kind_phys - dum1)
        shdlfice(i, k) = dlf2(i, k)*(dum1)
        shdlft(i, k) = dlf2(i, k)*dum1*latice/cpair
      end if
    end do
    end do

    det_ice(:ncol) = det_ice(:ncol)/1000._kind_phys  ! divide by density of water

  end subroutine park_macrophysics_detrain_run

end module park_macrophysics_detrain
