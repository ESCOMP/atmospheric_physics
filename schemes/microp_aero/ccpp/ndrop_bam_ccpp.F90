! CCPP wrapper for BAM droplet activation (ndrop_bam).
! Looks up BAM aerosol objects from aerosol_instances and delegates to the
! portable ndrop_bam driver.
module ndrop_bam_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: ndrop_bam_ccpp_init
  public :: ndrop_bam_ccpp_run

  ! smallest mixing ratio considered in microphysics
  real(kind_phys), parameter :: qsmall = 1.e-18_kind_phys

contains

!> \section arg_table_ndrop_bam_ccpp_init Argument Table
!! \htmlinclude ndrop_bam_ccpp_init.html
  subroutine ndrop_bam_ccpp_init(amIRoot, iulog, mwh2o, r_universal, tmelt, rhoh2o, &
    errmsg, errflg)

    use ndrop_bam, only: ndrop_bam_init

    logical,         intent(in) :: amIRoot
    integer,         intent(in) :: iulog
    real(kind_phys), intent(in) :: mwh2o
    real(kind_phys), intent(in) :: r_universal
    real(kind_phys), intent(in) :: tmelt
    real(kind_phys), intent(in) :: rhoh2o

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call ndrop_bam_init(amIRoot, iulog, mwh2o, r_universal, tmelt, rhoh2o)

  end subroutine ndrop_bam_ccpp_init

!> \section arg_table_ndrop_bam_ccpp_run Argument Table
!! \htmlinclude ndrop_bam_ccpp_run.html
  subroutine ndrop_bam_ccpp_run( &
    ncol, pver, top_lev, &
    gravit, rair, tmelt, cpair, rh2o, rhoh2o, latvap, &
    rho, tair, wsub, qcld, ast, numliq, deltatin, &
    npccn, nacon, ccn, naer2_diag, &
    errmsg, errflg)

    use ndrop_bam,              only: ndrop_bam_calc, naer_all, psat
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_state, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod,      only: aerosol_state

    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    integer,         intent(in)  :: top_lev
    real(kind_phys), intent(in)  :: gravit
    real(kind_phys), intent(in)  :: rair
    real(kind_phys), intent(in)  :: tmelt
    real(kind_phys), intent(in)  :: cpair
    real(kind_phys), intent(in)  :: rh2o
    real(kind_phys), intent(in)  :: rhoh2o
    real(kind_phys), intent(in)  :: latvap
    real(kind_phys), intent(in)  :: rho(:, :)
    real(kind_phys), intent(in)  :: t(:, :)
    real(kind_phys), intent(in)  :: wsub(:, :)
    real(kind_phys), intent(in)  :: qcld(:, :)
    real(kind_phys), intent(in)  :: ast(:, :)
    real(kind_phys), intent(in)  :: numliq(:, :)
    real(kind_phys), intent(in)  :: deltatin
    real(kind_phys), intent(out) :: npccn(:, :)
    real(kind_phys), intent(out) :: nacon(:, :, :)
    real(kind_phys), intent(out) :: ccn(:, :, :)
    real(kind_phys), intent(out) :: naer2_diag(:, :, :)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: iaermod
    class(aerosol_properties), pointer :: aero_props_bam
    class(aerosol_state),      pointer :: aero_state_bam

    errmsg = ''
    errflg = 0

    ! Find BAM properties and state from aerosol instances
    aero_props_bam => null()
    aero_state_bam => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props_bam => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props_bam)) then
        if (aero_props_bam%model_is('BAM')) then
          aero_state_bam => aerosol_instances_get_state(iaermod, list_idx=0)
          exit
        end if
      end if
      aero_props_bam => null()
    end do

    ! No BAM model found but microp_aero is active (i.e., aquaplanet).
    ! Zero outputs and return.
    if (.not. associated(aero_props_bam) .or. &
        .not. associated(aero_state_bam)) then
      npccn(:, :)      = 0._kind_phys
      nacon(:, :, :)   = 0._kind_phys
      ccn(:, :, :)     = 0._kind_phys
      naer2_diag(:, :, :) = 0._kind_phys
      return
    end if

    ! TODO for Claude: newly added to pass rho to ndrop_bam_calc
    ! need to add pmid to args, rho to local vars (ncol,pver)
    do k = top_lev, pver
      do i = 1, ncol
        rho(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
    end do

    call ndrop_bam_calc( &
      aero_state = aero_state_bam, &
      aero_props = aero_props_bam, &
      ncol       = ncol,           &
      pver       = pver,           &
      top_lev    = top_lev,        &
      gravit     = gravit,         &
      rair       = rair,           &
      tmelt      = tmelt,          &
      cpair      = cpair,          &
      rh2o       = rh2o,           &
      rhoh2o     = rhoh2o,         &
      latvap     = latvap,         &
      rho        = rho,            &
      tair       = t,              &
      wsub       = wsub,           &
      qcld       = qcld,           &
      qsmall_in  = qsmall,         &
      ast        = ast,            &
      numliq     = numliq,         &
      deltatin   = deltatin,       & ! below output:
      npccn      = npccn,          &
      nacon      = nacon,          &
      ccn        = ccn,            &
      naer2_diag = naer2_diag,     &
      errmsg     = errmsg,         &
      errflg     = errflg)

  end subroutine ndrop_bam_ccpp_run

end module ndrop_bam_ccpp
