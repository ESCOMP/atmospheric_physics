! Compute cloud fractions.
! This was used in combination with the two-moment MG microphysics in CAM,
! hence the "two moment" distinction. Formerly called "cldfrc2m".
!
! Provides liquid stratus (PDF and RHU methods) and ice stratus cloud fraction.
! Original authors: Sungsu Park, Andrew Gettelman
module compute_cloud_fraction_two_moment
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! CCPP-compliant subroutines:
  public :: compute_cloud_fraction_two_moment_init

  ! Underlying subroutines that are not CCPP schemes
  ! but are used by other schemes and internally use the parameters
  ! read from the above init phase.
  public :: astG_PDF_single, astG_PDF
  public :: astG_RHU_single, astG_RHU
  public :: aist_single, aist_vector

  !REMOVECAM: this public parameter is used by CAM code and this public declaration can be removed
  ! after CAM is retired and all code that uses this parameter is migrated to receive this via
  ! the CCPP framework.
  public :: CAMstfrac
  !REMOVECAM_END

  ! Hardcoded parameters
  logical, parameter :: CAMstfrac = .false.  ! If .true. (.false.),
  ! use Slingo (triangular PDF-based) liquid stratus fraction
  logical, parameter :: freeze_dry = .false. ! If .true., use 'freeze dry' in liquid stratus fraction formula

  ! Module variables set by init (used internally by subroutines)
  real(kind_phys) :: premit            ! Top pressure bound for mid-level clouds [Pa]
  real(kind_phys) :: premib            ! Bottom pressure bound for mid-level clouds [Pa]

  ! Ice cloud closure option
  ! 1 = Wang and Sassen
  ! 2 = Schiller (iciwc)
  ! 3 = Wood and Field
  ! 4 = Wilson (based on Smith)
  ! 5 = modified Slingo (ssat & empty cloud)
  integer         :: iceopt

  ! Ice cloud fraction fit parameters (used by aist_single and aist_vector)
  ! Wang and Sassen IWC parameters (iceopt=1)
  ! DOI: 10.1175/1520-0469(2002)059<2291:CCMPRU>2.0.CO;2
  real(kind_phys), parameter :: wang_sassen_a  = 26.87_kind_phys
  real(kind_phys), parameter :: wang_sassen_b  = 0.569_kind_phys
  real(kind_phys), parameter :: wang_sassen_c  = 0.002892_kind_phys
  ! Schiller parameters (iceopt=2)
  ! DOI: 10.1029/2008JD010342
  real(kind_phys), parameter :: schiller_a     = -68.4202_kind_phys
  real(kind_phys), parameter :: schiller_b     = 0.983917_kind_phys
  real(kind_phys), parameter :: schiller_c     = 2.81795_kind_phys
  ! Wood and Field parameters (iceopt=3)
  ! DOI: 10.1175/1520-0469(2000)057<1888:RBTWCW>2.0.CO;2
  real(kind_phys), parameter :: wood_field_Kc  = 75._kind_phys
  ! Minimum grid box avg ice for having a 'cloud'
  real(kind_phys), parameter :: minice         = 1.e-12_kind_phys
  ! Minimum ice cloud fraction threshold
  real(kind_phys), parameter :: mincld         = 1.e-4_kind_phys

  real(kind_phys) :: icecrit           ! Critical RH for ice clouds in Wilson & Ballard closure
                                       ! (smaller = more ice clouds)
  real(kind_phys) :: qist_min          ! Minimum in-stratus IWC constraint [kg kg-1]
  real(kind_phys) :: qist_max          ! Maximum in-stratus IWC constraint [kg kg-1]
  logical         :: do_subgrid_growth
  logical         :: do_avg_aist_algs
  real(kind_phys) :: rair              ! Gas constant of dry air [J kg-1 K-1]

contains

  pure elemental logical function is_land(landfrac)
    real(kind_phys), intent(in) :: landfrac
    is_land = nint(landfrac) == 1
  end function is_land

!> \section arg_table_compute_cloud_fraction_two_moment_init Argument Table
!! \htmlinclude compute_cloud_fraction_two_moment_init.html
  subroutine compute_cloud_fraction_two_moment_init( &
    amIRoot, iulog, &
    premit_in, premib_in, iceopt_in, icecrit_in, &
    qist_min_in, qist_max_in, do_subgrid_growth_in, do_avg_aist_algs_in, &
    rair_in, &
    errmsg, errflg)

    ! Input arguments
    logical,         intent(in) :: amIRoot
    integer,         intent(in) :: iulog
    real(kind_phys), intent(in) :: premit_in
    real(kind_phys), intent(in) :: premib_in
    integer,         intent(in) :: iceopt_in
    real(kind_phys), intent(in) :: icecrit_in
    real(kind_phys), intent(in) :: qist_min_in
    real(kind_phys), intent(in) :: qist_max_in
    logical,         intent(in) :: do_subgrid_growth_in
    logical,         intent(in) :: do_avg_aist_algs_in
    real(kind_phys), intent(in) :: rair_in

    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errflg = 0
    errmsg = ''

    ! Set module variables from input arguments
    premit = premit_in
    premib = premib_in
    iceopt = iceopt_in
    icecrit = icecrit_in
    qist_min = qist_min_in
    qist_max = qist_max_in
    do_subgrid_growth = do_subgrid_growth_in
    do_avg_aist_algs = do_avg_aist_algs_in
    rair = rair_in

    if (amIRoot) then
      write (iulog, *) 'compute_cloud_fraction_two_moment_init parameters:'
      write (iulog, *) '  premit        = ', premit
      write (iulog, *) '  premib        = ', premib
      write (iulog, *) '  iceopt        = ', iceopt
      write (iulog, *) '  icecrit        = ', icecrit
      write (iulog, *) '  qist_min       = ', qist_min
      write (iulog, *) '  qist_max       = ', qist_max
      write (iulog, *) '  do_subgrid_growth = ', do_subgrid_growth
      write (iulog, *) '  do_avg_aist_algs  = ', do_avg_aist_algs
    end if

  end subroutine compute_cloud_fraction_two_moment_init

  pure subroutine astG_PDF_single(U, p, qv, landfrac, snowh, a, Ga, rhminl, rhminl_adj_land, rhminh, orhmin)
    ! ---------------------------------------------------------
    ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the
    ! analytical formulation of triangular PDF.
    ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)',
    ! so using constant 'dV' assume that width is proportional
    ! to the saturation specific humidity.
    !   dV ~ 0.1.
    !   cldrh : RH of in-stratus( = 1 if no supersaturation)
    ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is
    ! G is discontinuous across U = 1.  In fact, it does not
    ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that
    ! they will produce the same results.
    ! ---------------------------------------------------------

    real(kind_phys), intent(in)  :: U               ! Relative humidity
    real(kind_phys), intent(in)  :: p               ! Pressure [Pa]
    real(kind_phys), intent(in)  :: qv              ! Grid-mean water vapor specific humidity [kg kg-1]
    real(kind_phys), intent(in)  :: landfrac        ! Land fraction
    real(kind_phys), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

    real(kind_phys), intent(out) :: a               ! Stratus fraction
    real(kind_phys), intent(out) :: Ga              ! dU/da

    real(kind_phys), intent(in)  :: rhminl          ! Critical relative humidity for low-level  liquid stratus
    real(kind_phys), intent(in)  :: rhminl_adj_land ! Adjustment drop of rhminl over the land
    real(kind_phys), intent(in)  :: rhminh          ! Critical relative humidity for high-level liquid stratus

    real(kind_phys), optional, intent(out) :: orhmin ! Critical RH

    ! Local variables
    real(kind_phys), parameter :: cldrh = 1.0_kind_phys  ! RH of stratus cloud

    real(kind_phys) :: dV     ! Width of triangular PDF
    real(kind_phys) :: rhmin  ! Critical RH
    real(kind_phys) :: rhwght


    ! ---------------- !
    ! Main computation !
    ! ---------------- !

    if (p >= premib) then
      if (is_land(landfrac) .and. (snowh <= 0.000001_kind_phys)) then
        rhmin = rhminl - rhminl_adj_land
      else
        rhmin = rhminl
      end if

      dV = cldrh - rhmin

      if (U >= 1._kind_phys) then
        a = 1._kind_phys
        Ga = 1.e10_kind_phys
      elseif (U > (cldrh - dV/6._kind_phys) .and. U < 1._kind_phys) then
        a = 1._kind_phys - (-3._kind_phys/sqrt(2._kind_phys)*(U - cldrh)/dV)**(2._kind_phys/3._kind_phys)
        Ga = dV/sqrt(2._kind_phys)*sqrt(1._kind_phys - a)
      elseif (U > (cldrh - dV) .and. U <= (cldrh - dV/6._kind_phys)) then
        a = 4._kind_phys*(cos((1._kind_phys/3._kind_phys)*&
            (acos((3._kind_phys/2._kind_phys/sqrt(2._kind_phys))* &
            (1._kind_phys + (U - cldrh)/dV)) - 2._kind_phys*3.141592_kind_phys)))**2._kind_phys
        Ga = dV/sqrt(2._kind_phys)*(1._kind_phys/sqrt(a) - sqrt(a))
      elseif (U <= (cldrh - dV)) then
        a = 0._kind_phys
        Ga = 1.e10_kind_phys
      end if

      if (freeze_dry) then
        a = a*max(0.15_kind_phys, min(1.0_kind_phys, qv/0.0030_kind_phys))
        Ga = Ga/max(0.15_kind_phys, min(1.0_kind_phys, qv/0.0030_kind_phys))
      end if
    elseif (p < premit) then
      rhmin = rhminh
      dV = cldrh - rhmin

      if (U >= 1._kind_phys) then
        a = 1._kind_phys
        Ga = 1.e10_kind_phys
      elseif (U > (cldrh - dV/6._kind_phys) .and. U < 1._kind_phys) then
        a = 1._kind_phys - (-3._kind_phys/sqrt(2._kind_phys)*(U - cldrh)/dV)**(2._kind_phys/3._kind_phys)
        Ga = dV/sqrt(2._kind_phys)*sqrt(1._kind_phys - a)
      elseif (U > (cldrh - dV) .and. U <= (cldrh - dV/6._kind_phys)) then
        a = 4._kind_phys*(cos((1._kind_phys/3._kind_phys)*&
            (acos((3._kind_phys/2._kind_phys/sqrt(2._kind_phys))* &
            (1._kind_phys + (U - cldrh)/dV)) - 2._kind_phys*3.141592_kind_phys)))**2._kind_phys
        Ga = dV/sqrt(2._kind_phys)*(1._kind_phys/sqrt(a) - sqrt(a))
      elseif (U <= (cldrh - dV)) then
        a = 0._kind_phys
        Ga = 1.e10_kind_phys
      end if
    else
      rhwght = (premib - (max(p, premit)))/(premib - premit)
      rhmin = rhminh*rhwght + rhminl*(1.0_kind_phys - rhwght)

      dV = cldrh - rhmin

      if (U >= 1._kind_phys) then
        a = 1._kind_phys
        Ga = 1.e10_kind_phys
      elseif (U > (cldrh - dV/6._kind_phys) .and. U < 1._kind_phys) then
        a = 1._kind_phys - (-3._kind_phys/sqrt(2._kind_phys)*(U - cldrh)/dV)**(2._kind_phys/3._kind_phys)
        Ga = dV/sqrt(2._kind_phys)*sqrt(1._kind_phys - a)
      elseif (U > (cldrh - dV) .and. U <= (cldrh - dV/6._kind_phys)) then
        a = 4._kind_phys*(cos((1._kind_phys/3._kind_phys)*&
            (acos((3._kind_phys/2._kind_phys/sqrt(2._kind_phys))* &
            (1._kind_phys + (U - cldrh)/dV)) - 2._kind_phys*3.141592_kind_phys)))**2._kind_phys
        Ga = dV/sqrt(2._kind_phys)*(1._kind_phys/sqrt(a) - sqrt(a))
      elseif (U <= (cldrh - dV)) then
        a = 0._kind_phys
        Ga = 1.e10_kind_phys
      end if
    end if

    if (present(orhmin)) orhmin = rhmin

  end subroutine astG_PDF_single

  pure subroutine astG_PDF(U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol, &
                           rhminl_in, rhminl_adj_land_in, rhminh_in)
    ! ---------------------------------------------------------
    ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the
    ! analytical formulation of triangular PDF.
    ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)',
    ! so using constant 'dV' assume that width is proportional
    ! to the saturation specific humidity.
    !   dV ~ 0.1.
    !   cldrh : RH of in-stratus( = 1 if no supersaturation)
    ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is
    ! G is discontinuous across U = 1.  In fact, it does not
    ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that
    ! they will produce the same results.
    ! ---------------------------------------------------------

    real(kind_phys), intent(in)  :: U_in(:)               ! Relative humidity
    real(kind_phys), intent(in)  :: p_in(:)               ! Pressure [Pa]
    real(kind_phys), intent(in)  :: qv_in(:)              ! Grid-mean water vapor specific humidity [kg kg-1]
    real(kind_phys), intent(in)  :: landfrac_in(:)        ! Land fraction
    real(kind_phys), intent(in)  :: snowh_in(:)           ! Snow depth (liquid water equivalent)

    real(kind_phys), intent(out) :: a_out(:)              ! Stratus fraction
    real(kind_phys), intent(out) :: Ga_out(:)             ! dU/da
    integer,         intent(in)  :: ncol

    real(kind_phys), intent(in)  :: rhminl_in(:)          ! Critical relative humidity for low-level  liquid stratus
    real(kind_phys), intent(in)  :: rhminl_adj_land_in(:) ! Adjustment drop of rhminl over the land
    real(kind_phys), intent(in)  :: rhminh_in(:)          ! Critical relative humidity for high-level liquid stratus

    real(kind_phys) :: rhminl          ! Critical relative humidity for low-level  liquid stratus
    real(kind_phys) :: rhminl_adj_land ! Adjustment drop of rhminl over the land
    real(kind_phys) :: rhminh          ! Critical relative humidity for high-level liquid stratus

    real(kind_phys) :: U               ! Relative humidity
    real(kind_phys) :: p               ! Pressure [Pa]
    real(kind_phys) :: qv              ! Grid-mean water vapor specific humidity [kg kg-1]
    real(kind_phys) :: landfrac        ! Land fraction
    real(kind_phys) :: snowh           ! Snow depth (liquid water equivalent)

    real(kind_phys) :: a               ! Stratus fraction
    real(kind_phys) :: Ga              ! dU/da

    ! Local variables
    real(kind_phys), parameter :: cldrh = 1.0_kind_phys          ! RH of stratus cloud
    integer         :: i               ! Loop indexes
    real(kind_phys) :: dV              ! Width of triangular PDF
    real(kind_phys) :: rhmin           ! Critical RH
    real(kind_phys) :: rhwght


    ! ---------------- !
    ! Main computation !
    ! ---------------- !

    a_out(:) = 0._kind_phys
    Ga_out(:) = 0._kind_phys

    do i = 1, ncol
      U = U_in(i)
      p = p_in(i)
      qv = qv_in(i)
      landfrac = landfrac_in(i)
      snowh = snowh_in(i)

      rhminl = rhminl_in(i)
      rhminl_adj_land = rhminl_adj_land_in(i)
      rhminh = rhminh_in(i)

      if (p >= premib) then
        if (is_land(landfrac_in(i)) .and. (snowh <= 0.000001_kind_phys)) then
          rhmin = rhminl - rhminl_adj_land
        else
          rhmin = rhminl
        end if

        dV = cldrh - rhmin

        if (U >= 1._kind_phys) then
          a = 1._kind_phys
          Ga = 1.e10_kind_phys
        elseif (U > (cldrh - dV/6._kind_phys) .and. U < 1._kind_phys) then
          a = 1._kind_phys - (-3._kind_phys/sqrt(2._kind_phys)*(U - cldrh)/dV)**(2._kind_phys/3._kind_phys)
          Ga = dV/sqrt(2._kind_phys)*sqrt(1._kind_phys - a)
        elseif (U > (cldrh - dV) .and. U <= (cldrh - dV/6._kind_phys)) then
          a = 4._kind_phys*(cos((1._kind_phys/3._kind_phys)*&
              (acos((3._kind_phys/2._kind_phys/sqrt(2._kind_phys))* &
               (1._kind_phys + (U - cldrh)/dV)) - 2._kind_phys*3.141592_kind_phys)))**2._kind_phys
          Ga = dV/sqrt(2._kind_phys)*(1._kind_phys/sqrt(a) - sqrt(a))
        elseif (U <= (cldrh - dV)) then
          a = 0._kind_phys
          Ga = 1.e10_kind_phys
        end if

        if (freeze_dry) then
          a = a*max(0.15_kind_phys, min(1.0_kind_phys, qv/0.0030_kind_phys))
          Ga = Ga/max(0.15_kind_phys, min(1.0_kind_phys, qv/0.0030_kind_phys))
        end if
      elseif (p < premit) then
        rhmin = rhminh
        dV = cldrh - rhmin

        if (U >= 1._kind_phys) then
          a = 1._kind_phys
          Ga = 1.e10_kind_phys
        elseif (U > (cldrh - dV/6._kind_phys) .and. U < 1._kind_phys) then
          a = 1._kind_phys - (-3._kind_phys/sqrt(2._kind_phys)*(U - cldrh)/dV)**(2._kind_phys/3._kind_phys)
          Ga = dV/sqrt(2._kind_phys)*sqrt(1._kind_phys - a)
        elseif (U > (cldrh - dV) .and. U <= (cldrh - dV/6._kind_phys)) then
          a = 4._kind_phys*(cos((1._kind_phys/3._kind_phys)*&
             (acos((3._kind_phys/2._kind_phys/sqrt(2._kind_phys))* &
              (1._kind_phys + (U - cldrh)/dV)) - 2._kind_phys*3.141592_kind_phys)))**2._kind_phys
          Ga = dV/sqrt(2._kind_phys)*(1._kind_phys/sqrt(a) - sqrt(a))
        elseif (U <= (cldrh - dV)) then
          a = 0._kind_phys
          Ga = 1.e10_kind_phys
        end if
      else
        rhwght = (premib - (max(p, premit)))/(premib - premit)
        rhmin = rhminh*rhwght + rhminl*(1.0_kind_phys - rhwght)
        dV = cldrh - rhmin

        if (U >= 1._kind_phys) then
          a = 1._kind_phys
          Ga = 1.e10_kind_phys
        elseif (U > (cldrh - dV/6._kind_phys) .and. U < 1._kind_phys) then
          a = 1._kind_phys - (-3._kind_phys/sqrt(2._kind_phys)*(U - cldrh)/dV)**(2._kind_phys/3._kind_phys)
          Ga = dV/sqrt(2._kind_phys)*sqrt(1._kind_phys - a)
        elseif (U > (cldrh - dV) .and. U <= (cldrh - dV/6._kind_phys)) then
          a = 4._kind_phys*(cos((1._kind_phys/3._kind_phys)*&
              (acos((3._kind_phys/2._kind_phys/sqrt(2._kind_phys))* &
               (1._kind_phys + (U - cldrh)/dV)) - 2._kind_phys*3.141592_kind_phys)))**2._kind_phys
          Ga = dV/sqrt(2._kind_phys)*(1._kind_phys/sqrt(a) - sqrt(a))
        elseif (U <= (cldrh - dV)) then
          a = 0._kind_phys
          Ga = 1.e10_kind_phys
        end if
      end if

      a_out(i) = a
      Ga_out(i) = Ga
    end do

  end subroutine astG_PDF

  pure subroutine astG_RHU_single(U, p, qv, landfrac, snowh, a, Ga, rhminl, rhminl_adj_land, rhminh, orhmin)

    ! ---------------------------------------------------------
    ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the
    ! CAM3.5 cloud fraction formula.
    ! Below is valid only for CAMUW at 1.9x2.5 fv dynamics core
    ! For the other cases, I should re-define 'rhminl,rhminh' &
    ! 'premib,premit'.
    ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is
    ! G is discontinuous across U = 1.
    ! ---------------------------------------------------------

    real(kind_phys), intent(in)  :: U               ! Relative humidity
    real(kind_phys), intent(in)  :: p               ! Pressure [Pa]
    real(kind_phys), intent(in)  :: qv              ! Grid-mean water vapor specific humidity [kg kg-1]
    real(kind_phys), intent(in)  :: landfrac        ! Land fraction
    real(kind_phys), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

    real(kind_phys), intent(out) :: a               ! Stratus fraction
    real(kind_phys), intent(out) :: Ga              ! dU/da

    real(kind_phys), intent(in)  :: rhminl          ! Critical relative humidity for low-level  liquid stratus
    real(kind_phys), intent(in)  :: rhminl_adj_land ! Adjustment drop of rhminl over the land
    real(kind_phys), intent(in)  :: rhminh          ! Critical relative humidity for high-level liquid stratus

    real(kind_phys), optional, intent(out) :: orhmin ! Critical RH

    ! Local variables
    real(kind_phys) rhmin                      ! Critical RH
    real(kind_phys) rhdif                      ! Factor for stratus fraction
    real(kind_phys) rhwght


    ! ---------------- !
    ! Main computation !
    ! ---------------- !

    if (p >= premib) then

      if (is_land(landfrac) .and. (snowh <= 0.000001_kind_phys)) then
        rhmin = rhminl - rhminl_adj_land
      else
        rhmin = rhminl
      end if
      rhdif = (U - rhmin)/(1.0_kind_phys - rhmin)
      a = min(1._kind_phys, (max(rhdif, 0.0_kind_phys))**2)
      if ((U >= 1._kind_phys) .or. (U <= rhmin)) then
        Ga = 1.e20_kind_phys
      else
        Ga = 0.5_kind_phys*(1._kind_phys - rhmin)*((1._kind_phys - rhmin)/(U - rhmin))
      end if
      if (freeze_dry) then
        a = a*max(0.15_kind_phys, min(1.0_kind_phys, qv/0.0030_kind_phys))
        Ga = Ga/max(0.15_kind_phys, min(1.0_kind_phys, qv/0.0030_kind_phys))
      end if

    else if (p < premit) then
      rhmin = rhminh
      rhdif = (U - rhmin)/(1.0_kind_phys - rhmin)
      a = min(1._kind_phys, (max(rhdif, 0._kind_phys))**2)
      if ((U >= 1._kind_phys) .or. (U <= rhmin)) then
        Ga = 1.e20_kind_phys
      else
        Ga = 0.5_kind_phys*(1._kind_phys - rhmin)*((1._kind_phys - rhmin)/(U - rhmin))
      end if
    else
      rhwght = (premib - (max(p, premit)))/(premib - premit)
      rhmin = rhminh*rhwght + rhminl*(1.0_kind_phys - rhwght)

      rhdif = (U - rhmin)/(1.0_kind_phys - rhmin)
      a = min(1._kind_phys, (max(rhdif, 0._kind_phys))**2)
      if ((U >= 1._kind_phys) .or. (U <= rhmin)) then
        Ga = 1.e10_kind_phys
      else
        Ga = 0.5_kind_phys*(1._kind_phys - rhmin)*((1._kind_phys - rhmin)/(U - rhmin))
      end if
    end if

    if (present(orhmin)) orhmin = rhmin

  end subroutine astG_RHU_single

  pure subroutine astG_RHU(U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol, &
                           rhminl_in, rhminl_adj_land_in, rhminh_in)

    ! ---------------------------------------------------------
    ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the
    ! CAM3.5 cloud fraction formula.
    ! Below is valid only for CAMUW at 1.9x2.5 fv dynamics core
    ! For the other cases, I should re-define 'rhminl,rhminh' &
    ! 'premib,premit'.
    ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is
    ! G is discontinuous across U = 1.
    ! ---------------------------------------------------------

    real(kind_phys), intent(in)  :: U_in(:)               ! Relative humidity
    real(kind_phys), intent(in)  :: p_in(:)               ! Pressure [Pa]
    real(kind_phys), intent(in)  :: qv_in(:)              ! Grid-mean water vapor specific humidity [kg kg-1]
    real(kind_phys), intent(in)  :: landfrac_in(:)        ! Land fraction
    real(kind_phys), intent(in)  :: snowh_in(:)           ! Snow depth (liquid water equivalent)

    real(kind_phys), intent(out) :: a_out(:)              ! Stratus fraction
    real(kind_phys), intent(out) :: Ga_out(:)             ! dU/da
    integer,         intent(in)  :: ncol

    real(kind_phys), intent(in)  :: rhminl_in(:)          ! Critical relative humidity for low-level  liquid stratus
    real(kind_phys), intent(in)  :: rhminl_adj_land_in(:) ! Adjustment drop of rhminl over the land
    real(kind_phys), intent(in)  :: rhminh_in(:)          ! Critical relative humidity for high-level liquid stratus

    real(kind_phys) :: U               ! Relative humidity
    real(kind_phys) :: p               ! Pressure [Pa]
    real(kind_phys) :: qv              ! Grid-mean water vapor specific humidity [kg kg-1]
    real(kind_phys) :: landfrac        ! Land fraction
    real(kind_phys) :: snowh           ! Snow depth (liquid water equivalent)

    real(kind_phys) :: rhminl          ! Critical relative humidity for low-level  liquid stratus
    real(kind_phys) :: rhminl_adj_land ! Adjustment drop of rhminl over the land
    real(kind_phys) :: rhminh          ! Critical relative humidity for high-level liquid stratus

    real(kind_phys) :: a               ! Stratus fraction
    real(kind_phys) :: Ga              ! dU/da

    ! Local variables
    integer :: i
    real(kind_phys) :: rhmin                      ! Critical RH
    real(kind_phys) :: rhdif                      ! Factor for stratus fraction
    real(kind_phys) :: rhwght


    ! ---------------- !
    ! Main computation !
    ! ---------------- !

    a_out(:) = 0._kind_phys
    Ga_out(:) = 0._kind_phys

    do i = 1, ncol
      U = U_in(i)
      p = p_in(i)
      qv = qv_in(i)
      landfrac = landfrac_in(i)
      snowh = snowh_in(i)

      rhminl = rhminl_in(i)
      rhminl_adj_land = rhminl_adj_land_in(i)
      rhminh = rhminh_in(i)

      if (p >= premib) then
        if (is_land(landfrac_in(i)) .and. (snowh <= 0.000001_kind_phys)) then
          rhmin = rhminl - rhminl_adj_land
        else
          rhmin = rhminl
        end if
        rhdif = (U - rhmin)/(1.0_kind_phys - rhmin)
        a = min(1._kind_phys, (max(rhdif, 0.0_kind_phys))**2)
        if ((U >= 1._kind_phys) .or. (U <= rhmin)) then
          Ga = 1.e20_kind_phys
        else
          Ga = 0.5_kind_phys*(1._kind_phys - rhmin)*((1._kind_phys - rhmin)/(U - rhmin))
        end if
        if (freeze_dry) then
          a = a*max(0.15_kind_phys, min(1.0_kind_phys, qv/0.0030_kind_phys))
          Ga = Ga/max(0.15_kind_phys, min(1.0_kind_phys, qv/0.0030_kind_phys))
        end if
      else if (p < premit) then
        rhmin = rhminh
        rhdif = (U - rhmin)/(1.0_kind_phys - rhmin)
        a = min(1._kind_phys, (max(rhdif, 0._kind_phys))**2)
        if ((U >= 1._kind_phys) .or. (U <= rhmin)) then
          Ga = 1.e20_kind_phys
        else
          Ga = 0.5_kind_phys*(1._kind_phys - rhmin)*((1._kind_phys - rhmin)/(U - rhmin))
        end if
      else
        rhwght = (premib - (max(p, premit)))/(premib - premit)
        rhmin = rhminh*rhwght + rhminl*(1.0_kind_phys - rhwght)

        rhdif = (U - rhmin)/(1.0_kind_phys - rhmin)
        a = min(1._kind_phys, (max(rhdif, 0._kind_phys))**2)
        if ((U >= 1._kind_phys) .or. (U <= rhmin)) then
          Ga = 1.e10_kind_phys
        else
          Ga = 0.5_kind_phys*(1._kind_phys - rhmin)*((1._kind_phys - rhmin)/(U - rhmin))
        end if
      end if

      a_out(i) = a
      Ga_out(i) = Ga
    end do

  end subroutine astG_RHU

  subroutine aist_single(qv, T, p, qi, landfrac, snowh, aist, &
                         rhmaxi, rhmini, rhminl, rhminl_adj_land, rhminh, &
                         qsatfac_out)

    ! ------------------------------------------
    ! Compute non-physical ice stratus fraction
    ! ------------------------------------------

    use wv_saturation, only: qsat_water, svp_water, svp_ice

    real(kind_phys), intent(in)  :: qv              ! Grid-mean water vapor[kg kg-1]
    real(kind_phys), intent(in)  :: T               ! Temperature
    real(kind_phys), intent(in)  :: p               ! Pressure [Pa]
    real(kind_phys), intent(in)  :: qi              ! Grid-mean ice water content [kg kg-1]
    real(kind_phys), intent(in)  :: landfrac        ! Land fraction
    real(kind_phys), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

    real(kind_phys), intent(out) :: aist            ! Non-physical ice stratus fraction (0<=aist<=1)

    real(kind_phys), intent(in)  :: rhmaxi
    real(kind_phys), intent(in)  :: rhmini          ! Critical relative humidity for ice stratus
    real(kind_phys), intent(in)  :: rhminl          ! Critical relative humidity for low-level  liquid stratus
    real(kind_phys), intent(in)  :: rhminl_adj_land ! Adjustment drop of rhminl over the land
    real(kind_phys), intent(in)  :: rhminh          ! Critical relative humidity for high-level liquid stratus
    real(kind_phys), optional, intent(out) :: qsatfac_out ! Subgrid scaling factor for qsat

    ! Local variables
    real(kind_phys) :: rhmin               ! Critical RH
    real(kind_phys) :: rhwght

    real(kind_phys) :: ttmp                ! Limited temperature
    real(kind_phys) :: icicval             ! Empirical IWC value [ kg kg-1 ]
    real(kind_phys) :: rho                 ! Local air density
    real(kind_phys) :: esl                 ! Liq sat vapor pressure
    real(kind_phys) :: esi                 ! Ice sat vapor pressure
    real(kind_phys) :: ncf, phi            ! Wilson and Ballard parameters
    real(kind_phys) :: es, qs

    real(kind_phys) :: rhi                 ! grid box averaged relative humidity over ice
    real(kind_phys) :: icimr               ! in cloud ice mixing ratio
    real(kind_phys) :: rhdif               ! working variable for slingo scheme


    if (present(qsatfac_out)) qsatfac_out = 1.0_kind_phys

    ! ---------------- !
    ! Main computation !
    ! ---------------- !

    call qsat_water(T, p, es, qs)
    esl = svp_water(T)
    esi = svp_ice(T)

    if (iceopt < 3) then
      if (iceopt == 1) then
        ttmp = max(195._kind_phys, min(T, 253._kind_phys)) - 273.16_kind_phys
        icicval = wang_sassen_a + wang_sassen_b*ttmp + wang_sassen_c*ttmp**2._kind_phys
        rho = p/(rair*T)
        icicval = icicval*1.e-6_kind_phys/rho
      else
        ttmp = max(190._kind_phys, min(T, 273.16_kind_phys))
        icicval = 10._kind_phys**(schiller_a*schiller_b**ttmp + schiller_c)
        icicval = icicval*1.e-6_kind_phys*18._kind_phys/28.97_kind_phys
      end if
      aist = max(0._kind_phys, min(qi/icicval, 1._kind_phys))
    elseif (iceopt == 3) then
      aist = 1._kind_phys - exp(-wood_field_Kc*qi/(qs*(esi/esl)))
      aist = max(0._kind_phys, min(aist, 1._kind_phys))
    elseif (iceopt == 4) then
      if (p >= premib) then
        if (is_land(landfrac) .and. (snowh <= 0.000001_kind_phys)) then
          rhmin = rhminl - rhminl_adj_land
        else
          rhmin = rhminl
        end if
      else if (p < premit) then
        rhmin = rhminh
      else
        rhwght = (premib - (max(p, premit)))/(premib - premit)
        rhmin = rhminh*rhwght + rhminl*(1.0_kind_phys - rhwght)
      end if
      ncf = qi/((1._kind_phys - icecrit)*qs)
      if (ncf <= 0._kind_phys) then
        aist = 0._kind_phys
      elseif (ncf > 0._kind_phys .and. ncf <= 1._kind_phys/6._kind_phys) then
        aist = 0.5_kind_phys*(6._kind_phys*ncf)**(2._kind_phys/3._kind_phys)
      elseif (ncf > 1._kind_phys/6._kind_phys .and. ncf < 1._kind_phys) then
        phi = (acos(3._kind_phys*(1._kind_phys - ncf)/2._kind_phys**(3._kind_phys/2._kind_phys)) + &
               4._kind_phys*3.1415927_kind_phys)/3._kind_phys
        aist = (1._kind_phys - 4._kind_phys*cos(phi)*cos(phi))
      else
        aist = 1._kind_phys
      end if
      aist = max(0._kind_phys, min(aist, 1._kind_phys))
    else if (iceopt == 5) then
      ! set rh ice cloud fraction
      rhi = (qv + qi)/qs*(esl/esi)
      if (rhmaxi == rhmini) then
        if (rhi > rhmini) then
          rhdif = 1._kind_phys
        else
          rhdif = 0._kind_phys
        end if
      else
        rhdif = (rhi - rhmini)/(rhmaxi - rhmini)
      end if
      aist = min(1.0_kind_phys, max(rhdif, 0._kind_phys)**2)

      ! Similar to alpha in Wilson & Ballard (1999), determine a
      ! scaling factor for saturation vapor pressure that reflects
      ! the cloud fraction, rhmini, and rhmaxi.
      !
      ! NOTE: Limit qsatfac so that adjusted RHliq would be 1. or less.
      if (present(qsatfac_out) .and. do_subgrid_growth) then
        qsatfac_out = max(min(qv/qs, 1._kind_phys), (1._kind_phys - aist)*rhmini + aist*rhmaxi)
      end if

      ! limiter to remove empty cloud and ice with no cloud
      ! and set icecld fraction to mincld if ice exists

      if (qi < minice) then
        aist = 0._kind_phys
      else
        aist = max(mincld, aist)
      end if

      ! enforce limits on icimr
      if (qi >= minice) then
        icimr = qi/aist

        !minimum
        if (icimr < qist_min) then
          if (do_avg_aist_algs) then
            !
            ! Take the geometric mean of the iceopt=4 and iceopt=5 values.
            ! Mods developed by Thomas Toniazzo for NorESM.
            aist = max(0._kind_phys, min(1._kind_phys, sqrt(aist*qi/qist_min)))
          else
            !
            ! Default for iceopt=5
            aist = max(0._kind_phys, min(1._kind_phys, qi/qist_min))
          end if
        end if
        !maximum
        if (icimr > qist_max) then
          aist = max(0._kind_phys, min(1._kind_phys, qi/qist_max))
        end if
      end if
    end if

    ! 0.999_kind_phys is added to prevent infinite 'ql_st' at the end of instratus_condensate
    ! computed after updating 'qi_st'.
    aist = max(0._kind_phys, min(aist, 0.999_kind_phys))

  end subroutine aist_single

  subroutine aist_vector(qv_in, T_in, p_in, qi_in, ni_in, landfrac_in, snowh_in, aist_out, ncol, &
                         rhmaxi_in, rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in, &
                         qsatfac_out)

    ! ------------------------------------------
    ! Compute non-physical ice stratus fraction
    ! ------------------------------------------

    use wv_saturation, only: qsat_water, svp_water_vect, svp_ice_vect

    real(kind_phys), intent(in)  :: qv_in(:)              ! Grid-mean water vapor [kg kg-1]
    real(kind_phys), intent(in)  :: T_in(:)               ! Temperature
    real(kind_phys), intent(in)  :: p_in(:)               ! Pressure [Pa]
    real(kind_phys), intent(in)  :: qi_in(:)              ! Grid-mean ice water content [kg kg-1]
    real(kind_phys), intent(in)  :: ni_in(:)              ! Grid-mean ice water number concentration [# kg-1]
    real(kind_phys), intent(in)  :: landfrac_in(:)        ! Land fraction
    real(kind_phys), intent(in)  :: snowh_in(:)           ! Snow depth (liquid water equivalent)

    real(kind_phys), intent(out) :: aist_out(:)           ! Non-physical ice stratus fraction ( 0<= aist <= 1 )
    integer,         intent(in)  :: ncol

    real(kind_phys), intent(in)  :: rhmaxi_in(:)
    real(kind_phys), intent(in)  :: rhmini_in(:)          ! Critical relative humidity for ice stratus
    real(kind_phys), intent(in)  :: rhminl_in(:)          ! Critical relative humidity for low-level liquid stratus
    real(kind_phys), intent(in)  :: rhminl_adj_land_in(:) ! Adjustment drop of rhminl over the land
    real(kind_phys), intent(in)  :: rhminh_in(:)          ! Critical relative humidity for high-level liquid stratus
    real(kind_phys), optional, intent(out) :: qsatfac_out(:) ! Subgrid scaling factor for qsat

    ! Local variables

    real(kind_phys) :: qv                              ! Grid-mean water vapor [kg kg-1]
    real(kind_phys) :: T                               ! Temperature
    real(kind_phys) :: p                               ! Pressure [Pa]
    real(kind_phys) :: qi                              ! Grid-mean ice water content [kg kg-1]
    real(kind_phys) :: ni
    real(kind_phys) :: landfrac                        ! Land fraction
    real(kind_phys) :: snowh                           ! Snow depth (liquid water equivalent)

    real(kind_phys) :: rhmaxi                          ! Critical relative humidity for ice stratus
    real(kind_phys) :: rhmini                          ! Critical relative humidity for ice stratus
    real(kind_phys) :: rhminl                          ! Critical relative humidity for low-level  liquid stratus
    real(kind_phys) :: rhminl_adj_land                 ! Adjustment drop of rhminl over the land
    real(kind_phys) :: rhminh                          ! Critical relative humidity for high-level liquid stratus

    real(kind_phys) :: aist                            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

    real(kind_phys) :: rhmin                           ! Critical RH
    real(kind_phys) :: rhwght

    real(kind_phys) :: ttmp                            ! Limited temperature
    real(kind_phys) :: icicval                         ! Empirical IWC value [ kg kg-1 ]
    real(kind_phys) :: rho                             ! Local air density
    real(kind_phys) :: esl(ncol)                       ! Liq sat vapor pressure
    real(kind_phys) :: esi(ncol)                       ! Ice sat vapor pressure
    real(kind_phys) :: ncf, phi                        ! Wilson and Ballard parameters
    real(kind_phys) :: qs
    real(kind_phys) :: esat_in(ncol)
    real(kind_phys) :: qsat_in(ncol)

    real(kind_phys) :: rhi                             ! grid box averaged relative humidity over ice
    real(kind_phys) :: icimr                           ! in cloud ice mixing ratio
    real(kind_phys) :: rhdif                           ! working variable for slingo scheme
    real(kind_phys) :: nil                             ! Ice number concentration [#/L]

    ! Heymsfield/Boudala fit parameters (iceopt=6)
    ! Boudala et al. (2002), Heymsfield et al. (2012)
    real(kind_phys), parameter :: heymsfield_a = 6.73834e-08_kind_phys
    real(kind_phys), parameter :: heymsfield_b = 0.0533110_kind_phys
    real(kind_phys), parameter :: heymsfield_c = 0.3493813_kind_phys

    integer         :: i


    if (present(qsatfac_out)) qsatfac_out = 1.0_kind_phys

    ! ---------------- !
    ! Main computation !
    ! ---------------- !

    aist_out(:) = 0._kind_phys
    esat_in(:) = 0._kind_phys
    qsat_in(:) = 0._kind_phys

    call qsat_water(T_in(1:ncol), p_in(1:ncol), esat_in(1:ncol), qsat_in(1:ncol), ncol)
    call svp_water_vect(T_in(1:ncol), esl(1:ncol), ncol)
    call svp_ice_vect(T_in(1:ncol), esi(1:ncol), ncol)

    do i = 1, ncol

      landfrac = landfrac_in(i)
      snowh = snowh_in(i)
      T = T_in(i)
      qv = qv_in(i)
      p = p_in(i)
      qi = qi_in(i)
      ni = ni_in(i)
      qs = qsat_in(i)

      rhmaxi = rhmaxi_in(i)
      rhmini = rhmini_in(i)
      rhminl = rhminl_in(i)
      rhminl_adj_land = rhminl_adj_land_in(i)
      rhminh = rhminh_in(i)

      if (iceopt < 3) then
        if (iceopt == 1) then
          ttmp = max(195._kind_phys, min(T, 253._kind_phys)) - 273.16_kind_phys
          icicval = wang_sassen_a + wang_sassen_b*ttmp + wang_sassen_c*ttmp**2._kind_phys
          rho = p/(rair*T)
          icicval = icicval*1.e-6_kind_phys/rho
        else
          ttmp = max(190._kind_phys, min(T, 273.16_kind_phys))
          icicval = 10._kind_phys**(schiller_a*schiller_b**ttmp + schiller_c)
          icicval = icicval*1.e-6_kind_phys*18._kind_phys/28.97_kind_phys
        end if
        aist = max(0._kind_phys, min(qi/icicval, 1._kind_phys))
      elseif (iceopt == 3) then
        aist = 1._kind_phys - exp(-wood_field_Kc*qi/(qs*(esi(i)/esl(i))))
        aist = max(0._kind_phys, min(aist, 1._kind_phys))
      elseif (iceopt == 4) then
        if (p >= premib) then
          if (is_land(landfrac_in(i)) .and. (snowh <= 0.000001_kind_phys)) then
            rhmin = rhminl - rhminl_adj_land
          else
            rhmin = rhminl
          end if
        elseif (p < premit) then
          rhmin = rhminh
        else
          rhwght = (premib - (max(p, premit)))/(premib - premit)
          rhmin = rhminh*rhwght + rhminl*(1.0_kind_phys - rhwght)
        end if
        ncf = qi/((1._kind_phys - icecrit)*qs)
        if (ncf <= 0._kind_phys) then
          aist = 0._kind_phys
        elseif (ncf > 0._kind_phys .and. ncf <= 1._kind_phys/6._kind_phys) then
          aist = 0.5_kind_phys*(6._kind_phys*ncf)**(2._kind_phys/3._kind_phys)
        elseif (ncf > 1._kind_phys/6._kind_phys .and. ncf < 1._kind_phys) then
         phi = (acos(3._kind_phys*(1._kind_phys-ncf)/2._kind_phys**(3._kind_phys/2._kind_phys))+&
               4._kind_phys*3.1415927_kind_phys)/3._kind_phys
          aist = (1._kind_phys - 4._kind_phys*cos(phi)*cos(phi))
        else
          aist = 1._kind_phys
        end if
        aist = max(0._kind_phys, min(aist, 1._kind_phys))
      elseif (iceopt == 5) then
        ! set rh ice cloud fraction
        rhi = (qv + qi)/qs*(esl(i)/esi(i))
        if (rhmaxi == rhmini) then
          if (rhi > rhmini) then
            rhdif = 1._kind_phys
          else
            rhdif = 0._kind_phys
          end if
        else
          rhdif = (rhi - rhmini)/(rhmaxi - rhmini)
        end if
        aist = min(1.0_kind_phys, max(rhdif, 0._kind_phys)**2)

      elseif (iceopt == 6) then
        !----- ICE CLOUD OPTION 6: fit based on T and Number (Gettelman: based on Heymsfield obs)
        ! Use observations from
        ! Heymsfield et al. 2013 (https://doi.org/10.1175/JAS-D-12-0124.1) of IWC and Ni v. Temp
        ! Multivariate fit follows form of
        ! Boudala 2002 (https://doi.org/10.1002/joc.774): ICIWC = a * exp(b*T) * N^c
        ! a=6.73e-8, b=0.05, c=0.349
        ! N is #/L, so need to convert Ni_L=N*rhoa/1000.
        rho = p/(rair*T)
        nil = ni*rho/1000._kind_phys
        icicval = heymsfield_a*exp(heymsfield_b*T)*nil**heymsfield_c
        ! result is in g m-3, convert to kg H2O / kg air (icimr...)
        icicval = icicval/rho/1000._kind_phys
        aist = max(0._kind_phys, min(qi/icicval, 1._kind_phys))
        aist = min(aist, 1._kind_phys)

      end if

      if (iceopt == 5 .or. iceopt == 6) then
        ! Similar to alpha in Wilson & Ballard (1999), determine a
        ! scaling factor for saturation vapor pressure that reflects
        ! the cloud fraction, rhmini, and rhmaxi.
        ! https://doi.org/10.1002/qj.49712555707
        !
        ! NOTE: Limit qsatfac so that adjusted RHliq would be 1. or less.
        if (present(qsatfac_out) .and. do_subgrid_growth) then
          qsatfac_out(i) = max(min(qv/qs, 1._kind_phys), (1._kind_phys - aist)*rhmini + aist*rhmaxi)
        end if

        ! limiter to remove empty cloud and ice with no cloud
        ! and set icecld fraction to mincld if ice exists
        if (qi < minice) then
          aist = 0._kind_phys
        else
          aist = max(mincld, aist)
        end if

        ! enforce limits on icimr
        if (qi >= minice) then
          icimr = qi/aist

          ! minimum:
          if (icimr < qist_min) then
            if (do_avg_aist_algs) then
              ! Take the geometric mean of the iceopt=4 and iceopt=5 values.
              ! Mods developed by Thomas Toniazzo for NorESM.
              aist = max(0._kind_phys, min(1._kind_phys, sqrt(aist*qi/qist_min)))
            else
              ! Default for iceopt=5
              aist = max(0._kind_phys, min(1._kind_phys, qi/qist_min))
            end if
          end if
          !maximum
          if (icimr > qist_max) then
            aist = max(0._kind_phys, min(1._kind_phys, qi/qist_max))
          end if

        end if
      end if

      ! 0.999_kind_phys is added to prevent infinite 'ql_st' at the end of instratus_condensate
      ! computed after updating 'qi_st'.
      aist = max(0._kind_phys, min(aist, 0.999_kind_phys))
      aist_out(i) = aist
    end do

  end subroutine aist_vector

end module compute_cloud_fraction_two_moment
