!=============================================================================
! Common dust module
!   Portable: host constants are passed as arguments; errors are reported
!   through errmsg/errflg for the caller to handle.
!=============================================================================
module dust_common

  use ccpp_kinds,       only: kind_phys

  implicit none
  private

  public :: dust_set_params

  ! Dust aerosol material density used by the emissions number/mass flux
  ! conversion and the settling parameters below (value from CAM mo_constants).
  real(kind_phys), public, parameter :: dust_density = 2.5e+3_kind_phys  ! [kg m-3]

contains

  !=============================================================================
  !
  ! !DESCRIPTION:
  !
  ! Compute source efficiency factor from topography
  ! Initialize other variables used in subroutine Dust:
  ! ovr_src_snk_mss(m,n) and tmp1.
  ! Define particle diameter and density needed by atm model
  ! as well as by dry dep model
  ! Source: Paul Ginoux (for source efficiency factor)
  ! Modifications by C. Zender and later by S. Levis
  ! Rest of subroutine from C. Zender's dust model
  !=============================================================================
  subroutine dust_set_params(nbin, dmt_grd, dmt_vwr, stk_crc, pi, rair, gravit, errmsg, errflg)

    !
    ! !USES
    !
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)

    !
    ! !ARGUMENTS:
    !
    integer, intent(in)  :: nbin
    real(kind_phys),intent(in)  :: dmt_grd(:)
    real(kind_phys),intent(out) :: dmt_vwr(:)
    real(kind_phys),intent(out) :: stk_crc(:)
    real(kind_phys),intent(in)  :: pi          ! host model constants
    real(kind_phys),intent(in)  :: rair        ! gas constant for dry air (J/K/kg)
    real(kind_phys),intent(in)  :: gravit      ! gravitational acceleration (m/s2)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !
    ! !REVISION HISTORY
    ! Created by Samual Levis
    ! Revised for CAM by Natalie Mahowald
    !EOP
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    !Local Variables
    integer, parameter :: dst_src_nbr =3
    integer, parameter :: sz_nbr =200

    integer  :: m,n                     !indices
    real(kind_phys) :: dmt_min(nbin)      ![m] Size grid minimum
    real(kind_phys) :: dmt_max(nbin)      ![m] Size grid maximum
    real(kind_phys) :: dmt_ctr(nbin)      ![m] Diameter at bin center
    real(kind_phys) :: dmt_dlt(nbin)      ![m] Width of size bin
    real(kind_phys) :: slp_crc(nbin)      ![frc] Slip correction factor
    real(kind_phys) :: vlm_rsl(nbin)      ![m3 m-3] Volume concentration resolved
    real(kind_phys) :: vlc_stk(nbin)      ![m s-1] Stokes settling velocity
    real(kind_phys) :: vlc_grv(nbin)      ![m s-1] Settling velocity
    real(kind_phys) :: ryn_nbr_grv(nbin)  ![frc] Reynolds number at terminal velocity
    real(kind_phys) :: cff_drg_grv(nbin)  ![frc] Drag coefficient at terminal velocity
    real(kind_phys) :: tmp                     !temporary
    real(kind_phys) :: ln_gsd                  ![frc] ln(gsd)
    real(kind_phys) :: gsd_anl                 ![frc] Geometric standard deviation
    real(kind_phys) :: dmt_vma                 ![m] Mass median diameter analytic She84 p.75 Tabl.1
    real(kind_phys) :: dmt_nma                 ![m] Number median particle diameter
    real(kind_phys) :: lgn_dst                 !Lognormal distribution at sz_ctr
    real(kind_phys) :: eps_max                 ![frc] Relative accuracy for convergence
    real(kind_phys) :: eps_crr                 ![frc] Current relative accuracy
    real(kind_phys) :: itr_idx                 ![idx] Counting index
    real(kind_phys) :: dns_mdp                 ![kg m-3] Midlayer density
    real(kind_phys) :: mfp_atm                 ![m] Mean free path of air
    real(kind_phys) :: vsc_dyn_atm             ![kg m-1 s-1] Dynamic viscosity of air
    real(kind_phys) :: vsc_knm_atm             ![kg m-1 s-1] Kinematic viscosity of air
    real(kind_phys) :: vlc_grv_old             ![m s-1] Previous gravitational settling velocity
    real(kind_phys) :: series_ratio            !Factor for logarithmic grid
    real(kind_phys) :: lngsdsqrttwopi_rcp      !Factor in lognormal distribution
    real(kind_phys) :: sz_min(sz_nbr)          ![m] Size Bin minima
    real(kind_phys) :: sz_max(sz_nbr)          ![m] Size Bin maxima
    real(kind_phys) :: sz_ctr(sz_nbr)          ![m] Size Bin centers
    real(kind_phys) :: sz_dlt(sz_nbr)          ![m] Size Bin widths

    errmsg = ''
    errflg = 0

    stk_crc(:) = nan
    dmt_vwr(:) = nan

    ! Introducing particle diameter. Needed by atm model and by dry dep model.
    ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
    ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)

    ! Charlie allows logarithmic or linear option for size distribution
    ! however, he hardwires the distribution to logarithmic in his code
    ! therefore, I take his logarithmic code only
    ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
    ! he currently works with dst_nbr = 4, so I only take the relevant code
    ! if dust_number ever becomes different from 4, must add call grd_mk (dstpsd.F90)
    ! as done in subroutine dst_psd_ini
    ! note that here dust_number = dst_nbr

    ! Override automatic grid with preset grid if available
    do n = 1, nbin
       dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
       dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
       dmt_ctr(n) = 0.5_kind_phys * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
       dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
    end do

 ! sets dust_dmt_vwr ....

    ! Bin physical properties
    gsd_anl = 2.0_kind_phys      ! [frc] Geometric std dev PaG77 p. 2080 Table1
    ln_gsd = log(gsd_anl)

    ! Set a fundamental statistic for each bin
    dmt_vma = 2.524e-6_kind_phys ! [m] Mass median diameter analytic She84 p.75 Table1
    dmt_vma = 3.5e-6_kind_phys
    ! Compute analytic size statistics
    ! Convert mass median diameter to number median diameter (call vma2nma)
    dmt_nma = dmt_vma * exp(-3.0_kind_phys*ln_gsd*ln_gsd) ! [m]
    ! Compute resolved size statistics for each size distribution
    ! In C. Zender's code call dst_sz_rsl
    do n = 1, nbin
       series_ratio = (dmt_max(n)/dmt_min(n))**(1.0_kind_phys/sz_nbr)
       sz_min(1) = dmt_min(n)
       do m = 2, sz_nbr                            ! Loop starts at 2
          sz_min(m) = sz_min(m-1) * series_ratio
       end do

       ! Derived grid values
       do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
          sz_max(m) = sz_min(m+1)                  ! [m]
       end do
       sz_max(sz_nbr) = dmt_max(n)                 ! [m]

       ! Final derived grid values
       do m = 1, sz_nbr
          sz_ctr(m) = 0.5_kind_phys * (sz_min(m)+sz_max(m))
          sz_dlt(m) = sz_max(m)-sz_min(m)
       end do
       lngsdsqrttwopi_rcp = 1.0_kind_phys / (ln_gsd*sqrt(2.0_kind_phys*pi))
       dmt_vwr(n) = 0.0_kind_phys ! [m] Mass wgted diameter resolved
       vlm_rsl(n) = 0.0_kind_phys ! [m3 m-3] Volume concentration resolved
       do m = 1, sz_nbr
          ! Evaluate lognormal distribution for these sizes (call lgn_evl)
          tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
          lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_kind_phys*tmp*tmp) / sz_ctr(m)
          ! Integrate moments of size distribution
          dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
               pi / 6.0_kind_phys * (sz_ctr(m)**3.0_kind_phys) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
          vlm_rsl(n) = vlm_rsl(n) +                                &
               pi / 6.0_kind_phys * (sz_ctr(m)**3.0_kind_phys) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
       end do
       dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved
    end do

  ! sets stk_crc ...

    ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)
    eps_max = 1.0e-4_kind_phys
    dns_mdp = 100000._kind_phys / (295.0_kind_phys*rair) ![kg m-3] const prs_mdp & tpt_vrt
    ! Size-independent thermokinetic properties
    vsc_dyn_atm = 1.72e-5_kind_phys * ((295.0_kind_phys/273.0_kind_phys)**1.5_kind_phys) * 393.0_kind_phys / &
         (295.0_kind_phys+120.0_kind_phys)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
    mfp_atm = 2.0_kind_phys * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
         (100000._kind_phys*sqrt(8.0_kind_phys/(pi*rair*295.0_kind_phys)))
    vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

    do m = 1, nbin
       slp_crc(m) = 1.0_kind_phys + 2.0_kind_phys * mfp_atm *                      &
            (1.257_kind_phys+0.4_kind_phys*exp(-1.1_kind_phys*dmt_vwr(m)/(2.0_kind_phys*mfp_atm))) / &
            dmt_vwr(m)                      ! [frc] Slip correction factor SeP97 p.464
       vlc_stk(m) = (1.0_kind_phys/18.0_kind_phys) * dmt_vwr(m) * dmt_vwr(m) * dust_density * &
            gravit * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
    end do

    ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
    ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
    ! important and empirical drag coefficients must be employed
    ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
    ! Using Stokes' velocity rather than iterative solution with empirical
    ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

    ! Iterative solution for drag coefficient, Reynolds number, and terminal veloc
    do m = 1, nbin

       ! Initialize accuracy and counter
       eps_crr = eps_max + 1.0_kind_phys  ![frc] Current relative accuracy
       itr_idx = 0                 ![idx] Counting index

       ! Initial guess for vlc_grv is exact for Re < 0.1
       vlc_grv(m) = vlc_stk(m)     ![m s-1]
       eps_loop: do while(eps_crr > eps_max)

          ! Save terminal velocity for convergence test
          vlc_grv_old = vlc_grv(m) ![m s-1]
          ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

          ! Update drag coefficient based on new Reynolds number
          if (ryn_nbr_grv(m) < 0.1_kind_phys) then
             cff_drg_grv(m) = 24.0_kind_phys / ryn_nbr_grv(m) !Stokes' law Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0_kind_phys) then
             cff_drg_grv(m) = (24.0_kind_phys/ryn_nbr_grv(m)) *    &
                  (1.0_kind_phys + 3.0_kind_phys*ryn_nbr_grv(m)/16.0_kind_phys + &
                  9.0_kind_phys*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                  log(2.0_kind_phys*ryn_nbr_grv(m))/160.0_kind_phys)        !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 500.0_kind_phys) then
             cff_drg_grv(m) = (24.0_kind_phys/ryn_nbr_grv(m)) * &
                  (1.0_kind_phys + 0.15_kind_phys*ryn_nbr_grv(m)**0.687_kind_phys) !Sep97 p.463 (8.32)
          else if (ryn_nbr_grv(m) < 2.0e5_kind_phys) then
             cff_drg_grv(m) = 0.44_kind_phys                         !Sep97 p.463 (8.32)
          else
             write(errmsg,'(a,es9.2)') &
                  'Dustini error: Reynolds number too large in stk_crc_get(): ryn_nbr_grv(m) = ', ryn_nbr_grv(m)
             errflg = 1
             return
          end if

          ! Update terminal velocity based on new Reynolds number and drag coeff
          ! [m s-1] Terminal veloc SeP97 p.467 (8.44)
          vlc_grv(m) = sqrt(4.0_kind_phys * gravit * dmt_vwr(m) * slp_crc(m) * dust_density / &
               (3.0_kind_phys*cff_drg_grv(m)*dns_mdp))
          eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
          if (itr_idx == 12) then
             ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
             ! due to discontinuities in derivative of drag coefficient
             vlc_grv(m) = 0.5_kind_phys * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
          end if
          if (itr_idx > 20) then
             ! Dustini error: terminal velocity not converging in stk_crc_get();
             ! break the loop and accept the current value (was an iulog warning)
             exit eps_loop
          end if
          itr_idx = itr_idx + 1
       end do eps_loop  !end while
    end do   !end loop over size

    ! Compute factors to convert Stokes' settling velocities to
    ! actual settling velocities
    do m = 1, nbin
       stk_crc(m) = vlc_grv(m) / vlc_stk(m)
    end do

  end subroutine dust_set_params

end module dust_common
