module kessler

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: kessler_run ! Main routine
   public :: kessler_init ! init routine
   public :: kessler_timestep_init ! init timestep routine

   ! Private module data (constants set at initialization)
   real(kind_phys) :: rd    ! gas constant for dry air, J/(kgK)
   real(kind_phys) :: cp    ! heat capacity at constant pressure, J/(kgK)
   real(kind_phys) :: lv    ! latent heat of vaporization, J/kg
   real(kind_phys) :: psl   ! reference pressure at sea level, mb
   real(kind_phys) :: rhoqr ! density of liquid water, kg/m^3

CONTAINS

   !> \section arg_table_kessler_init  Argument Table
   !! \htmlinclude kessler_init.html
   subroutine kessler_init(rd_in, cp_in, lv_in, psl_in, rhoqr_in, errmsg, errflg)
      ! Set physical constants to be consistent with calling model
      real(kind_phys),    intent(in)  :: rd_in    ! gas constant for dry air, J/(kgK)
      real(kind_phys),    intent(in)  :: cp_in    ! heat capacity at constant pres., J/(kgK)
      real(kind_phys),    intent(in)  :: lv_in    ! latent heat of vaporization, J/kg
      real(kind_phys),    intent(in)  :: psl_in   ! reference pressure at sea level, mb
      real(kind_phys),    intent(in)  :: rhoqr_in ! density of liquid water, kg/m^3

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      rd    = rd_in
      cp    = cp_in
      lv    = lv_in
      psl   = psl_in/100._kind_phys
      rhoqr = rhoqr_in

   end subroutine kessler_init

   !> \section arg_table_kessler_timestep_init  Argument Table
   !! \htmlinclude kessler_timestep_init.html
   subroutine kessler_timestep_init(ncol, nz, pdel, pdeldry, qv, qc, qr, errmsg, errflg)
      use state_converters, only : wet_to_dry_run

      ! Dummy arguments
      integer,         intent(in)    :: ncol
      integer,         intent(in)    :: nz
      real(kind_phys), intent(in)    :: pdel(:,:)
      real(kind_phys), intent(in)    :: pdeldry(:,:)
      real(kind_phys), intent(inout) :: qv(:,:)
      real(kind_phys), intent(inout) :: qc(:,:)
      real(kind_phys), intent(inout) :: qr(:,:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer  :: k

      errflg = 0
      errmsg = ''

      call wet_to_dry_run(ncol, nz, pdel, pdeldry, qv, qc, qr, errmsg, errflg)

   end subroutine kessler_timestep_init

   !-----------------------------------------------------------------------
   !
   !  Version:  2.0
   !
   !  Date:  January 22nd, 2015
   !
   !  Change log:
   !  v2 - Added sub-cycling of rain sedimentation so as not to violate
   !       CFL condition.
   !
   !  The KESSLER subroutine implements the Kessler (1969) microphysics
   !  parameterization as described by Soong and Ogura (1973) and Klemp
   !  and Wilhelmson (1978, KW). KESSLER is called at the end of each
   !  time step and makes the final adjustments to the potential
   !  temperature and moisture variables due to microphysical processes
   !  occurring during that time step. KESSLER is called once for each
   !  vertical column of grid cells. Increments are computed and added
   !  into the respective variables. The Kessler scheme contains three
   !  moisture categories: water vapor, cloud water (liquid water that
   !  moves with the flow), and rain water (liquid water that falls
   !  relative to the surrounding air). There  are no ice categories.
   !  Variables in the column are ordered from the surface to the top.
   !
   !  SUBROUTINE KESSLER(theta, qv, qc, qr, rho, pk, dt, z, nz, rainnc)
   !
   !  Input variables:
   !     temp   - temperature (K)
   !     qv     - water vapor mixing ratio (gm/gm)
   !     qc     - cloud water mixing ratio (gm/gm)
   !     qr     - rain  water mixing ratio (gm/gm)
   !     rho    - dry air density (not mean state as in KW) (kg/m^3)
   !     pk     - Exner function  (not mean state as in KW) (p/p0)**(R/cp)
   !     dt     - time step (s)
   !     z      - heights of thermodynamic levels in the grid column (m)
   !     nz     - number of thermodynamic levels in the column
   !     precl  - Precipitation rate (m_water/s)
   !
   ! Output variables:
   !     Increments are added into t, qv, qc, qr, and rainnc which are
   !     returned to the routine from which KESSLER was called. To obtain
   !     the total precip qt, after calling the KESSLER routine, compute:
   !
   !       qt = sum over surface grid cells of (rainnc * cell area)  (kg)
   !       [here, the conversion to kg uses (10^3 kg/m^3)*(10^-3 m/mm) = 1]
   !
   !
   !  Authors: Paul Ullrich
   !           University of California, Davis
   !           Email: paullrich@ucdavis.edu
   !
   !           Based on a code by Joseph Klemp
   !           (National Center for Atmospheric Research)
   !
   !  Reference:
   !
   !    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
   !    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
   !    Radius Sphere. Journal of Advances in Modeling Earth Systems.
   !    doi:10.1002/2015MS000435
   !
   !=======================================================================

   !> \section arg_table_kessler_run  Argument Table
   !! \htmlinclude kessler_run.html
   subroutine kessler_run(ncol, nz, dt, rho, z, pk, theta,                    &
        qv, qc, qr, precl, errmsg, errflg)

      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer,          intent(in)    :: ncol       ! Number of columns
      integer,          intent(in)    :: nz         ! Number of vertical levels
      real(kind_phys),  intent(in)    :: dt         ! Time step (s)
      real(kind_phys),  intent(in)    :: rho(:,:)   ! Dry air density (kg/m^3)
      real(kind_phys),  intent(in)    :: z(:,:)     ! Heights of thermo. levels (m)
      real(kind_phys),  intent(in)    :: pk(:,:)    ! Exner function (p/p0)**(R/cp)

      real(kind_phys),  intent(inout) :: theta(:,:) ! temperature (K)
      real(kind_phys),  intent(inout) :: qv(:,:)    ! Water vapor mixing ratio (gm/gm)
      real(kind_phys),  intent(inout) :: qc(:,:)    ! Cloud water mixing ratio (gm/gm)
      real(kind_phys),  intent(inout) :: qr(:,:)    ! Rain  water mixing ratio (gm/gm)

      real(kind_phys),  intent(out)   :: precl(:)   ! Precipi tation rate (m_water / s)
      character(len=*), intent(out)   :: errmsg
      integer,          intent(out)   :: errflg

      !------------------------------------------------
      !   Local variables
      !------------------------------------------------
      real(kind_phys) :: r(nz), rhalf(nz), velqr(nz), sed(nz), pc(nz)
      real(kind_phys) :: f5, f2x, xk, ern, qrprod, prod, qvs, dt_max, dt0
      integer         :: col, klev, rainsplit, nt

      ! Initialize output variables
      precl = 0._kind_phys
      errmsg = ''
      errflg = 0

      ! Check inputs
      if (dt <= 0._kind_phys) then
         write(errmsg,*) 'KESSLER called with nonpositive dt'
         errflg = 1
         return
      end if

      !------------------------------------------------
      !   Begin calculation
      !------------------------------------------------
      f2x   = 17.27_kind_phys
      f5    = 237.3_kind_phys * f2x * lv / cp
      xk    = .2875_kind_phys  !  kappa (r/cp)
      ! Loop through columns
      do col = 1, ncol
         do klev = 1, nz
            r(klev)     = 0.001_kind_phys * rho(col, klev)
            rhalf(klev) = sqrt(rho(col, 1) / rho(col, klev))
            pc(klev)    = 3.8_kind_phys / (pk(col, klev)**(1._kind_phys/xk)*psl)
            !
            ! if qr is (round-off) negative then the computation of
            ! velqr triggers floating point exception error when running
            ! in debugging mode with NAG
            !
            qr(col,klev) = MAX(qr(col,klev),0.0_kind_phys)
            !
            ! Liquid water terminal velocity (m/s) following KW eq. 2.15
            velqr(klev)  = 36.34_kind_phys * rhalf(klev) *                    &
                 (qr(col, klev) * r(klev))**0.1364_kind_phys
         end do

         ! Maximum time step size in accordance with CFL condition
         dt_max = dt
         do klev = 1, nz - 1
            ! NB: Original test for velqr /= 0 numerically unstable
            if (abs(velqr(klev)) > 1.0E-12_kind_phys) then
               dt_max = min(dt_max, 0.8_kind_phys*(z(col, klev+1) -           &
                    z(col, klev)) / velqr(klev))
            end if
         end do

         ! Number of subcycles
         rainsplit = ceiling(dt / dt_max)
         if (rainsplit < 1) then
            write(errmsg, *) 'KESSLER: bad rainsplit ',dt,dt_max,rainsplit
            errflg = 1
            return
         end if
         dt0 = dt / real(rainsplit, kind_phys)

         ! Subcycle through rain process
         nt = 1
         do while (nt <= rainsplit)

            ! Precipitation rate (m/s)
            precl(col) = precl(col) + rho(col, 1) * qr(col, 1) *  &
                 velqr(1) / rhoqr

            ! Sedimentation term using upstream differencing
            do klev = 1, nz-1
               sed(klev) = dt0 *                                              &
                    ((r(klev+1) * qr(col, klev+1) * velqr(klev+1)) -          &
                     (r(klev) * qr(col, klev) * velqr(klev))) /               &
                    (r(klev) * (z(col, klev+1) - z(col, klev)))
            end do
            sed(nz) = -dt0 * qr(col, nz) * velqr(nz) / &
                 (0.5_kind_phys * (z(col, nz)-z(col, nz-1)))

            ! Adjustment terms
            do klev = 1, nz

               ! Autoconversion and accretion rates following KW eq. 2.13a,b
               qrprod = qc(col, klev) - (qc(col, klev) - dt0 * &
                    max(.001_kind_phys * (qc(col, klev)-.001_kind_phys), &
                        0._kind_phys)) / &
                        (1._kind_phys + dt0 * 2.2_kind_phys * &
                         qr(col, klev)**.875_kind_phys)
               qc(col, klev) = max(qc(col, klev) - qrprod, 0._kind_phys)
               qr(col, klev) = max(qr(col, klev) + qrprod + sed(klev), &
                                   0._kind_phys)

               ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
               qvs = pc(klev) * exp(f2x*(pk(col, klev)*theta(col, klev) - 273._kind_phys) / (pk(col, klev)*theta(col, klev) - 36._kind_phys))
               prod = (qv(col, klev) - qvs) / (1._kind_phys + qvs*f5 / (pk(col, klev)*theta(col, klev) - 36._kind_phys)**2)


               ! Evaporation rate following KW eq. 2.14a,b
               ern = min(dt0 * (((1.6_kind_phys + 124.9_kind_phys*(r(klev)*qr(col, klev))**.2046_kind_phys) * &
                    (r(klev) * qr(col, klev))**.525_kind_phys) /                                &
                    (2550000._kind_phys * pc(klev) / (3.8_kind_phys*qvs) + 540000._kind_phys)) *           &
                    (dim(qvs,qv(col, klev)) / (r(klev)*qvs)),                            &
                    max(-prod-qc(col, klev),0._kind_phys),qr(col, klev))

               ! Saturation adjustment following KW eq. 3.10
               theta(col, klev)= theta(col, klev) + (lv / (cp * pk(col, klev)) * (max(prod,-qc(col, klev)) - ern))
               qv(col, klev) = max(qv(col, klev) - max(prod, -qc(col, klev)) + ern, 0._kind_phys)
               qc(col, klev) = qc(col, klev) + max(prod, -qc(col, klev))
               qr(col, klev) = qr(col, klev) - ern
            end do

            ! Recalculate liquid water terminal velocity
            if (nt /= rainsplit) then
               do klev = 1, nz
                  velqr(klev)  = 36.34_kind_phys * rhalf(klev) * (qr(col, klev)*r(klev))**0.1364_kind_phys
               end do
               !
               ! recompute rainsplit since velqr has changed
               !
               do klev = 1, nz - 1
                  if (abs(velqr(klev)) > 1.0E-12_kind_phys) then
                     dt_max = min(dt_max, 0.8_kind_phys*(z(col, klev+1) - z(col, klev)) / velqr(klev))
                  end if
               end do
               ! Number of subcycles
               rainsplit = ceiling(dt / dt_max)
               if (rainsplit < 1) then
                  write(errmsg, *) 'KESSLER: bad rainsplit ',dt,dt_max,rainsplit
                  errflg = 1
                  return
               end if
               dt0 = dt / real(rainsplit, kind_phys)
            end if
            nt=nt+1
         end do

         precl(col) = precl(col) / real(rainsplit, kind_phys)

      end do ! column loop

   end subroutine kessler_run

   !=======================================================================

end module kessler
