! Contains routines to compute tendencies from sedimentation of dust
! Author: Phil Rasch
module dust_sediment_mod
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: dust_sediment_tend

  real(kind_phys), parameter :: mxsedfac = 0.99_kind_phys       ! maximum sedimentation flux factor

contains

  ! Apply Particle Gravitational Sedimentation
  subroutine dust_sediment_tend( &
    ncol, dtime, pint, pdel, &
    dustmr, pvdust, dusttend, sfdust, &
    pver, gravit, errmsg, errflg)

    integer, intent(in)  :: ncol
    real(kind_phys), intent(in)  :: dtime                     ! time step
    real(kind_phys), intent(in)  :: pint(:, :)               ! interfaces pressure (Pa)
    real(kind_phys), intent(in)  :: pdel(:, :)               ! pressure diff across layer (Pa)
    real(kind_phys), intent(in)  :: dustmr(:, :)               ! dust (kg/kg)
    real(kind_phys), intent(in)  :: pvdust(:, :)              ! vertical velocity of dust drops  (Pa/s)
    ! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))
    real(kind_phys), intent(out) :: dusttend(:, :)             ! dust tend
    real(kind_phys), intent(out) :: sfdust(:)               ! surface flux of dust (rain, kg/m/s)

    integer, intent(in)  :: pver                      ! number of vertical levels
    real(kind_phys), intent(in)  :: gravit                    ! gravitational acceleration (m/s2)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Local variables
    real(kind_phys) :: fxdust(ncol, pver + 1)                     ! fluxes at the interfaces, dust (positive = down)
    integer :: i, k

    errmsg = ''
    errflg = 0

    ! initialize variables
    fxdust(:ncol, :) = 0._kind_phys ! flux at interfaces (dust)
    dusttend(:ncol, :) = 0._kind_phys ! tend (dust)
    sfdust(:ncol) = 0._kind_phys ! sedimentation flux out bot of column (dust)

    ! fluxes at interior points
    call getflx(ncol, pint, dustmr, pvdust, dtime, fxdust, pver, errmsg, errflg)
    if (errflg /= 0) return

    ! calculate fluxes at boundaries
    do i = 1, ncol
      fxdust(i, 1) = 0
      ! surface flux by upstream scheme
      fxdust(i, pver + 1) = dustmr(i, pver)*pvdust(i, pver + 1)*dtime
    end do

    ! filter out any negative fluxes from the getflx routine
    do k = 2, pver
      fxdust(:ncol, k) = max(0._kind_phys, fxdust(:ncol, k))
    end do

    ! Limit the flux out of the bottom of each cell to the water content in each phase.
    ! Apply mxsedfac to prevent generating very small negative cloud water/ice
    ! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
    ! ***Should we include the flux in the top, to allow for thin surface layers?
    ! ***Requires simple treatment of cloud overlap, already included below.
    do k = 1, pver
      do i = 1, ncol
        fxdust(i, k + 1) = min(fxdust(i, k + 1), mxsedfac*dustmr(i, k)*pdel(i, k))
      end do
    end do

    ! Now calculate the tendencies
    do k = 1, pver
      do i = 1, ncol
        ! net flux into cloud changes cloud dust/ice (all flux is out of cloud)
        dusttend(i, k) = (fxdust(i, k) - fxdust(i, k + 1))/(dtime*pdel(i, k))
      end do
    end do

    ! convert flux out the bottom to mass units Pa -> kg/m2/s
    sfdust(:ncol) = fxdust(:ncol, pver + 1)/(dtime*gravit)

  end subroutine dust_sediment_tend


  subroutine getflx(ncol, xw, phi, vel, deltat, flux, pver, errmsg, errflg)
    !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
    !....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    !.........phi1......phi2.......phi3.....phi4.......phi5.......
    integer :: ncol                      ! number of colums to process

    integer :: i
    integer :: k

    real(kind_phys) :: vel(:, :)
    real(kind_phys) :: flux(:, :)
    real(kind_phys) :: xw(:, :)
    real(kind_phys) :: phi(:, :)
    real(kind_phys) :: deltat

    integer, intent(in) :: pver       ! number of vertical levels
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    real(kind_phys) :: psi(ncol, pver + 1)
    real(kind_phys) :: fdot(ncol, pver + 1)
    real(kind_phys) :: xx(ncol)
    real(kind_phys) :: fxdot(ncol)
    real(kind_phys) :: fxdd(ncol)

    real(kind_phys) :: psistar(ncol)

    real(kind_phys) :: xxk(ncol, pver)

    errmsg = ''
    errflg = 0

    do i = 1, ncol
!        integral of phi
      psi(i, 1) = 0._kind_phys
!        fluxes at boundaries
      flux(i, 1) = 0
      flux(i, pver + 1) = 0._kind_phys
    end do

!     integral function
    do k = 2, pver + 1
      do i = 1, ncol
        psi(i, k) = phi(i, k - 1)*(xw(i, k) - xw(i, k - 1)) + psi(i, k - 1)
      end do
    end do

!     calculate the derivatives for the interpolating polynomial
    call cfdotmc_pro(ncol, xw, psi, fdot, pver)

!  NEW WAY
!     calculate fluxes at interior pts
    do k = 2, pver
      do i = 1, ncol
        xxk(i, k) = xw(i, k) - vel(i, k)*deltat
      end do
    end do
    do k = 2, pver
      call cfint2(ncol, xw, psi, fdot, xxk(:, k), fxdot, fxdd, psistar, pver, errmsg, errflg)
      if (errflg /= 0) return
      do i = 1, ncol
        flux(i, k) = (psi(i, k) - psistar(i))
      end do
    end do

  end subroutine getflx


  subroutine cfint2(ncol, x, f, fdot, xin, fxdot, fxdd, psistar, pver, errmsg, errflg)
    ! input
    integer :: ncol                      ! number of colums to process

    real(kind_phys) :: x(:, :)
    real(kind_phys) :: f(:, :)
    real(kind_phys) :: fdot(:, :)
    real(kind_phys) :: xin(:)

    integer, intent(in) :: pver       ! number of vertical levels

    ! output
    real(kind_phys) :: fxdot(:)
    real(kind_phys) :: fxdd(:)
    real(kind_phys) :: psistar(:)

    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: i
    integer :: k
    integer :: intz(ncol)
    real(kind_phys) :: dx
    real(kind_phys) :: s
    real(kind_phys) :: c2
    real(kind_phys) :: c3
    real(kind_phys) :: xx
    real(kind_phys) :: xinf
    real(kind_phys) :: psi1, psi2, psi3, psim
    real(kind_phys) :: cfint
    real(kind_phys) :: cfnew
    real(kind_phys) :: xins(ncol)

!     the minmod function
    real(kind_phys) :: a, b, c
    real(kind_phys) :: minmod
    real(kind_phys) :: medan
    minmod(a, b) = 0.5_kind_phys*(sign(1._kind_phys, a) + sign(1._kind_phys, b))*min(abs(a), abs(b))
    medan(a, b, c) = a + minmod(b - a, c - a)

    errmsg = ''
    errflg = 0

    do i = 1, ncol
      xins(i) = medan(x(i, 1), xin(i), x(i, pver + 1))
      intz(i) = 0
    end do

! first find the interval
    do k = 1, pver
      do i = 1, ncol
        if ((xins(i) - x(i, k))*(x(i, k + 1) - xins(i)) >= 0._kind_phys) then
          intz(i) = k
        end if
      end do
    end do

    do i = 1, ncol
      if (intz(i) == 0) then
        write (errmsg, *) 'DUST_SEDIMENT_MOD:cfint2 -- interval was not found for col i ', i
        errflg = 1
        return
      end if
    end do

! now interpolate
    do i = 1, ncol
      k = intz(i)
      dx = (x(i, k + 1) - x(i, k))
      s = (f(i, k + 1) - f(i, k))/dx
      c2 = (3*s - 2*fdot(i, k) - fdot(i, k + 1))/dx
      c3 = (fdot(i, k) + fdot(i, k + 1) - 2*s)/dx**2
      xx = (xins(i) - x(i, k))
      fxdot(i) = (3*c3*xx + 2*c2)*xx + fdot(i, k)
      fxdd(i) = 6*c3*xx + 2*c2
      cfint = ((c3*xx + c2)*xx + fdot(i, k))*xx + f(i, k)

!        limit the interpolant
      psi1 = f(i, k) + (f(i, k + 1) - f(i, k))*xx/dx
      if (k == 1) then
        psi2 = f(i, 1)
      else
        psi2 = f(i, k) + (f(i, k) - f(i, k - 1))*xx/(x(i, k) - x(i, k - 1))
      end if
      if (k + 1 == pver + 1) then
        psi3 = f(i, pver + 1)
      else
        psi3 = f(i, k + 1) - (f(i, k + 2) - f(i, k + 1))*(dx - xx)/(x(i, k + 2) - x(i, k + 1))
      end if
      psim = medan(psi1, psi2, psi3)
      cfnew = medan(cfint, psi1, psim)
      if (abs(cfnew - cfint)/(abs(cfnew) + abs(cfint) + 1.e-36_kind_phys) > .03_kind_phys) then
!     CHANGE THIS BACK LATER!!!
!     $        .gt..1) then

!     UNCOMMENT THIS LATER!!!
!            write(iulog,*) ' cfint2 limiting important ', cfint, cfnew

      end if
      psistar(i) = cfnew
    end do

  end subroutine cfint2

  subroutine cfdotmc_pro(ncol, x, f, fdot, pver)

    !     prototype version; eventually replace with final SPITFIRE scheme

    !     calculate the derivative for the interpolating polynomial
    !     multi column version

    ! input
    integer :: ncol                      ! number of colums to process

    real(kind_phys) :: x(:, :)
    real(kind_phys) :: f(:, :)
    integer, intent(in) :: pver       ! number of vertical levels
    ! output
    real(kind_phys) :: fdot(:, :)          ! derivative at nodes

    ! assumed variable distribution
    !     x1.......x2.......x3.......x4.......x5.......x6     1,pverp points
    !     f1.......f2.......f3.......f4.......f5.......f6     1,pverp points
    !     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
    !     .........d2.......d3.......d4.......d5.........     2,pver points
    !     .........s2.......s3.......s4.......s5.........     2,pver points
    !     .............dh2......dh3......dh4.............     2,pver-1 points
    !     .............eh2......eh3......eh4.............     2,pver-1 points
    !     ..................e3.......e4..................     3,pver-1 points
    !     .................ppl3......ppl4................     3,pver-1 points
    !     .................ppr3......ppr4................     3,pver-1 points
    !     .................t3........t4..................     3,pver-1 points
    !     ................fdot3.....fdot4................     3,pver-1 points

! work variables

    integer :: i
    integer :: k

    real(kind_phys) :: a                    ! work var
    real(kind_phys) :: b                    ! work var
    real(kind_phys) :: c                    ! work var
    real(kind_phys) :: s(ncol, pver + 1)             ! first divided differences at nodes
    real(kind_phys) :: sh(ncol, pver + 1)            ! first divided differences between nodes
    real(kind_phys) :: d(ncol, pver + 1)             ! second divided differences at nodes
    real(kind_phys) :: dh(ncol, pver + 1)            ! second divided differences between nodes
    real(kind_phys) :: e(ncol, pver + 1)             ! third divided differences at nodes
    real(kind_phys) :: eh(ncol, pver + 1)            ! third divided differences between nodes
    real(kind_phys) :: pp                   ! p prime
    real(kind_phys) :: ppl(ncol, pver + 1)           ! p prime on left
    real(kind_phys) :: ppr(ncol, pver + 1)           ! p prime on right
    real(kind_phys) :: qpl
    real(kind_phys) :: qpr
    real(kind_phys) :: ttt
    real(kind_phys) :: t
    real(kind_phys) :: tmin
    real(kind_phys) :: tmax
    real(kind_phys) :: delxh(ncol, pver + 1)

!     the minmod function
    real(kind_phys) :: minmod
    real(kind_phys) :: medan
    minmod(a, b) = 0.5_kind_phys*(sign(1._kind_phys, a) + sign(1._kind_phys, b))*min(abs(a), abs(b))
    medan(a, b, c) = a + minmod(b - a, c - a)

    do k = 1, pver

!        first divided differences between nodes
      do i = 1, ncol
        delxh(i, k) = (x(i, k + 1) - x(i, k))
        sh(i, k) = (f(i, k + 1) - f(i, k))/delxh(i, k)
      end do

!        first and second divided differences at nodes
      if (k >= 2) then
        do i = 1, ncol
          d(i, k) = (sh(i, k) - sh(i, k - 1))/(x(i, k + 1) - x(i, k - 1))
          s(i, k) = minmod(sh(i, k), sh(i, k - 1))
        end do
      end if
    end do

!     second and third divided diffs between nodes
    do k = 2, pver - 1
      do i = 1, ncol
        eh(i, k) = (d(i, k + 1) - d(i, k))/(x(i, k + 2) - x(i, k - 1))
        dh(i, k) = minmod(d(i, k), d(i, k + 1))
      end do
    end do

!     treat the boundaries
    do i = 1, ncol
      e(i, 2) = eh(i, 2)
      e(i, pver) = eh(i, pver - 1)
!        outside level
      fdot(i, 1) = sh(i, 1) - d(i, 2)*delxh(i, 1) &
                   - eh(i, 2)*delxh(i, 1)*(x(i, 1) - x(i, 3))
      fdot(i, 1) = minmod(fdot(i, 1), 3*sh(i, 1))
      fdot(i, pver + 1) = sh(i, pver) + d(i, pver)*delxh(i, pver) &
                          + eh(i, pver - 1)*delxh(i, pver)*(x(i, pver + 1) - x(i, pver - 1))
      fdot(i, pver + 1) = minmod(fdot(i, pver + 1), 3*sh(i, pver))
!        one in from boundary
      fdot(i, 2) = sh(i, 1) + d(i, 2)*delxh(i, 1) - eh(i, 2)*delxh(i, 1)*delxh(i, 2)
      fdot(i, 2) = minmod(fdot(i, 2), 3*s(i, 2))
      fdot(i, pver) = sh(i, pver) - d(i, pver)*delxh(i, pver) &
                      - eh(i, pver - 1)*delxh(i, pver)*delxh(i, pver - 1)
      fdot(i, pver) = minmod(fdot(i, pver), 3*s(i, pver))
    end do

    do k = 3, pver - 1
      do i = 1, ncol
        e(i, k) = minmod(eh(i, k), eh(i, k - 1))
      end do
    end do

    do k = 3, pver - 1

      do i = 1, ncol

!           p prime at k-0.5
        ppl(i, k) = sh(i, k - 1) + dh(i, k - 1)*delxh(i, k - 1)
!           p prime at k+0.5
        ppr(i, k) = sh(i, k) - dh(i, k)*delxh(i, k)

        t = minmod(ppl(i, k), ppr(i, k))

!           derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
        pp = sh(i, k - 1) + d(i, k)*delxh(i, k - 1)

!           quartic estimate of fdot
        fdot(i, k) = pp &
                     - delxh(i, k - 1)*delxh(i, k) &
                     *(eh(i, k - 1)*(x(i, k + 2) - x(i, k)) &
                       + eh(i, k)*(x(i, k) - x(i, k - 2)) &
                       )/(x(i, k + 2) - x(i, k - 2))

!           now limit it
        qpl = sh(i, k - 1) &
              + delxh(i, k - 1)*minmod(d(i, k - 1) + e(i, k - 1)*(x(i, k) - x(i, k - 2)), &
                                       d(i, k) - e(i, k)*delxh(i, k))
        qpr = sh(i, k) &
              + delxh(i, k)*minmod(d(i, k) + e(i, k)*delxh(i, k - 1), &
                                   d(i, k + 1) + e(i, k + 1)*(x(i, k) - x(i, k + 2)))

        fdot(i, k) = medan(fdot(i, k), qpl, qpr)

        ttt = minmod(qpl, qpr)
        tmin = min(0._kind_phys, 3*s(i, k), 1.5_kind_phys*t, ttt)
        tmax = max(0._kind_phys, 3*s(i, k), 1.5_kind_phys*t, ttt)

        fdot(i, k) = fdot(i, k) + minmod(tmin - fdot(i, k), tmax - fdot(i, k))

      end do

    end do

  end subroutine cfdotmc_pro
end module dust_sediment_mod
