!===============================================================================
! used to compute sea salt surface emissions for modal and sectional aerosol models
!===============================================================================
module sslt_sections
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: sslt_sections_init
  public :: fluxes
  public :: nsections
  public :: Dg
  public :: rdry

  integer, parameter :: nsections = 31

  ! only use up to ~20um
  real(kind_phys), parameter :: Dg(nsections) = [ &
                                0.0020e-5_kind_phys, 0.0025e-5_kind_phys, 0.0032e-5_kind_phys, &
                                0.0040e-5_kind_phys, 0.0051e-5_kind_phys, 0.0065e-5_kind_phys, &
                                0.0082e-5_kind_phys, 0.0104e-5_kind_phys, 0.0132e-5_kind_phys, &
                                0.0167e-5_kind_phys, 0.0211e-5_kind_phys, 0.0267e-5_kind_phys, &
                                0.0338e-5_kind_phys, 0.0428e-5_kind_phys, 0.0541e-5_kind_phys, &
                                0.0685e-5_kind_phys, 0.0867e-5_kind_phys, 0.1098e-5_kind_phys, &
                                0.1389e-5_kind_phys, 0.1759e-5_kind_phys, 0.2226e-5_kind_phys, &
                                0.2818e-5_kind_phys, 0.3571e-5_kind_phys, 0.4526e-5_kind_phys, &
                                0.5735e-5_kind_phys, 0.7267e-5_kind_phys, 0.9208e-5_kind_phys, &
                                1.1668e-5_kind_phys, 1.4786e-5_kind_phys, 1.8736e-5_kind_phys, &
                                2.3742e-5_kind_phys]

  real(kind_phys), dimension(nsections) :: bm, rdry, rm
  real(kind_phys), dimension(4, nsections) :: consta, constb  !constants for calculating emission polynomial

contains

  subroutine sslt_sections_init()

    integer :: m

    ! use Ekman's ss
    rdry(:) = Dg(:)/2._kind_phys   ! meter
    ! multiply rm with 1.814 because it should be RH=80% and not dry particles
    ! for the parameterization
    rm(:) = 1.814_kind_phys*rdry(:)*1.e6_kind_phys   ! um
    bm(:) = (0.380_kind_phys - log10(rm(:)))/0.65_kind_phys  ! use in Manahan

    ! calculate constants form emission polynomials
    do m = 1, nsections
      if ((m) <= 9) then
        consta(1, m) = (-2.576_kind_phys)*10._kind_phys**35*Dg(m)**4 + 5.932_kind_phys*10._kind_phys**28 &
                       *Dg(m)**3 + (-2.867_kind_phys)*10._kind_phys**21*Dg(m)**2 + (-3.003_kind_phys) &
                       *10._kind_phys**13*Dg(m) + (-2.881_kind_phys)*10._kind_phys**6
        constb(1, m) = 7.188_kind_phys*10._kind_phys**37 &
                       *Dg(m)**4 + (-1.616_kind_phys)*10._kind_phys**31*Dg(m)**3 + 6.791_kind_phys*10._kind_phys**23 &
                       *Dg(m)**2 + 1.829_kind_phys*10._kind_phys**16*Dg(m) + 7.609_kind_phys*10._kind_phys**8
      else if ((m) >= 10 .and. (m) <= 13) then
        consta(2, m) = (-2.452_kind_phys)*10._kind_phys**33*Dg(m)**4 + 2.404_kind_phys*10._kind_phys**27 &
                       *Dg(m)**3 + (-8.148_kind_phys)*10._kind_phys**20*Dg(m)**2 + (1.183_kind_phys)*10._kind_phys**14 &
                       *Dg(m) + (-6.743_kind_phys)*10._kind_phys**6
        constb(2, m) = 7.368_kind_phys*10._kind_phys**35 &
                       *Dg(m)**4 + (-7.310_kind_phys)*10._kind_phys**29*Dg(m)**3 + 2.528_kind_phys*10._kind_phys**23 &
                       *Dg(m)**2 + (-3.787_kind_phys)*10._kind_phys**16*Dg(m) + 2.279_kind_phys*10._kind_phys**9
      else if ((m) >= 14 .and. (m) < 22) then
        consta(3, m) = (1.085_kind_phys)*10._kind_phys**29*Dg(m)**4 + (-9.841_kind_phys)*10._kind_phys**23 &
                       *Dg(m)**3 + (3.132_kind_phys)*10._kind_phys**18*Dg(m)**2 + (-4.165_kind_phys)*10._kind_phys**12 &
                       *Dg(m) + (2.181_kind_phys)*10._kind_phys**6
        constb(3, m) = (-2.859_kind_phys)*10._kind_phys**31 &
                       *Dg(m)**4 + (2.601_kind_phys)*10._kind_phys**26*Dg(m)**3 + (-8.297_kind_phys)*10._kind_phys**20 &
                       *Dg(m)**2 + (1.105_kind_phys)*10._kind_phys**15*Dg(m) + (-5.800_kind_phys)*10._kind_phys**8
      else if (m >= 22 .and. m <= 40) then
        ! use monahan
        consta(4, m) = (1.373_kind_phys*rm(m)**(-3)*(1 + 0.057_kind_phys*rm(m)**1.05_kind_phys) &
                        *10**(1.19_kind_phys*exp(-bm(m)**2))) &
                       *(rm(m) - rm(m - 1))
      end if
    end do
  end subroutine sslt_sections_init

  function fluxes(sst, u10cubed, ncol) result(fi)

    real(kind_phys), intent(in) :: sst(:)
    real(kind_phys), intent(in) :: u10cubed(:)
    integer, intent(in) :: ncol

    real(kind_phys) :: fi(ncol, nsections)

    integer :: m
    real(kind_phys) :: W(ncol)

    ! Calculations of source strength and size distribution
    ! NB the 0.1 is the dlogDp we have to multiplie with to get the flux, but the value dependence
    ! of course on what dlogDp you have. You will also have to change the sections of Dg if you use
    ! a different number of size bins with different intervals.

    W(:ncol) = 3.84e-6_kind_phys*u10cubed(:ncol)*0.1_kind_phys ! whitecap area

    ! calculate number flux fi (#/m2/s)
    fi(:, :) = 0._kind_phys
    do m = 1, nsections
      if (m <= 9) then
        fi(:ncol, m) = W(:ncol)*((sst(:ncol))*consta(1, m) + constb(1, m))
      else if (m >= 10 .and. m <= 13) then
        fi(:ncol, m) = W(:ncol)*((sst(:ncol))*consta(2, m) + constb(2, m))
      else if (m >= 14 .and. m < 22) then
        fi(:ncol, m) = W(:ncol)*((sst(:ncol))*consta(3, m) + constb(3, m))
      else if (m >= 22 .and. m <= 40) then
        ! use Monahan
        fi(:ncol, m) = consta(4, m)*u10cubed(:ncol)
      end if
    end do

  end function fluxes

end module sslt_sections
