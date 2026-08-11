!----------------------------------------------------------------------------------
! low level utility module for cloud aerosols
!
! Created by Francis Vitt
!
! Portable (CCPP-ready): array sizes are runtime arguments and host constants
! are passed in; no CAM infrastructure dependencies.
!----------------------------------------------------------------------------------
module cldaero_mod

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: cldaero_uptakerate
  public :: cldaero_conc_t
  public :: cldaero_allocate
  public :: cldaero_deallocate

  type cldaero_conc_t
    real(kind_phys), pointer :: so4c(:, :)
    real(kind_phys), pointer :: nh4c(:, :)
    real(kind_phys), pointer :: no3c(:, :)
    real(kind_phys), pointer :: xlwc(:, :)
    real(kind_phys) :: so4_fact
  end type cldaero_conc_t

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  function cldaero_allocate(ncol, pver) result(cldconc)
    integer, intent(in) :: ncol  ! number of columns in chunk
    integer, intent(in) :: pver  ! number of vertical levels

    type(cldaero_conc_t), pointer :: cldconc

    allocate (cldconc)
    allocate (cldconc%so4c(ncol, pver))
    allocate (cldconc%nh4c(ncol, pver))
    allocate (cldconc%no3c(ncol, pver))
    allocate (cldconc%xlwc(ncol, pver))

    cldconc%so4c(:, :) = 0._kind_phys
    cldconc%nh4c(:, :) = 0._kind_phys
    cldconc%no3c(:, :) = 0._kind_phys
    cldconc%xlwc(:, :) = 0._kind_phys
    cldconc%so4_fact = 2._kind_phys

  end function cldaero_allocate

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  subroutine cldaero_deallocate(cldconc)
    type(cldaero_conc_t), pointer :: cldconc

    if (associated(cldconc%so4c)) then
      deallocate (cldconc%so4c)
      nullify (cldconc%so4c)
    end if

    if (associated(cldconc%nh4c)) then
      deallocate (cldconc%nh4c)
      nullify (cldconc%nh4c)
    end if

    if (associated(cldconc%no3c)) then
      deallocate (cldconc%no3c)
      nullify (cldconc%no3c)
    end if

    if (associated(cldconc%xlwc)) then
      deallocate (cldconc%xlwc)
      nullify (cldconc%xlwc)
    end if

    deallocate (cldconc)
    nullify (cldconc)

  end subroutine cldaero_deallocate

!----------------------------------------------------------------------------------
! utility function for cloud-borne aerosols
!----------------------------------------------------------------------------------

  function cldaero_uptakerate(xl, cldnum, cfact, cldfrc, tfld, press, pi) result(uptkrate)

    real(kind_phys), intent(in) :: xl, cldnum, cfact, cldfrc, tfld, press
    real(kind_phys), intent(in) :: pi   ! host value of pi (passed for bit-for-bit consistency)

    real(kind_phys) :: uptkrate

    real(kind_phys) :: &
      rad_cd, radxnum_cd, num_cd, &
      gasdiffus, gasspeed, knudsen, &
      fuchs_sutugin, volx34pi_cd

!-----------------------------------------------------------------------
! compute uptake of h2so4 and msa to cloud water
!
! first-order uptake rate is
! 4*pi*(drop radius)*(drop number conc)
! *(gas diffusivity)*(fuchs sutugin correction)

! num_cd = (drop number conc in 1/cm^3)
    num_cd = 1.0e-3_kind_phys*cldnum*cfact/cldfrc
    num_cd = max(num_cd, 0.0_kind_phys)

! rad_cd = (drop radius in cm), computed from liquid water and drop number,
! then bounded by 0.5 and 50.0 micrometers
! radxnum_cd = (drop radius)*(drop number conc)
! volx34pi_cd = (3/4*pi) * (liquid water volume in cm^3/cm^3)

    volx34pi_cd = xl*0.75_kind_phys/pi

! following holds because volx34pi_cd = num_cd*(rad_cd**3)
    radxnum_cd = (volx34pi_cd*num_cd*num_cd)**0.3333333_kind_phys

! apply bounds to rad_cd to avoid the occasional unphysical value
    if (radxnum_cd <= volx34pi_cd*4.0e4_kind_phys) then
      radxnum_cd = volx34pi_cd*4.0e4_kind_phys
      rad_cd = 50.0e-4_kind_phys
    else if (radxnum_cd >= volx34pi_cd*4.0e8_kind_phys) then
      radxnum_cd = volx34pi_cd*4.0e8_kind_phys
      rad_cd = 0.5e-4_kind_phys
    else
      rad_cd = radxnum_cd/num_cd
    end if

! gasdiffus = h2so4 gas diffusivity from mosaic code (cm^2/s)
! (pmid must be Pa)
    gasdiffus = 0.557_kind_phys*(tfld**1.75_kind_phys)/press

! gasspeed = h2so4 gas mean molecular speed from mosaic code (cm/s)
    gasspeed = 1.455e4_kind_phys*sqrt(tfld/98.0_kind_phys)

! knudsen number
    knudsen = 3.0_kind_phys*gasdiffus/(gasspeed*rad_cd)

! following assumes accomodation coefficient = 0.65
! (Adams & Seinfeld, 2002, JGR, and references therein)
! fuchs_sutugin = (0.75*accom*(1. + knudsen)) /
! (knudsen*(1.0 + knudsen + 0.283*accom) + 0.75*accom)
    fuchs_sutugin = (0.4875_kind_phys*(1._kind_phys + knudsen))/ &
                    (knudsen*(1.184_kind_phys + knudsen) + 0.4875_kind_phys)

! instantaneous uptake rate
    uptkrate = 12.56637_kind_phys*radxnum_cd*gasdiffus*fuchs_sutugin

  end function cldaero_uptakerate

end module cldaero_mod
