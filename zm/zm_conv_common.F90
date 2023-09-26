module zm_conv_common

  use ccpp_kinds, only:  kind_phys


  implicit none

  save
!  CACNote ------------  private                         ! Make default type private to the module?
  public


!
! Private data
!
   real(kind_phys) rl         ! wg latent heat of vaporization.
   real(kind_phys) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(kind_phys) :: capelmt ! namelist configurable:
                       ! threshold value for cape for deep convection.
   real(kind_phys) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
   real(kind_phys) :: ke_lnd
   real(kind_phys) :: c0_lnd       ! set from namelist input zmconv_c0_lnd
   real(kind_phys) :: c0_ocn       ! set from namelist input zmconv_c0_ocn
   integer  :: num_cin      ! set from namelist input zmconv_num_cin
                            ! The number of negative buoyancy regions that are allowed
                            ! before the convection top and CAPE calculations are completed.
   logical  :: zm_org
   real(kind_phys) tau   ! convective time scale
! CACNOTE - unused?   real(kind_phys),parameter :: c1 = 6.112_kind_phys
! CACNOTE - unused?   real(kind_phys),parameter :: c2 = 17.67_kind_phys
! CACNOTE - unused?   real(kind_phys),parameter :: c3 = 243.5_kind_phys
   real(kind_phys) :: tfreez
   real(kind_phys) :: eps1
   real(kind_phys) :: momcu
   real(kind_phys) :: momcd

   logical :: no_deep_pbl ! default = .false.
                          ! no_deep_pbl = .true. eliminates deep convection entirely within PBL


!moved from moistconvection.F90
   real(kind_phys) :: rgrav       ! reciprocal of grav
   real(kind_phys) :: rgas        ! gas constant for dry air
   real(kind_phys) :: grav        ! = gravit
   real(kind_phys) :: cp          ! = cpres = cpair

   integer  limcnv       ! top interface level limit for convection

   logical :: lparcel_pbl     ! Switch to turn on mixing of parcel MSE air, and picking launch level to be the top of the PBL.


   real(kind_phys) :: tiedke_add      ! namelist configurable
   real(kind_phys) :: dmpdz_param     ! namelist configurable

end module zm_conv_common
