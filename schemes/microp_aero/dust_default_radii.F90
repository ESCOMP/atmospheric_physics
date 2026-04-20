! Set default dust radii for 4 size bins used in contact freezing.
! For modal aerosols, bin 3 is overwritten later with wet diameter.
module dust_default_radii

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: dust_default_radii_timestep_init

  real(kind_phys), parameter :: rn_dst1 = 0.258e-6_kind_phys
  real(kind_phys), parameter :: rn_dst2 = 0.717e-6_kind_phys
  real(kind_phys), parameter :: rn_dst3 = 1.576e-6_kind_phys
  real(kind_phys), parameter :: rn_dst4 = 3.026e-6_kind_phys

contains

!> \section arg_table_dust_default_radii_timestep_init Argument Table
!! \htmlinclude dust_default_radii_timestep_init.html
  subroutine dust_default_radii_timestep_init(ncol, pver, ndust, rndst, errmsg, errflg)
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    integer,            intent(in)  :: ndust

    real(kind_phys),    intent(out) :: rndst(:,:,:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    rndst(:ncol, :pver, 1) = rn_dst1
    rndst(:ncol, :pver, 2) = rn_dst2
    rndst(:ncol, :pver, 3) = rn_dst3
    rndst(:ncol, :pver, 4) = rn_dst4

  end subroutine dust_default_radii_timestep_init

end module dust_default_radii
