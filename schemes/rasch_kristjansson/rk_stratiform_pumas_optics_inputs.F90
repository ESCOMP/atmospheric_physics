! Derive the PUMAS-convention cloud optics inputs (dei, pgam, lamc)
! for RRTMGP gammadist liquid / Mitchell ice cloud optics
! from the RK (CAM4 cldefr) climatological effective radii rel/rei.
!
! RK microphysics carries no size distribution information, so the gamma
! distribution shape parameter pgam is fixed and the slope lamc is chosen to
! reproduce the RK liquid effective radius
!                           re = (pgam+3)/(2*lamc)
! Points without condensate are left at zero and are skipped by the optics.
module rk_stratiform_pumas_optics_inputs
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: rk_stratiform_pumas_optics_inputs_run

  ! Fixed gamma distribution shape parameter. The PUMAS (Rotstayn & Liu 2003)
  ! https://doi.org/10.1175/1520-0442(2003)016<3476:SOTFIA>2.0.CO;2
  ! pgam fit floors at 2 for droplet numbers > 60 cm-3 so we use this value
  ! as it is used for the optics tables by PUMAS for nearly all clouds.
  ! Optical properties at fixed effective radius are only weakly shape dependent.
  real(kind_phys), parameter :: pgam_fixed = 2._kind_phys

  ! Bulk ice densities used by PUMAS to convert effective radius to the
  ! generalized effective diameter the ice optics tables were built for
  ! (dei = rei * rhoi/rhows * 2, as in PUMAS post-interstitial).
  real(kind_phys), parameter :: rhoi  = 500._kind_phys ! bulk density ice [kg m-3]
  real(kind_phys), parameter :: rhows = 917._kind_phys ! bulk density water solid [kg m-3]

contains

!> \section arg_table_rk_stratiform_pumas_optics_inputs_run Argument Table
!! \htmlinclude arg_table_rk_stratiform_pumas_optics_inputs_run.html
  subroutine rk_stratiform_pumas_optics_inputs_run( &
    ncol, pver, &
    rel, rei, &
    iclwp, iciwp, &
    dei, pgam, lamc, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver

    real(kind_phys),    intent(in)  :: rel(:,:)   ! effective_radius_of_stratiform_cloud_liquid_water_droplet [um]
    real(kind_phys),    intent(in)  :: rei(:,:)   ! effective_radius_of_stratiform_cloud_ice_crystal [um]
    real(kind_phys),    intent(in)  :: iclwp(:,:) ! in_cloud_liquid_water_path_for_radiation [kg m-2]
    real(kind_phys),    intent(in)  :: iciwp(:,:) ! in_cloud_ice_water_path [kg m-2]

    ! Output arguments
    real(kind_phys),    intent(out) :: dei(:,:)   ! effective_diameter_of_stratiform_cloud_ice_crystal_for_radiation [um]
    real(kind_phys),    intent(out) :: pgam(:,:)  ! size_distribution_shape_parameter_for_microphysics [1]
    real(kind_phys),    intent(out) :: lamc(:,:)  ! slope_of_droplet_distribution_for_optics (in m-1) [1]
    character(len=*),   intent(out) :: errmsg     ! error message
    integer,            intent(out) :: errflg     ! error flag

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    dei(:ncol,:)  = 0._kind_phys
    pgam(:ncol,:) = 0._kind_phys
    lamc(:ncol,:) = 0._kind_phys

    do k = 1, pver
      do i = 1, ncol
        ! lamc > 0 flags cloud presence to the liquid optics.
        if (iclwp(i,k) > 0._kind_phys .and. rel(i,k) > 0._kind_phys) then
          pgam(i,k) = pgam_fixed
          lamc(i,k) = (pgam_fixed + 3._kind_phys)/(2._kind_phys*rel(i,k)*1.e-6_kind_phys)
        end if
        ! The ice optics skip points with dei == 0.
        if (iciwp(i,k) > 0._kind_phys) then
          dei(i,k) = rei(i,k)*rhoi/rhows*2._kind_phys
        end if
      end do
    end do

  end subroutine rk_stratiform_pumas_optics_inputs_run

end module rk_stratiform_pumas_optics_inputs
