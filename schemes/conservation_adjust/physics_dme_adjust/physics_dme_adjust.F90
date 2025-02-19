module physics_dme_adjust

   use ccpp_kinds, only:  kind_phys

   implicit none

contains
!===============================================================================
!> \section arg_table_physics_dme_adjust_run Argument Table
!! \htmlinclude physics_dme_adjust_run.html
!!
  subroutine physics_dme_adjust_run(ncol, pver, pcnst, ps, pint, pdel, lnpint, rpdel, const_props, const_array, qini, liqini, iceini, &
                            errmsg, errflg)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Adjust the dry mass in each layer back to the value of physics input state
    !
    ! Method: Conserve the integrated mass, momentum and total energy in each layer
    !         by scaling the specific mass of consituents, specific momentum (velocity)
    !         and specific total energy by the relative change in layer mass. Solve for
    !         the new temperature by subtracting the new kinetic energy from total energy
    !         and inverting the hydrostatic equation
    !
    !         The mass in each layer is modified, changing the relationship of the layer
    !         interfaces and midpoints to the surface pressure. The result is no longer in
    !         the original hybrid coordinate.
    !
    !         This procedure cannot be applied to the "eul" or "sld" dycores because they
    !         require the hybrid coordinate.
    !
    ! Author: Byron Boville

    ! !REVISION HISTORY:
    !   03.03.28  Boville    Created, partly from code by Lin in p_d_adjust
    !
    !-----------------------------------------------------------------------

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    implicit none
    !
    ! Arguments
    !
    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: pcnst
    real(kind_phys),                   intent(inout) :: ps(:)
    real(kind_phys),                   intent(inout) :: pint(:,:)
    real(kind_phys),                   intent(inout) :: pdel(:,:)
    real(kind_phys),                   intent(inout) :: lnpint(:,:)
    real(kind_phys),                   intent(inout) :: rpdel(:,:)
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)
    real(kind_phys),                   intent(inout) :: const_array(:,:,:)

    real(kind_phys),                   intent(in)    :: qini(:,:)    ! initial specific humidity
    real(kind_phys),                   intent(in)    :: liqini(:,:)  ! initial total liquid
    real(kind_phys),                   intent(in)    :: iceini(:,:)  ! initial total ice
    character(len=512),                intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: k,m           ! Longitude, level indices
    real(kind_phys) :: fdq(ncol)    ! mass adjustment factor

    real(kind_phys) :: tot_water (ncol,2)  ! total water (initial, present)
    real(kind_phys) :: tot_water_chg(ncol) ! total water change


    integer :: m_cnst
    logical :: water_flag

    errmsg = ' '
    errflg = 0

    !
    !-----------------------------------------------------------------------
    ! adjust dry mass in each layer back to input value, while conserving
    ! constituents, momentum, and total energy
    ps(:) = pint(:,1)

    do k = 1, pver
      tot_water(:,1) = qini(:,k) +liqini(:,k)+iceini(:,k) !initial total H2O
      tot_water(:,2) = 0.0_kind_phys
      do m_cnst=1,pcnst
        call const_props(m_cnst)%is_water_species(water_flag)
        if(water_flag) then
          tot_water(:,2) = tot_water(:,2)+const_array(:,k,m_cnst)
        end if
      end do
      fdq(:) = 1._kind_phys + tot_water(:,2) - tot_water(:,1)
      ! adjust constituents to conserve mass in each layer
      do m = 1, pcnst
        const_array(:,k,m) = const_array(:,k,m) / fdq(:)
      end do
      ! compute new total pressure variables
      pdel  (:,k  ) = pdel(:,k  ) * fdq(:)
      ps(:)         = ps(:)       + pdel(:,k)
      pint  (:,k+1) = pint(:,k  ) + pdel(:,k)
      lnpint(:,k+1) = log(pint(:,k+1))
      rpdel (:,k  ) = 1._kind_phys/ pdel(:,k  )
      !note that mid-level variables (e.g. pmid) are not recomputed
    end do

  end subroutine physics_dme_adjust_run

end module physics_dme_adjust
