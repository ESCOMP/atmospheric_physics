module geopotential

!---------------------------------------------------------------------------------
! Compute geopotential from temperature or
! compute geopotential and temperature from dry static energy.
!
! The hydrostatic matrix elements must be consistent with the dynamics algorithm.
! The diagonal element is the itegration weight from interface k minus vert_step to midpoint k.
! The offdiagonal element is the weight between interfaces.
! 
! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
!---------------------------------------------------------------------------------

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  public geopotential_dse
  public geopotential_t

contains
!===============================================================================
  subroutine geopotential_dse(pver, pverp, fvdyn,  vert_surf, vert_toa,  &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  ,  &
       dse    , q      , phis   , rair   , gravit , cpair  ,  &
       zvir   , t      , zi     , zm     , ncol             )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the temperature  and geopotential height (above the surface) at the
! midpoints and interfaces from the input dry static energy and pressures.
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!
! Input arguments
    integer, intent(in) :: ncol                  ! Number of longitudes

    ! rair, and cpair are passed in as slices of rank 3 arrays allocated
    ! at runtime. Don't specify size to avoid temporary copy.
    integer,         intent(in) :: pver
    integer,         intent(in) :: pverp
    integer,         intent(in) :: vert_surf
    integer,         intent(in) :: vert_toa
    logical,         intent(in) :: fvdyn         ! finite volume dynamics
    real(kind_phys), intent(in) :: piln (:,:)    ! (pcols,pverp) - Log interface pressures
    real(kind_phys), intent(in) :: pmln (:,:)    ! (pcols,pver)  - Log midpoint pressures
    real(kind_phys), intent(in) :: pint (:,:)    ! (pcols,pverp) - Interface pressures
    real(kind_phys), intent(in) :: pmid (:,:)    ! (pcols,pver)  - Midpoint pressures
    real(kind_phys), intent(in) :: pdel (:,:)    ! (pcols,pver)  - layer thickness
    real(kind_phys), intent(in) :: rpdel(:,:)    ! (pcols,pver)  - inverse of layer thickness
    real(kind_phys), intent(in) :: dse  (:,:)    ! (pcols,pver)  - dry static energy
    real(kind_phys), intent(in) :: q    (:,:)    ! (pcols,pver)  - specific humidity
    real(kind_phys), intent(in) :: phis (:)      ! (pcols)       - surface geopotential
    real(kind_phys), intent(in) :: rair (:,:)    !               - Gas constant for dry air
    real(kind_phys), intent(in) :: gravit        !               - Acceleration of gravity
    real(kind_phys), intent(in) :: cpair(:,:)    !               - specific heat at constant p for dry air
    real(kind_phys), intent(in) :: zvir (:,:)    ! (pcols,pver)  - rh2o/rair - 1

! Output arguments

    real(kind_phys), intent(out) :: t(:,:)       ! (pcols,pver)  - temperature
    real(kind_phys), intent(out) :: zi(:,:)      ! (pcols,pverp) - Height above surface at interfaces
    real(kind_phys), intent(out) :: zm(:,:)      ! (pcols,pver)  - Geopotential height at mid level
!
!---------------------------Local variables-----------------------------------------
!
    integer  :: i,k,vert_step,indx_surf        ! Lon, level, level indices, vertical_step, index of surface level
    integer  :: indx_level, indx_above         ! Index for level, index for above
    real(kind_phys) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(kind_phys) :: hkl(ncol)               ! off-diagonal element
    real(kind_phys) :: rog(ncol,pver)          ! Rair / gravit
    real(kind_phys) :: tv                      ! virtual temperature
    real(kind_phys) :: tvfac                   ! Tv/T
!
!----------------------------------------------------------------------------------
    rog(:ncol,:) = rair(:ncol,:) / gravit

    if (vert_surf > vert_toa) then 
       vert_step = -1
       indx_surf = vert_toa +1
       indx_level = 1
       indx_above = 0
    else
       vert_step = 1
       indx_surf = vert_surf
       indx_level = 0
       indx_above = 1
    end if

! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,indx_surf) = 0.0_kind_phys
    end do

! Compute the virtual temperature, zi, zm from bottom up
! Note, zi(i,k) is the interface above zm(i,k) when ordered top to bot
    do k = vert_surf, vert_toa, vert_step

! First set hydrostatic elements consistent with dynamics
       if (fvdyn) then
          do i = 1,ncol
             hkl(i) = piln(i,k+indx_level) - piln(i,k+indx_above)
             hkk(i) = 1._kind_phys - pint(i,k+indx_above) * hkl(i) * rpdel(i,k)
          end do
       else
          do i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_kind_phys * hkl(i)
          end do
       end if

! Now compute tv, t, zm, zi
       do i = 1,ncol
          tvfac   = 1._kind_phys + zvir(i,k) * q(i,k)
          tv      = (dse(i,k) - phis(i) - gravit*zi(i,k+indx_level)) / ((cpair(i,k) / tvfac) + &
	                                                               rair(i,k)*hkk(i))

          t (i,k) = tv / tvfac

          zm(i,k)            = zi(i,k+indx_level) + rog(i,k) * tv * hkk(i)
          zi(i,k+indx_above) = zi(i,k+indx_level) + rog(i,k) * tv * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_dse

!===============================================================================
  subroutine geopotential_t(pver, pverp, fvdyn, vert_surf, vert_toa,              &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the geopotential height (above the surface) at the midpoints and 
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------

use ppgrid, only : pcols

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer,         intent(in) :: pver
    integer,         intent(in) :: pverp
    integer,         intent(in) :: vert_surf
    integer,         intent(in) :: vert_toa
    logical,         intent(in) :: fvdyn         ! finite volume dynamics
    integer, intent(in) :: ncol                  ! Number of longitudes

    real(kind_phys), intent(in) :: piln (:,:)    ! (pcols,pverp) - Log interface pressures
    real(kind_phys), intent(in) :: pmln (:,:)    ! (pcols,pver)  - Log midpoint pressures
    real(kind_phys), intent(in) :: pint (:,:)    ! (pcols,pverp) - Interface pressures
    real(kind_phys), intent(in) :: pmid (:,:)    ! (pcols,pver)  - Midpoint pressures
    real(kind_phys), intent(in) :: pdel (:,:)    ! (pcols,pver)  - layer thickness
    real(kind_phys), intent(in) :: rpdel(:,:)    ! (pcols,pver)  - inverse of layer thickness
    real(kind_phys), intent(in) :: t    (:,:)    ! (pcols,pver)  - temperature
    real(kind_phys), intent(in) :: q    (:,:)    ! (pcols,pver)  - specific humidity
    real(kind_phys), intent(in) :: rair (:,:)    ! (pcols,pver)  - Gas constant for dry air
    real(kind_phys), intent(in) :: gravit        !               - Acceleration of gravity
    real(kind_phys), intent(in) :: zvir (:,:)    ! (pcols,pver)  - rh2o/rair - 1

! Output arguments

    real(kind_phys), intent(out) :: zi(:,:)      ! (pcols,pverp) - Height above surface at interfaces
    real(kind_phys), intent(out) :: zm(:,:)      ! (pcols,pver)  - Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
    integer  :: i,k,vert_step,indx_surf        ! Lon, level indices, vertical step, index of surface
    integer  :: indx_level, indx_above         ! Index for level, index for above
    real(kind_phys) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(kind_phys) :: hkl(ncol)               ! off-diagonal element
    real(kind_phys) :: rog(ncol,pver)          ! Rair / gravit
    real(kind_phys) :: tv                      ! virtual temperature
    real(kind_phys) :: tvfac                   ! Tv/T
!
!-----------------------------------------------------------------------
!
    rog(:ncol,:) = rair(:ncol,:) / gravit

    if (vert_surf > vert_toa) then 
       vert_step = -1
       indx_surf = vert_toa +1
       indx_level = 1
       indx_above = 0
    else
       vert_step = 1
       indx_surf = vert_surf
       indx_level = 0
       indx_above = 1
    end if


! The surface height is zero by definition. 
    do i = 1,ncol
       zi(i,indx_surf) = 0.0_kind_phys
    end do

! Compute zi, zm from bottom up. 
! Note, zi(i,k) is the interface above zm(i,k) when ordered top to bot

    do k = vert_surf, vert_toa, vert_step

! First set hydrostatic elements consistent with dynamics

       if (fvdyn) then
          do i = 1,ncol
             hkl(i) = piln(i,k+indx_level) - piln(i,k+indx_above)
             hkk(i) = 1._kind_phys - pint(i,k+indx_above) * hkl(i) * rpdel(i,k)
          end do
       else
          do i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_kind_phys * hkl(i)
          end do
       end if

! Now compute tv, zm, zi

       do i = 1,ncol
          tvfac   = 1._kind_phys + zvir(i,k) * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k)            = zi(i,k+indx_level) + rog(i,k) * tv * hkk(i)
          zi(i,k+indx_above) = zi(i,k+indx_level) + rog(i,k) * tv * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_t
end module geopotential
