module geopotential_t

   !---------------------------------------------------------------------------
   ! Compute geopotential from temperature
   !
   ! The hydrostatic matrix elements must be consistent with the dynamics
   !      algorithm.
   ! The diagonal element is the itegration weight from interface k minus
   !      vert_step to midpoint k.
   ! The offdiagonal element is the weight between interfaces.
   !
   ! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
   !---------------------------------------------------------------------------

   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public geopotential_t_run

CONTAINS
   !===========================================================================

   !> \section arg_table_geopotential_t_run  Argument Table
   !! \htmlinclude geopotential_t_run.html
   subroutine geopotential_t_run(pver, pverp, vert_surf, vert_toa,            &
        piln, pmln, pint, pmid, pdel, rpdel,                                  &
        t, q, rair, gravit, zvir,                                             &
        zi, zm, ncol, errflg, errmsg)
      use dyn_comp, only: dycore_is

      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! Compute the geopotential height (above the surface) at the midpoints and
      ! interfaces using the input temperatures and pressures.
      !
      !-----------------------------------------------------------------------

      !------------------------------Arguments--------------------------------
      !
      ! Input arguments
      !
      integer,         intent(in)  :: pver
      integer,         intent(in)  :: pverp
      integer,         intent(in)  :: vert_surf
      integer,         intent(in)  :: vert_toa
      integer,         intent(in)  :: ncol       ! Number of longitudes

      real(kind_phys), intent(in)  :: piln (:,:) ! (ncol,pverp) - Log interface pressures
      real(kind_phys), intent(in)  :: pmln (:,:) ! (ncol,pver)  - Log midpoint pressures
      real(kind_phys), intent(in)  :: pint (:,:) ! (ncol,pverp) - Interface pressures
      real(kind_phys), intent(in)  :: pmid (:,:) ! (ncol,pver)  - Midpoint pressures
      real(kind_phys), intent(in)  :: pdel (:,:) ! (ncol,pver)  - layer thickness
      real(kind_phys), intent(in)  :: rpdel(:,:) ! (ncol,pver)  - inverse of layer thickness
      real(kind_phys), intent(in)  :: t    (:,:) ! (ncol,pver)  - temperature
      real(kind_phys), intent(in)  :: q    (:,:) ! (ncol,pver)  - specific humidity
      real(kind_phys), intent(in)  :: rair (:,:) ! (ncol,pver)  - Gas constant for dry air
      real(kind_phys), intent(in)  :: gravit     !               - Acceleration of gravity
      real(kind_phys), intent(in)  :: zvir (:,:) ! (ncol,pver)  - rh2o/rair - 1

      ! Output arguments

      real(kind_phys), intent(out) :: zi(:,:)    ! (ncol,pverp) - Height above surface at interfaces
      real(kind_phys), intent(out) :: zm(:,:)    ! (ncol,pver)  - Geopotential height at mid level
      integer,            intent(out) :: errflg
      character(len=512), intent(out) :: errmsg
      !
      !---------------------------Local variables-----------------------------
      !
      integer                      :: i,k,vert_step,indx_surf ! Lon, level indices, vertical step, index of surface
      integer                      :: indx_level, indx_above  ! Index for level, index for above
      real(kind_phys)              :: hkk(ncol)               ! diagonal element of hydrostatic matrix
      real(kind_phys)              :: hkl(ncol)               ! off-diagonal element
      real(kind_phys)              :: rog(ncol,pver)          ! Rair / gravit
      real(kind_phys)              :: tv                      ! virtual temperature
      real(kind_phys)              :: tvfac                   ! Tv/T
      logical                      :: fvdyn
      !
      !-----------------------------------------------------------------------
      !

      errmsg = ''
      errflg = 0

      fvdyn = dycore_is('LR')

      rog(:ncol,:) = rair(:ncol,:) / gravit

      !!XXgoldyXX ==> cacraigucar: Please explain index_xxx
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
         zi(i, vert_surf) = 0.0_kind_phys
      end do

      ! Compute zi, zm from bottom up.
      ! Note, zi(i,k) is the interface above zm(i,k) when ordered top to bot

      do k = vert_surf, vert_toa, vert_step

         ! First set hydrostatic elements consistent with dynamics

         if (fvdyn) then
            do i = 1,ncol
               hkl(i) = piln(i,k-vert_step) - piln(i,k)
               hkk(i) = 1._kind_phys - pint(i,k) * hkl(i) * rpdel(i,k)
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

            zm(i,k)            = zi(i,k-vert_step) + rog(i,k) * tv * hkk(i)
            zi(i,k+indx_above) = zi(i,k-vert_step) + rog(i,k) * tv * hkl(i)
         end do
      end do

   end subroutine geopotential_t_run

end module geopotential_t
