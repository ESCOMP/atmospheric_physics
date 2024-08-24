module zm_conv_convtran

  use ccpp_kinds, only:  kind_phys

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zm_conv_convtran_run     ! convective transport


contains

!===============================================================================
!> \section arg_table_zm_conv_convtran_run Argument Table
!! \htmlinclude zm_conv_convtran_run.html
!!
subroutine zm_conv_convtran_run(ncol, pver, &
                    doconvtran,q       ,ncnst   ,mu      ,md      , &
                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,fracis  ,dqdt    ,dpdry)
!-----------------------------------------------------------------------
!
! Purpose:
! Convective transport of trace species
!
! Mixing ratios may be with respect to either dry or moist air
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: P. Rasch
!
!-----------------------------------------------------------------------
! CACNOTE - replace with CCPP constituents
   use constituents,    only: cnst_get_type_byind

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol
   integer, intent(in) :: pver
   integer, intent(in) :: ncnst          ! number of tracers to transport
   logical, intent(in) :: doconvtran(:)  ! flag for doing convective transport      (ncnst)
   real(kind_phys), intent(in) :: q(:,:,:)      ! Tracer array including moisture          (ncol,pver,ncnst)
   real(kind_phys), intent(in) :: mu(:,:)       ! Mass flux up                             (ncol,pver)
   real(kind_phys), intent(in) :: md(:,:)       ! Mass flux down                           (ncol,pver)
   real(kind_phys), intent(in) :: du(:,:)       ! Mass detraining from updraft             (ncol,pver)
   real(kind_phys), intent(in) :: eu(:,:)       ! Mass entraining from updraft             (ncol,pver)
   real(kind_phys), intent(in) :: ed(:,:)       ! Mass entraining from downdraft           (ncol,pver)
   real(kind_phys), intent(in) :: dp(:,:)       ! Delta pressure between interfaces        (ncol,pver)
   real(kind_phys), intent(in) :: dsubcld(:)    ! Delta pressure from cloud base to sfc    (ncol)
   real(kind_phys), intent(in) :: fracis(:,:,:) ! fraction of tracer that is insoluble     (ncol,pver,ncnst)

   integer, intent(in) :: jt(:)          ! Index of cloud top for each column       (ncol)
   integer, intent(in) :: mx(:)          ! Index of cloud top for each column       (ncol)
   integer, intent(in) :: ideep(:)       ! Gathering array                          (ncol)
   integer, intent(in) :: il1g           ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g           ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep          ! Time step index

   real(kind_phys), intent(in) :: dpdry(:,:)    ! Delta pressure between interfaces        (ncol,pver)

! input/output

   real(kind_phys), intent(out) :: dqdt(:,:,:)  ! Tracer tendency array  (ncol,pver,ncnst)

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index

   real(kind_phys) cabv                 ! Mix ratio of constituent above
   real(kind_phys) cbel                 ! Mix ratio of constituent below
   real(kind_phys) cdifr                ! Normalized diff between cabv and cbel
   real(kind_phys) chat(ncol,pver)     ! Mix ratio in env at interfaces
   real(kind_phys) cond(ncol,pver)     ! Mix ratio in downdraft at interfaces
   real(kind_phys) const(ncol,pver)    ! Gathered tracer array
   real(kind_phys) fisg(ncol,pver)     ! gathered insoluble fraction of tracer
   real(kind_phys) conu(ncol,pver)     ! Mix ratio in updraft at interfaces
   real(kind_phys) dcondt(ncol,pver)   ! Gathered tend array
   real(kind_phys) small                ! A small number
   real(kind_phys) mbsth                ! Threshold for mass fluxes
   real(kind_phys) mupdudp              ! A work variable
   real(kind_phys) minc                 ! A work variable
   real(kind_phys) maxc                 ! A work variable
   real(kind_phys) fluxin               ! A work variable
   real(kind_phys) fluxout              ! A work variable
   real(kind_phys) netflux              ! A work variable

   real(kind_phys) dutmp(ncol,pver)       ! Mass detraining from updraft
   real(kind_phys) eutmp(ncol,pver)       ! Mass entraining from updraft
   real(kind_phys) edtmp(ncol,pver)       ! Mass entraining from downdraft
   real(kind_phys) dptmp(ncol,pver)    ! Delta pressure between interfaces
   real(kind_phys) total(ncol)
   real(kind_phys) negadt,qtmp

!-----------------------------------------------------------------------
!
   small = 1.e-36_kind_phys
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_kind_phys

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each constituent
   do m = 2, ncnst
      if (doconvtran(m)) then

         if (cnst_get_type_byind(m).eq.'dry') then
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dpdry(i,k)
                  dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               end do
            end do
         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dp(i,k)
                  dutmp(i,k) = du(i,k)
                  eutmp(i,k) = eu(i,k)
                  edtmp(i,k) = ed(i,k)
               end do
            end do
         endif
!        dptmp = dp

! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc < 0) then
                  cdifr = 0._kind_phys
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr > 1.E-6_kind_phys) then
                  cabv = max(const(i,km1),maxc*1.e-12_kind_phys)
                  cbel = max(const(i,k),maxc*1.e-12_kind_phys)
                  chat(i,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5_kind_phys* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._kind_phys

            end do
         end do

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               conu(i,kk) = (+eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1))/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eutmp(i,kk)*fisg(i,kk)* &
                                  const(i,kk)*dptmp(i,kk) )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-edtmp(i,km1)*fisg(i,km1)*const(i,km1) &
                                  *dptmp(i,km1) )/md(i,k)
               endif
            end do
         end do


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

! limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*min(chat(i,k),const(i,km1)) &
                         -(md(i,k)  *cond(i,k) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*min(chat(i,kp1),const(i,k)) &
                         -(md(i,kp1)*cond(i,kp1) + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux) < max(fluxin,fluxout)*1.e-12_kind_phys) then
                  netflux = 0._kind_phys
               endif
               dcondt(i,k) = netflux/dptmp(i,k)
            end do
         end do
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then

                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_kind_phys) then
                     netflux = 0._kind_phys
                  endif
                  dcondt(i,k) = netflux/dptmp(i,k)
               else if (k > mx(i)) then
                  dcondt(i,k) = 0._kind_phys
               end if
            end do
         end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._kind_phys
         do k = 1,pver
            kp1 = min(pver,k+1)
            do i = il1g,il2g
               dqdt(ideep(i),k,m) = dcondt(i,k)
            end do
         end do

      end if      ! for doconvtran

   end do

   return
end subroutine zm_conv_convtran_run


end module zm_conv_convtran
