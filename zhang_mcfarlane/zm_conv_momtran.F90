module zm_conv_momtran

  use ccpp_kinds, only:  kind_phys

  implicit none

  save
  private                         ! Make default type private to the module
  public zm_conv_momtran_run      ! convective momentum transport
  integer, parameter, private :: num_winds=2  ! Number of wind directions (for historical purposes)


contains

!===============================================================================
!> \section arg_table_zm_conv_momtran_run Argument Table
!! \htmlinclude zm_conv_momtran_run.html
!!
subroutine zm_conv_momtran_run(ncol, pver, pverp, &
                    domomtran,windu, windv, mu, md, &
                    momcu, momcd, &
                    du, eu, ed, dp, dsubcld , &
                    jt, mx, ideep , il1g, il2g, &
                    nstep, windu_tend, windv_tend, pguallu, pguallv, pgdallu, pgdallv, &
                    icwuu, icwuv, icwdu, icwdv, dt, seten)
!-----------------------------------------------------------------------
!
! Purpose:
! Convective transport of momentum
!
! Mixing ratios may be with respect to either dry or moist air
!
! Method:
! Based on the convtran subroutine by P. Rasch
! <Also include any applicable external references.>
!
! Author: J. Richter and P. Rasch
!
!-----------------------------------------------------------------------

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: pver, pverp
   logical, intent(in) :: domomtran      ! flag for doing convective transport    (num_winds)
   real(kind_phys), intent(in) :: windu(:,:)  ! U Wind array                                    (ncol,pver)
   real(kind_phys), intent(in) :: windv(:,:)  ! V Wind array                                    (ncol,pver)
   real(kind_phys), intent(in) :: mu(:,:)       ! Mass flux up                              (ncol,pver)
   real(kind_phys), intent(in) :: md(:,:)       ! Mass flux down                            (ncol,pver)
   real(kind_phys), intent(in) :: momcu
   real(kind_phys), intent(in) :: momcd
   real(kind_phys), intent(in) :: du(:,:)       ! Mass detraining from updraft              (ncol,pver)
   real(kind_phys), intent(in) :: eu(:,:)       ! Mass entraining from updraft              (ncol,pver)
   real(kind_phys), intent(in) :: ed(:,:)       ! Mass entraining from downdraft            (ncol,pver)
   real(kind_phys), intent(in) :: dp(:,:)       ! Delta pressure between interfaces         (ncol,pver)
   real(kind_phys), intent(in) :: dsubcld(:)       ! Delta pressure from cloud base to sfc  (ncol)
   real(kind_phys), intent(in) :: dt    ! time step in seconds

   integer, intent(in) :: jt(:)         ! Index of cloud top for each column         (ncol)
   integer, intent(in) :: mx(:)         ! Index of cloud top for each column         (ncol)
   integer, intent(in) :: ideep(:)      ! Gathering array                            (ncol)
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index



! input/output

   real(kind_phys), intent(out) :: windu_tend(:,:)  ! U wind tendency
   real(kind_phys), intent(out) :: windv_tend(:,:)  ! V wind tendency

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer kkm1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index
   integer ii                 ! Work index

   real(kind_phys) cabv                 ! Mix ratio of constituent above
   real(kind_phys) cbel                 ! Mix ratio of constituent below
   real(kind_phys) cdifr                ! Normalized diff between cabv and cbel
   real(kind_phys) chat(ncol,pver)     ! Mix ratio in env at interfaces
   real(kind_phys) cond(ncol,pver)     ! Mix ratio in downdraft at interfaces
   real(kind_phys) const(ncol,pver)    ! Gathered wind array
   real(kind_phys) conu(ncol,pver)     ! Mix ratio in updraft at interfaces
   real(kind_phys) dcondt(ncol,pver)   ! Gathered tend array
   real(kind_phys) mbsth                ! Threshold for mass fluxes
   real(kind_phys) mupdudp              ! A work variable
   real(kind_phys) minc                 ! A work variable
   real(kind_phys) maxc                 ! A work variable
   real(kind_phys) fluxin               ! A work variable
   real(kind_phys) fluxout              ! A work variable
   real(kind_phys) netflux              ! A work variable


   real(kind_phys) sum                  ! sum
   real(kind_phys) sum2                  ! sum2

   real(kind_phys) mududp(ncol,pver) ! working variable
   real(kind_phys) mddudp(ncol,pver)     ! working variable

   real(kind_phys) pgu(ncol,pver)      ! Pressure gradient term for updraft
   real(kind_phys) pgd(ncol,pver)      ! Pressure gradient term for downdraft

   real(kind_phys),intent(out) ::  pguallu(:,:)      ! Apparent force from  updraft PG on U winds  ! (ncol,pver)
   real(kind_phys),intent(out) ::  pguallv(:,:)      ! Apparent force from  updraft PG on V winds  ! (ncol,pver)
   real(kind_phys),intent(out) ::  pgdallu(:,:)      ! Apparent force from  downdraft PG on U winds! (ncol,pver)
   real(kind_phys),intent(out) ::  pgdallv(:,:)      ! Apparent force from  downdraft PG on V winds! (ncol,pver)

   real(kind_phys),intent(out) ::  icwuu(:,:)      ! In-cloud U winds in updraft           ! (ncol,pver)
   real(kind_phys),intent(out) ::  icwuv(:,:)      ! In-cloud V winds in updraft           ! (ncol,pver)
   real(kind_phys),intent(out) ::  icwdu(:,:)      ! In-cloud U winds in downdraft         ! (ncol,pver)
   real(kind_phys),intent(out) ::  icwdv(:,:)      ! In-cloud V winds in downdraft         ! (ncol,pver)

   real(kind_phys),intent(out) ::  seten(:,:) ! Dry static energy tendency                ! (ncol,pver)
   real(kind_phys)                 gseten(ncol,pver) ! Gathered dry static energy tendency

   real(kind_phys) :: winds(ncol,pver,num_winds)       ! combined winds array
   real(kind_phys) :: wind_tends(ncol,pver,num_winds)  ! combined tendency array
   real(kind_phys) :: pguall(ncol,pver,num_winds)      ! Combined apparent force from  updraft PG on U winds
   real(kind_phys) :: pgdall(ncol,pver,num_winds)      ! Combined apparent force from  downdraft PG on U winds
   real(kind_phys) :: icwu(ncol,pver,num_winds)        ! Combined In-cloud winds in updraft
   real(kind_phys) :: icwd(ncol,pver,num_winds)        ! Combined In-cloud winds in downdraft

   real(kind_phys)  mflux(ncol,pverp,num_winds)   ! Gathered momentum flux

   real(kind_phys)  wind0(ncol,pver,num_winds)       !  gathered  wind before time step
   real(kind_phys)  windf(ncol,pver,num_winds)       !  gathered  wind after time step
   real(kind_phys) fkeb, fket, ketend_cons, ketend, utop, ubot, vtop, vbot, gset2


!-----------------------------------------------------------------------
!
! Combine winds in single array
   winds(:,:,1) = windu(:,:)
   winds(:,:,2) = windv(:,:)

! Initialize outgoing fields
   pguall(:,:,:)     = 0.0_kind_phys
   pgdall(:,:,:)     = 0.0_kind_phys
! Initialize in-cloud winds to environmental wind
   icwu(:ncol,:,:)       = winds(:ncol,:,:)
   icwd(:ncol,:,:)       = winds(:ncol,:,:)

! Initialize momentum flux and  final winds
   mflux(:,:,:)       = 0.0_kind_phys
   wind0(:,:,:)         = 0.0_kind_phys
   windf(:,:,:)         = 0.0_kind_phys

! Initialize dry static energy

   seten(:,:)         = 0.0_kind_phys
   gseten(:,:)         = 0.0_kind_phys

! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_kind_phys

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each wind component
   do m = 1, num_winds                    !start at m = 1 to transport momentum
      if (domomtran) then

! Gather up the winds and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = winds(ideep(i),k,m)
                wind0(i,k,m) = const(i,k)
            end do
         end do


! From now on work only with gathered data

! Interpolate winds to interfaces

         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g

               ! use arithmetic mean
               chat(i,k) = 0.5_kind_phys* (const(i,k)+const(i,km1))

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._kind_phys

            end do
         end do


!
! Pressure Perturbation Term
!

      !Top boundary:  assume mu is zero

         k=1
         pgu(:il2g,k) = 0.0_kind_phys
         pgd(:il2g,k) = 0.0_kind_phys

         do k=2,pver-1
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

               !interior points

               mududp(i,k) =  ( mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  mu(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))

               pgu(i,k) = - momcu * 0.5_kind_phys * mududp(i,k)


               mddudp(i,k) =  ( md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  md(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))

               pgd(i,k) = - momcd * 0.5_kind_phys * mddudp(i,k)


            end do
         end do

       ! bottom boundary
       k = pver
       km1 = max(1,k-1)
       do i=il1g,il2g

          mududp(i,k) =   mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1)
          pgu(i,k) = - momcu *  mududp(i,k)

          mddudp(i,k) =   md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1)

          pgd(i,k) = - momcd * mddudp(i,k)

       end do


!
! In-cloud velocity calculations
!

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         kkm1 = max(1,kk-1)
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
            if (mupdudp > mbsth) then

               conu(i,kk) = (+eu(i,kk)*const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-ed(i,km1)*const(i,km1)*dp(i,km1))-pgd(i,km1)*dp(i,km1)/md(i,k)
            endif


         end do



! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkm1 = max(1,kk-1)
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
               if (mupdudp > mbsth) then

                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eu(i,kk)* &
                                  const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
               endif
            end do

         end do


! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then

                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-ed(i,km1)*const(i,km1) &
                                  *dp(i,km1)-pgd(i,km1)*dp(i,km1) )/md(i,k)

               endif
            end do
         end do


         sum = 0._kind_phys
         sum2 = 0._kind_phys


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
               ii = ideep(i)

! version 1 hard to check for roundoff errors
               dcondt(i,k) =  &
                           +(mu(i,kp1)* (conu(i,kp1)-chat(i,kp1)) &
                           -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                           +md(i,kp1)* (cond(i,kp1)-chat(i,kp1)) &
                           -md(i,k)*   (cond(i,k)-chat(i,k)) &
                          )/dp(i,k)

            end do
         end do

  ! dcont for bottom layer
          !
          do k = kbm,pver
             km1 = max(1,k-1)
             do i = il1g,il2g
                if (k == mx(i)) then

                   ! version 1
                   dcondt(i,k) = (1._kind_phys/dp(i,k))*   &
                        (-mu(i,k)*(conu(i,k)-chat(i,k)) &
                        -md(i,k)*(cond(i,k)-chat(i,k)) &
                        )
                end if
             end do
          end do

! Initialize to zero everywhere, then scatter tendency back to full array
         wind_tends(:,:,m) = 0._kind_phys

         do k = 1,pver
            do i = il1g,il2g
               ii = ideep(i)
               wind_tends(ii,k,m) = dcondt(i,k)
    ! Output apparent force on the mean flow from pressure gradient
               pguall(ii,k,m) = -pgu(i,k)
               pgdall(ii,k,m) = -pgd(i,k)
               icwu(ii,k,m)   =  conu(i,k)
               icwd(ii,k,m)   =  cond(i,k)
            end do
         end do

          ! Calculate momentum flux in units of mb*m/s2

          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                mflux(i,k,m) = &
                     -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                     -md(i,k)*   (cond(i,k)-chat(i,k))
             end do
          end do


          ! Calculate winds at the end of the time step

          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                km1 = max(1,k-1)
                kp1 = k+1
                windf(i,k,m) = const(i,k)    -   (mflux(i,kp1,m) - mflux(i,k,m)) * dt /dp(i,k)

             end do
          end do

       end if      ! for domomtran
   end do

 ! Need to add an energy fix to account for the dissipation of kinetic energy
    ! Formulation follows from Boville and Bretherton (2003)
    ! formulation by PJR

    do k = ktm,pver
       km1 = max(1,k-1)
       kp1 = min(pver,k+1)
       do i = il1g,il2g

          ii = ideep(i)

          ! calculate the KE fluxes at top and bot of layer
          ! based on a discrete approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
          utop = (wind0(i,k,1)+wind0(i,km1,1))/2._kind_phys
          vtop = (wind0(i,k,2)+wind0(i,km1,2))/2._kind_phys
          ubot = (wind0(i,kp1,1)+wind0(i,k,1))/2._kind_phys
          vbot = (wind0(i,kp1,2)+wind0(i,k,2))/2._kind_phys
          fket = utop*mflux(i,k,1)   + vtop*mflux(i,k,2)    ! top of layer
          fkeb = ubot*mflux(i,k+1,1) + vbot*mflux(i,k+1,2)  ! bot of layer

          ! divergence of these fluxes should give a conservative redistribution of KE
          ketend_cons = (fket-fkeb)/dp(i,k)

          ! tendency in kinetic energy resulting from the momentum transport
          ketend = ((windf(i,k,1)**2 + windf(i,k,2)**2) - (wind0(i,k,1)**2 + wind0(i,k,2)**2))/dt

          ! the difference should be the dissipation
          gset2 = ketend_cons - ketend
          gseten(i,k) = gset2

       end do

    end do

    ! Scatter dry static energy to full array
    do k = 1,pver
       do i = il1g,il2g
          ii = ideep(i)
          seten(ii,k) = gseten(i,k)

       end do
    end do

! Split out the wind tendencies
   windu_tend(:,:) = wind_tends(:,:,1)
   windv_tend(:,:) = wind_tends(:,:,2)

   pguallu(:,:)     = pguall(:,:,1)
   pguallv(:,:)     = pguall(:,:,2)
   pgdallu(:,:)     = pgdall(:,:,1)
   pgdallv(:,:)     = pgdall(:,:,2)
   icwuu(:ncol,:)       = icwu(:,:,1)
   icwuv(:ncol,:)       = icwu(:,:,2)
   icwdu(:ncol,:)       = icwd(:,:,1)
   icwdv(:ncol,:)       = icwd(:,:,2)

   return
end subroutine zm_conv_momtran_run


end module zm_conv_momtran
