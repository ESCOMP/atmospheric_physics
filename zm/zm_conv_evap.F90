module zm_conv_evap

  use ccpp_kinds, only:  kind_phys

! CACNOTE - Need to ccpp'ize cloud_fraction
  use cloud_fraction,  only: cldfrc_fice

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zm_conv_evap_run         ! evaporation of precip from ZM schemea

contains


!===============================================================================
!> \section arg_table_zm_conv_evap_run Argument Table
!! \htmlinclude zm_conv_evap_run.html
!!
subroutine zm_conv_evap_run(ncol, pver, pverp, &
     gravit, latice, latvap, tmelt, &
     cpres, ke, ke_lnd, zm_org, &
     t,pmid,pdel,q, &
     landfrac, &
     tend_s, tend_s_snwprd, tend_s_snwevmlt, tend_q, &
     prdprec, cldfrc, deltat,  &
     prec, snow, ntprprd, ntsnprd, flxprec, flxsnow, prdsnow)


!-----------------------------------------------------------------------
! Compute tendencies due to evaporation of rain from ZM scheme
!--
! Compute the total precipitation and snow fluxes at the surface.
! Add in the latent heat of fusion for snow formation and melt, since it not dealt with
! in the Zhang-MacFarlane parameterization.
! Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
!-----------------------------------------------------------------------

!CACNOTE - Not sure what to do about qsat_water
    use wv_saturation,  only: qsat

!------------------------------Arguments--------------------------------
    integer,intent(in) :: ncol                               ! number of columns
    integer,intent(in) :: pver, pverp
    real(kind_phys),intent(in) :: gravit                     ! gravitational acceleration (m s-2)
    real(kind_phys),intent(in) :: latice                     ! Latent heat of fusion (J kg-1)
    real(kind_phys),intent(in) :: latvap                     ! Latent heat of vaporization (J kg-1)
    real(kind_phys),intent(in) :: tmelt                      ! Freezing point of water (K)
    real(kind_phys), intent(in) :: cpres      ! specific heat at constant pressure in j/kg-degk.
    real(kind_phys), intent(in) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
    real(kind_phys), intent(in) :: ke_lnd
    logical, intent(in)         :: zm_org
    real(kind_phys),intent(in), dimension(:,:) :: t          ! temperature (K)                              (ncol,pver)
    real(kind_phys),intent(in), dimension(:,:) :: pmid       ! midpoint pressure (Pa)                       (ncol,pver)
    real(kind_phys),intent(in), dimension(:,:) :: pdel       ! layer thickness (Pa)                         (ncol,pver)
    real(kind_phys),intent(in), dimension(:,:) :: q          ! water vapor (kg/kg)                          (ncol,pver)
    real(kind_phys),intent(in), dimension(:) :: landfrac     ! land fraction                                (ncol)
    real(kind_phys),intent(inout), dimension(:,:) :: tend_s     ! heating rate (J/kg/s)                     (ncol,pver)
    real(kind_phys),intent(inout), dimension(:,:) :: tend_q     ! water vapor tendency (kg/kg/s)            (ncol,pver)
    real(kind_phys),intent(out), dimension(:,:) :: tend_s_snwprd ! Heating rate of snow production        (ncol,pver)
    real(kind_phys),intent(out), dimension(:,:) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow (ncol,pver)



    real(kind_phys), intent(in   ) :: prdprec(:,:)! precipitation production (kg/ks/s)                      (ncol,pver)
    real(kind_phys), intent(in   ) :: cldfrc(:,:) ! cloud fraction                                          (ncol,pver)
    real(kind_phys), intent(in   ) :: deltat             ! time step

    real(kind_phys), intent(inout) :: prec(:)        ! Convective-scale preciptn rate                       (ncol)
    real(kind_phys), intent(out)   :: snow(:)        ! Convective-scale snowfall rate                       (ncol)

    real(kind_phys), optional, intent(in), allocatable  :: prdsnow(:,:) ! snow production (kg/ks/s)

!
!---------------------------Local storage-------------------------------

    real(kind_phys) :: es    (ncol,pver)    ! Saturation vapor pressure
    real(kind_phys) :: fice   (ncol,pver)    ! ice fraction in precip production
    real(kind_phys) :: fsnow_conv(ncol,pver) ! snow fraction in precip production
    real(kind_phys) :: qs   (ncol,pver)    ! saturation specific humidity
    real(kind_phys),intent(out) :: flxprec(:,:)   ! Convective-scale flux of precip at interfaces (kg/m2/s) ! (ncol,pverp)
    real(kind_phys),intent(out) :: flxsnow(:,:)   ! Convective-scale flux of snow   at interfaces (kg/m2/s) ! (ncol,pverp)
    real(kind_phys),intent(out) :: ntprprd(:,:)   ! net precip production in layer                          ! (ncol,pver)
    real(kind_phys),intent(out) :: ntsnprd(:,:)   ! net snow production in layer                            ! (ncol,pver)
    real(kind_phys) :: work1                  ! temp variable (pjr)
    real(kind_phys) :: work2                  ! temp variable (pjr)

    real(kind_phys) :: evpvint(ncol)         ! vertical integral of evaporation
    real(kind_phys) :: evpprec(ncol)         ! evaporation of precipitation (kg/kg/s)
    real(kind_phys) :: evpsnow(ncol)         ! evaporation of snowfall (kg/kg/s)
    real(kind_phys) :: snowmlt(ncol)         ! snow melt tendency in layer
    real(kind_phys) :: flxsntm(ncol)         ! flux of snow into layer, after melting

    real(kind_phys) :: kemask
    real(kind_phys) :: evplimit               ! temp variable for evaporation limits
    real(kind_phys) :: rlat(ncol)
    real(kind_phys) :: dum
    real(kind_phys) :: omsm

    integer :: i,k                     ! longitude,level indices
    logical :: old_snow


!-----------------------------------------------------------------------

    ! If prdsnow is passed in and allocated, then use it in the calculation, otherwise
    ! use the old snow calculation
    old_snow=.true.
    if (present(prdsnow)) then
       if (allocated(prdsnow)) then
          old_snow=.false.
       end if
    end if

! convert input precip to kg/m2/s
    prec(:ncol) = prec(:ncol)*1000._kind_phys

! determine saturation vapor pressure
    do k = 1,pver
       call qsat(t(1:ncol,k), pmid(1:ncol,k), es(1:ncol,k), qs(1:ncol,k), ncol)
    end do
! determine ice fraction in rain production (use cloud water parameterization fraction at present)
!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
    fice(:,:) = 0._kind_phys
    fsnow_conv(:,:) = 0._kind_phys
!REMOVECAM_END
    call cldfrc_fice(ncol, t(1:ncol,:), fice(1:ncol,:), fsnow_conv(1:ncol,:))

! zero the flux integrals on the top boundary
    flxprec(:ncol,1) = 0._kind_phys
    flxsnow(:ncol,1) = 0._kind_phys
    evpvint(:ncol)   = 0._kind_phys
    omsm=0.9999_kind_phys

    do k = 1, pver
       do i = 1, ncol

! Melt snow falling into layer, if necessary.
        if( old_snow ) then
          if (t(i,k) > tmelt) then
             flxsntm(i) = 0._kind_phys
             snowmlt(i) = flxsnow(i,k) * gravit/ pdel(i,k)
          else
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0._kind_phys
          end if
        else
          ! make sure melting snow doesn't reduce temperature below threshold
          if (t(i,k) > tmelt) then
              dum = -latice/cpres*flxsnow(i,k)*gravit/pdel(i,k)*deltat
              if (t(i,k) + dum .le. tmelt) then
                dum = (t(i,k)-tmelt)*cpres/latice/deltat
                dum = dum/(flxsnow(i,k)*gravit/pdel(i,k))
                dum = max(0._kind_phys,dum)
                dum = min(1._kind_phys,dum)
              else
                dum = 1._kind_phys
              end if
              dum = dum*omsm
              flxsntm(i) = flxsnow(i,k)*(1.0_kind_phys-dum)
              snowmlt(i) = dum*flxsnow(i,k)*gravit/ pdel(i,k)
          else
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0._kind_phys
          end if
        end if

! relative humidity depression must be > 0 for evaporation
          evplimit = max(1._kind_phys - q(i,k)/qs(i,k), 0._kind_phys)

          if (zm_org) then
             kemask = ke * (1._kind_phys - landfrac(i)) + ke_lnd * landfrac(i)
          else
             kemask = ke
          endif

! total evaporation depends on flux in the top of the layer
! flux prec is the net production above layer minus evaporation into environmet
          evpprec(i) = kemask * (1._kind_phys - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k))

! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
! Currently does not include heating/cooling change to qs
          evplimit   = max(0._kind_phys, (qs(i,k)-q(i,k)) / deltat)

! Don't evaporate more than is falling into the layer - do not evaporate rain formed
! in this layer but if precip production is negative, remove from the available precip
! Negative precip production occurs because of evaporation in downdrafts.
          evplimit   = min(evplimit, flxprec(i,k) * gravit / pdel(i,k))

! Total evaporation cannot exceed input precipitation
          evplimit   = min(evplimit, (prec(i) - evpvint(i)) * gravit / pdel(i,k))

          evpprec(i) = min(evplimit, evpprec(i))
          if( .not.old_snow ) then
            evpprec(i) = max(0._kind_phys, evpprec(i))
            evpprec(i) = evpprec(i)*omsm
          end if


! evaporation of snow depends on snow fraction of total precipitation in the top after melting
          if (flxprec(i,k) > 0._kind_phys) then
!            prevent roundoff problems
             work1 = min(max(0._kind_phys,flxsntm(i)/flxprec(i,k)),1._kind_phys)
             evpsnow(i) = evpprec(i) * work1
          else
             evpsnow(i) = 0._kind_phys
          end if

! vertically integrated evaporation
          evpvint(i) = evpvint(i) + evpprec(i) * pdel(i,k)/gravit

! net precip production is production - evaporation
          ntprprd(i,k) = prdprec(i,k) - evpprec(i)
! net snow production is precip production * ice fraction - evaporation - melting
! the small amount added to flxprec in the work1 expression has been increased from
! 1e-36 to 8.64e-11 (1e-5 mm/day).  This causes the temperature based partitioning
! scheme to be used for small flxprec amounts.  This is to address error growth problems.

      if( old_snow ) then
          if (flxprec(i,k).gt.0._kind_phys) then
             work1 = min(max(0._kind_phys,flxsnow(i,k)/flxprec(i,k)),1._kind_phys)
          else
             work1 = 0._kind_phys
          endif

          work2 = max(fsnow_conv(i,k), work1)
          if (snowmlt(i).gt.0._kind_phys) work2 = 0._kind_phys
          ntsnprd(i,k) = prdprec(i,k)*work2 - evpsnow(i) - snowmlt(i)
          tend_s_snwprd  (i,k) = prdprec(i,k)*work2*latice
          tend_s_snwevmlt(i,k) = - ( evpsnow(i) + snowmlt(i) )*latice
      else
          ntsnprd(i,k) = prdsnow(i,k) - min(flxsnow(i,k)*gravit/pdel(i,k), evpsnow(i)+snowmlt(i))
          tend_s_snwprd  (i,k) = prdsnow(i,k)*latice
          tend_s_snwevmlt(i,k) = -min(flxsnow(i,k)*gravit/pdel(i,k), evpsnow(i)+snowmlt(i) )*latice
      end if

! precipitation fluxes
          flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * pdel(i,k)/gravit
          flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * pdel(i,k)/gravit

! protect against rounding error
          flxprec(i,k+1) = max(flxprec(i,k+1), 0._kind_phys)
          flxsnow(i,k+1) = max(flxsnow(i,k+1), 0._kind_phys)

! heating (cooling) and moistening due to evaporation
! - latent heat of vaporization for precip production has already been accounted for
! - snow is contained in prec
          if( old_snow ) then
             tend_s(i,k)   =-evpprec(i)*latvap + ntsnprd(i,k)*latice
          else
             tend_s(i,k)   =-evpprec(i)*latvap + tend_s_snwevmlt(i,k)
          end if
          tend_q(i,k) = evpprec(i)
       end do
    end do

! set output precipitation rates (m/s)
    prec(:ncol) = flxprec(:ncol,pver+1) / 1000._kind_phys
    snow(:ncol) = flxsnow(:ncol,pver+1) / 1000._kind_phys

  end subroutine zm_conv_evap_run


end module zm_conv_evap
