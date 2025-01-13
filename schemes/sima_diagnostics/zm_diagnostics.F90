module zm_diagnostics
   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: zm_diagnostics_init ! init routine
   public :: zm_diagnostics_run  ! main routine

CONTAINS

 !> \section arg_table_zm_diagnostics_init  Argument Table
 !! \htmlinclude zm_diagnostics_init.html
 subroutine zm_diagnostics_init(errmsg, errflg)
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables:

    errmsg = ''
    errflg = 0

    call history_add_field ('PRECZ', 'lwe_precipitation_rate_at_surface_due_to_deep_convection_due_to_Zhang-McFarlane', &
                            horiz_only, 'avg', 'm s-1')
    call history_add_field ('ZMFLXPRC', 'precipitation_flux_at_interface_due_to_deep_convection', 'ilev', 'avg', 'kg m-2 s-1')
    call history_add_field ('ZMFLXSNW', 'frozen_precipitation_flux_at_interface_due_to_deep_convection', &
                            'ilev', 'avg', 'kg m-2 s-1')
    call history_add_field ('ZMNTPRPD', 'tendency_of_precipitation_wrt_moist_air_and_condensed_water_due_to_deep_convection', &
                            'lev', 'avg', 'kg kg-1 s-1')
    call history_add_field ('ZMNTSNPD', &
                            'tendency_of_frozen_precipitation_wrt_moist_air_and_condensed_water_due_to_deep_convection', &
                            'lev', 'avg', 'kg kg-1 s-1')

    call history_add_field ('CMFMC_DP', 'Convection mass flux from ZM deep ', 'ilev', 'avg', 'kg m-2 s-1')
    call history_add_field ('PRECCDZM', 'lwe_precipitation_rate_at_surface_due_to_deep_convection', horiz_only, 'avg', 'm s-1')

    call history_add_field ('PCONVB', 'convection base pressure',  horiz_only ,  'avg', 'Pa'    )
    call history_add_field ('PCONVT', 'convection top  pressure',  horiz_only ,  'avg', 'Pa'    )

    call history_add_field ('CAPE',   'zhang_mcfarlane_convective_available_potential_energy',  horiz_only,   'avg', 'J kg-1')
    call history_add_field ('FREQZM', 'Fractional occurance of ZM convection',  horiz_only  , 'avg', 'fraction')

    call history_add_field ('ZMMU',   'ZM convection updraft mass flux',  'lev',  'avg', 'kg m-2 s-1')
    call history_add_field ('ZMMD',   'ZM convection downdraft mass flux', 'lev',  'avg', 'kg m-2 s-1')

    call history_add_field ('ZMUPGU',&
                            'tendency_of_eastward_wind_due_to_zhang_mcfarlane_deep_convective_updraft_pressure_gradient_term',&
                            'lev',  'avg', 'm s-2')
    call history_add_field ('ZMUPGD',&
                            'tendency_of_eastward_wind_due_to_zhang_mcfarlane_deep_convective_downdraft_pressure_gradient_term',&
                            'lev',  'avg', 'm s-2')
    call history_add_field ('ZMVPGU',&
                            'tendency_of_northward_wind_due_to_zhang_mcfarlane_deep_convective_updraft_pressure_gradient_term',&
                            'lev',  'avg', 'm s-2')
    call history_add_field ('ZMVPGD',&
                            'tendency_of_northward_wind_due_to_zhang_mcfarlane_deep_convective_downdraft_pressure_gradient_term',&
                            'lev',  'avg', 'm s-2')

    call history_add_field ('ZMICUU', 'in_cloud_eastward_wind_in_updraft_due_to_deep_convection',  'lev',  'avg', 'm s-1')
    call history_add_field ('ZMICUD', 'in_cloud_eastward_wind_in_downdraft_due_to_deep_convection',  'lev',  'avg', 'm s-1')
    call history_add_field ('ZMICVU', 'in_cloud_northward_wind_in_updraft_due_to_deep_convection',  'lev',  'avg', 'm s-1')
    call history_add_field ('ZMICVD', 'in_cloud_northward_wind_in_downdraft_due_to_deep_convection', 'lev',  'avg', 'm s-1')

    call history_add_field ('DLFZM',   'detrainment_of_cloud_liquid_due_to_deep_convection', 'lev', 'avg','kg kg-1 s-1 ')


   end subroutine zm_diagnostics_init

   !> \section arg_table_zm_diagnostics_run  Argument Table
   !! \htmlinclude zm_diagnostics_run.html
   subroutine zm_diagnostics_run(ncol, pver, pverp, ideep, cpair, prec, cape, gravit, mu, md, &
!      dif, dlf, ps, pap, maxg, jt, flxprec, flxsnow, ntprprd, ntsnprd, pguallu, pguallv, &
      dlf, ps, pap, maxg, jt, flxprec, flxsnow, ntprprd, ntsnprd, pguallu, pguallv, &
      pgdallu, pgdallv, icwuu, icwuv, icwdu, icwdv, mcon, errmsg, errflg)

      use cam_history, only: history_out_field

      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: pver
      integer, intent(in) :: pverp
!      real(kind_phys), intent(in) :: const_array(:,:,:)
      integer, intent(in) :: ideep(:)

      real(kind_phys), intent(in) :: cpair
      real(kind_phys), intent(in) :: prec(:)
      real(kind_phys), intent(in) :: cape(:)
      real(kind_phys), intent(in) :: gravit
      real(kind_phys), intent(in) :: mu(:,:)
      real(kind_phys), intent(in) :: md(:,:)
!      real(kind_phys), intent(in) :: dif(:,:)
      real(kind_phys), intent(in) :: dlf(:,:)
      real(kind_phys), intent(in) :: ps(:)
      real(kind_phys), intent(in) :: pap(:,:)
      integer, intent(in) :: maxg(ncol)
      integer, intent(in) :: jt(ncol)
      real(kind_phys),intent(in) :: flxprec(:,:)
      real(kind_phys),intent(in) :: flxsnow(:,:)
      real(kind_phys),intent(in) :: ntprprd(:,:)
      real(kind_phys),intent(in) :: ntsnprd(:,:)
      real(kind_phys),intent(in) :: pguallu(:,:)
      real(kind_phys),intent(in) :: pguallv(:,:)
      real(kind_phys),intent(in) :: pgdallu(:,:)
      real(kind_phys),intent(in) :: pgdallv(:,:)
      real(kind_phys),intent(in) :: icwuu(:,:)
      real(kind_phys),intent(in) :: icwuv(:,:)
      real(kind_phys),intent(in) :: icwdu(:,:)
      real(kind_phys),intent(in) :: icwdv(:,:)

      real(kind_phys),intent(inout) :: mcon(:,:)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: lengath  ! number of columns with deep convection

      real(kind_phys) :: freqzm(ncol)
      real(kind_phys) :: pcont(ncol)
      real(kind_phys) :: pconb(ncol)
      real(kind_phys) :: ftem(ncol,pver)
      real(kind_phys) :: mconzm(ncol,pverp)
      real(kind_phys) :: mu_out(ncol,pver)
      real(kind_phys) :: md_out(ncol,pver)

      integer :: index_cldliq
      integer i, ii, k

      errmsg = ''
      errflg = 0

      mu_out(:,:) = 0._kind_phys
      md_out(:,:) = 0._kind_phys

      lengath = count(ideep > 0)
      if (lengath > ncol) lengath = ncol  ! should not happen, but force it to not be larger than ncol for safety sake

      call history_out_field('CAPE', cape)

      freqzm(:) = 0._kind_phys
      do i = 1,lengath
         freqzm(ideep(i)) = 1.0_kind_phys
      end do
      call history_out_field('FREQZM  ',freqzm)

      mcon(:ncol,:pverp) = mcon(:ncol,:pverp) * 100._kind_phys/gravit
      mconzm(:ncol,:pverp) = mcon(:ncol,:pverp)

      call history_out_field('CMFMC_DP', mconzm)


   ! Store upward and downward mass fluxes in un-gathered arrays
   ! + convert from mb/s to kg/m^2/s
   do i=1,lengath
      do k=1,pver
         ii = ideep(i)
         mu_out(ii,k) = mu(i,k) * 100._kind_phys/gravit
         md_out(ii,k) = md(i,k) * 100._kind_phys/gravit
      end do
   end do

   call history_out_field('ZMMU', mu_out)
   call history_out_field('ZMMD', md_out)

!   call history_out_field('DIFZM'   ,dif)
   call history_out_field('DLFZM'   ,dlf)

   pcont(:ncol) = ps(:ncol)
   pconb(:ncol) = ps(:ncol)
   do i = 1,lengath
       if (maxg(i).gt.jt(i)) then
          pcont(ideep(i)) = pap(ideep(i),jt(i))  ! gathered array (or jctop ungathered)
          pconb(ideep(i)) = pap(ideep(i),maxg(i))! gathered array
       endif
       !     write(iulog,*) ' pcont, pconb ', pcont(i), pconb(i), cnt(i), cnb(i)
    end do
    call history_out_field('PCONVT  ',pcont)
    call history_out_field('PCONVB  ',pconb)

   call history_out_field('ZMFLXPRC', flxprec)
   call history_out_field('ZMFLXSNW', flxsnow)
   call history_out_field('ZMNTPRPD', ntprprd)
   call history_out_field('ZMNTSNPD', ntsnprd)

!CAM was outputting the exact same quantity to both fields
   call history_out_field('PRECCDZM   ',prec)
   call history_out_field('PRECZ   ', prec)

   ! Output apparent force from  pressure gradient
   call history_out_field('ZMUPGU', pguallu)
   call history_out_field('ZMUPGD', pgdallu)
   call history_out_field('ZMVPGU', pguallv)
   call history_out_field('ZMVPGD', pgdallv)

   ! Output in-cloud winds
   call history_out_field('ZMICUU', icwuu)
   call history_out_field('ZMICUD', icwdu)
   call history_out_field('ZMICVU', icwuv)
   call history_out_field('ZMICVD', icwdv)

   end subroutine zm_diagnostics_run

   !=======================================================================

end module zm_diagnostics
