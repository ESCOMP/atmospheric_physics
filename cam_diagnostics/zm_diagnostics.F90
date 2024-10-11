module zm_diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! THIS IS A TEMPLATE
!   1. copy this file to a new file with the correct name
!        (SCHEME_diagnostics.F90)
!   2. do a search and replace for "SCHEME" in this file and
!        replace with your scheme name
!   3. Add desired history_add_field calls to the init phase
!   4. Add all fields that are being output as inputs to the run phase
!   5. Add desired history_out_field calls to the run phase
!   6. Run $ccpp_framework/scripts/ccpp_fortran_to_metadata.py on this .F90
!        file to generate the metadata
!   7. Complete the metadata (fill out standard names, units, dimensions)
!   8. Add this scheme to the SDF file for your suite (likely will be at end)
!   9. Delete this header section
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    call history_add_field ('ZM_ORG     ', 'zhang_mcfarlane_organization_parameter_of_deep_convection', 'lev', 'avg', '1')
    call history_add_field ('ZM_ORG2D   ', 'zhang_mcfarlane_organization_parameter_of_deep_convection_copied_to_whole_column',&
                            'lev', 'avg', '1')

    call history_add_field ('PRECZ', 'lwe_precipitation_rate_at_surface_due_to_deep_convection_due_to_Zhang-McFarlane', &
                            horiz_only, 'avg', 'm s-1')
    call history_add_field ('ZMFLXPRC', 'Flux of precipitation from ZM convection', 'ilev', 'avg', 'kg m-2 s-1')
    call history_add_field ('ZMFLXSNW', 'Flux of snow from ZM convection', 'ilev', 'avg', 'kg m-2 s-1')
    call history_add_field ('ZMNTPRPD', 'Net precipitation production from ZM convection', 'lev', 'avg', 'kg kg-1 s-1')
    call history_add_field ('ZMNTSNPD', 'Net snow production from ZM convection', 'lev', 'avg', 'kg kg-1 s-1')
    call history_add_field ('ZMEIHEAT', 'Heating by ice and evaporation in ZM convection', 'lev', 'avg', 'W kg-1')

    call history_add_field ('CMFMC_DP', 'Convection mass flux from ZM deep ', 'ilev', 'avg', 'kg m-2 s-1')
    call history_add_field ('PRECCDZM', 'Convective precipitation rate from ZM deep', horiz_only, 'avg', 'm s-1')

    call history_add_field ('PCONVB', 'convection base pressure',  horiz_only ,  'avg', 'Pa'    )
    call history_add_field ('PCONVT', 'convection top  pressure',  horiz_only ,  'avg', 'Pa'    )

    call history_add_field ('CAPE',   'zhang_mcfarlane_convective_available_potential_energycap',  horiz_only,   'avg', 'J kg-1')
    call history_add_field ('FREQZM', 'Fractional occurance of ZM convection',  horiz_only  , 'avg', 'fraction')

    call history_add_field ('ZMMU',   'ZM convection updraft mass flux',  'lev',  'avg', 'kg m-2 s-1')
    call history_add_field ('ZMMD',   'ZM convection downdraft mass flux', 'lev',  'avg', 'kg m-2 s-1')

    call history_add_field ('ZMUPGU', 'zonal force from ZM updraft pressure gradient term',  'lev',  'avg', 'm s-2')
    call history_add_field ('ZMUPGD', 'zonal force from ZM downdraft pressure gradient term',  'lev',  'avg', 'm s-2')
    call history_add_field ('ZMVPGU', 'meridional force from ZM updraft pressure gradient term', 'lev',  'avg', 'm s-2')
    call history_add_field ('ZMVPGD', 'merdional force from ZM downdraft pressure gradient term', 'lev',  'avg', 'm s-2')

    call history_add_field ('ZMICUU', 'ZM in-cloud U updrafts',  'lev',  'avg', 'm/s')
    call history_add_field ('ZMICUD', 'ZM in-cloud U downdrafts',  'lev',  'avg', 'm/s')
    call history_add_field ('ZMICVU', 'ZM in-cloud V updrafts',  'lev',  'avg', 'm/s')
    call history_add_field ('ZMICVD', 'ZM in-cloud V downdrafts', 'lev',  'avg', 'm/s')

    call history_add_field ('DLFZM',   'Detrained liquid water from ZM convection', 'lev', 'avg','kg kg-1 s-1 ')


   end subroutine zm_diagnostics_init

   !> \section arg_table_zm_diagnostics_run  Argument Table
   !! \htmlinclude zm_diagnostics_run.html
   subroutine zm_diagnostics_run(ncol, pver, pverp, ideep, cpair, heat, prec, cape, errmsg, errflg)

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
      real(kind_phys), intent(in) :: heat(:,:)
      real(kind_phys), intent(in) :: prec(:)
      real(kind_phys), intent(in) :: cape(:)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: lengath  ! number of columns with deep convection

      real(kind_phys) :: ftem(ncol,pver)
      real(kind_phys) :: mcon(ncol,pverp)
      real(kind_phys) :: mconzm(ncol,pverp)

      integer :: index_cldliq

      errmsg = ''
      errflg = 0

      lengath = count(ideep > 0)
      if (lengath > ncol) lengath = ncol  ! should not happen, but force it to not be larger than ncol for safety sake

      call history_out_field('PRECZ', prec)

!      ftem(:,:) = 0._kind_phys
!      ftem(:ncol,:pver) = heat(:ncol,:pver)/cpair
!      call history_out_field('ZMDT', ftem)

!      call history_out_field('CAPE', cape)

!      freqzm(:) = 0._r8
!      do i = 1,lengath
!         freqzm(ideep(i)) = 1.0_r8
!      end do
!      call history_out_field('FREQZM  ',freqzm)

!      mcon(:ncol,:pverp) = mcon(:ncol,:pverp) * 100._r8/gravit
!      mconzm(:ncol,:pverp) = mcon(:ncol,:pverp)

!      call outfld('CMFMC_DP', mconzm)


   end subroutine zm_diagnostics_run

   !=======================================================================

end module zm_diagnostics
