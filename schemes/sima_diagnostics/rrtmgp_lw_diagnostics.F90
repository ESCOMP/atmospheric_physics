module rrtmgp_lw_diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! THIS IS A TEMPLATE
!   1. copy this file to a new file with the correct name
!        (rrtmgp_lw_diagnostics.F90)
!   2. do a search and replace for "rrtmgp_lw" in this file and
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

   public :: rrtmgp_lw_diagnostics_init ! init routine
   public :: rrtmgp_lw_diagnostics_run  ! main routine

   character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ',&
                                      '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

CONTAINS

   !> \section arg_table_rrtmgp_lw_diagnostics_init  Argument Table
   !! \htmlinclude rrtmgp_lw_diagnostics_init.html
   subroutine rrtmgp_lw_diagnostics_init(num_diag_subcycles, active_calls, errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      integer,             intent(in) :: num_diag_subcycles  ! Number of diagnostic subcycles
      logical,             intent(in) :: active_calls(:)     ! Logical array of flags for whether a specified subcycle is active

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:
      integer :: icall

      errmsg = ''
      errflg = 0

      ! Loop over number of diagnostics subcycles
      !  and add the relevant fields for each cycle if it's active
      do icall = 1, num_diag_subcycles
         if (active_calls(icall)) then
            call history_add_field('QRL'//diag(icall-1),     'Longwave heating rate',                            'lev',      'avg',  'K s-1')
            call history_add_field('QRLC'//diag(icall-1),    'Clearsky longwave heating rate',                   'lev',      'avg',  'K s-1')
            call history_add_field('FLNT'//diag(icall-1),    'Net longwave flux at top of model',                horiz_only, 'avg',  'W m-2')
            call history_add_field('FLNTC'//diag(icall-1),   'Clearky net longwave flux at top of model',        horiz_only, 'avg',  'W m-2')
            call history_add_field('FLUT'//diag(icall-1),    'Upwelling longwave flux at top of model',          horiz_only, 'avg',  'W m-2')
            call history_add_field('FLUTC'//diag(icall-1),   'Clearsky upwelling longwave flux at top of model', horiz_only, 'avg',  'W m-2')
            call history_add_field('LWCF'//diag(icall-1),    'Longwave cloud forcing',                           horiz_only, 'avg',  'W m-2')
            call history_add_field('FLN200'//diag(icall-1),  'Net longwave flux at 200 mb',                      horiz_only, 'avg',  'W m-2')
            call history_add_field('FLN200C'//diag(icall-1), 'Clearsky net longwave flux at 200 mb',             horiz_only, 'avg',  'W m-2')
            call history_add_field('FLNR'//diag(icall-1),    'Net longwave flux at tropopause',                  horiz_only, 'avg',  'W m-2')
            call history_add_field('FLNS'//diag(icall-1),    'Net longwave flux at surface',                     horiz_only, 'avg',  'W m-2')
            call history_add_field('FLNSC'//diag(icall-1),   'Clearsky net longwave flux at surface',            horiz_only, 'avg',  'W m-2')
            call history_add_field('FLDS'//diag(icall-1),    'Downwelling longwave flux at surface',             horiz_only, 'avg',  'W m-2')
            call history_add_field('FLDSC'//diag(icall-1),   'Clearky Downwelling longwave flux at surface',     horiz_only, 'avg',  'W m-2')

            ! Fluxes on CAM grid
            call history_add_field('FUL'//diag(icall-1),     'Longwave upward flux',                             'ilev',     'inst', 'W m-2')
            call history_add_field('FDL'//diag(icall-1),     'Longwave downward flux',                           'ilev',     'inst', 'W m-2')
            call history_add_field('FULC'//diag(icall-1),    'Longwave clear-sky upward flux',                   'ilev',     'inst', 'W m-2')
            call history_add_field('FDLC'//diag(icall-1),    'Longwave clear-sky downward flux',                 'ilev',     'inst', 'W m-2')
         end if
      end do

      call history_add_field('EMIS', 'Cloud longwave emissivity',                       'lev', 'avg', '1')

      ! Heating rate needed for d(theta)/dt computation
      call history_add_field('HR',   'Heating rate needed for d(theat)/dt computation', 'lev', 'avg', 'K s-1')

   end subroutine rrtmgp_lw_diagnostics_init

   !> \section arg_table_rrtmgp_lw_diagnostics_run  Argument Table
   !! \htmlinclude rrtmgp_lw_diagnostics_run.html
   subroutine rrtmgp_lw_diagnostics_run(num_diag_subcycles, icall, active_calls, flw, flwc, rpdel, ncol, &
                  nlay, pver, pverp, pint, gravit, p_trop, ktopcam, ktoprad, write_output, errmsg, errflg)

      use cam_history,        only: history_out_field
      use ccpp_fluxes,        only: ty_fluxes_broadband_ccpp
      use ccpp_fluxes_byband, only: ty_fluxes_byband_ccpp
      use interpolate_data,   only: vertinterp
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer,                        intent(in) :: num_diag_subcycles  ! Number of diagnostics subcycles
      integer,                        intent(in) :: icall               ! Current diagnostic subcycle
      integer,                        intent(in) :: ncol                ! Number of horizontal points
      integer,                        intent(in) :: nlay                ! Number of vertical layers used in radiation calcluation
      integer,                        intent(in) :: pver                ! Number of vertical layers
      integer,                        intent(in) :: pverp               ! Number of vertical layer interfaces
      integer,                        intent(in) :: ktopcam             ! Index in CAM arrays of top level (layer or interface) at which RRTMGP is active
      integer,                        intent(in) :: ktoprad             ! Index in RRTMGP array corresponding to top layer or interface of CAM arrays
      logical,                        intent(in) :: active_calls(:)     ! Logical array of flags for whether a specified subcycle is active
      logical,                        intent(in) :: write_output        ! Flag to write output for radiation
      real(kind_phys),                intent(in) :: gravit              ! Standard gravitiational acceleration
      real(kind_phys),                intent(in) :: pint(:,:)           ! Air pressure at layer interfaces [Pa]
      real(kind_phys),                intent(in) :: p_trop(:)           ! Tropopause air pressure [Pa]
      real(kind_phys),                intent(in) :: rpdel(:,:)          ! Reciprocal of layer thickness [Pa-1]
      type(ty_fluxes_byband_ccpp),    intent(in) :: flw                 ! Longwave all-sky flux object
      type(ty_fluxes_broadband_ccpp), intent(in) :: flwc                ! Longwave clear-sky flux object
      
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: diag_index, idx
      real(kind_phys) :: fnl(ncol, pverp)
      real(kind_phys) :: fcnl(ncol, pverp)
      real(kind_phys) :: qrl(ncol, pver)
      real(kind_phys) :: qrlc(ncol, pver)
      real(kind_phys) :: fln200(ncol)
      real(kind_phys) :: fln200c(ncol)
      real(kind_phys) :: flnr(ncol)

      errmsg = ''
      errflg = 0

      ! Diagnostic indices are reversed
      diag_index = num_diag_subcycles - icall

      ! Don't do anything if this subcycle is inactive or we're not configured to write radiation output
      if ((.not. active_calls(diag_index)) .or. (.not. write_output)) then
         return
      end if

      fnl = 0.0_kind_phys
      fcnl = 0.0_kind_phys

      ! RTE-RRTMGP convention for net is (down - up) **CAM assumes (up - down) !!
      fnl( :,ktopcam:) = -1._kind_phys * flw%fluxes%flux_net(    :, ktoprad:)
      fcnl(:,ktopcam:) = -1._kind_phys * flwc%fluxes%flux_net(   :, ktoprad:)

      call heating_rate(ncol, ktopcam, pver, fnl, gravit, rpdel, qrl)
      call heating_rate(ncol, ktopcam, pver, fcnl, gravit, rpdel, qrlc)

      ! History out field calls
      call history_out_field('QRL'//diag(icall),     qrl(:ncol,:)/cpair)
      call history_out_field('QRLC'//diag(icall),    qrlc(:ncol,:)/cpair)

      call history_out_field('FLNT'//diag(icall),    fnl(:,ktopcam))
      call history_out_field('FLNTC'//diag(icall),   fcnl(:,ktopcam))

      call history_out_field('FLUT'//diag(icall),    flw%fluxes%flux_up(:, ktoprad))
      call history_out_field('FLUTC'//diag(icall),   flwc%fluxes%flux_up(:, ktoprad))

      ftem(:) = flwc%fluxes%flux_up(:, ktoprad) - flw%fluxes%flux_up(:, ktoprad)
      call history_out_field('LWCF'//diag(icall),    ftem)

      ! Output fluxes at 200 mb
      call vertinterp(ncol, ncol, pverp, pint, 20000._r8, fnl,  fln200)
      call vertinterp(ncol, ncol, pverp, pint, 20000._r8, fcnl, fln200c)
      call history_out_field('FLN200'//diag(icall),  fln200)
      call history_out_field('FLN200C'//diag(icall), fln200c)

      do idx = 1,ncol
         call vertinterp(1, 1, pverp, pint(idx,:), p_trop(idx), fnl(idx,:), flnr(idx))
      end do
      call history_out_field('FLNR'//diag(icall),    flnr)

      call history_out_field('FLNS'//diag(icall),    fnl(:,pverp))
      call history_out_field('FLNSC'//diag(icall),   fcnl(:,pverp))

      call history_out_field('FLDS'//diag(icall),    flw%fluxes%flux_dn(:, nlay+1))
      call history_out_field('FLDSC'//diag(icall),   flwc%fluxes%flux_dn(:, nlay+1))

      ! Fluxes on the CAM grid
      call history_out_field('FDL'//diag(icall),     flw%fluxes%flux_dn( :, ktoprad:))
      call history_out_field('FDLC'//diag(icall),    flwc%fluxes%flux_dn(:, ktoprad:))
      call history_out_field('FUL'//diag(icall),     flw%fluxes%flux_up( :, ktoprad:))
      call history_out_field('FULC'//diag(icall),    flwc%fluxes%flux_up(:, ktoprad:))

   end subroutine rrtmgp_lw_diagnostics_run

   !=======================================================================

   subroutine lw_heating_rate(ncol, ktopcam, pver, flux_net, gravit, rpdel, hrate)
      ! Compute heating rate as a dry static energy tendency

      ! arguments
      integer,          intent(in) :: ncol
      real(kind_phys),  intent(in) :: flux_net(:,:)   ! W m-2
      real(kind_phys),  intent(in) :: gravit          ! m s-2
      real(kind_phys),  intent(in) :: rpdel(:,:)      ! Pa
      real(kind_phys), intent(out) :: hrate(:,:)      ! J kg-1 s-1

      ! local vars
      integer :: kdx

      do kdx = ktopcam, pver
         ! (flux divergence as bottom-MINUS-top) * g/dp
         hrate(:,kdx) = (flux_net(:,kdx+1) - flux_net(:,kdx)) * &
                       gravit * state%rpdel(:,kdx)
      end do
   end subroutine lw_heating_rate

end module rrtmgp_lw_diagnostics
