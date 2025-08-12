module rrtmgp_sw_diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! THIS IS A TEMPLATE
!   1. copy this file to a new file with the correct name
!        (rrtmgp_sw_diagnostics.F90)
!   2. do a search and replace for "rrtmgp_sw" in this file and
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

   public :: rrtmgp_sw_diagnostics_init ! init routine
   public :: rrtmgp_sw_diagnostics_run  ! main routine

   integer, parameter :: N_DIAG=10
   character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ',&
                                      '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

CONTAINS

   !> \section arg_table_rrtmgp_sw_diagnostics_init  Argument Table
   !! \htmlinclude rrtmgp_sw_diagnostics_init.html
   subroutine rrtmgp_sw_diagnostics_init(num_diag_subcycles, active_calls, errmsg, errflg)
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
            call history_add_field('QRS'//diag(icall-1),     'Shortwave heating rate',                            'lev',      'avg',  'K s-1')
            call history_add_field('QRSC'//diag(icall-1),    'Clearsky shortwave heating rate',                   'lev',      'avg',  'K s-1')
            call history_add_field('FSNT'//diag(icall-1),    'Net shortwave flux at top of model',                horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNTC'//diag(icall-1),   'Clearky net shortwave flux at top of model',        horiz_only, 'avg',  'W m-2')
            call history_add_field('FSUT'//diag(icall-1),    'Upwelling shortwave flux at top of model',          horiz_only, 'avg',  'W m-2')
            call history_add_field('FSUTC'//diag(icall-1),   'Clearsky upwelling shortwave flux at top of model', horiz_only, 'avg',  'W m-2')
            call history_add_field('SWCF'//diag(icall-1),    'Shortwave cloud forcing',                           horiz_only, 'avg',  'W m-2')
            call history_add_field('FSN200'//diag(icall-1),  'Net shortwave flux at 200 mb',                      horiz_only, 'avg',  'W m-2')
            call history_add_field('FSN200C'//diag(icall-1), 'Clearsky net shortwave flux at 200 mb',             horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNR'//diag(icall-1),    'Net shortwave flux at tropopause',                  horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNS'//diag(icall-1),    'Net shortwave flux at surface',                     horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNSC'//diag(icall-1),   'Clearsky net shortwave flux at surface',            horiz_only, 'avg',  'W m-2')
            call history_add_field('FSDS'//diag(icall-1),    'Downwelling shortwave flux at surface',             horiz_only, 'avg',  'W m-2')
            call history_add_field('FSDSC'//diag(icall-1),   'Clearky Downwelling shortwave flux at surface',     horiz_only, 'avg',  'W m-2')

            ! Fluxes on CAM grid
            call history_add_field('FUS'//diag(icall-1),     'Shortwave upward flux',                             'ilev',     'inst', 'W m-2')
            call history_add_field('FDS'//diag(icall-1),     'Shortwave downward flux',                           'ilev',     'inst', 'W m-2')
            call history_add_field('FUSC'//diag(icall-1),    'Shortwave clear-sky upward flux',                   'ilev',     'inst', 'W m-2')
            call history_add_field('FDSC'//diag(icall-1),    'Shortwave clear-sky downward flux',                 'ilev',     'inst', 'W m-2')
         end if
      end do

   end subroutine rrtmgp_sw_diagnostics_init

   !> \section arg_table_rrtmgp_sw_diagnostics_run  Argument Table
   !! \htmlinclude rrtmgp_sw_diagnostics_run.html
   subroutine rrtmgp_sw_diagnostics_run(num_diag_subcycles, icall, active_calls, fsw, fswc, rpdel, ncol, &
                  nlay, pver, pverp, pint, gravit, cpair, p_trop, ktopcam, ktoprad, write_output, errmsg, errflg)

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
      real(kind_phys),                intent(in) :: cpair
      real(kind_phys),                intent(in) :: pint(:,:)           ! Air pressure at layer interfaces [Pa]
      real(kind_phys),                intent(in) :: p_trop(:)           ! Tropopause air pressure [Pa]
      real(kind_phys),                intent(in) :: rpdel(:,:)          ! Reciprocal of layer thickness [Pa-1]
      type(ty_fluxes_byband_ccpp),    intent(in) :: fsw                 ! Shortwave all-sky flux object
      type(ty_fluxes_broadband_ccpp), intent(in) :: fswc                ! Shortwave clear-sky flux object
      
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
      real(kind_phys) :: ftem(ncol)

      errmsg = ''
      errflg = 0

      ! Diagnostic indices are reversed
      diag_index = num_diag_subcycles - icall

      ! Don't do anything if this subcycle is inactive or we're not configured to write radiation output
      if ((.not. active_calls(diag_index+1)) .or. (.not. write_output)) then
         return
      end if

      fnl = 0.0_kind_phys
      fcnl = 0.0_kind_phys

      ! RTE-RRTMGP convention for net is (down - up) **CAM assumes (up - down) !!
      fnl( :,ktopcam:) = -1._kind_phys * fsw%fluxes%flux_net(    :, ktoprad:)
      fcnl(:,ktopcam:) = -1._kind_phys * fswc%fluxes%flux_net(   :, ktoprad:)

      call sw_heating_rate(ncol, ktopcam, pver, fnl, gravit, rpdel, qrl)
      call sw_heating_rate(ncol, ktopcam, pver, fcnl, gravit, rpdel, qrlc)

      ! History out field calls
      call history_out_field('QRS'//diag(diag_index),     qrl(:,:)/cpair)
      call history_out_field('QRSC'//diag(diag_index),    qrlc(:,:)/cpair)

      call history_out_field('FSNT'//diag(diag_index),    fnl(:,ktopcam))
      call history_out_field('FSNTC'//diag(diag_index),   fcnl(:,ktopcam))

      call history_out_field('FSUT'//diag(diag_index),    fsw%fluxes%flux_up(:, ktoprad))
      call history_out_field('FSUTC'//diag(diag_index),   fswc%fluxes%flux_up(:, ktoprad))

      ftem(:) = fswc%fluxes%flux_up(:, ktoprad) - fsw%fluxes%flux_up(:, ktoprad)
      call history_out_field('SWCF'//diag(diag_index),    ftem)

      ! Output fluxes at 200 mb
      call vertinterp(ncol, ncol, pverp, pint, 20000._kind_phys, fnl,  fln200)
      call vertinterp(ncol, ncol, pverp, pint, 20000._kind_phys, fcnl, fln200c)
      call history_out_field('FSN200'//diag(diag_index),  fln200)
      call history_out_field('FSN200C'//diag(diag_index), fln200c)

      do idx = 1,ncol
         call vertinterp(1, 1, pverp, pint(idx,:), p_trop(idx), fnl(idx,:), flnr(idx))
      end do
      call history_out_field('FSNR'//diag(diag_index),    flnr)

      call history_out_field('FSNS'//diag(diag_index),    fnl(:,pverp))
      call history_out_field('FSNSC'//diag(diag_index),   fcnl(:,pverp))

      call history_out_field('FSDS'//diag(diag_index),    fsw%fluxes%flux_dn(:, nlay+1))
      call history_out_field('FSDSC'//diag(diag_index),   fswc%fluxes%flux_dn(:, nlay+1))

      ! Fluxes on the CAM grid
      call history_out_field('FDS'//diag(diag_index),     fsw%fluxes%flux_dn( :, ktoprad:))
      call history_out_field('FDSC'//diag(diag_index),    fswc%fluxes%flux_dn(:, ktoprad:))
      call history_out_field('FUS'//diag(diag_index),     fsw%fluxes%flux_up( :, ktoprad:))
      call history_out_field('FUSC'//diag(diag_index),    fswc%fluxes%flux_up(:, ktoprad:))

   end subroutine rrtmgp_sw_diagnostics_run

   !=======================================================================

   subroutine sw_heating_rate(ncol, ktopcam, pver, flux_net, gravit, rpdel, hrate)
      ! Compute heating rate as a dry static energy tendency

      ! arguments
      integer,          intent(in) :: ncol
      integer,          intent(in) :: ktopcam
      integer,          intent(in) :: pver
      real(kind_phys),  intent(in) :: flux_net(:,:)   ! W m-2
      real(kind_phys),  intent(in) :: gravit          ! m s-2
      real(kind_phys),  intent(in) :: rpdel(:,:)      ! Pa
      real(kind_phys), intent(out) :: hrate(:,:)      ! J kg-1 s-1

      ! local vars
      integer :: kdx

      do kdx = ktopcam, pver
         ! (flux divergence as bottom-MINUS-top) * g/dp
         hrate(:,kdx) = (flux_net(:,kdx+1) - flux_net(:,kdx)) * &
                       gravit * rpdel(:,kdx)
      end do
   end subroutine sw_heating_rate

end module rrtmgp_sw_diagnostics
