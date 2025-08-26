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
            call history_add_field('SOLIN'//diag(icall-1),   'Solar isolation',                                   horiz_only, 'avg',  'W m-2')
            call history_add_field('QRS'//diag(icall-1),     'Solar heating rate',                                'lev',      'avg',  'K s-1')
            call history_add_field('QRSC'//diag(icall-1),    'Clearsky solar heating rate',                       'lev',      'avg',  'K s-1')
            call history_add_field('FSNT'//diag(icall-1),    'Net solar flux at top of model',                    horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNTC'//diag(icall-1),   'Clearky net solar flux at top of model',            horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNTOA'//diag(icall-1),  'Net solar flux at top of atmosphere',               horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNTOAC'//diag(icall-1), 'Clearsky net solar flux at top of atmosphere',      horiz_only, 'avg',  'W m-2')
            call history_add_field('SWCF'//diag(icall-1),    'Shortwave cloud forcing',                           horiz_only, 'avg',  'W m-2')
            call history_add_field('FSUTOA'//diag(icall-1),  'Upwelling solar flux at top of atmospehre',         horiz_only, 'avg',  'W m-2')
            call history_add_field('FSN200'//diag(icall-1),  'Net shortwave flux at 200 mb',                      horiz_only, 'avg',  'W m-2')
            call history_add_field('FSN200C'//diag(icall-1), 'Clearsky net shortwave flux at 200 mb',             horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNR'//diag(icall-1),    'Net solar flux at tropopause',                      horiz_only, 'avg',  'W m-2')
            call history_add_field('SOLL'//diag(icall-1),    'Solar downward near infrared direct to surface',    horiz_only, 'avg',  'W m-2')
            call history_add_field('SOLS'//diag(icall-1),    'Solar downward visible direct to surface',          horiz_only, 'avg',  'W m-2')
            call history_add_field('SOLLD'//diag(icall-1),   'Solar downward near infrared diffuse to surface',   horiz_only, 'avg',  'W m-2')
            call history_add_field('SOLSD'//diag(icall-1),   'Solar downward visible diffuse to surface',         horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNS'//diag(icall-1),    'Net solar flux at surface',                         horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNSC'//diag(icall-1),   'Clearsky net solar flux at surface',                horiz_only, 'avg',  'W m-2')
            call history_add_field('FSDS'//diag(icall-1),    'Downwelling solar flux at surface',                 horiz_only, 'avg',  'W m-2')
            call history_add_field('FSDSC'//diag(icall-1),   'Clearky downwelling solar flux at surface',         horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNIRTOA'//diag(icall-1),'Net near-infrared flux (Nimbus-7 WFOV at top of atmosphere',          horiz_only, 'avg', 'W m-2')
            call history_add_field('FSNRTOAC'//diag(icall-1),'Clearsky net near-infrared flux (Nimbus-7 WFOV at top of atmosphere', horiz_only, 'avg', 'W m-2')
            call history_add_field('FSNRTOAS'//diag(icall-1),'Net near-infrared flux (>= 0.7 microns) at top of atmosphere',        horiz_only, 'avg', 'W m-2')

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
   subroutine rrtmgp_sw_diagnostics_run(num_diag_subcycles, icall, active_calls, fsw, fswc, rpdel, ncol, nday, idxday, &
                  nlay, pver, pverp, pint, gravit, cpair, p_trop, fns, fcns, qrs, qrsc, fsnt, fsns, sols, soll, solsd, &
                  solld, ktopcam, ktoprad, write_output, errmsg, errflg)

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
      integer,                        intent(in) :: nday                ! Number of daytime points
      integer,                        intent(in) :: idxday(:)           ! Daytime points
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
      real(kind_phys),                intent(in) :: fns(:,:)            ! Net shortwave all-sky flux [W m-2]
      real(kind_phys),                intent(in) :: fcns(:,:)           ! Net shortwave clear-sky flux [W m-2]
      real(kind_phys),                intent(in) :: qrs(:,:)            ! Heating rate (all-sky, shortwave) [J kg-1 s-1]
      real(kind_phys),                intent(in) :: qrsc(:,:)           ! Heating rate (clear-sky, shortwave) [J kg-1 s-1]
      real(kind_phys),                intent(in) :: fsnt(:)             ! Shortwave outgoing flux at model top [W m-2]
      real(kind_phys),                intent(in) :: fsns(:)             ! Shortwave upward flux at surface [W m-2]
      real(kind_phys),                intent(in) :: sols(:)
      real(kind_phys),                intent(in) :: soll(:)
      real(kind_phys),                intent(in) :: solsd(:)
      real(kind_phys),                intent(in) :: solld(:)
      type(ty_fluxes_byband_ccpp),    intent(in) :: fsw                 ! Shortwave all-sky flux object
      type(ty_fluxes_broadband_ccpp), intent(in) :: fswc                ! Shortwave clear-sky flux object
      
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: diag_index, idx
      real(kind_phys) :: solin(ncol)
      real(kind_phys) :: fcns(ncol)
      real(kind_phys) :: fsntoa(ncol)
      real(kind_phys) :: fsntoac(ncol)
      real(kind_phys) :: fsutoa(ncol)
      real(kind_phys) :: fsdsc(ncol)
      real(kind_phys) :: flux_sw_up(ncol,pver)
      real(kind_phys) :: flux_sw_dn(ncol,pver)
      real(kind_phys) :: flux_sw_clr_up(ncol,pver)
      real(kind_phys) :: flux_sw_clr_dn(ncol,pver)
      real(kind_phys) :: fsntc(ncol)
      real(kind_phys) :: fsnsc(ncol)

      errmsg = ''
      errflg = 0

      ! Diagnostic indices are reversed
      diag_index = num_diag_subcycles - icall

      ! Don't do anything if this subcycle is inactive or we're not configured to write radiation output
      if ((.not. active_calls(diag_index+1)) .or. (.not. write_output)) then
         return
      end if

      ! Initialize to provide 0.0 values for night columns.
      solin = 0.0_kind_phys
      fcns = 0.0_kind_phys
      fsntoa = 0.0_kind_phys
      fsntoac = 0.0_kind_phys
      fsutoa = 0.0_kind_phys
      fsdsc = 0.0_kind_phys
      flux_sw_up = 0.0_kind_phys
      flux_sw_dn = 0.0_kind_phys
      flux_sw_clr_up = 0.0_kind_phys
      flux_sw_clr_dn = 0.0_kind_phys

      ! Load up diagnostic arrays
      do idx = 1, nday
         solin(idxday(idx) = fswc%fluxes%flux_dn(idx, 1)
         fsntoa(idxday(idx)) = fsw%fluxes%flux_net(idx, 1)
         fsntoac(idxday(idx)) = fswc%fluxes%flux_net(idx, 1)
         fsutoa(idxday(idx)) = fsw%fluxes%flux_up(idx, 1)
         fsdsc(idxday(idx))  = fswc%fluxes%flux_dn(idx, nlay+1)
         flux_sw_up(idxday(idx),ktopcam:) = fsw%fluxes%flux_up(idx,ktoprad:)
         flux_sw_dn(idxday(idx),ktopcam:) = fsw%fluxes%flux_dn(idx,ktoprad:)
         flux_sw_clr_up(idxday(idx),ktopcam:) = fswc%fluxes%flux_up(idx,ktoprad:)
         flux_sw_clr_dn(idxday(idx),ktopcam:) = fswc%fluxes%flux_dn(idx,ktoprad:)
      end do

      fsntc(:) = fcns(:, pverp)   ! net sw clearsky flux at top
      fsnsc(:) = fcns(:, pverp)   ! net sw clearsky flux at surface

      ! History out field calls
      call history_out_field('SOLIN'//diag(diag_index),   solin)
      call history_out_field('QRS'//diag(diag_index),     qrs(:,:)/cpair)
      call history_out_field('QRSC'//diag(diag_index),    qrsc(:,:)/cpair)

      call history_out_field('FSNT'//diag(diag_index),    fsnt)
      call history_out_field('FSNTC'//diag(diag_index),   fsntc)
      call history_out_field('FSNTOA'//diag(diag_index),  fsntoa)
      call history_out_field('FSNTOAC'//diag(diag_index), fsntoac)

      call history_out_field('SWCF'//diag(diag_index),    fsntoa - fsntoac)

      call history_out_field('FSUTOA'//diag(diag_index),  fsutoa)

      ! Output fluxes at 200 mb
      call vertinterp(ncol, ncol, pverp, pint, 20000._kind_phys, fns,  fsn200)
      call vertinterp(ncol, ncol, pverp, pint, 20000._kind_phys, fcns, fsn200c)
      call history_out_field('FSN200'//diag(diag_index),  fsn200)
      call history_out_field('FSN200C'//diag(diag_index), fsn200c)

      do idx = 1,ncol
         call vertinterp(1, 1, pverp, pint(idx,:), p_trop(idx), fns(idx,:), fsnr(idx))
      end do
      call history_out_field('FSNR'//diag(diag_index),    fsnr)

      call history_out_field('SOLS'//diag(diag_index),    sols)
      call history_out_field('SOLL'//diag(diag_index),    soll)
      call history_out_field('SOLSD'//diag(diag_index),   solsd)
      call history_out_field('SOLLD'//diag(diag_index),   solld)

      call history_out_field('FSNS'//diag(diag_index),    fsns)
      call history_out_field('FSNSC'//diag(diag_index),   fsnsc)

      call history_out_field('FSDS'//diag(diag_index),    fsds)
      call history_out_field('FSDSC'//diag(diag_index),   fsdsc)

      ! Fluxes on the CAM grid
      call history_out_field('FDS'//diag(diag_index),     flux_sw_dn)
      call history_out_field('FDSC'//diag(diag_index),    flux_sw_clr_dn)
      call history_out_field('FUS'//diag(diag_index),     flux_sw_up)
      call history_out_field('FUSC'//diag(diag_index),    flux_sw_clr_up)

   end subroutine rrtmgp_sw_diagnostics_run

end module rrtmgp_sw_diagnostics
