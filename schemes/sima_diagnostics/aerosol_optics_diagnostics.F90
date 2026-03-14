! Aerosol optics diagnostics scheme for CAM-SIMA.
! Outputs aggregate SW/LW optical property diagnostics and per-species
! AODs/burdens from the aerosol_optics CCPP scheme.
!
! Follows rrtmgp_lw_diagnostics.F90 pattern. Day-only fields use fillvalue
! for night columns; day+night (dn) fields carry data for all columns.
!
! Author: Haipeng Lin, NSF-NCAR/CGD/AMP, March 2026
module aerosol_optics_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: aerosol_optics_diagnostics_init
   public :: aerosol_optics_diagnostics_run

   ! ODV_ aerosol diagnostics for bulk aerosol
   integer :: num_odv = 0
   character(len=64), allocatable :: odv_names(:)

contains

   !> \section arg_table_aerosol_optics_diagnostics_init  Argument Table
   !! \htmlinclude aerosol_optics_diagnostics_init.html
   subroutine aerosol_optics_diagnostics_init(errmsg, errflg)
      use cam_history,                  only: history_add_field
      use cam_history_support,          only: horiz_only
      use radiative_aerosol_definitions, only: bulk_aerosol_list

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: i

      errmsg = ''
      errflg = 0

      !-----------------------------------------------------------------
      ! SW column diagnostics (horiz_only, avg)
      !-----------------------------------------------------------------
      call history_add_field('AODVIS',      'Aerosol optical depth at 550 nm',                    horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODVISdn',    'Aerosol optical depth at 550 nm, day and night',      horiz_only, 'avg', '1')
      call history_add_field('AODUV',       'Aerosol optical depth in UV band',                    horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODUVdn',     'Aerosol optical depth in UV band, day and night',     horiz_only, 'avg', '1')
      call history_add_field('AODNIR',      'Aerosol optical depth in near-IR band',               horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODNIRdn',    'Aerosol optical depth in near-IR band, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODABS',      'Aerosol absorption optical depth at 550 nm',          horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODABSdn',    'Aerosol absorption optical depth at 550 nm, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODTOT',      'Total aerosol optical depth (all SW bands)',          horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODTOTdn',    'Total aerosol optical depth (all SW bands), day and night', horiz_only, 'avg', '1')
      call history_add_field('SSAVIS',      'Aerosol single scatter albedo at 550 nm',             horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('SSAVISdn',    'Aerosol single scatter albedo at 550 nm, day and night', horiz_only, 'avg', '1')

      ! Stratospheric column AODs
      call history_add_field('AODVISst',    'Stratospheric aerosol optical depth at 550 nm',       horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODVISstdn',  'Stratospheric aerosol optical depth at 550 nm, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODUVst',     'Stratospheric aerosol optical depth in UV band',      horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODUVstdn',   'Stratospheric aerosol optical depth in UV band, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODNIRst',    'Stratospheric aerosol optical depth in near-IR band', horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODNIRstdn',  'Stratospheric aerosol optical depth in near-IR band, day and night', horiz_only, 'avg', '1')

      !-----------------------------------------------------------------
      ! SW per-level diagnostics (lev, avg)
      !-----------------------------------------------------------------
      call history_add_field('EXTINCT',      'Aerosol extinction at 550 nm',                       'lev', 'avg', 'm-1', flag_xyfill=.true.)
      call history_add_field('EXTINCTdn',    'Aerosol extinction at 550 nm, day and night',         'lev', 'avg', 'm-1')
      call history_add_field('EXTINCTUV',    'Aerosol extinction in UV band',                       'lev', 'avg', 'm-1', flag_xyfill=.true.)
      call history_add_field('EXTINCTUVdn',  'Aerosol extinction in UV band, day and night',        'lev', 'avg', 'm-1')
      call history_add_field('EXTINCTNIR',   'Aerosol extinction in near-IR band',                  'lev', 'avg', 'm-1', flag_xyfill=.true.)
      call history_add_field('EXTINCTNIRdn', 'Aerosol extinction in near-IR band, day and night',   'lev', 'avg', 'm-1')
      call history_add_field('ABSORB',       'Aerosol absorption at 550 nm',                        'lev', 'avg', 'm-1', flag_xyfill=.true.)
      call history_add_field('ABSORBdn',     'Aerosol absorption at 550 nm, day and night',         'lev', 'avg', 'm-1')

      !-----------------------------------------------------------------
      ! Per-species column diagnostics (horiz_only, avg)
      !-----------------------------------------------------------------
      call history_add_field('AODDUST',        'Dust aerosol optical depth',          horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODDUSTdn',      'Dust aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODSO4',         'Sulfate aerosol optical depth',        horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODSO4dn',       'Sulfate aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODBC',          'Black carbon aerosol optical depth',   horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODBCdn',        'Black carbon aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODPOM',         'POM aerosol optical depth',            horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODPOMdn',       'POM aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODSOA',         'SOA aerosol optical depth',            horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODSOAdn',       'SOA aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODSS',          'Seasalt aerosol optical depth',        horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODSSdn',        'Seasalt aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODABSBC',       'Black carbon absorption optical depth', horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODABSBCdn',     'Black carbon absorption optical depth, day and night', horiz_only, 'avg', '1')

      call history_add_field('BURDENDUST',     'Dust column burden',                   horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENDUSTdn',   'Dust column burden, day and night',    horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENSO4',      'Sulfate column burden',                horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENSO4dn',    'Sulfate column burden, day and night',  horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENBC',       'Black carbon column burden',           horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENBCdn',     'Black carbon column burden, day and night', horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENPOM',      'POM column burden',                    horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENPOMdn',    'POM column burden, day and night',      horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENSOA',      'SOA column burden',                    horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENSOAdn',    'SOA column burden, day and night',      horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENSEASALT',  'Seasalt column burden',                horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENSEASALTdn','Seasalt column burden, day and night',  horiz_only, 'avg', 'kg m-2')

      !-----------------------------------------------------------------
      ! LW per-level diagnostics (lev, avg)
      !-----------------------------------------------------------------
      call history_add_field('TOTABSLW',  'Total LW aerosol absorption optical depth', 'lev', 'avg', '1')
      call history_add_field('AODABSLW',  'LW aerosol absorption optical depth at diagnostic band', 'lev', 'avg', '1')

      !-----------------------------------------------------------------
      ! Additional CAM-standard column diagnostics (aliases)
      !-----------------------------------------------------------------
      call history_add_field('AEROD_v',  'Total aerosol optical depth in visible band', horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODvstrt', 'Stratospheric aerosol optical depth in visible band', horiz_only, 'avg', '1', flag_xyfill=.true.)

      !-----------------------------------------------------------------
      ! Per-aerosol visible OD diagnostics (ODV_<name>) for bulk aerosols
      !-----------------------------------------------------------------
      num_odv = bulk_aerosol_list(0)%numaerosols
      if (num_odv > 0) then
        allocate(odv_names(num_odv))
        do i = 1, num_odv
          odv_names(i) = 'ODV_' // trim(bulk_aerosol_list(0)%aer(i)%camname)
          call history_add_field(trim(odv_names(i)), &
               trim(bulk_aerosol_list(0)%aer(i)%camname)//' optical depth in visible band', &
               horiz_only, 'avg', '1', flag_xyfill=.true.)
        end do
      end if

   end subroutine aerosol_optics_diagnostics_init

   !> \section arg_table_aerosol_optics_diagnostics_run  Argument Table
   !! \htmlinclude aerosol_optics_diagnostics_run.html
   subroutine aerosol_optics_diagnostics_run( &
      ncol, pver, nswbands, nlwbands, &
      aer_tau, aer_tau_w, aer_lw_abs, &
      idx_sw_diag, idx_uv_diag, idx_nir_diag, idx_lw_diag, &
      pdeldry, pmid, t, rga, rair, &
      nnite, idxnite, troplev, &
      dustaod, sulfaod, bcaod, pomaod, soaaod, ssltaod, aodabsbc, &
      burdendust, burdenso4, burdenbc, burdenpom, burdensoa, burdenseasalt, &
      ssavis, aodvis, &
      num_bulk_aer, odv_col_aod, &
      errmsg, errflg)

      use cam_history,         only: history_out_field
      use cam_history_support, only: fillvalue

      !-----------------------------------------------------------------
      ! Input arguments
      !-----------------------------------------------------------------
      integer,         intent(in) :: ncol
      integer,         intent(in) :: pver
      integer,         intent(in) :: nswbands
      integer,         intent(in) :: nlwbands
      real(kind_phys), intent(in) :: aer_tau(:,:,:)        ! SW extinction OD (ncol,pver,nswbands)
      real(kind_phys), intent(in) :: aer_tau_w(:,:,:)      ! SW ssa*tau (ncol,pver,nswbands)
      real(kind_phys), intent(in) :: aer_lw_abs(:,:,:)     ! LW absorption OD (ncol,pver,nlwbands)
      integer,         intent(in) :: idx_sw_diag           ! Index of SW diagnostic (vis) band
      integer,         intent(in) :: idx_uv_diag           ! Index of UV band
      integer,         intent(in) :: idx_nir_diag          ! Index of near-IR band
      integer,         intent(in) :: idx_lw_diag           ! Index of LW diagnostic band
      real(kind_phys), intent(in) :: pdeldry(:,:)          ! Dry air pressure thickness (ncol,pver) [Pa]
      real(kind_phys), intent(in) :: pmid(:,:)             ! Air pressure (ncol,pver) [Pa]
      real(kind_phys), intent(in) :: t(:,:)                ! Air temperature (ncol,pver) [K]
      real(kind_phys), intent(in) :: rga                   ! 1/g [s2 m-1]
      real(kind_phys), intent(in) :: rair                  ! Gas constant of dry air [J kg-1 K-1]
      integer,         intent(in) :: nnite                 ! Number of night columns
      integer,         intent(in) :: idxnite(:)            ! Indices of night columns
      integer,         intent(in) :: troplev(:)            ! Tropopause level index (ncol)

      ! Per-species inputs from aerosol_optics
      real(kind_phys), intent(in) :: dustaod(:)
      real(kind_phys), intent(in) :: sulfaod(:)
      real(kind_phys), intent(in) :: bcaod(:)
      real(kind_phys), intent(in) :: pomaod(:)
      real(kind_phys), intent(in) :: soaaod(:)
      real(kind_phys), intent(in) :: ssltaod(:)
      real(kind_phys), intent(in) :: aodabsbc(:)
      real(kind_phys), intent(in) :: burdendust(:)
      real(kind_phys), intent(in) :: burdenso4(:)
      real(kind_phys), intent(in) :: burdenbc(:)
      real(kind_phys), intent(in) :: burdenpom(:)
      real(kind_phys), intent(in) :: burdensoa(:)
      real(kind_phys), intent(in) :: burdenseasalt(:)
      real(kind_phys), intent(in) :: ssavis(:)
      real(kind_phys), intent(in) :: aodvis(:)

      ! Per-constituent visible OD inputs from aerosol_optics
      integer,         intent(in) :: num_bulk_aer
      real(kind_phys), intent(in) :: odv_col_aod(:,:)      ! (ncol, num_bulk_aer)

      ! CCPP error handling
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      integer :: icol, i, iaer

      ! Derived column quantities
      real(kind_phys) :: aodvis_l(ncol)
      real(kind_phys) :: aoduv(ncol)
      real(kind_phys) :: aodnir(ncol)
      real(kind_phys) :: aodabs(ncol)
      real(kind_phys) :: aodtot(ncol)
      real(kind_phys) :: ssavis_l(ncol)

      ! Stratospheric column AODs
      real(kind_phys) :: aodvisst(ncol)
      real(kind_phys) :: aoduvst(ncol)
      real(kind_phys) :: aodnirst(ncol)

      ! Per-level quantities
      real(kind_phys) :: mass(ncol, pver)
      real(kind_phys) :: air_density(ncol, pver)
      real(kind_phys) :: extinct(ncol, pver)
      real(kind_phys) :: extinctuv(ncol, pver)
      real(kind_phys) :: extinctnir(ncol, pver)
      real(kind_phys) :: absorb(ncol, pver)

      ! LW diagnostics
      real(kind_phys) :: totabslw(ncol, pver)
      real(kind_phys) :: aodabslw(ncol, pver)

      ! Local copy for ODV_ output with fillvalue
      real(kind_phys) :: odv_tmp(ncol)

      ! Local copies of per-species fields (for fillvalue modification)
      real(kind_phys) :: dustaod_l(ncol), sulfaod_l(ncol), bcaod_l(ncol)
      real(kind_phys) :: pomaod_l(ncol), soaaod_l(ncol), ssltaod_l(ncol)
      real(kind_phys) :: aodabsbc_l(ncol)
      real(kind_phys) :: burdendust_l(ncol), burdenso4_l(ncol), burdenbc_l(ncol)
      real(kind_phys) :: burdenpom_l(ncol), burdensoa_l(ncol), burdenseasalt_l(ncol)

      errmsg = ''
      errflg = 0

      !-----------------------------------------------------------------
      ! Derived quantities
      !-----------------------------------------------------------------
      mass(:ncol, :) = pdeldry(:ncol, :) * rga
      air_density(:ncol, :) = pmid(:ncol, :) / (rair * t(:ncol, :))

      !-----------------------------------------------------------------
      ! Column AODs (sum over levels)
      !-----------------------------------------------------------------
      ! aodvis and ssavis are already computed in aerosol_optics_run
      aodvis_l(:ncol) = aodvis(:ncol)
      ssavis_l(:ncol) = ssavis(:ncol)

      do icol = 1, ncol
        aoduv(icol)     = sum(aer_tau(icol, :, idx_uv_diag))
        aodnir(icol)    = sum(aer_tau(icol, :, idx_nir_diag))
        aodabs(icol)    = sum(aer_tau(icol, :, idx_sw_diag) - aer_tau_w(icol, :, idx_sw_diag))
        aodtot(icol)    = sum(aer_tau(icol, :, :))
      end do

      !-----------------------------------------------------------------
      ! Stratospheric column AODs (levels above tropopause)
      !-----------------------------------------------------------------
      do icol = 1, ncol
        aodvisst(icol) = sum(aer_tau(icol, 1:troplev(icol), idx_sw_diag))
        aoduvst(icol)  = sum(aer_tau(icol, 1:troplev(icol), idx_uv_diag))
        aodnirst(icol) = sum(aer_tau(icol, 1:troplev(icol), idx_nir_diag))
      end do

      !-----------------------------------------------------------------
      ! Per-level extinction [m-1] and absorption [m-1]
      !-----------------------------------------------------------------
      extinct(:ncol, :)    = aer_tau(:ncol, :, idx_sw_diag) * air_density(:ncol, :) / mass(:ncol, :)
      extinctuv(:ncol, :)  = aer_tau(:ncol, :, idx_uv_diag) * air_density(:ncol, :) / mass(:ncol, :)
      extinctnir(:ncol, :) = aer_tau(:ncol, :, idx_nir_diag) * air_density(:ncol, :) / mass(:ncol, :)
      absorb(:ncol, :)     = (aer_tau(:ncol, :, idx_sw_diag) - aer_tau_w(:ncol, :, idx_sw_diag)) &
                             * air_density(:ncol, :) / mass(:ncol, :)

      !-----------------------------------------------------------------
      ! LW diagnostics
      !-----------------------------------------------------------------
      totabslw(:ncol, :) = sum(aer_lw_abs(:ncol, :, :), dim=3)
      aodabslw(:ncol, :) = aer_lw_abs(:ncol, :, idx_lw_diag)

      !-----------------------------------------------------------------
      ! Copy per-species fields for fillvalue modification
      !-----------------------------------------------------------------
      dustaod_l(:ncol)       = dustaod(:ncol)
      sulfaod_l(:ncol)       = sulfaod(:ncol)
      bcaod_l(:ncol)         = bcaod(:ncol)
      pomaod_l(:ncol)        = pomaod(:ncol)
      soaaod_l(:ncol)        = soaaod(:ncol)
      ssltaod_l(:ncol)       = ssltaod(:ncol)
      aodabsbc_l(:ncol)      = aodabsbc(:ncol)
      burdendust_l(:ncol)    = burdendust(:ncol)
      burdenso4_l(:ncol)     = burdenso4(:ncol)
      burdenbc_l(:ncol)      = burdenbc(:ncol)
      burdenpom_l(:ncol)     = burdenpom(:ncol)
      burdensoa_l(:ncol)     = burdensoa(:ncol)
      burdenseasalt_l(:ncol) = burdenseasalt(:ncol)

      !-----------------------------------------------------------------
      ! Output day+night (dn) variants first (all columns valid)
      !-----------------------------------------------------------------

      ! SW column
      call history_out_field('AODVISdn',   aodvis_l)
      call history_out_field('AODUVdn',    aoduv)
      call history_out_field('AODNIRdn',   aodnir)
      call history_out_field('AODABSdn',   aodabs)
      call history_out_field('AODTOTdn',   aodtot)
      call history_out_field('SSAVISdn',   ssavis_l)

      ! Stratospheric
      call history_out_field('AODVISstdn', aodvisst)
      call history_out_field('AODUVstdn',  aoduvst)
      call history_out_field('AODNIRstdn', aodnirst)

      ! SW per-level
      call history_out_field('EXTINCTdn',    extinct)
      call history_out_field('EXTINCTUVdn',  extinctuv)
      call history_out_field('EXTINCTNIRdn', extinctnir)
      call history_out_field('ABSORBdn',     absorb)

      ! Per-species column (dn)
      call history_out_field('AODDUSTdn',       dustaod_l)
      call history_out_field('AODSO4dn',        sulfaod_l)
      call history_out_field('AODBCdn',         bcaod_l)
      call history_out_field('AODPOMdn',        pomaod_l)
      call history_out_field('AODSOAdn',        soaaod_l)
      call history_out_field('AODSSdn',         ssltaod_l)
      call history_out_field('AODABSBCdn',      aodabsbc_l)

      call history_out_field('BURDENDUSTdn',    burdendust_l)
      call history_out_field('BURDENSO4dn',     burdenso4_l)
      call history_out_field('BURDENBCdn',      burdenbc_l)
      call history_out_field('BURDENPOMdn',     burdenpom_l)
      call history_out_field('BURDENSOAdn',     burdensoa_l)
      call history_out_field('BURDENSEASALTdn', burdenseasalt_l)

      ! LW per-level (no day/night distinction)
      call history_out_field('TOTABSLW',  totabslw)
      call history_out_field('AODABSLW',  aodabslw)

      !-----------------------------------------------------------------
      ! Apply fillvalue to night columns for day-only variants
      !-----------------------------------------------------------------
      do i = 1, nnite
        aodvis_l(idxnite(i))  = fillvalue
        aoduv(idxnite(i))     = fillvalue
        aodnir(idxnite(i))    = fillvalue
        aodabs(idxnite(i))    = fillvalue
        aodtot(idxnite(i))    = fillvalue
        ssavis_l(idxnite(i))  = fillvalue
        aodvisst(idxnite(i))  = fillvalue
        aoduvst(idxnite(i))   = fillvalue
        aodnirst(idxnite(i))  = fillvalue

        extinct(idxnite(i), :)    = fillvalue
        extinctuv(idxnite(i), :)  = fillvalue
        extinctnir(idxnite(i), :) = fillvalue
        absorb(idxnite(i), :)     = fillvalue

        dustaod_l(idxnite(i))       = fillvalue
        sulfaod_l(idxnite(i))       = fillvalue
        bcaod_l(idxnite(i))         = fillvalue
        pomaod_l(idxnite(i))        = fillvalue
        soaaod_l(idxnite(i))        = fillvalue
        ssltaod_l(idxnite(i))       = fillvalue
        aodabsbc_l(idxnite(i))      = fillvalue

        burdendust_l(idxnite(i))    = fillvalue
        burdenso4_l(idxnite(i))     = fillvalue
        burdenbc_l(idxnite(i))      = fillvalue
        burdenpom_l(idxnite(i))     = fillvalue
        burdensoa_l(idxnite(i))     = fillvalue
        burdenseasalt_l(idxnite(i)) = fillvalue
      end do

      !-----------------------------------------------------------------
      ! Output day-only variants (night columns have fillvalue)
      !-----------------------------------------------------------------

      ! SW column
      call history_out_field('AODVIS',   aodvis_l)
      call history_out_field('AODUV',    aoduv)
      call history_out_field('AODNIR',   aodnir)
      call history_out_field('AODABS',   aodabs)
      call history_out_field('AODTOT',   aodtot)
      call history_out_field('SSAVIS',   ssavis_l)

      ! Stratospheric
      call history_out_field('AODVISst', aodvisst)
      call history_out_field('AODUVst',  aoduvst)
      call history_out_field('AODNIRst', aodnirst)

      ! SW per-level
      call history_out_field('EXTINCT',    extinct)
      call history_out_field('EXTINCTUV',  extinctuv)
      call history_out_field('EXTINCTNIR', extinctnir)
      call history_out_field('ABSORB',     absorb)

      ! Per-species column
      call history_out_field('AODDUST',       dustaod_l)
      call history_out_field('AODSO4',        sulfaod_l)
      call history_out_field('AODBC',         bcaod_l)
      call history_out_field('AODPOM',        pomaod_l)
      call history_out_field('AODSOA',        soaaod_l)
      call history_out_field('AODSS',         ssltaod_l)
      call history_out_field('AODABSBC',      aodabsbc_l)

      call history_out_field('BURDENDUST',    burdendust_l)
      call history_out_field('BURDENSO4',     burdenso4_l)
      call history_out_field('BURDENBC',      burdenbc_l)
      call history_out_field('BURDENPOM',     burdenpom_l)
      call history_out_field('BURDENSOA',     burdensoa_l)
      call history_out_field('BURDENSEASALT', burdenseasalt_l)

      ! AEROD_v and AODvstrt (aliases for AODVIS/AODVISst, day-only)
      call history_out_field('AEROD_v',  aodvis_l)
      call history_out_field('AODvstrt', aodvisst)

      ! Per-aerosol visible OD (ODV_<name>), day-only with fillvalue
      do iaer = 1, min(num_bulk_aer, num_odv)
        odv_tmp(:ncol) = odv_col_aod(:ncol, iaer)
        do i = 1, nnite
          odv_tmp(idxnite(i)) = fillvalue
        end do
        call history_out_field(trim(odv_names(iaer)), odv_tmp)
      end do

   end subroutine aerosol_optics_diagnostics_run

end module aerosol_optics_diagnostics
