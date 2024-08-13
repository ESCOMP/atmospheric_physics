! This module is used to diagnose the location of the tropopause. Multiple
! algorithms are provided, some of which may not be able to identify a
! tropopause in all situations. To handle these cases, an analytic
! definition and a climatology are provided that can be used to fill in
! when the original algorithm fails. The tropopause temperature and
! pressure are determined and can be output to the history file.
!
! These routines are based upon code in the WACCM chemistry module
! including mo_tropoause.F90 and llnl_set_chem_trop.F90. The code
! for the Reichler et al. [2003] algorithm is from:
!
!   http://www.gfdl.noaa.gov/~tjr/TROPO/tropocode.htm
!
! Author: Charles Bardeen
! Created: April, 2009

module tropopause_find
  !---------------------------------------------------------------
  ! ... variables for the tropopause module
  !---------------------------------------------------------------

  use ccpp_kinds,           only : kind_phys
  use shr_const_mod,        only : pi => shr_const_pi
  use ppgrid,               only : pcols, pver, begchunk, endchunk
  use cam_abortutils,       only : endrun
  use cam_logfile,          only : iulog
  use cam_history_support,  only : fillvalue
  use physics_types,        only : physics_state
  use spmd_utils,           only : masterproc

  implicit none

  private

  ! CCPP-compliant subroutines
  public  :: tropopause_find_init
  public  :: tropopause_find_run
  
  public  :: tropopause_readnl, tropopause_output
  public  :: tropopause_findChemTrop
  public  :: TROP_ALG_NONE, TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  public  :: TROP_ALG_STOBIE, TROP_ALG_HYBSTOB, TROP_ALG_TWMO, TROP_ALG_WMO
  public  :: TROP_ALG_CPP
  public  :: NOTFOUND

  save

  ! These parameters define and enumeration to be used to define the primary
  ! and backup algorithms to be used with the tropopause_find() method. The
  ! backup algorithm is meant to provide a solution when the primary algorithm
  ! fail. The algorithms that can't fail are: TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  ! and TROP_ALG_STOBIE.
  integer, parameter    :: TROP_ALG_NONE      = 1    ! Don't evaluate
  integer, parameter    :: TROP_ALG_ANALYTIC  = 2    ! Analytic Expression
  integer, parameter    :: TROP_ALG_CLIMATE   = 3    ! Climatology
  integer, parameter    :: TROP_ALG_STOBIE    = 4    ! Stobie Algorithm
  integer, parameter    :: TROP_ALG_TWMO      = 5    ! WMO Definition, Reichler et al. [2003]
  integer, parameter    :: TROP_ALG_WMO       = 6    ! WMO Definition
  integer, parameter    :: TROP_ALG_HYBSTOB   = 7    ! Hybrid Stobie Algorithm
  integer, parameter    :: TROP_ALG_CPP       = 8    ! Cold Point Parabolic
  
  integer, parameter    :: TROP_NALG          = 8    ! Number of Algorithms  
  character,parameter   :: TROP_LETTER(TROP_NALG) = (/ ' ', 'A', 'C', 'S', 'T', 'W', 'H', 'F' /)
                                                     ! unique identifier for output, don't use P

  ! These variables should probably be controlled by namelist entries.
  logical ,parameter    :: output_all         = .False.              ! output tropopause info from all algorithms
  integer ,parameter    :: default_primary    = TROP_ALG_TWMO        ! default primary algorithm
  integer ,parameter    :: default_backup     = TROP_ALG_CLIMATE     ! default backup algorithm

  ! Namelist variables
  character(len=256)    :: tropopause_climo_file = 'trop_climo'      ! absolute filepath of climatology file

  ! These variables are used to store the climatology data.
  real(kind_phys)              :: days(12)                                  ! days in the climatology
  real(kind_phys), pointer     :: tropp_p_loc(:,:,:)                        ! climatological tropopause pressures

  integer, parameter :: NOTFOUND = -1

  real(kind_phys),parameter :: ALPHA  = 0.03_kind_phys
    
  ! physical constants
  ! These constants are set in module variables rather than as parameters 
  ! to support the aquaplanet mode in which the constants have values determined
  ! by the experiment protocol
  real(kind_phys) :: cnst_kap     ! = cappa
  real(kind_phys) :: cnst_faktor  ! = -gravit/rair
  real(kind_phys) :: cnst_rga     ! = 1/gravit
  real(kind_phys) :: cnst_ka1     ! = cnst_kap - 1._kind_phys

!================================================================================================
contains
!================================================================================================

   ! Read namelist variables.
   subroutine tropopause_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'tropopause_readnl'

      namelist /tropopause_nl/ tropopause_climo_file
      !-----------------------------------------------------------------------------

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'tropopause_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, tropopause_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(tropopause_climo_file, len(tropopause_climo_file), mpichar, 0, mpicom)
#endif

   end subroutine tropopause_readnl


  ! This routine is called during intialization and must be called before the
  ! other methods in this module can be used. Its main tasks are to read in the
  ! climatology from a file and to define the output fields. Much of this code
  ! is taken from mo_tropopause.
  !> \section arg_table_tropopause_find_init Argument Table
  !! \htmlinclude tropopause_find_init.html
  subroutine tropopause_find_init(cappa, rair, gravit, errmsg, errflg)
  
    !use cam_history,   only: addfld, horiz_only
    real(kind_phys), intent(in)         :: cappa   ! R/Cp
    real(kind_phys), intent(in)         :: rair    ! Dry air gas constant (J K-1 kg-1)
    real(kind_phys), intent(in)         :: gravit  ! Gravitational acceleration (m s-2)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! define physical constants
    cnst_kap    = cappa
    cnst_faktor = -gravit/rair
    cnst_rga    = 1._kind_phys/gravit              ! Reciprocal of gravit (s2 m-1)
    cnst_ka1    = cnst_kap - 1._kind_phys

    ! Define the output fields.
    ! call addfld('TROP_P',          horiz_only,  'A', 'Pa',          'Tropopause Pressure',              flag_xyfill=.True.)
    ! call addfld('TROP_T',          horiz_only,  'A', 'K',           'Tropopause Temperature',           flag_xyfill=.True.)
    ! call addfld('TROP_Z',          horiz_only,  'A', 'm',           'Tropopause Height',                flag_xyfill=.True.)
    ! call addfld('TROP_DZ',         (/ 'lev' /), 'A', 'm',           'Relative Tropopause Height')
    ! call addfld('TROP_PD',         (/ 'lev' /), 'A', 'probability', 'Tropopause Probabilty')
    ! call addfld('TROP_FD',         horiz_only,  'A', 'probability', 'Tropopause Found')
    
    ! call addfld('TROPP_P',         horiz_only,  'A', 'Pa',          'Tropopause Pressure (primary)',    flag_xyfill=.True.)
    ! call addfld('TROPP_T',         horiz_only,  'A', 'K',           'Tropopause Temperature (primary)', flag_xyfill=.True.)
    ! call addfld('TROPP_Z',         horiz_only,  'A', 'm',           'Tropopause Height (primary)',      flag_xyfill=.True.)
    ! call addfld('TROPP_DZ',        (/ 'lev' /), 'A', 'm',           'Relative Tropopause Height (primary)')
    ! call addfld('TROPP_PD',        (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (primary)')
    ! call addfld('TROPP_FD',        horiz_only,  'A', 'probability', 'Tropopause Found (primary)')
    
    ! call addfld('TROPF_P',         horiz_only,  'A',  'Pa',         'Tropopause Pressure (cold point)',    flag_xyfill=.True.)
    ! call addfld('TROPF_T',         horiz_only,  'A',  'K',          'Tropopause Temperature (cold point)', flag_xyfill=.True.)
    ! call addfld('TROPF_Z',         horiz_only,  'A',  'm',          'Tropopause Height (cold point)',      flag_xyfill=.True.)
    ! call addfld('TROPF_DZ',        (/ 'lev' /),  'A', 'm',          'Relative Tropopause Height (cold point)', flag_xyfill=.True.)
    ! call addfld('TROPF_PD',        (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (cold point)')
    ! call addfld('TROPF_FD',        horiz_only,  'A', 'probability', 'Tropopause Found (cold point)')

    ! call addfld( 'hstobie_trop',   (/ 'lev' /), 'I',  'fraction of model time', 'Lowest level with stratospheric chemsitry')
    ! call addfld( 'hstobie_linoz',  (/ 'lev' /), 'I',  'fraction of model time', 'Lowest possible Linoz level')
    ! call addfld( 'hstobie_tropop', (/ 'lev' /), 'I', 'fraction of model time', &
    !      'Troposphere boundary calculated in chemistry' )

    ! ! If requested, be prepared to output results from all of the methods.
    ! if (output_all) then
    !   call addfld('TROPA_P',  horiz_only,  'A',  'Pa',          'Tropopause Pressure (analytic)',        flag_xyfill=.True.)
    !   call addfld('TROPA_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (analytic)',      flag_xyfill=.True.)
    !   call addfld('TROPA_Z',  horiz_only,  'A',  'm',          'Tropopause Height (analytic)',           flag_xyfill=.True.)
    !   call addfld('TROPA_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (analytic)')
    !   call addfld('TROPA_FD', horiz_only,  'A', 'probability', 'Tropopause Found (analytic)')

    !   call addfld('TROPC_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (climatology)',      flag_xyfill=.True.)
    !   call addfld('TROPC_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (climatology)',   flag_xyfill=.True.)
    !   call addfld('TROPC_Z',  horiz_only,  'A',  'm',          'Tropopause Height (climatology)',        flag_xyfill=.True.)
    !   call addfld('TROPC_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (climatology)')
    !   call addfld('TROPC_FD', horiz_only,  'A', 'probability', 'Tropopause Found (climatology)')

    !   call addfld('TROPS_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (stobie)',           flag_xyfill=.True.)
    !   call addfld('TROPS_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (stobie)',        flag_xyfill=.True.)
    !   call addfld('TROPS_Z',  horiz_only,  'A',  'm',          'Tropopause Height (stobie)',             flag_xyfill=.True.)
    !   call addfld('TROPS_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (stobie)')
    !   call addfld('TROPS_FD', horiz_only,  'A', 'probability', 'Tropopause Found (stobie)')

    !   call addfld('TROPT_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (twmo)',             flag_xyfill=.True.)
    !   call addfld('TROPT_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (twmo)',          flag_xyfill=.True.)
    !   call addfld('TROPT_Z',  horiz_only,  'A',  'm',          'Tropopause Height (twmo)',               flag_xyfill=.True.)
    !   call addfld('TROPT_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (twmo)')
    !   call addfld('TROPT_FD', horiz_only,  'A', 'probability', 'Tropopause Found (twmo)')

    !   call addfld('TROPW_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (WMO)',              flag_xyfill=.True.)
    !   call addfld('TROPW_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (WMO)',           flag_xyfill=.True.)
    !   call addfld('TROPW_Z',  horiz_only,  'A',  'm',          'Tropopause Height (WMO)',                flag_xyfill=.True.)
    !   call addfld('TROPW_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (WMO)')
    !   call addfld('TROPW_FD', horiz_only,  'A', 'probability', 'Tropopause Found (WMO)')

    !   call addfld('TROPH_P',  horiz_only,  'A',  'Pa',         'Tropopause Pressure (Hybrid Stobie)',    flag_xyfill=.True.)
    !   call addfld('TROPH_T',  horiz_only,  'A',  'K',          'Tropopause Temperature (Hybrid Stobie)', flag_xyfill=.True.)
    !   call addfld('TROPH_Z',  horiz_only,  'A',  'm',          'Tropopause Height (Hybrid Stobie)',      flag_xyfill=.True.)
    !   call addfld('TROPH_PD', (/ 'lev' /), 'A', 'probability', 'Tropopause Distribution (Hybrid Stobie)')
    !   call addfld('TROPH_FD', horiz_only,  'A', 'probability', 'Tropopause Found (Hybrid Stobie)')
    ! end if


    call tropopause_read_file()


  end subroutine tropopause_find_init


  ! Searches all the columns in the chunk and attempts to identify the tropopause.
  ! Two routines can be specifed, a primary routine which is tried first and a
  ! backup routine which will be tried only if the first routine fails. If the
  ! tropopause can not be identified by either routine, then a NOTFOUND is returned
  ! for the tropopause level, temperature and pressure.
  !> \section arg_table_tropopause_find_run Argument Table
  !! \htmlinclude tropopause_find_run.html
  subroutine tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
    tropLev, tropP, tropT, tropZ, primary, backup, &
    errmsg, errflg)

    implicit none

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    integer, optional, intent(in)       :: primary                   ! primary detection algorithm
    integer, optional, intent(in)       :: backup                    ! backup detection algorithm
    integer,            intent(out)     :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(out)     :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(out)     :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(out)     :: tropZ(:)              ! tropopause height (m)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! Local Variable
    integer       :: primAlg            ! Primary algorithm
    integer       :: backAlg            ! Backup algorithm

    ! Initialize the results to a missing value, so that the algorithms will
    ! attempt to find the tropopause for all of them.
    tropLev(:) = NOTFOUND
    if (present(tropP)) tropP(:) = fillvalue
    if (present(tropT)) tropT(:) = fillvalue
    if (present(tropZ)) tropZ(:) = fillvalue

    ! Set the algorithms to be used, either the ones provided or the defaults.
    if (present(primary)) then
      primAlg = primary
    else
      primAlg = default_primary
    end if

    if (present(backup)) then
      backAlg = backup
    else
      backAlg = default_backup
    end if

    ! Try to find the tropopause using the primary algorithm.
    if (primAlg /= TROP_ALG_NONE) then
      call tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                primAlg, tropLev, tropP, tropT, tropZ)
    end if

    if ((backAlg /= TROP_ALG_NONE) .and. any(tropLev(:) == NOTFOUND)) then
      call tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                backAlg, tropLev, tropP, tropT, tropZ)
    end if

    return
  end subroutine tropopause_find_run

  subroutine tropopause_read_file
    !------------------------------------------------------------------
    ! ... initialize upper boundary values
    !------------------------------------------------------------------
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish
    use dyn_grid,     only : get_dyn_grid_parm
    use phys_grid,    only : get_ncols_p, get_rlat_all_p, get_rlon_all_p	
    use ioFileMod,    only : getfil
    use time_manager, only : get_calday
    use physconst,    only : pi
    use cam_pio_utils, only: cam_pio_openfile
    use pio,          only : file_desc_t, var_desc_t, pio_inq_dimid, pio_inq_dimlen, &
         pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite

    !------------------------------------------------------------------
    ! ... local variables
    !------------------------------------------------------------------
    integer :: i, j, n
    integer :: ierr
    type(file_desc_t) :: pio_id
    integer :: dimid
    type(var_desc_t) :: vid
    integer :: nlon, nlat, ntimes
    integer :: start(3)
    integer :: count(3)
    integer, parameter :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
         716, 816, 915, 1016, 1115, 1216 /)
    integer :: plon, plat
    type(interp_type) :: lon_wgts, lat_wgts
    real(kind_phys), allocatable :: tropp_p_in(:,:,:)
    real(kind_phys), allocatable :: lat(:)
    real(kind_phys), allocatable :: lon(:)
    real(kind_phys) :: to_lats(pcols), to_lons(pcols)
    real(kind_phys), parameter :: d2r=pi/180._kind_phys, zero=0._kind_phys, twopi=pi*2._kind_phys
    character(len=256) :: locfn
    integer  :: c, ncols


    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')


    !-----------------------------------------------------------------------
    !       ... open netcdf file
    !-----------------------------------------------------------------------
    call getfil (tropopause_climo_file, locfn, 0)
    call cam_pio_openfile(pio_id, trim(locfn), PIO_NOWRITE)

    !-----------------------------------------------------------------------
    !       ... get time dimension
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'time', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, ntimes )
    if( ntimes /= 12 )then
       write(iulog,*) 'tropopause_find_init: number of months = ',ntimes,'; expecting 12'
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... get latitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lat', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlat )
    allocate( lat(nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_find_init: lat allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( pio_id, 'lat', vid )
    ierr = pio_get_var( pio_id, vid, lat )
    lat(:nlat) = lat(:nlat) * d2r
    !-----------------------------------------------------------------------
    !       ... get longitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lon', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlon )
    allocate( lon(nlon), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_find_init: lon allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( pio_id, 'lon', vid )
    ierr = pio_get_var( pio_id, vid, lon )
    lon(:nlon) = lon(:nlon) * d2r

    !------------------------------------------------------------------
    !  ... allocate arrays
    !------------------------------------------------------------------
    allocate( tropp_p_in(nlon,nlat,ntimes), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_find_init: tropp_p_in allocation error = ',ierr
       call endrun
    end if
    !------------------------------------------------------------------
    !  ... read in the tropopause pressure
    !------------------------------------------------------------------
    ierr = pio_inq_varid( pio_id, 'trop_p', vid )
    start = (/ 1, 1, 1 /)
    count = (/ nlon, nlat, ntimes /)
    ierr = pio_get_var( pio_id, vid, start, count, tropp_p_in )

    !------------------------------------------------------------------
    !  ... close the netcdf file
    !------------------------------------------------------------------
    call pio_closefile( pio_id )

    !--------------------------------------------------------------------
    !  ... regrid
    !--------------------------------------------------------------------

    allocate( tropp_p_loc(pcols,begchunk:endchunk,ntimes), stat=ierr )

    if( ierr /= 0 ) then
      write(iulog,*) 'tropopause_find_init: tropp_p_loc allocation error = ',ierr
      call endrun
    end if

    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       call lininterp_init(lon, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(lat, nlat, to_lats, ncols, 1, lat_wgts)
       do n=1,ntimes
          call lininterp(tropp_p_in(:,:,n), nlon, nlat, tropp_p_loc(1:ncols,c,n), ncols, lon_wgts, lat_wgts)    
       end do
       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)
    end do
    deallocate(lon)
    deallocate(lat)
    deallocate(tropp_p_in)

    !--------------------------------------------------------
    ! ... initialize the monthly day of year times
    !--------------------------------------------------------

    do n = 1,12
       days(n) = get_calday( dates(n), 0 )
    end do
    if (masterproc) then
       write(iulog,*) 'tropopause_find_init : days'
       write(iulog,'(1p,5g15.8)') days(:)
    endif

  end subroutine tropopause_read_file
  

  ! This analytic expression closely matches the mean tropopause determined
  ! by the NCEP reanalysis and has been used by the radiation code.
  subroutine tropopause_analytic(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                 tropLev, tropP, tropT, tropZ)

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)   :: tropLev(:)             ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)
    
    ! Local Variables
    integer       :: i
    integer       :: k
    real(kind_phys)      :: tP                       ! tropopause pressure (Pa)

    ! Iterate over all of the columns.
    do i = 1, ncol
     
      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then
      
        ! Calculate the pressure of the tropopause.
        tP = (25000.0_kind_phys - 15000.0_kind_phys * (cos(lat(i)))**2)
      
        ! Find the level that contains the tropopause.
        do k = pver, 2, -1
          if (tP >= pint(i, k)) then
            tropLev(i) = k
            exit
          end if
        end do
        
        ! Return the optional outputs
        if (present(tropP)) tropP(i) = tP
        
        if (present(tropT)) then
          tropT(i) = tropopause_interpolateT(pmid, t, i, tropLev(i), tP)
        end if

        if (present(tropZ)) then
          tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
        end if
      end if
    end do
  end subroutine tropopause_analytic


  ! Read the tropopause pressure in from a file containging a climatology. The
  ! data is interpolated to the current dat of year and latitude.
  !
  ! NOTE: The data is read in during tropopause_init and stored in the module
  ! variable trop
  subroutine tropopause_climate(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                tropLev, tropP, tropT, tropZ)
     use time_manager, only : get_curr_calday

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)
 
    ! Local Variables
    integer       :: i
    integer       :: k
    integer       :: m
    integer       :: lchnk                    ! chunk identifier
    real(kind_phys)      :: tP                       ! tropopause pressure (Pa)
    real(kind_phys)      :: calday                   ! day of year including fraction
    real(kind_phys)      :: dels
    integer       :: last
    integer       :: next

    ! Information about the chunk.  
    lchnk = pstate%lchnk
    
    ! If any columns remain to be indentified, the nget the current
    ! day from the calendar.
    
    if (any(tropLev == NOTFOUND)) then
    
      ! Determine the calendar day.
      calday = get_curr_calday()
      
      !--------------------------------------------------------
      ! ... setup the time interpolation
      !--------------------------------------------------------
      if( calday < days(1) ) then
        next = 1
        last = 12
        dels = (365._kind_phys + calday - days(12)) / (365._kind_phys + days(1) - days(12))
      else if( calday >= days(12) ) then
        next = 1
        last = 12
        dels = (calday - days(12)) / (365._kind_phys + days(1) - days(12))
      else
        do m = 11,1,-1
           if( calday >= days(m) ) then
              exit
           end if
        end do
        last = m
        next = m + 1
        dels = (calday - days(m)) / (days(m+1) - days(m))
      end if
      
      dels = max( min( 1._kind_phys,dels ),0._kind_phys )
        

      ! Iterate over all of the columns.
      do i = 1, ncol
       
        ! Skip column in which the tropopause has already been found.
        if (tropLev(i) == NOTFOUND) then
        
        !--------------------------------------------------------
        ! ... get tropopause level from climatology
        !--------------------------------------------------------
          ! Interpolate the tropopause pressure.
          tP = tropp_p_loc(i,lchnk,last) &
            + dels * (tropp_p_loc(i,lchnk,next) - tropp_p_loc(i,lchnk,last))
                
          ! Find the associated level.
          do k = pver, 2, -1
            if (tP >= pint(i, k)) then
              tropLev(i) = k
              exit
            end if
          end do
      
          ! Return the optional outputs
          if (present(tropP)) tropP(i) = tP
          
          if (present(tropT)) then
            tropT(i) = tropopause_interpolateT(pmid, t, i, tropLev(i), tP)
          end if

          if (present(tropZ)) then
            tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
          end if
        end if
      end do
    end if        

    return    
  end subroutine tropopause_climate
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine tropopause_hybridstobie(ncol, pver, pmid, t, zm, &
                                     tropLev, tropP, tropT, tropZ)
    !use cam_history,  only : outfld
  
    !-----------------------------------------------------------------------
    ! Originally written by Philip Cameron-Smith, LLNL
    !
    !   Stobie-Linoz hybrid: the highest altitude of 
    !          a) Stobie algorithm, or 
    !          b) minimum Linoz pressure.
    !
    ! NOTE: the ltrop(i) gridbox itself is assumed to be a STRATOSPHERIC gridbox.
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !        ... Local variables
    !-----------------------------------------------------------------------
    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levelserp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver

    integer,            intent(inout)   :: tropLev(:)             ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)
    
    real(kind_phys),parameter  ::  min_Stobie_Pressure= 40.E2_kind_phys !For case 2 & 4.  [Pa]
    real(kind_phys),parameter  ::  max_Linoz_Pressure =208.E2_kind_phys !For case     4.  [Pa]

    integer      :: i, k, ncol
    real(kind_phys)     :: stobie_min, shybrid_temp      !temporary variable for case 2 & 3.
    integer      :: ltrop_linoz(pcols)            !Lowest possible Linoz vertical level
    integer      :: ltrop_trop(pcols)             !Tropopause level for hybrid case.
    logical      :: ltrop_linoz_set               !Flag that lowest linoz level already found.
    real(kind_phys)     :: trop_output(pcols,pver)        !For output purposes only.
    real(kind_phys)     :: trop_linoz_output(pcols,pver)  !For output purposes only.
    real(kind_phys)     :: trop_trop_output(pcols,pver)   !For output purposes only.

    !    write(iulog,*) 'In set_chem_trop, o3_ndx =',o3_ndx
    ltrop_linoz(:) = 1  ! Initialize to default value.
    ltrop_trop(:) = 1   ! Initialize to default value.

    LOOP_COL4: do i=1,ncol

       ! Skip column in which the tropopause has already been found.
       not_found: if (tropLev(i) == NOTFOUND) then

          stobie_min = 1.e10_kind_phys    ! An impossibly large number
          ltrop_linoz_set = .FALSE.
          LOOP_LEV: do k=pver,1,-1
             IF (pmid(i,k) < min_stobie_pressure) cycle
             shybrid_temp = ALPHA * t(i,k) - Log10(pmid(i,k))
             !PJC_NOTE: the units of pmid won't matter, because it is just an additive offset.
             IF (shybrid_temp<stobie_min) then 
                ltrop_trop(i)=k     
                stobie_min = shybrid_temp
             ENDIF
             IF (pmid(i,k) < max_Linoz_pressure .AND. .NOT. ltrop_linoz_set) THEN
                ltrop_linoz(i) = k
                ltrop_linoz_set = .TRUE.
             ENDIF
          enddo LOOP_LEV

          tropLev(i) = MIN(ltrop_trop(i),ltrop_linoz(i))

          if (present(tropP)) then
             tropP(i) = pmid(i,tropLev(i))
          endif
          if (present(tropT)) then
             tropT(i) = t(i,tropLev(i))
          endif
          if (present(tropZ)) then
             tropZ(i) = zm(i,tropLev(i))
          endif

       endif not_found

    enddo LOOP_COL4

    trop_output(:,:)=0._kind_phys
    trop_linoz_output(:,:)=0._kind_phys
    trop_trop_output(:,:)=0._kind_phys
    do i=1,ncol
       if (tropLev(i)>0) then
          trop_output(i,tropLev(i))=1._kind_phys
          trop_linoz_output(i,ltrop_linoz(i))=1._kind_phys
          trop_trop_output(i,ltrop_trop(i))=1._kind_phys
       endif
    enddo

    !call outfld( 'hstobie_trop',   trop_output(:ncol,:),       ncol, pstate%lchnk )
    !call outfld( 'hstobie_linoz',  trop_linoz_output(:ncol,:), ncol, pstate%lchnk )
    !call outfld( 'hstobie_tropop', trop_trop_output(:ncol,:),  ncol, pstate%lchnk )

  end subroutine tropopause_hybridstobie
  
  ! This routine originates with Stobie at NASA Goddard, but does not have a
  ! known reference. It was supplied by Philip Cameron-Smith of LLNL.
  !
  subroutine tropopause_stobie(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                               tropLev, tropP, tropT, tropZ)

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    integer,            intent(inout)   :: tropLev(:)             ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)
    
    ! Local Variables
    integer       :: i
    integer       :: k
    integer       :: tLev                     ! tropopause level
    real(kind_phys)      :: tP                       ! tropopause pressure (Pa)
    real(kind_phys)      :: stobie(pver)             ! stobie weighted temperature
    real(kind_phys)      :: sTrop                    ! stobie value at the tropopause

    ! Iterate over all of the columns.
    do i = 1, ncol
     
      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then
      
        ! Caclulate a pressure weighted temperature.
        stobie(:) = ALPHA * t(i,:) - log10(pmid(i, :))

        ! Search from the bottom up, looking for the first minimum.
        tLev  = -1
  
        do k = pver-1, 1, -1
    
          if (pmid(i, k) <= 4000._kind_phys) then
            exit
          end if
    
          if (pmid(i, k) >= 55000._kind_phys) then
            cycle
          end if
          
          if ((tLev == -1) .or. (stobie(k) < sTrop)) then
            tLev  = k
            tP    = pmid(i, k)
            sTrop = stobie(k)
          end if
        end do
        
        if (tLev /= -1) then
          tropLev(i) = tLev
        
          ! Return the optional outputs
          if (present(tropP)) tropP(i) = tP
          
          if (present(tropT)) then
            tropT(i) = tropopause_interpolateT(pmid, t, i, tropLev(i), tP)
          end if

          if (present(tropZ)) then
            tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
          end if
        end if
      end if
    end do
    
    return
  end subroutine tropopause_stobie


  ! This routine is an implementation of Reichler et al. [2003] done by
  ! Reichler and downloaded from his web site. Minimal modifications were
  ! made to have the routine work within the CAM framework (i.e. using
  ! CAM constants and types).
  !
  ! NOTE: I am not a big fan of the goto's and multiple returns in this
  ! code, but for the moment I have left them to preserve as much of the
  ! original and presumably well tested code as possible.
  ! UPDATE: The most "obvious" substitutions have been made to replace
  ! goto/return statements with cycle/exit. The structure is still
  ! somewhat tangled.
  ! UPDATE 2: "gamma" renamed to "gam" in order to avoid confusion
  ! with the Fortran 2008 intrinsic. "level" argument removed because
  ! a physics column is not contiguous, so using explicit dimensions
  ! will cause the data to be needlessly copied.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! determination of tropopause height from gridded temperature data
  !
  ! reference: Reichler, T., M. Dameris, and R. Sausen (2003)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine twmo(t, p, plimu, pliml, gam, trp)

    real(kind_phys), intent(in), dimension(:)      :: t, p
    real(kind_phys), intent(in)                    :: plimu, pliml, gam
    real(kind_phys), intent(out)                   :: trp
    
    real(kind_phys), parameter                     :: deltaz = 2000.0_kind_phys

    real(kind_phys)                                :: pmk, pm, a, b, tm, dtdp, dtdz
    real(kind_phys)                                :: ag, bg, ptph
    real(kind_phys)                                :: pm0, pmk0, dtdz0
    real(kind_phys)                                :: p2km, asum, aquer
    real(kind_phys)                                :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
    integer                                 :: level
    integer                                 :: icount, jj
    integer                                 :: j
    

    trp=-99.0_kind_phys                           ! negative means not valid
    
    ! initialize start level
    ! dt/dz
    level = size(t)
    pmk= .5_kind_phys * (p(level-1)**cnst_kap+p(level)**cnst_kap)
    pm = pmk**(1/cnst_kap)               
    a = (t(level-1)-t(level))/(p(level-1)**cnst_kap-p(level)**cnst_kap)
    b = t(level)-(a*p(level)**cnst_kap)
    tm = a * pmk + b               
    dtdp = a * cnst_kap * (pm**cnst_ka1)
    dtdz = cnst_faktor*dtdp*pm/tm

    main_loop: do j=level-1,2,-1
      pm0 = pm
      pmk0 = pmk
      dtdz0  = dtdz
    
      ! dt/dz
      pmk= .5_kind_phys * (p(j-1)**cnst_kap+p(j)**cnst_kap)
      pm = pmk**(1/cnst_kap)               
      a = (t(j-1)-t(j))/(p(j-1)**cnst_kap-p(j)**cnst_kap)
      b = t(j)-(a*p(j)**cnst_kap)
      tm = a * pmk + b               
      dtdp = a * cnst_kap * (pm**cnst_ka1)
      dtdz = cnst_faktor*dtdp*pm/tm
      ! dt/dz valid?
      if (dtdz.le.gam) cycle main_loop    ! no, dt/dz < -2 K/km
      if (pm.gt.plimu)   cycle main_loop    ! no, too low
  
      ! dtdz is valid, calculate tropopause pressure
      if (dtdz0.lt.gam) then
        ag = (dtdz-dtdz0) / (pmk-pmk0)     
        bg = dtdz0 - (ag * pmk0)          
        ptph = exp(log((gam-bg)/ag)/cnst_kap)
      else
        ptph = pm
      endif
  
      if (ptph.lt.pliml) cycle main_loop
      if (ptph.gt.plimu) cycle main_loop
  
      ! 2nd test: dtdz above 2 km must not exceed gam
      p2km = ptph + deltaz*(pm/tm)*cnst_faktor     ! p at ptph + 2km
      asum = 0.0_kind_phys                                ! dtdz above
      icount = 0                                   ! number of levels above
  
      ! test until apm < p2km
      in_loop: do jj=j,2,-1
    
        pmk2 = .5_kind_phys * (p(jj-1)**cnst_kap+p(jj)**cnst_kap) ! p mean ^kappa
        pm2 = pmk2**(1/cnst_kap)                           ! p mean
        if(pm2.gt.ptph) cycle in_loop            ! doesn't happen
        if(pm2.lt.p2km) exit in_loop             ! ptropo is valid

        a2 = (t(jj-1)-t(jj))                     ! a
        a2 = a2/(p(jj-1)**cnst_kap-p(jj)**cnst_kap)
        b2 = t(jj)-(a2*p(jj)**cnst_kap)          ! b
        tm2 = a2 * pmk2 + b2                     ! T mean
        dtdp2 = a2 * cnst_kap * (pm2**(cnst_kap-1))  ! dt/dp
        dtdz2 = cnst_faktor*dtdp2*pm2/tm2
        asum = asum+dtdz2
        icount = icount+1
        aquer = asum/float(icount)               ! dt/dz mean
   
        ! discard ptropo ?
        if (aquer.le.gam) cycle main_loop      ! dt/dz above < gam
    
      enddo in_loop  ! test next level
    
      trp = ptph
      exit main_loop
    enddo main_loop
    
  end subroutine twmo
  

  ! This routine uses an implementation of Reichler et al. [2003] done by
  ! Reichler and downloaded from his web site. This is similar to the WMO
  ! routines, but is designed for GCMs with a coarse vertical grid.
  subroutine tropopause_twmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                             tropLev, tropP, tropT, tropZ)

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables 
    real(kind_phys), parameter     :: gam    = -0.002_kind_phys         ! K/m
    real(kind_phys), parameter     :: plimu    = 45000._kind_phys         ! Pa
    real(kind_phys), parameter     :: pliml    = 7500._kind_phys          ! Pa
     
    integer                 :: i
    integer                 :: k
    real(kind_phys)                :: tP                       ! tropopause pressure (Pa)

    ! Iterate over all of the columns.
    do i = 1, ncol
     
      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        ! Use the routine from Reichler.
        call twmo(t(i, :), pmid(i, :), plimu, pliml, gam, tP)
     
        ! if successful, store of the results and find the level and temperature.
        if (tP > 0) then
        
          ! Find the associated level.
          do k = pver, 2, -1
            if (tP >= pint(i, k)) then
              tropLev(i) = k
              exit
            end if
          end do
          
          ! Return the optional outputs
          if (present(tropP)) tropP(i) = tP
          
          if (present(tropT)) then
            tropT(i) = tropopause_interpolateT(pmid, t, i, tropLev(i), tP)
          end if

          if (present(tropZ)) then
            tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
          end if
        end if
      end if
    end do
    
    return
  end subroutine tropopause_twmo
  
  ! This routine implements the WMO definition of the tropopause (WMO, 1957; Seidel and Randel, 2006).
  ! This requires that the lapse rate be less than 2 K/km for an altitude range
  ! of 2 km. The search starts at the surface and stops the first time this
  ! criteria is met.
  !
  ! NOTE: This code was modeled after the code in mo_tropopause; however, the
  ! requirement that dt be greater than 0 was removed and the check to make
  ! sure that the lapse rate is maintained for 2 km was added.
  subroutine tropopause_wmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                            tropLev, tropP, tropT, tropZ)

    type(physics_state), intent(in)    :: pstate 
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables 
    real(kind_phys), parameter    :: ztrop_low   = 5000._kind_phys        ! lowest tropopause level allowed (m)
    real(kind_phys), parameter    :: ztrop_high  = 20000._kind_phys       ! highest tropopause level allowed (m)
    real(kind_phys), parameter    :: max_dtdz    = 0.002_kind_phys        ! max dt/dz for tropopause level (K/m)
    real(kind_phys), parameter    :: min_trop_dz = 2000._kind_phys        ! min tropopause thickness (m)

    integer                 :: i
    integer                 :: k
    integer                 :: k2
    integer                 :: ncol                         ! number of columns in the chunk
    real(kind_phys)                :: tP                           ! tropopause pressure (Pa)
    real(kind_phys)                :: dt

    ! Iterate over all of the columns.
    do i = 1, ncol
     
      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        kloop: do k = pver-1, 2, -1
         
          ! Skip levels below the minimum and stop if nothing is found
          ! before the maximum.
          if (zm(i, k) < ztrop_low) then
            cycle kloop
          else if (zm(i, k) > ztrop_high) then
            exit kloop
          end if
          
          ! Compare the actual lapse rate to the threshold
          dt = t(i, k) - t(i, k-1)
            
          if (dt <= (max_dtdz * (zm(i, k-1) - zm(i, k)))) then
            
            ! Make sure that the lapse rate stays below the threshold for the
            ! specified range.
            k2loop: do k2 = k-1, 2, -1
              if ((zm(i, k2) - zm(i, k)) >= min_trop_dz) then
                tP = pmid(i, k)
                tropLev(i) = k
                exit k2loop
              end if
              
              dt = t(i, k) - t(i, k2)
              if (dt > (max_dtdz * (zm(i, k2) - zm(i, k)))) then
                exit k2loop
              end if
           end do k2loop

           if (tropLev(i) == NOTFOUND) then
              cycle kloop
           else 

              ! Return the optional outputs
              if (present(tropP)) tropP(i) = tP
              
              if (present(tropT)) then
                tropT(i) = tropopause_interpolateT(pmid, t, i, tropLev(i), tP)
              end if
  
              if (present(tropZ)) then
                tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
              end if

              exit kloop
            end if
          end if
        end do kloop
      end if
    end do
    
    return
  end subroutine tropopause_wmo
  
  
  ! This routine searches for the cold point tropopause, and uses a parabolic
  ! fit of the coldest point and two adjacent points to interpolate the cold point
  ! between model levels.
  subroutine tropopause_cpp(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                            tropLev, tropP, tropT, tropZ)

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables 
    real(kind_phys), parameter    :: ztrop_low   = 5000._kind_phys        ! lowest tropopause level allowed (m)
    real(kind_phys), parameter    :: ztrop_high  = 25000._kind_phys       ! highest tropopause level allowed (m)

    integer                 :: i
    integer                 :: k, firstk, lastk
    integer                 :: k2
    integer                 :: ncol                         ! number of columns in the chunk
    real(kind_phys)                :: tZ                           ! tropopause height (m)
    real(kind_phys)                :: tmin
    real(kind_phys)                :: f0, f1, f2
    real(kind_phys)                :: x0, x1, x2
    real(kind_phys)                :: c0, c1, c2
    real(kind_phys)                :: a, b, c

    ! Iterate over all of the columns.
    do i = 1, ncol
     
      firstk = 0
      lastk  = pver+1

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then
        tmin = 1e6_kind_phys

        kloop: do k = pver-1, 2, -1
         
          ! Skip levels below the minimum and stop if nothing is found
          ! before the maximum.
          if (zm(i, k) < ztrop_low) then
            firstk = k
            cycle kloop
          else if (zm(i, k) > ztrop_high) then
            lastk = k
            exit kloop
          end if
          
          ! Find the coldest point
          if (t(i, k) < tmin) then
            tropLev(i) = k
            tmin = t(i,k)
          end if
        end do kloop

        
        ! If the minimum is at the edge of the search range, then don't
        ! consider this to be a minima
        if ((tropLev(i) >= (firstk-1)) .or. (tropLev(i) <= (lastk+1))) then
          tropLev(i) = NOTFOUND
        else

          ! If returning P, Z, or T, then do a parabolic fit using the
          ! cold point and it its 2 surrounding points to interpolate
          ! between model levels.
          if (present(tropP) .or. present(tropZ) .or. present(tropT)) then
            f0 = t(i, tropLev(i)-1)
            f1 = t(i, tropLev(i))
            f2 = t(i, tropLev(i)+1)

            x0 = zm(i, tropLev(i)-1)
            x1 = zm(i, tropLev(i))
            x2 = zm(i, tropLev(i)+1)

            c0 = (x0-x1)*(x0-x2)
            c1 = (x1-x0)*(x1-x2)
            c2 = (x2-x0)*(x2-x1)

            ! Determine the quadratic coefficients of:
            !   T = a * z^2 - b*z + c
            a = (f0/c0 + f1/c1 + f2/c2)
            b = (f0/c0*(x1+x2) + f1/c1*(x0+x2) + f2/c2*(x0+x1))
            c = f0/c0*x1*x2 + f1/c1*x0*x2 + f2/c2*x0*x1

            ! Find the altitude of the minimum temperature
            tZ = 0.5_kind_phys * b / a
            
            ! The fit should be between the upper and lower points,
            ! so skip the point if the fit fails.
            if ((tZ >= x0) .or. (tZ <= x2)) then
              tropLev(i) = NOTFOUND
            else

              ! Return the optional outputs
              if (present(tropP)) then
                tropP(i) = tropopause_interpolateP(pmid, zm, i, tropLev(i), tZ)
              end if
            
              if (present(tropT)) then
                tropT(i) = a * tZ*tZ - b*tZ + c
              end if

              if (present(tropZ)) then
                tropZ(i) = tZ
              end if
            end if
          end if
        end if
      end if
    end do
    
    return
  end subroutine tropopause_cpp
  
  ! Searches all the columns in the chunk and attempts to identify the "chemical"
  ! tropopause. This is the lapse rate tropopause, backed up by the climatology
  ! if the lapse rate fails to find the tropopause at pressures higher than a certain
  ! threshold. This pressure threshold depends on latitude. Between 50S and 50N, 
  ! the climatology is used if the lapse rate tropopause is not found at P > 75 hPa. 
  ! At high latitude (poleward of 50), the threshold is increased to 125 hPa to 
  ! eliminate false events that are sometimes detected in the cold polar stratosphere.
  !
  ! NOTE: This routine was adapted from code in chemistry.F90 and mo_gasphase_chemdr.F90.
  subroutine tropopause_findChemTrop(ncol, pver, lat, pint, pmid, t, zm, phis, &
                                     tropLev, primary, backup)

    implicit none

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    integer, optional, intent(in)       :: primary                   ! primary detection algorithm
    integer, optional, intent(in)       :: backup                    ! backup detection algorithm
    integer,            intent(out)     :: tropLev(:)            ! tropopause level index

    ! Local Variable
    real(kind_phys), parameter :: rad2deg = 180._kind_phys/pi                      ! radians to degrees conversion factor
    real(kind_phys)            :: dlats(pcols)
    integer             :: i
    integer             :: backAlg

    ! First use the lapse rate tropopause.
    call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zm, phis, tropLev, primary=primary, backup=TROP_ALG_NONE)
   
    ! Now check high latitudes (poleward of 50) and set the level to the
    ! climatology if the level was not found or is at P <= 125 hPa.
    dlats(:ncol) = lat(:ncol) * rad2deg ! convert to degrees

    if (present(backup)) then
      backAlg = backup
    else
      backAlg = default_backup
    end if
    
    do i = 1, ncol
      if (abs(dlats(i)) > 50._kind_phys) then
        if (tropLev(i) .ne. NOTFOUND) then
          if (pmid(i, tropLev(i)) <= 12500._kind_phys) then
            tropLev(i) = NOTFOUND
          end if
        end if
      end if
    end do
        
    ! Now use the backup algorithm
    if ((backAlg /= TROP_ALG_NONE) .and. any(tropLev(:) == NOTFOUND)) then
      call tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zm, phis, backAlg, tropLev)
    end if
    
    return
  end subroutine tropopause_findChemTrop
  
  
  ! Call the appropriate tropopause detection routine based upon the algorithm
  ! specifed.
  !
  ! NOTE: It is assumed that the output fields have been initialized by the
  ! caller, and only output values set to fillvalue will be detected.
  subroutine tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                  algorithm, tropLev, tropP, tropT, tropZ)

    implicit none

    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    integer,            intent(in)      :: algorithm                 ! detection algorithm
    integer,            intent(inout)   :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)              ! tropopause height (m)

    ! Dispatch the request to the appropriate routine.
    select case(algorithm)
      case(TROP_ALG_ANALYTIC)
        call tropopause_analytic(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_CLIMATE)
        call tropopause_climate(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_STOBIE)
        call tropopause_stobie(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_HYBSTOB)
        call tropopause_hybridstobie(ncol, pver, pmid, t, zm, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_TWMO)
        call tropopause_twmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_WMO)
        call tropopause_wmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_CPP)
        call tropopause_cpp(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case default
        write(iulog, *) 'tropopause: Invalid detection algorithm (',  algorithm, ') specified.'
        call endrun
    end select
    
    return
  end subroutine tropopause_findUsing


  ! This routine interpolates the pressures in the physics state to
  ! find the pressure at the specified tropopause altitude.
  function tropopause_interpolateP(pmid, zm, icol, tropLev, tropZ)
 
    implicit none

    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index   
    real(kind_phys), optional, intent(in)      :: tropZ              ! tropopause pressure (m)
    real(kind_phys)                            :: tropopause_interpolateP
    
    ! Local Variables
    real(kind_phys)   :: tropP              ! tropopause pressure (Pa)
    real(kind_phys)   :: dlogPdZ            ! dlog(p)/dZ
    
    ! Interpolate the temperature linearly against log(P)
    
    ! Is the tropopause at the midpoint?
    if (tropZ == zm(icol, tropLev)) then
      tropP = pmid(icol, tropLev)
    
    else if (tropZ > zm(icol, tropLev)) then
    
      ! It is above the midpoint? Make sure we aren't at the top.
      if (tropLev > 1) then
        dlogPdZ = (log(pmid(icol, tropLev)) - log(pmid(icol, tropLev - 1))) / &
          (zm(icol, tropLev) - zm(icol, tropLev - 1))
        tropP = pmid(icol, tropLev) + exp((tropZ - zm(icol, tropLev)) * dlogPdZ)
      end if
    else
      
      ! It is below the midpoint. Make sure we aren't at the bottom.
      if (tropLev < pver) then
        dlogPdZ =  (log(pmid(icol, tropLev + 1)) - log(pmid(icol, tropLev))) / &
          (zm(icol, tropLev + 1) - zm(icol, tropLev))
        tropP = pmid(icol, tropLev) + exp((tropZ - zm(icol, tropLev)) * dlogPdZ)
      end if
    end if
    
    tropopause_interpolateP = tropP
  end function tropopause_interpolateP

  
  ! This routine interpolates the temperatures in the physics state to
  ! find the temperature at the specified tropopause pressure.
  function tropopause_interpolateT(pmid, t, icol, tropLev, tropP)
 
    implicit none

    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index   
    real(kind_phys), optional, intent(in)      :: tropP              ! tropopause pressure (Pa)
    real(kind_phys)                            :: tropopause_interpolateT
    
    ! Local Variables
    real(kind_phys)   :: tropT              ! tropopause temperature (K)
    real(kind_phys)   :: dTdlogP            ! dT/dlog(P)
    
    ! Interpolate the temperature linearly against log(P)
    
    ! Is the tropopause at the midpoint?
    if (tropP == pmid(icol, tropLev)) then
      tropT = t(icol, tropLev)
    
    else if (tropP < pmid(icol, tropLev)) then
    
      ! It is above the midpoint? Make sure we aren't at the top.
      if (tropLev > 1) then
        dTdlogP = (t(icol, tropLev) - t(icol, tropLev - 1)) / &
          (log(pmid(icol, tropLev)) - log(pmid(icol, tropLev - 1)))
        tropT = t(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dTdlogP
      end if
    else
      
      ! It is below the midpoint. Make sure we aren't at the bottom.
      if (tropLev < pver) then
        dTdlogP = (t(icol, tropLev + 1) - t(icol, tropLev)) / &
          (log(pmid(icol, tropLev + 1)) - log(pmid(icol, tropLev)))
        tropT = t(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dTdlogP
      end if
    end if
    
    tropopause_interpolateT = tropT
  end function tropopause_interpolateT

  
  ! This routine interpolates the geopotential height in the physics state to
  ! find the geopotential height at the specified tropopause pressure.
  function tropopause_interpolateZ(pint, pmid, zi, zm, phis, icol, tropLev, tropP)

    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index   
    real(kind_phys), optional, intent(in)      :: tropP              ! tropopause pressure (Pa)
    real(kind_phys)                            :: tropopause_interpolateZ
    
    ! Local Variables
    real(kind_phys)   :: tropZ              ! tropopause geopotential height (m)
    real(kind_phys)   :: dZdlogP            ! dZ/dlog(P)
    
    ! Interpolate the geopotential height linearly against log(P)
    
    ! Is the tropoause at the midpoint?
    if (tropP == pmid(icol, tropLev)) then
      tropZ = zm(icol, tropLev)
    
    else if (tropP < pmid(icol, tropLev)) then
    
      ! It is above the midpoint? Make sure we aren't at the top.
      dZdlogP = (zm(icol, tropLev) - zi(icol, tropLev)) / &
        (log(pmid(icol, tropLev)) - log(pint(icol, tropLev)))
      tropZ = zm(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dZdlogP
    else
      
      ! It is below the midpoint. Make sure we aren't at the bottom.
      dZdlogP = (zm(icol, tropLev) - zi(icol, tropLev+1)) / &
        (log(pmid(icol, tropLev)) - log(pint(icol, tropLev+1)))
      tropZ = zm(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dZdlogP
    end if
    
    tropopause_interpolateZ = tropZ + phis(icol)*cnst_rga
  end function tropopause_interpolateZ

  
  ! Output the tropopause pressure and temperature to the history files. Two sets
  ! of output will be generated, one for the default algorithm and another one
  ! using the default routine, but backed by a climatology when the default
  ! algorithm fails.
  subroutine tropopause_output(ncol, pver, lat, pint, pmid, t, zi, zm, phis)
    !use cam_history,  only : outfld
    
    real(kind_phys), intent(in)         :: ncol          ! Number of atmospheric columns
    real(kind_phys), intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:,:)      ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m), pverp
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m), pver
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
  
    ! Local Variables
    integer       :: i
    integer       :: alg
    integer       :: ncol                     ! number of cloumns in the chunk
    !integer       :: lchnk                    ! chunk identifier
    integer       :: tropLev(pcols)           ! tropopause level index   
    real(kind_phys)      :: tropP(pcols)             ! tropopause pressure (Pa)
    real(kind_phys)      :: tropT(pcols)             ! tropopause temperature (K)
    real(kind_phys)      :: tropZ(pcols)             ! tropopause height (m)
    real(kind_phys)      :: tropFound(pcols)         ! tropopause found
    real(kind_phys)      :: tropDZ(pcols, pver)      ! relative tropopause height (m)
    real(kind_phys)      :: tropPdf(pcols, pver)     ! tropopause probability distribution

    ! Information about the chunk.  
    !lchnk = pstate%lchnk

    ! Find the tropopause using the default algorithm backed by the climatology.
    call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ)
    
    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue 
    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ(i)
      end if
    end do

    ! call outfld('TROP_P',   tropP(:ncol),      ncol, lchnk)
    ! call outfld('TROP_T',   tropT(:ncol),      ncol, lchnk)
    ! call outfld('TROP_Z',   tropZ(:ncol),      ncol, lchnk)
    ! call outfld('TROP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    ! call outfld('TROP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    ! call outfld('TROP_FD',  tropFound(:ncol),  ncol, lchnk)
    
    
    ! Find the tropopause using just the primary algorithm.
    call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, backup=TROP_ALG_NONE)

    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue 
    
    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ(i)
      end if
    end do

    ! call outfld('TROPP_P',   tropP(:ncol),      ncol, lchnk)
    ! call outfld('TROPP_T',   tropT(:ncol),      ncol, lchnk)
    ! call outfld('TROPP_Z',   tropZ(:ncol),      ncol, lchnk)
    ! call outfld('TROPP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    ! call outfld('TROPP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    ! call outfld('TROPP_FD',  tropFound(:ncol),  ncol, lchnk)


    ! Find the tropopause using just the cold point algorithm.
    call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, primary=TROP_ALG_CPP, backup=TROP_ALG_NONE)

    tropPdf(:,:) = 0._kind_phys
    tropFound(:) = 0._kind_phys
    tropDZ(:,:) = fillvalue 
    
    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._kind_phys
        tropFound(i) = 1._kind_phys
        tropDZ(i,:) = zm(i,:) - tropZ(i)
      end if
    end do

    ! call outfld('TROPF_P',   tropP(:ncol),      ncol, lchnk)
    ! call outfld('TROPF_T',   tropT(:ncol),      ncol, lchnk)
    ! call outfld('TROPF_Z',   tropZ(:ncol),      ncol, lchnk)
    ! call outfld('TROPF_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    ! call outfld('TROPF_PD',  tropPdf(:ncol, :), ncol, lchnk)
    ! call outfld('TROPF_FD',  tropFound(:ncol),  ncol, lchnk)
    
    
    ! If requested, do all of the algorithms.
    if (output_all) then
    
      do alg = 2, TROP_NALG
    
        ! Find the tropopause using just the analytic algorithm.
        call tropopause_find_run(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, primary=alg, backup=TROP_ALG_NONE)
  
        tropPdf(:,:) = 0._kind_phys
        tropFound(:) = 0._kind_phys
      
        do i = 1, ncol
          if (tropLev(i) /= NOTFOUND) then
            tropPdf(i, tropLev(i)) = 1._kind_phys
            tropFound(i) = 1._kind_phys
          end if
        end do
  
        ! call outfld('TROP' // TROP_LETTER(alg) // '_P',   tropP(:ncol),      ncol, lchnk)
        ! call outfld('TROP' // TROP_LETTER(alg) // '_T',   tropT(:ncol),      ncol, lchnk)
        ! call outfld('TROP' // TROP_LETTER(alg) // '_Z',   tropZ(:ncol),      ncol, lchnk)
        ! call outfld('TROP' // TROP_LETTER(alg) // '_PD',  tropPdf(:ncol, :), ncol, lchnk)
        ! call outfld('TROP' // TROP_LETTER(alg) // '_FD',  tropFound(:ncol),  ncol, lchnk)
      end do
    end if
    
    return
  end subroutine tropopause_output
end module tropopause_find
