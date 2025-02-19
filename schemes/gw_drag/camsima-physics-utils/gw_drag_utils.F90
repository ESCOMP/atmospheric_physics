module gw_drag_utils
  !-------------------------------------------------------------------
  ! Support module for CAM-SIMA to read in gravity wave data.
  !
  ! Remarks: This module is not CCPP-ized, but has been cleaned up
  ! for use within CAM-SIMA, particularly removal of chunk support.
  !-------------------------------------------------------------------

  use ccpp_kinds,     only: kind_phys
  use runtime_obj,    only: unset_str
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc

  implicit none
  private

  public :: gw_readnl
  public :: gw_read_rdg_file
  public :: gw_read_movmtn_file
  public :: gw_read_beres_dp_file
  public :: gw_read_beres_sh_file
  public :: gw_read_topo_file

!> \section arg_table_gw_drag_utils  Argument Table
!! \htmlinclude gw_drag_utils.html
  ! days in year for climatological tropopause pressure data
  integer,         public, parameter   :: gw_slices = 12

  ! climatological tropopause pressures (pcols,ntimes)
  real(kind_phys), public, allocatable :: gw_p_loc(:,:)

  ! monthly day-of-year times corresponding to climatological data (12)
  real(kind_phys), public, allocatable :: gw_days(:)

  ! Private module data
  character(len=256)                   :: gw_rdg_file = unset_str

contains
  ! Read namelist variable gw_rdg_file.
  ! Containing this within CAM as otherwise the climo file has to be passed from physics -> here
  ! then the data from here -> physics. Keeping it at CAM level for now. (hplin, 8/16/24)
  subroutine gw_readnl(nlfile)
    use shr_nl_mod,      only: find_group_name => shr_nl_find_group_name
    use mpi,             only: mpi_character
    use spmd_utils,      only: mpicom

    ! filepath for file containing namelist input
    character(len=*), intent(in) :: nlfile

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'gw_readnl'

    namelist /gw_nl/ gw_rdg_file

    if (masterproc) then
       open(newunit=unitn, file=trim(nlfile), status='old')
       call find_group_name(unitn, 'gw_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, gw_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(gw_rdg_file, len(gw_rdg_file), mpi_character, 0, mpicom, ierr)

    ! Print out namelist variables
    if (masterproc) then
      write(iulog,*) subname, ' options:'
      write(iulog,*) '  Tropopause climatology file will be read: ', trim(gw_rdg_file)
    endif

  end subroutine gw_readnl

  subroutine gw_read_rdg_file()
    !------------------------------------------------------------------
    ! ... initialize upper boundary values
    !------------------------------------------------------------------
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish
    use physics_grid, only : get_rlat_all_p, get_rlon_all_p
    use physics_grid, only : pcols => columns_on_task
    use physconst,    only : pi
    use time_manager, only : get_calday
    use ioFileMod,    only : cam_get_file
    use cam_pio_utils, only: cam_pio_openfile
    use pio,          only : file_desc_t, var_desc_t, pio_inq_dimid, pio_inq_dimlen, &
         pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite

    ! character(len=256), intent(in) :: gw_rdg_file  ! absolute path of climatology file

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
    type(interp_type) :: lon_wgts, lat_wgts
    real(kind_phys), allocatable :: gw_p_in(:,:,:)
    real(kind_phys), allocatable :: lat(:)
    real(kind_phys), allocatable :: lon(:)
    real(kind_phys) :: to_lats(pcols), to_lons(pcols)
    real(kind_phys), parameter :: d2r=pi/180._kind_phys, zero=0._kind_phys, twopi=pi*2._kind_phys
    character(len=256) :: locfn

    !-----------------------------------------------------------------------
    !       ... open netcdf file
    !-----------------------------------------------------------------------
    call cam_get_file (gw_rdg_file, locfn, allow_fail=.false.)
    call cam_pio_openfile(pio_id, trim(locfn), PIO_NOWRITE)

    !-----------------------------------------------------------------------
    !       ... get time dimension
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'time', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, ntimes )
    if( ntimes /= 12 )then
       write(iulog,*) 'gw_read_rdg_file: number of months = ',ntimes,'; expecting 12'
       call endrun('gw_read_rdg_file: number of months in file incorrect')
    end if
    !-----------------------------------------------------------------------
    !       ... get latitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lat', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlat )
    allocate( lat(nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'gw_read_rdg_file: lat allocation error = ',ierr
       call endrun('gw_read_rdg_file: failed to allocate lat')
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
       write(iulog,*) 'gw_read_rdg_file: lon allocation error = ',ierr
       call endrun('gw_read_rdg_file: failed to allocate lon')
    end if
    ierr = pio_inq_varid( pio_id, 'lon', vid )
    ierr = pio_get_var( pio_id, vid, lon )
    lon(:nlon) = lon(:nlon) * d2r

    !------------------------------------------------------------------
    !  ... allocate arrays
    !------------------------------------------------------------------
    allocate( gw_p_in(nlon,nlat,ntimes), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'gw_read_rdg_file: gw_p_in allocation error = ',ierr
       call endrun('gw_read_rdg_file: failed to allocate gw_p_in')
    end if
    !------------------------------------------------------------------
    !  ... read in the tropopause pressure
    !------------------------------------------------------------------
    ierr = pio_inq_varid( pio_id, 'trop_p', vid )
    start = (/ 1, 1, 1 /)
    count = (/ nlon, nlat, ntimes /)
    ierr = pio_get_var( pio_id, vid, start, count, gw_p_in )

    !------------------------------------------------------------------
    !  ... close the netcdf file
    !------------------------------------------------------------------
    call pio_closefile( pio_id )

    !--------------------------------------------------------------------
    !  ... regrid
    !--------------------------------------------------------------------

    ! gw_p_loc is allocated with dimensions (pcols, begchunk:endchunk, ntimes) in CAM.
    ! in CAM-SIMA, the chunk dimension is collapsed as it is unused.
    allocate( gw_p_loc(pcols,ntimes), stat=ierr )

    if( ierr /= 0 ) then
      write(iulog,*) 'gw_read_rdg_file: gw_p_loc allocation error = ',ierr
      call endrun('gw_read_rdg_file: failed to allocate gw_p_loc')
    end if

    call get_rlat_all_p(pcols, to_lats)
    call get_rlon_all_p(pcols, to_lons)
    call lininterp_init(lon, nlon, to_lons, pcols, 2, lon_wgts, zero, twopi)
    call lininterp_init(lat, nlat, to_lats, pcols, 1, lat_wgts)
    do n=1,ntimes
       call lininterp(gw_p_in(:,:,n), nlon, nlat, gw_p_loc(1:pcols,n), pcols, lon_wgts, lat_wgts)
    end do
    call lininterp_finish(lon_wgts)
    call lininterp_finish(lat_wgts)
    deallocate(lon)
    deallocate(lat)
    deallocate(gw_p_in)

    !--------------------------------------------------------
    ! ... initialize the monthly day of year times
    !--------------------------------------------------------

    allocate( gw_days(gw_slices), stat=ierr )
    if( ierr /= 0 ) then
      write(iulog,*) 'gw_read_rdg_file: gw_days allocation error = ',ierr
      call endrun('gw_read_rdg_file: failed to allocate gw_days')
    end if

    do n = 1,12
       gw_days(n) = get_calday( dates(n), 0 )
    end do
    if (masterproc) then
       write(iulog,*) 'gw_read_rdg_file : gw_days'
       write(iulog,'(1p,5g15.8)') gw_days(:)
    endif

  end subroutine gw_read_rdg_file
end module gw_drag_utils
