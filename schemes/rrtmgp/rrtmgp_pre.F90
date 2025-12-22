module rrtmgp_pre
 implicit none
 private

 public :: rrtmgp_pre_init
 public :: rrtmgp_pre_timestep_init
 public :: rrtmgp_pre_run
 public :: radiation_do_ccpp ! Public because it needs to be accessed elsewhere in CAM

CONTAINS

!> \section arg_table_rrtmgp_pre_init Argument Table
!! \htmlinclude rrtmgp_pre_init.html
!!
  subroutine rrtmgp_pre_init(nradgas, available_gases, gaslist, gaslist_lc, errmsg, errflg)
     use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
     use atmos_phys_string_utils, only: to_lower
     integer,                    intent(in) :: nradgas          ! Number of radiatively active gases
     character(len=5),           intent(in) :: gaslist(:)       ! Radiatively active gas list
     type(ty_gas_concs_ccpp),   intent(out) :: available_gases  ! Gas concentrations object
     character(len=5),          intent(out) :: gaslist_lc(:)    ! Lowercase verison of radiatively active gas list
     character(len=*),          intent(out) :: errmsg
     integer,                   intent(out) :: errflg

     ! Local variables
     integer :: idx

     ! Set error variables
     errmsg = ''
     errflg = 0

     ! Create lowercase version of the gaslist for RRTMGP.  The ty_gas_concs_ccpp objects
     ! work with CAM's uppercase names, but other objects that get input from the gas
     ! concs objects don't work.
     do idx = 1, nradgas
        gaslist_lc(idx) = to_lower(gaslist(idx))
     end do

     errmsg = available_gases%gas_concs%init(gaslist_lc)
     if (len_trim(errmsg) /= 0) then
        errflg = 1
     end if

  end subroutine rrtmgp_pre_init

!> \section arg_table_rrtmgp_pre_timestep_init Argument Table
!! \htmlinclude rrtmgp_pre_timestep_init.html
!!
  subroutine rrtmgp_pre_timestep_init(ncol, coszrs, nstep, dtime, iradsw, irad_always, offset, &
                  idxday, nday, idxnite, nnite, errmsg, errflg)
     use ccpp_kinds, only: kind_phys
     real(kind_phys),       intent(in)  :: coszrs(:)     ! Cosine solar zenith angle
     integer,               intent(in)  :: nstep        ! Current timestep number
     integer,               intent(in)  :: ncol         ! Number of horizontal columns
     real(kind_phys),       intent(in)  :: dtime        ! Timestep size
     integer,               intent(in)  :: iradsw       ! Freq. of shortwave radiation calc in time steps (positive) or hours (negative)
     integer,               intent(in)  :: irad_always  ! Number of time steps to execute radiation continuously
     integer,               intent(out) :: nday         ! Number of daylight columns
     integer,               intent(out) :: nnite        ! Number of nighttime columns
     integer,               intent(out) :: idxday(:)    ! Indices of daylight columns
     integer,               intent(out) :: idxnite(:)   ! Indices of nighttime columns
     integer,               intent(out) :: offset       ! Offset for next SW radiation timestep
     integer,               intent(out) :: errflg
     character(len=*),      intent(out) :: errmsg

     logical :: dosw_next
     integer :: nstepsw_next, idx

     ! Get timestep of next radiation calculation
     dosw_next = .false.
     nstepsw_next = nstep
     offset = 0
     do while (.not. dosw_next)
        nstepsw_next = nstepsw_next + 1
        offset = offset + dtime
        call radiation_do_ccpp('sw', nstepsw_next, iradsw, irad_always, dosw_next, errmsg, errflg)
        if (errflg /= 0) then
           return
        end if
     end do
     ! Gather night/day column indices.
     nday = 0
     nnite = 0
     idxday = 0
     idxnite = 0
     do idx = 1, ncol
        if ( coszrs(idx) > 0.0_kind_phys ) then
           nday = nday + 1
           idxday(nday) = idx
        else
           nnite = nnite + 1
           idxnite(nnite) = idx
        end if
     end do

  end subroutine rrtmgp_pre_timestep_init

!> \section arg_table_rrtmgp_pre_run Argument Table
!! \htmlinclude rrtmgp_pre_run.html
!!
  subroutine rrtmgp_pre_run(coszrs, nstep, dtime, iradsw, iradlw, irad_always, ncol, &
                  next_cday, idxday, nday, idxnite, nnite, nlay, nlwbands, nswbands, &
                  spectralflux, nextsw_cday, dosw, dolw, dosw_heat, dolw_heat, fsw,  &
                  fswc, flw, flwc, errmsg, errflg)
     use ccpp_kinds,              only: kind_phys
     use ccpp_fluxes,             only: ty_fluxes_broadband_ccpp
     use ccpp_fluxes_byband,      only: ty_fluxes_byband_ccpp
     ! Inputs
     real(kind_phys), dimension(:),    intent(in) :: coszrs        ! Cosine solar zenith angle
     real(kind_phys),                  intent(in) :: next_cday     ! The calendar day of the next timestep
     real(kind_phys),                  intent(in) :: dtime         ! Timestep size [s]
     integer,                          intent(in) :: nstep         ! Timestep number
     integer,                          intent(in) :: iradsw        ! Freq. of shortwave radiation calc in time steps (positive) or hours (negative)
     integer,                          intent(in) :: iradlw        ! Freq. of longwave radiation calc in time steps (positive) or hours (negative)
     integer,                          intent(in) :: irad_always   ! Number of time steps to execute radiation continuously
     integer,                          intent(in) :: ncol          ! Number of columns
     integer,                          intent(in) :: nlay          ! Number of vertical layers
     integer,                          intent(in) :: nlwbands      ! Number of longwave bands
     integer,                          intent(in) :: nswbands      ! Number of shortwave bands
     integer,                          intent(in) :: nday          ! Number of daylight columns
     integer,                          intent(in) :: nnite         ! Number of nighttime columns
     integer, dimension(:),            intent(in) :: idxday        ! Indices of daylight columns
     integer, dimension(:),            intent(in) :: idxnite       ! Indices of nighttime columns
     logical,                          intent(in) :: spectralflux  ! Flag to calculate fluxes (up and down) per band
     ! Outputs
     real(kind_phys),               intent(inout) :: nextsw_cday   ! The next calendar day during which calculation will be performed
     class(ty_fluxes_broadband_ccpp), intent(out) :: fswc          ! Clear-sky shortwave flux object
     class(ty_fluxes_byband_ccpp),    intent(out) :: fsw           ! All-sky shortwave flux object
     class(ty_fluxes_broadband_ccpp), intent(out) :: flwc          ! Clear-sky longwave flux object
     class(ty_fluxes_byband_ccpp),    intent(out) :: flw           ! All-sky longwave flux object
     logical,                         intent(out) :: dosw          ! Flag to do shortwave calculation
     logical,                         intent(out) :: dolw          ! Flag to do longwave calculation
     logical,                         intent(out) :: dosw_heat     ! Flag to calculate net shortwave heating
     logical,                         intent(out) :: dolw_heat     ! Flag to calculate net longwave heating
     character(len=*),                intent(out) :: errmsg
     integer,                         intent(out) :: errflg

     ! Local variables
     integer :: idx

     ! Set error variables
     errflg = 0
     errmsg = ''

     ! Determine if we're going to do longwave and/or shortwave this timestep
     call radiation_do_ccpp('sw', nstep, iradsw, irad_always, dosw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if
     call radiation_do_ccpp('lw', nstep, iradlw, irad_always, dolw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     dosw_heat = (.not. dosw)
     dolw_heat = (.not. dolw)

     ! determine if next radiation time-step not equal to next time-step
     if (nstep >= 1) then
        if (next_cday /= nextsw_cday) nextsw_cday = -1._kind_phys
     end if

     ! Allocate the flux arrays and init to zero.
     call initialize_rrtmgp_fluxes_byband(nday, nlay+1, nswbands, nswbands, spectralflux, fsw, errmsg, errflg, do_direct=.true.)
     if (errflg /= 0) then
        return
     end if
     call initialize_rrtmgp_fluxes_broadband(nday, nlay+1, nswbands, nswbands, spectralflux, fswc, errmsg, errflg, do_direct=.true.)
     if (errflg /= 0) then
        return
     end if
     call initialize_rrtmgp_fluxes_byband(ncol, nlay+1, nlwbands, nswbands, spectralflux, flw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if
     call initialize_rrtmgp_fluxes_broadband(ncol, nlay+1, nlwbands, nswbands, spectralflux, flwc, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

  end subroutine rrtmgp_pre_run

!================================================================================================

subroutine radiation_do_ccpp(op, nstep, irad, irad_always, radiation_do, errmsg, errflg)

   ! Return radiation_do set to .true. if the specified operation is done this timestep.

   character(len=*),  intent(in) :: op             ! name of operation
   integer,           intent(in) :: nstep
   integer,           intent(in) :: irad
   integer,           intent(in) :: irad_always
   integer,          intent(out) :: errflg
   character(len=*), intent(out) :: errmsg
   logical,          intent(out) :: radiation_do   ! return value

   !-----------------------------------------------------------------------

   ! Set error variables
   errflg = 0
   errmsg = ''

   select case (op)
      case ('sw') ! do a shortwave heating calc this timestep?
         radiation_do = nstep == 0  .or.  irad == 1                     &
                     .or. (mod(nstep-1,irad) == 0  .and.  nstep /= 1)   &
                     .or. nstep <= irad_always
      case ('lw') ! do a longwave heating calc this timestep?
         radiation_do = nstep == 0  .or.  irad == 1                     &
                     .or. (mod(nstep-1,irad) == 0  .and.  nstep /= 1)   &
                     .or. nstep <= irad_always
      case default
         errflg = 1
         errmsg = 'radiation_do_ccpp: unknown operation:'//op
   end select

end subroutine radiation_do_ccpp

!=========================================================================================

subroutine initialize_rrtmgp_fluxes_broadband(ncol, nlevels, nbands, nswbands, spectralflux, fluxes, errmsg, errflg, do_direct)
   use ccpp_fluxes,             only: ty_fluxes_broadband_ccpp

   ! Allocate flux arrays and set values to zero.

   ! Arguments
   integer,                    intent(in)    :: ncol, nlevels, nbands, nswbands
   logical,                    intent(in)    :: spectralflux
   class(ty_fluxes_broadband_ccpp), intent(inout) :: fluxes
   logical, optional,          intent(in)    :: do_direct
   character(len=*),           intent(out)   :: errmsg
   integer,                    intent(out)   :: errflg

   ! Local variables
   logical :: do_direct_local
   character(len=256) :: alloc_errmsg
   character(len=*), parameter :: sub = 'initialize_rrtmgp_fluxes'
   !----------------------------------------------------------------------------

   if (present(do_direct)) then
      do_direct_local = .true.
   else
      do_direct_local = .false.
   end if

   ! Broadband fluxes
   allocate(fluxes%fluxes%flux_up(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_up". Message: ', &
              alloc_errmsg
      return
   end if
   allocate(fluxes%fluxes%flux_dn(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_dn". Message: ', &
              alloc_errmsg
      return
   end if
   allocate(fluxes%fluxes%flux_net(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_net". Message: ', &
              alloc_errmsg
      return
   end if
   if (do_direct_local) then
      allocate(fluxes%fluxes%flux_dn_dir(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_dn_dir". Message: ', &
                 alloc_errmsg
         return
      end if
   end if

   ! Initialize
   call reset_fluxes_broadband(fluxes)

end subroutine initialize_rrtmgp_fluxes_broadband

!=========================================================================================

subroutine initialize_rrtmgp_fluxes_byband(ncol, nlevels, nbands, nswbands, spectralflux, fluxes, errmsg, errflg, do_direct)
   use ccpp_fluxes_byband,      only: ty_fluxes_byband_ccpp

   ! Allocate flux arrays and set values to zero.

   ! Arguments
   integer,                    intent(in)    :: ncol, nlevels, nbands, nswbands
   logical,                    intent(in)    :: spectralflux
   class(ty_fluxes_byband_ccpp), intent(inout) :: fluxes
   logical, optional,          intent(in)    :: do_direct
   character(len=*),           intent(out)   :: errmsg
   integer,                    intent(out)   :: errflg

   ! Local variables
   logical :: do_direct_local
   character(len=256) :: alloc_errmsg
   character(len=*), parameter :: sub = 'initialize_rrtmgp_fluxes_byband'
   !----------------------------------------------------------------------------

   if (present(do_direct)) then
      do_direct_local = .true.
   else
      do_direct_local = .false.
   end if

   ! Broadband fluxes
   allocate(fluxes%fluxes%flux_up(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_up". Message: ', &
              alloc_errmsg
      return
   end if
   allocate(fluxes%fluxes%flux_dn(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_dn". Message: ', &
              alloc_errmsg
      return
   end if
   allocate(fluxes%fluxes%flux_net(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_net". Message: ', &
              alloc_errmsg
      return
   end if
   if (do_direct_local) then
      allocate(fluxes%fluxes%flux_dn_dir(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_dn_dir". Message: ', &
                 alloc_errmsg
         return
      end if
   end if

   ! Fluxes by band always needed for SW.  Only allocate for LW
   ! when spectralflux is true.
   if (nbands == nswbands .or. spectralflux) then
      allocate(fluxes%fluxes%bnd_flux_up(ncol, nlevels, nbands), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%bnd_flux_up". Message: ', &
                alloc_errmsg
         return
      end if
      allocate(fluxes%fluxes%bnd_flux_dn(ncol, nlevels, nbands), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%bnd_flux_dn". Message: ', &
                 alloc_errmsg
         return
      end if
      allocate(fluxes%fluxes%bnd_flux_net(ncol, nlevels, nbands), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%bnd_flux_net". Message: ', &
                 alloc_errmsg
         return
      end if
      if (do_direct_local) then
         allocate(fluxes%fluxes%bnd_flux_dn_dir(ncol, nlevels, nbands), stat=errflg, errmsg=alloc_errmsg)
         if (errflg /= 0) then
            write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%bnd_flux_dn_dir". Message: ', &
                    alloc_errmsg
            return
         end if
      end if
   end if

   ! Initialize
   call reset_fluxes_byband(fluxes)

end subroutine initialize_rrtmgp_fluxes_byband

!=========================================================================================

subroutine reset_fluxes_broadband(fluxes)
   use ccpp_kinds,              only: kind_phys
   use ccpp_fluxes,             only: ty_fluxes_broadband_ccpp

   ! Reset flux arrays to zero.

   class(ty_fluxes_broadband_ccpp), intent(inout) :: fluxes
   !----------------------------------------------------------------------------

   ! Reset broadband fluxes
   fluxes%fluxes%flux_up(:,:) = 0._kind_phys
   fluxes%fluxes%flux_dn(:,:) = 0._kind_phys
   fluxes%fluxes%flux_net(:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%flux_dn_dir)) fluxes%fluxes%flux_dn_dir(:,:) = 0._kind_phys

end subroutine reset_fluxes_broadband

!=========================================================================================

subroutine reset_fluxes_byband(fluxes)
   use ccpp_kinds,              only: kind_phys
   use ccpp_fluxes_byband,      only: ty_fluxes_byband_ccpp

   ! Reset flux arrays to zero.

   class(ty_fluxes_byband_ccpp), intent(inout) :: fluxes
   !----------------------------------------------------------------------------

   ! Reset broadband fluxes
   fluxes%fluxes%flux_up(:,:) = 0._kind_phys
   fluxes%fluxes%flux_dn(:,:) = 0._kind_phys
   fluxes%fluxes%flux_net(:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%flux_dn_dir)) fluxes%fluxes%flux_dn_dir(:,:) = 0._kind_phys

   ! Reset band-by-band fluxes
   if (associated(fluxes%fluxes%bnd_flux_up)) fluxes%fluxes%bnd_flux_up(:,:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%bnd_flux_dn)) fluxes%fluxes%bnd_flux_dn(:,:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%bnd_flux_net)) fluxes%fluxes%bnd_flux_net(:,:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%bnd_flux_dn_dir)) fluxes%fluxes%bnd_flux_dn_dir(:,:,:) = 0._kind_phys

end subroutine reset_fluxes_byband

!=========================================================================================

end module rrtmgp_pre
