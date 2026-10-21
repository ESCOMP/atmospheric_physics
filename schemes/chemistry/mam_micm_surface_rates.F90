! MAM aerosol surface area density -> MICM surface reaction rate parameters
!
! Every MICM surface reaction carries two rate parameters per grid cell,
! '<label>.effective radius [m]' and '<label>.particle number concentration
! [# m-3]'. This scheme fills them from one MAM mode's surface area density
! and wet number-mode diameter (aero_state%surf_area_dens, the sfc/dm_aer pair
! CAM's mo_usrrxt hetrxtrate consumes), so MICM's rate equals CAM's per-mode
! heterogeneous rate.
!
! CAM sums hetrxtrate over modes; a MICM surface reaction has one (r, N) pair.
! A mechanism reproduces CAM with one surface reaction per mode, and the mode
! is taken from the label's '_a<mm>' suffix (CAM's per-mode naming):
! 'SURF.N2O5_hydrolysis_a1' reads mode 1. A surface reaction without that
! suffix, or naming a mode MAM does not have, is an error at init.
!
! It must run after musica_ccpp_photolysis (which resets every rate parameter)
! and before musica_ccpp_chemistry. mam_surfarea_diagnostics, placed between
! this scheme and the solve, checks the result.
module mam_micm_surface_rates
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: mam_micm_surface_rates_init
   public :: mam_micm_surface_rates_run

   ! CAM's default sad_chem_spec_types for MAM (aerosol_nl namelist in CAM),
   ! the same list mam_surfarea_diagnostics uses.
   character(len=9), parameter :: sad_chem_spec_types(5) = &
        [ character(len=9) :: 'sulfate', 's-organic', 'p-organic', 'black-c', 'ammonium' ]

   integer :: nmodes = 0
   character(len=:), allocatable :: surf_label(:)
   integer,          allocatable :: surf_iradius(:) ! slot of the effective radius [m]
   integer,          allocatable :: surf_inumber(:) ! slot of the particle number concentration [# m-3]
   integer,          allocatable :: surf_mode(:)    ! MAM mode feeding the reaction

contains

!> \section arg_table_mam_micm_surface_rates_init Argument Table
!! \htmlinclude mam_micm_surface_rates_init.html
   subroutine mam_micm_surface_rates_init(errmsg, errflg)
      use mam_mode_metadata,        only: ntot_amode_val
      use musica_ccpp_micm,         only: rate_parameters_ordering
      use musica_ccpp_micm_surface, only: find_surface_reactions_in_ordering, mode_from_label

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: j

      errmsg = ''
      errflg = 0

      nmodes = ntot_amode_val

      call find_surface_reactions_in_ordering(rate_parameters_ordering, surf_label, &
           surf_iradius, surf_inumber, errmsg, errflg)
      if (errflg /= 0) return

      allocate(surf_mode(size(surf_label)))
      do j = 1, size(surf_label)
         surf_mode(j) = mode_from_label(surf_label(j))
         if (surf_mode(j) < 1 .or. surf_mode(j) > nmodes) then
            errflg = 1
            write(errmsg, '(a,i0,a)') 'mam_micm_surface_rates_init: MICM surface reaction "'// &
                 trim(surf_label(j))//'" does not name a MAM mode 1..', nmodes, &
                 ' with an _a<mm> label suffix'
            return
         end if
      end do

   end subroutine mam_micm_surface_rates_init

!> \section arg_table_mam_micm_surface_rates_run Argument Table
!! \htmlinclude mam_micm_surface_rates_run.html
   subroutine mam_micm_surface_rates_run(ncol, pver, troplev, t, pmid, pi, rate_parameters, &
        errmsg, errflg)
      use aerosol_instances_mod,    only: aerosol_instances_get_props, &
                                          aerosol_instances_get_state, &
                                          aerosol_instances_get_num_models
      use aerosol_properties_mod,   only: aerosol_properties
      use aerosol_state_mod,        only: aerosol_state
      use musica_ccpp_micm_surface, only: surface_reaction_parameters_from_mode

      integer,          intent(in)    :: ncol
      integer,          intent(in)    :: pver
      integer,          intent(in)    :: troplev(:)  ! (ncol) tropopause vertical layer index
      real(kind_phys),  intent(in)    :: t(:,:)      ! (ncol,pver) air temperature [K]
      real(kind_phys),  intent(in)    :: pmid(:,:)   ! (ncol,pver) air pressure [Pa]
      real(kind_phys),  intent(in)    :: pi
      real(kind_phys),  intent(inout) :: rate_parameters(:,:,:)  ! (ncol,pver,nparams) MICM rate parameters
      character(len=*), intent(out)   :: errmsg
      integer,          intent(out)   :: errflg

      integer :: iaermod, j
      class(aerosol_properties), pointer :: aero_props
      class(aerosol_state),      pointer :: aero_state_obj

      integer         :: beglev(ncol), endlev(ncol)
      real(kind_phys) :: relhum(ncol, pver)
      real(kind_phys) :: sad(ncol, pver)             ! surface area density [cm2/cm3]
      real(kind_phys) :: reff(ncol, pver)            ! effective radius [cm]
      real(kind_phys) :: sfc(ncol, pver, nmodes)     ! per-mode surface area density [cm2/cm3]
      real(kind_phys) :: dm_aer(ncol, pver, nmodes)  ! per-mode wet number-mode diameter [cm]

      errmsg = ''
      errflg = 0

      if (size(surf_label) == 0) return

      ! Find MAM properties and state from aerosol instances
      aero_props => null()
      aero_state_obj => null()
      do iaermod = 1, aerosol_instances_get_num_models()
        aero_props => aerosol_instances_get_props(iaermod, 0)
        if (associated(aero_props)) then
          if (aero_props%model_is('MAM')) then
            aero_state_obj => aerosol_instances_get_state(iaermod, list_idx=0)
            exit
          end if
        end if
        aero_props => null()
      end do

      if (.not. associated(aero_props) .or. &
          .not. associated(aero_state_obj)) then
         errflg = 1
         errmsg = 'mam_micm_surface_rates_run: MAM aerosol model not found in aerosol_instances'
         return
      end if

      ! Tropospheric levels only, as CAM's aero_model_surfarea:
      ! from the level below the tropopause down to the surface.
      beglev(:ncol) = troplev(:ncol) + 1
      endlev(:ncol) = pver

      ! The modal implementation of surf_area_dens ignores relhum (it uses the
      ! wet number-mode diameter from water uptake).
      relhum(:,:) = 0.0_kind_phys
      sfc(:,:,:) = 0.0_kind_phys
      dm_aer(:,:,:) = 0.0_kind_phys

      call aero_state_obj%surf_area_dens(aero_props, sad_chem_spec_types, ncol, pver, &
           beglev, endlev, relhum, pmid, t, pi, sad, reff, sfc, dm_aer)

      do j = 1, size(surf_label)
         call surface_reaction_parameters_from_mode(sfc(:,:,surf_mode(j)), dm_aer(:,:,surf_mode(j)), &
              troplev(:ncol), rate_parameters(:ncol,:,surf_iradius(j)), &
              rate_parameters(:ncol,:,surf_inumber(j)), errmsg, errflg)
         if (errflg /= 0) then
            errmsg = 'mam_micm_surface_rates_run: MICM surface reaction "'// &
                 trim(surf_label(j))//'": '//trim(errmsg)
            return
         end if
      end do

   end subroutine mam_micm_surface_rates_run

end module mam_micm_surface_rates
