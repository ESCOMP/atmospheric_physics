module rrtmgp_sw_mcica_subcol_gen

implicit none
private
save

public :: rrtmgp_sw_mcica_subcol_gen_run

!==================================================================================================
contains
!==================================================================================================

!>
!> \section arg_table_rrtmgp_sw_mcica_subcol_gen_run Argument Table
!! \htmlinclude rrtmgp_sw_mcica_subcol_gen_run.html
subroutine rrtmgp_sw_mcica_subcol_gen_run(dosw, kdist_sw, nswbands, nswgpts, nday, nlay, &
              pver, tiny, idxday, ktopcam, ktoprad, cldfprime, c_cld_tau,       &
              c_cld_tau_w, c_cld_tau_w_g, cloud_sw, pmid_day, errmsg, errflg)
   use ccpp_kinds,              only: kind_phys
   use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
   use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_optical_props,      only: ty_optical_props_2str_ccpp
   use shr_RandNum_mod,         only: ShrKissRandGen ! SIMA-specific randum number generator

   ! Compute combined cloud optical properties.
   ! Create MCICA stochastic arrays for cloud SW optical properties.
   ! Initialize optical properties object (cloud_sw) and load with MCICA columns.

   ! arguments
   class(ty_gas_optics_rrtmgp_ccpp), intent(in)  :: kdist_sw             ! shortwave gas optics object
   integer,                          intent(in)  :: nswbands             ! number of shortwave bands
   integer,                          intent(in)  :: nswgpts              ! number of shortwave g-points
   integer,                          intent(in)  :: nlay                 ! number of layers in radiation calculation (may include "extra layer")
   integer,                          intent(in)  :: nday                 ! number of daylight columns
   integer,                          intent(in)  :: pver                 ! total number of vertical layers
   integer,                          intent(in)  :: ktopcam              ! index in CAM arrays of top level (layer or interface) at which RRTMGP is active
   integer,                          intent(in)  :: ktoprad              ! index in RRTMGP array corresponding to top layer or interface of CAM arrays
   integer,                          intent(in)  :: idxday(:)            ! indices of daylight columns in the chunk
   real(kind_phys),                  intent(in)  :: tiny                 ! definition of tiny in RRTMGP
   real(kind_phys),                  intent(in)  :: c_cld_tau(:,:,:)     ! combined cloud extinction optical depth
   real(kind_phys),                  intent(in)  :: c_cld_tau_w(:,:,:)   ! combined cloud single scattering albedo * tau
   real(kind_phys),                  intent(in)  :: c_cld_tau_w_g(:,:,:) ! combined cloud asymmetry parameter * w * tau
   real(kind_phys),                  intent(in)  :: cldfprime(:,:)       ! combined cloud fraction
   real(kind_phys),                  intent(in)  :: pmid_day(:,:)        ! air ressure at mid-points [Pa]
   logical,                          intent(in)  :: dosw                 ! Flag to do shortwave radiation this timestep

   type(ty_optical_props_2str_ccpp), intent(out) :: cloud_sw             ! SW cloud optical properties object
   character(len=512),               intent(out) :: errmsg
   integer,                          intent(out) :: errflg

   ! Local variables

   integer :: i, k, n
   integer :: igpt, nver, isubcol
   integer :: istat
   integer, parameter :: changeseed = 1
   character(len=256) :: alloc_errmsg
   type(ShrKissRandGen) :: kiss_gen  ! KISS RNG object
   integer  :: kiss_seed(nday,4)
   real(kind_phys) :: rand_num_1d(nday,1)   ! random number (kissvec)
   real(kind_phys) :: rand_num(nday,pver-ktopcam+1)   ! random number (kissvec)
   logical  :: iscloudy(nswgpts,nday,pver-ktopcam+1)   ! flag that says whether a gridbox is cloudy
   real(kind_phys), parameter :: cldmin = 1.0e-80_kind_phys  ! min cloud fraction
   real(kind_phys) :: cdf(nswgpts,nday,pver-ktopcam+1)
   real(kind_phys) :: cldfrac(nday,pver-ktopcam+1) ! Cloud fraction clipped to cldmin

   ! Arrays for converting from CAM chunks to RRTMGP inputs.
   real(kind_phys), allocatable :: cldf(:,:)
   real(kind_phys), allocatable :: tauc(:,:,:)
   real(kind_phys), allocatable :: ssac(:,:,:)
   real(kind_phys), allocatable :: asmc(:,:,:)
   real(kind_phys), allocatable :: taucmcl(:,:,:)
   real(kind_phys), allocatable :: ssacmcl(:,:,:)
   real(kind_phys), allocatable :: asmcmcl(:,:,:)
   real(kind_phys), allocatable :: day_cld_tau(:,:,:)
   real(kind_phys), allocatable :: day_cld_tau_w(:,:,:)
   real(kind_phys), allocatable :: day_cld_tau_w_g(:,:,:)
   !--------------------------------------------------------------------------------

   ! if no daylight columns the cloud_sw object isn't initialized
   if (nday > 0 .and. dosw) then

      ! Initialize object for SW cloud optical properties.
      errmsg = cloud_sw%optical_props%alloc_2str(nday, nlay, kdist_sw%gas_props)
      if (len_trim(errmsg) > 0) then
         errflg = 1
         return
      end if

      ! number of CAM's layers in radiation calculation.  Does not include the "extra layer".
      nver = pver - ktopcam + 1

      allocate( &
         cldf(nday,nver),                      &
         day_cld_tau(nswbands,nday,nver),      &
         day_cld_tau_w(nswbands,nday,nver),    &
         day_cld_tau_w_g(nswbands,nday,nver),  &
         tauc(nswbands,nday,nver), taucmcl(nswgpts,nday,nver), &
         ssac(nswbands,nday,nver), ssacmcl(nswgpts,nday,nver), &
         asmc(nswbands,nday,nver), asmcmcl(nswgpts,nday,nver), stat=istat, errmsg=alloc_errmsg)
      if (istat /= 0) then
         errflg = 1
         write(errmsg,*) 'rrtmgp_sw_mcica_subcol_gen_run: failed to allocate variable(s) - message: ', alloc_errmsg
         return
      end if

      ! Subset data so just the daylight columns, and the number of CAM layers in the
      ! radiation calculation are used by MCICA to produce subcolumns.
      cldf            = cldfprime(       idxday(1:nday), ktopcam:)
      day_cld_tau     = c_cld_tau(    :, idxday(1:nday), ktopcam:)
      day_cld_tau_w   = c_cld_tau_w(  :, idxday(1:nday), ktopcam:)
      day_cld_tau_w_g = c_cld_tau_w_g(:, idxday(1:nday), ktopcam:)

      ! Compute the optical properties needed for the 2-stream calculations.  These calculations
      ! are the same as the RRTMG version.

      ! set cloud optical depth, clip @ zero
      tauc = merge(day_cld_tau, 0.0_kind_phys, day_cld_tau > 0.0_kind_phys)
      ! set value of asymmetry
      asmc = merge(day_cld_tau_w_g / max(day_cld_tau_w, tiny), 0.0_kind_phys, day_cld_tau_w > 0.0_kind_phys)
      ! set value of single scattering albedo
      ssac = merge(max(day_cld_tau_w, tiny) / max(tauc, tiny), 1.0_kind_phys , tauc > 0.0_kind_phys)
      ! set asymmetry to zero when tauc = 0
      asmc = merge(asmc, 0.0_kind_phys, tauc > 0.0_kind_phys)

      ! clip cloud fraction
      cldfrac(:,:) = cldf(:nday,:)
      where (cldfrac(:,:) < cldmin)
         cldfrac(:,:) = 0._kind_phys
      end where

      ! Create a seed that depends on the state of the columns.
      ! Use pmid_day from bottom four layers.
      do i = 1, nday
         kiss_seed(i,1) = (pmid_day(i,nlay)   - int(pmid_day(i,nlay)))    * 1000000000
         kiss_seed(i,2) = (pmid_day(i,nlay-1) - int(pmid_day(i,nlay-1)))  * 1000000000
         kiss_seed(i,3) = (pmid_day(i,nlay-2) - int(pmid_day(i,nlay-2)))  * 1000000000
         kiss_seed(i,4) = (pmid_day(i,nlay-3) - int(pmid_day(i,nlay-3)))  * 1000000000
      end do

      ! create the RNG object
      kiss_gen = ShrKissRandGen(kiss_seed)

      ! Advance randum number generator by changeseed values
      do i = 1, changeSeed
         call kiss_gen%random(rand_num_1d)
      end do

      ! Generate random numbers in each subcolumn at every level
      do isubcol = 1,nswgpts
         call kiss_gen%random(rand_num)
         cdf(isubcol,:,:) = rand_num(:,:)
      enddo

      ! Maximum-Random overlap
      ! i) pick a random number for top layer.
      ! ii) walk down the column:
      !    - if the layer above is cloudy, use the same random number as in the layer above
      !    - if the layer above is clear, use a new random number

      do k = 2, nver
         do i = 1, nday
            do isubcol = 1, nswgpts
               if (cdf(isubcol,i,k-1) > 1._kind_phys - cldfrac(i,k-1) ) then
                  cdf(isubcol,i,k) = cdf(isubcol,i,k-1)
               else
                  cdf(isubcol,i,k) = cdf(isubcol,i,k) * (1._kind_phys - cldfrac(i,k-1))
               end if
            end do
         end do
      end do

      do k = 1, nver
         iscloudy(:,:,k) = (cdf(:,:,k) >= 1._kind_phys - spread(cldfrac(:,k), dim=1, nCopies=nswgpts) )
      end do

      ! -- generate subcolumns for homogeneous clouds -----
      ! where there is a cloud, set the subcolumn cloud properties;
      ! incoming tauc should be in-cloud quantites and not grid-averaged quantities
      do k = 1,nver
         do i = 1,nday
            do isubcol = 1,nswgpts
               if (iscloudy(isubcol,i,k) .and. (cldfrac(i,k) > 0._kind_phys) ) then
                  n = kdist_sw%gas_props%convert_gpt2band(isubcol)
                  taucmcl(isubcol,i,k) = tauc(n,i,k)
                  ssacmcl(isubcol,i,k) = ssac(n,i,k)
                  asmcmcl(isubcol,i,k) = asmc(n,i,k)
               else
                  taucmcl(isubcol,i,k) = 0._kind_phys
                  ssacmcl(isubcol,i,k) = 1._kind_phys
                  asmcmcl(isubcol,i,k) = 0._kind_phys
               end if
            end do
         end do
      end do

      call kiss_gen%finalize()

      ! If there is an extra layer in the radiation then this initialization
      ! will provide the optical properties there.
      cloud_sw%optical_props%tau = 0.0_kind_phys
      cloud_sw%optical_props%ssa = 1.0_kind_phys
      cloud_sw%optical_props%g   = 0.0_kind_phys

      ! Set the properties on g-points.
      do igpt = 1,nswgpts
         cloud_sw%optical_props%g  (:,ktoprad:,igpt) = asmcmcl(igpt,:,:)
         cloud_sw%optical_props%ssa(:,ktoprad:,igpt) = ssacmcl(igpt,:,:)
         cloud_sw%optical_props%tau(:,ktoprad:,igpt) = taucmcl(igpt,:,:)
      end do

      ! validate checks that: tau > 0, ssa is in range [0,1], and g is in range [-1,1].
      errmsg = cloud_sw%optical_props%validate()
      if (len_trim(errmsg) > 0) then
         errflg = 1
         return
      end if

      ! delta scaling adjusts for forward scattering
      errmsg = cloud_sw%optical_props%delta_scale()
      if (len_trim(errmsg) > 0) then
         errflg = 1
         return
      end if

      ! All information is in cloud_sw, now deallocate local vars.
      deallocate( &
         cldf, tauc, ssac, asmc, &
         taucmcl, ssacmcl, asmcmcl,&
         day_cld_tau, day_cld_tau_w, day_cld_tau_w_g )

   end if

end subroutine rrtmgp_sw_mcica_subcol_gen_run

end module rrtmgp_sw_mcica_subcol_gen
