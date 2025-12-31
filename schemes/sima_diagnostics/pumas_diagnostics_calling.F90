module pumas_diagnostics_calling

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: pumas_diagnostics_init ! init routine
   public :: pumas_diagnostics_run  ! main routine

CONTAINS



   !> \section arg_table_pumas_diagnostics_init  Argument Table
   !! \htmlinclude pumas_diagnostics_init.html
   subroutine pumas_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0



   end subroutine pumas_diagnostics_init

   !> \section arg_table_pumas_diagnostics_run  Argument Table
   !! \htmlinclude pumas_diagnostics_run.html
      subroutine pumas_diagnostics_run(proc_rates,qcsinksum_rate1ord,naai,npccn,rndst,nacon,tlat,qvlat,qctend,qitend,nctend,       &
                                       nitend,qrtend,qstend,                                                                       &
                                       nrtend,nstend,qgtend,ngtend,effc,effc_fn,effi,sadice,sadsnow,prect,preci,nevapr,am_evp_st,  &
                                       prain,cmeout,deffi,pgamrad,lamcrad,qsout,dsout,qgout,ngout,dgout,lflx,iflx,gflx,rflx,sflx,  &
                                       qrout,reff_rain,reff_snow,reff_grau,nrout,nsout,refl,arefl,areflz,frefl,csrfl,acsrfl,fcsrfl,&
                                       refl10cm, reflz10cm,rercld,ncai,ncal,qrout2,qsout2,nrout2,nsout2,drout2,dsout2,qgout2,      &
                                       ngout2,dgout2,                                                                              &
                                       freqg,freqs,freqr,nfice,qcrat,prer_evap,                                                    &
                                       errmsg, errflg)

      use cam_history, only: history_in_field

      type (proc_rates_type), intent(inout)  :: proc_rates

      real(kind_phys), intent(in) :: qcsinksum_rate1ord(:,:) ! 1st order rate for direct cw to precip conversion
      real(kind_phys), intent(in) :: tlat(:,:)               ! latent heating rate       (W/kg)
      real(kind_phys), intent(in) :: qvlat(:,:)              ! microphysical tendency qv (1/s)
      real(kind_phys), intent(in) :: qctend(:,:)             ! microphysical tendency qc (1/s)
      real(kind_phys), intent(in) :: qitend(:,:)             ! microphysical tendency qi (1/s)
      real(kind_phys), intent(in) :: nctend(:,:)             ! microphysical tendency nc (1/(kg*s))
      real(kind_phys), intent(in) :: nitend(:,:)             ! microphysical tendency ni (1/(kg*s))

      real(kind_phys), intent(in) :: qrtend(:,:)             ! microphysical tendency qr (1/s)
      real(kind_phys), intent(in) :: qstend(:,:)             ! microphysical tendency qs (1/s)
      real(kind_phys), intent(in) :: nrtend(:,:)             ! microphysical tendency nr (1/(kg*s))
      real(kind_phys), intent(in) :: nstend(:,:)             ! microphysical tendency ns (1/(kg*s))
      real(kind_phys), intent(in) :: qgtend(:,:)             ! microphysical tendency qg (1/s)
      real(kind_phys), intent(in) :: ngtend(:,:)             ! microphysical tendency ng (1/(kg*s))

      real(kind_phys), intent(in) :: effc(:,:)               ! droplet effective radius (micron)
      real(kind_phys), intent(in) :: effc_fn(:,:)            ! droplet effective radius, assuming nc = 1.e8 kg-1
      real(kind_phys), intent(in) :: effi(:,:)               ! cloud ice effective radius (micron)
      real(kind_phys), intent(in) :: sadice(:,:)             ! cloud ice surface area density (cm2/cm3)
      real(kind_phys), intent(in) :: sadsnow(:,:)            ! cloud snow surface area density (cm2/cm3)
      real(kind_phys), intent(in) :: prect(:)                ! surface precip rate (m/s)
      real(kind_phys), intent(in) :: preci(:)                ! surface ice/snow precip rate (m/s)
      real(kind_phys), intent(in) :: nevapr(:,:)             ! evaporation rate of rain + snow (kg/kg/s)
      real(kind_phys), intent(in) :: am_evp_st(:,:)          ! stratiform evaporation area (frac)
      real(kind_phys), intent(in) :: prain(:,:)              ! production of rain + snow (kg/kg/s)
      real(kind_phys), intent(in) :: cmeout(:,:)             ! Rate of cond-evap of ice (kg/kg/s)
      real(kind_phys), intent(in) :: deffi(:,:)              ! ice effective diameter for optics (radiation) (micron)
      real(kind_phys), intent(in) :: pgamrad(:,:)            ! ice gamma parameter for optics (radiation) (no units)
      real(kind_phys), intent(in) :: lamcrad(:,:)            ! slope of droplet distribution for optics (radiation) (1/m)
      real(kind_phys), intent(in) :: qsout(:,:)              ! snow mixing ratio (kg/kg)
      real(kind_phys), intent(in) :: dsout(:,:)              ! snow diameter (m)
      real(kind_phys), intent(in) :: lflx(:,:)               ! grid-box average liquid condensate flux (kg m^-2 s^-1)
      real(kind_phys), intent(in) :: iflx(:,:)               ! grid-box average ice condensate flux (kg m^-2 s^-1)
      real(kind_phys), intent(in) :: rflx(:,:)               ! grid-box average rain flux (kg m^-2 s^-1)
      real(kind_phys), intent(in) :: sflx(:,:)               ! grid-box average snow flux (kg m^-2 s^-1)
      real(kind_phys), intent(in) :: gflx(:,:)               ! grid-box average graupel/hail flux (kg m^-2 s^-1)

      real(kind_phys), intent(in) :: qrout(:,:)              ! grid-box average rain mixing ratio (kg/kg)
      real(kind_phys), intent(in) :: reff_rain(:,:)          ! rain effective radius (micron)
      real(kind_phys), intent(in) :: reff_snow(:,:)          ! snow effective radius (micron)
      real(kind_phys), intent(in) :: reff_grau(:,:)          ! graupel effective radius (micron)

      real(kind_phys), intent(in) :: nrout(:,:)              ! rain number concentration (1/m3)
      real(kind_phys), intent(in) :: nsout(:,:)              ! snow number concentration (1/m3)
      real(kind_phys), intent(in) :: refl(:,:)               ! analytic radar reflectivity (94GHZ, cloud radar) (dBZ)
      real(kind_phys), intent(in) :: arefl(:,:)              ! average reflectivity will zero points inside valid range (dBZ)
      real(kind_phys), intent(in) :: areflz(:,:)             ! average reflectivity in z. (mm6 m-3)
      real(kind_phys), intent(in) :: frefl(:,:)              ! fractional occurrence of radar reflectivity
      real(kind_phys), intent(in) :: csrfl(:,:)              ! cloudsat reflectivity (dBZ)
      real(kind_phys), intent(in) :: acsrfl(:,:)             ! cloudsat average (dBZ)
      real(kind_phys), intent(in) :: fcsrfl(:,:)             ! cloudsat fractional occurrence of radar reflectivity
      real(kind_phys), intent(in) :: refl10cm(:,:)           ! 10cm (rain) analytic radar reflectivity (dBZ)
      real(kind_phys), intent(in) :: reflz10cm(:,:)          ! 10cm (rain) analytic radar reflectivity (mm6 m-3)
      real(kind_phys), intent(in) :: rercld(:,:)             ! effective radius calculation for rain + cloud
      real(kind_phys), intent(in) :: ncai(:,:)               ! input number conc of ice nuclei available (1/m3)
      real(kind_phys), intent(in) :: ncal(:,:)               ! input number conc of CCN (1/m3)
      real(kind_phys), intent(in) :: qrout2(:,:)             ! copy of qrin as used to compute drin2
      real(kind_phys), intent(in) :: qsout2(:,:)             ! copy of qsin as used to compute dsin2
      real(kind_phys), intent(in) :: nrout2(:,:)             ! copy of nrin as used to compute drin2
      real(kind_phys), intent(in) :: nsout2(:,:)             ! copy of nsin as used to compute dsin2
      real(kind_phys), intent(in) :: drout2(:,:)             ! mean rain particle diameter (m)
      real(kind_phys), intent(in) :: dsout2(:,:)             ! mean snow particle diameter (m)
      real(kind_phys), intent(in) :: freqs(:,:)              ! fractional occurrence of snow
      real(kind_phys), intent(in) :: freqr(:,:)              ! fractional occurrence of rain
      real(kind_phys), intent(in) :: nfice(:,:)              ! fraction of frozen water to total condensed water
      real(kind_phys), intent(in) :: qcrat(:,:)              ! limiter for qc process rates (1=no limit --> 0. no qc)
      real(kind_phys), intent(in) :: qgout(:,:)              ! graupel/hail mixing ratio (kg/kg)
      real(kind_phys), intent(in) :: dgout(:,:)              ! graupel/hail diameter (m)
      real(kind_phys), intent(in) :: ngout(:,:)              ! graupel/hail number concentration (1/m3)
      real(kind_phys), intent(in) :: qgout2(:,:)             ! copy of qgin as used to compute dgin2
      real(kind_phys), intent(in) :: ngout2(:,:)             ! copy of ngin as used to compute dgin2
      real(kind_phys), intent(in) :: dgout2(:,:)             ! mean graupel/hail particle diameter (m)
      real(kind_phys), intent(in) :: freqg(:,:)              ! fractional occurrence of graupel
      real(kind_phys), intent(in) :: prer_evap(:,:)          ! evaporation rate of rain (kg/kg/s)

      ! CCPP error handling variables
      character(len=512), intent(in) :: errmsg
      integer,            intent(in) :: errflg

      errmsg = ''
      errflg = 0

!!!!!!!! DIAGNOSTIC CODE GOES HERE   !!!!!!!!!

!! The calls for history_out_field are in the pumas_diagnostcis.F90 file, but also with all of the other pumas code.
!! Those calls need to be hooked up to the parameters in the call list here.
!! There is an issue with subcolumn/grid scale I'm not sure how to address the grid averaging that is done after the
!! pumas_tend function call.

      end subroutine pumas_diagnostics_run

end module pumas_diagnostics_calling
