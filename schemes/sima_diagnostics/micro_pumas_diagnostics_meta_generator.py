#!/usr/bin/env python3
"""
Generate micro_pumas_diagnostics.meta, the CCPP metadata for the PUMAS diagnostics
schemes in micro_pumas_diagnostics.F90 (micro_pumas_pre_diagnostics
and micro_pumas_diagnostics).

The metadata has well over 100 variables, most of which share a few dimension and
type patterns, so it is generated rather than edited by hand.

Usage (from this directory):
    python3 micro_pumas_diagnostics_meta_generator.py [output_file]

The output file defaults to micro_pumas_diagnostics.meta in this directory.

To add or change an argument:
    1. Add/update its entry in the variable table (V), which maps the Fortran
       argument name to (standard_name, units, dimensions, type)
    2. Add/move the argument name in the argument list for the routine
       (pre, mini, mrun), in the SAME ORDER as the Fortran argument list
    3. Rerun this script

Comment lines which need to appear in an entry (e.g. FOR_REVIEWER_FEEDBACK
notes) go in EXTRA.
"""
import os
import sys
L2 = "(horizontal_loop_extent, vertical_layer_dimension)"
LI = "(horizontal_loop_extent, vertical_interface_dimension)"
L1 = "(horizontal_loop_extent)"
S  = "()"
R  = "real | kind = kind_phys"
PR = "real | kind = pumas_r8"
LP = "(microphysics_horizontal_loop_extent, microphysics_vertical_layer_dimension)"
I  = "integer"
C  = "character | kind = len=*"
q  = "_mixing_ratio_wrt_moist_air_and_condensed_water"
n  = "mass_number_concentration_of_"
nq = "_in_moist_air_and_condensed_water"
# Constituents must be marked as advected so that the CCPP framework can find them
ADVECTED = {"q", "cldliq", "cldice", "numliq", "numice", "qr", "nr", "qs", "ns", "qg", "ng"}
INOUT = {"acprecl", "acgcme", "acnum"}
def v(name, sn, units, dims, typ=R, intent="in", extra=None):
    if name in INOUT: intent = "inout"
    s = f"[ {name} ]\n  standard_name = {sn}\n"
    if extra: s += extra
    s += f"  units = {units}\n  type = {typ}\n  dimensions = {dims}\n  intent = {intent}\n"
    if name in ADVECTED: s += "  advected = true\n"
    return s
def err():
    return [v("errmsg","ccpp_error_message","none",S,"character | kind = len=512","out"),
            v("errflg","ccpp_error_code","1",S,I,"out")]
def hdr(scheme):
    return (f"[ccpp-table-properties]\n  name = {scheme}\n  type = scheme\n"
            "  dependencies = ../pumas/pumas/pumas_kinds.F90,../pumas/pumas/micro_pumas_utils.F90\n"
            "  dependencies = ../pumas/pumas/pumas_gamma_function.F90,../pumas/pumas/micro_pumas_diags.F90\n"
            "  dependencies = ../pumas/pumas/pumas_stochastic_collect_tau.F90\n"
            "  dependencies = ../../to_be_ccppized/wv_saturation.F90\n"
            "  kind_spec = pumas_kinds:pumas_r8=>kind_r8\n\n")
def tbl(name):
    return f"[ccpp-arg-table]\n  name  = {name}\n  type  = scheme\n"
V = {
 "ncol": ("horizontal_loop_extent","count",S,I),
 "pver": ("vertical_layer_dimension","count",S,I),
 "pcols": ("horizontal_dimension","count",S,I),
 "top_lev": ("vertical_layer_index_of_troposphere_cloud_physics_top","index",S,I),
 "col_start": ("horizontal_loop_begin","count",S,I),
 "col_end": ("horizontal_loop_end","count",S,I),
 "dtime": ("timestep_for_physics","s",S,R),
 "gravit": ("standard_gravitational_acceleration","m s-2",S,R),
 "rair": ("gas_constant_of_dry_air","J kg-1 K-1",S,R),
 "tmelt": ("freezing_point_of_water","K",S,R),
 "rhoh2o": ("fresh_liquid_water_density_at_0c","kg m-3",S,R),
 "pi": ("pi_constant","1",S,R),
 "cpairv": ("composition_dependent_specific_heat_of_dry_air_at_constant_pressure","J kg-1 K-1",L2,R),
 "micro_mg_warm_rain": ("control_for_warm_rain_method","none",S,C),
 "pmid": ("air_pressure","Pa",L2,R),
 "pdel": ("air_pressure_thickness","Pa",L2,R),
 "temp": ("air_temperature","K",L2,R),
 "zm": ("geopotential_height_wrt_surface","m",L2,R),
 "q": ("water_vapor"+q,"kg kg-1",L2,R),
 "micro_mg_adjust_cpt": ("do_cold_point_tropopause_temperature_adjustment_for_microphysics","flag",S,"logical"),
 "tropLev_cpp": ("tropopause_vertical_layer_index_from_cold_point_method","index",L1,I),
 "tropT_cpp": ("tropopause_air_temperature_from_cold_point_method","K",L1,R),
 "tropZ_cpp": ("tropopause_geopotential_height_wrt_surface_from_cold_point_method","m",L1,R),
 "cldliq": ("cloud_liquid_water"+q,"kg kg-1",L2,R),
 "cldice": ("cloud_ice"+q,"kg kg-1",L2,R),
 "numliq": (n+"cloud_liquid_water_droplets"+nq,"kg-1",L2,R),
 "numice": (n+"cloud_ice_water_crystals"+nq,"kg-1",L2,R),
 "qr": ("rain"+q,"kg kg-1",L2,R),
 "nr": (n+"rain_drops"+nq,"kg-1",L2,R),
 "qs": ("snow"+q,"kg kg-1",L2,R),
 "ns": (n+"snow_crystals"+nq,"kg-1",L2,R),
 "qg": ("graupel_water"+q,"kg kg-1",L2,R),
 "ng": (n+"graupel_particles"+nq,"kg-1",L2,R),
 "ast": ("stratiform_cloud_area_fraction","fraction",L2,R),
 "cld": ("cloud_area_fraction","fraction",L2,R),
 "concld": ("convective_cloud_area_fraction","fraction",L2,R),
 "cmeliq": ("rate_of_condensation_minus_evaporation_for_cloud_liquid_water"+q,"kg kg-1 s-1",L2,R),
 "naai": ("tendency_of_mass_number_concentration_of_activated_ice_nuclei","kg-1 s-1",L2,R),
 "naai_hom": ("tendency_of_mass_number_concentration_of_activated_ice_nuclei_due_to_homogeneous_freezing","kg-1 s-1",L2,R),
 "prec_dp": ("lwe_precipitation_rate_at_surface_due_to_deep_convection","m s-1",L1,R),
 "prec_sh": ("lwe_precipitation_rate_at_surface_due_to_shallow_convection","m s-1",L1,R),
 "calday": ("fractional_calendar_days_on_end_of_current_timestep","days",S,R),
 "lat": ("latitude","rad",L1,R),
 "lon": ("longitude","rad",L1,R),
 "eccen": ("planet_orbital_eccentricity_factor","1",S,R),
 "obliqr": ("planet_obliquity","rad",S,R),
 "lambm0": ("mean_longitude_of_perihelion_at_vernal_equinox","rad",S,R),
 "mvelpp": ("moving_vernal_equinox_longitude_of_perihelion_plus_pi","rad",S,R),
 "proc_rates": ("microphysics_process_rates","none",S,"proc_rates_type"),
 "tlat": ("tendency_of_dry_air_enthalpy_at_constant_pressure","J kg-1 s-1",L2,R),
 "qvlat": ("tendency_of_pumas_water_vapor"+q,"kg kg-1 s-1",LP,PR),
 "qcten": ("tendency_of_pumas_cloud_liquid_water"+q,"kg kg-1 s-1",LP,PR),
 "qiten": ("tendency_of_pumas_cloud_ice"+q,"kg kg-1 s-1",LP,PR),
 "ncten": ("tendency_of_pumas_"+n+"cloud_liquid_water_droplets"+nq,"kg-1 s-1",LP,PR),
 "niten": ("tendency_of_pumas_"+n+"cloud_ice_water_crystals"+nq,"kg-1 s-1",LP,PR),
 "effi": ("effective_radius_of_stratiform_cloud_ice_crystal","um",L2,R),
 "rel": ("effective_radius_of_stratiform_cloud_liquid_water_droplet","um",L2,R),
 "des": ("effective_diameter_of_stratiform_snow_crystal_for_radiation","um",L2,R),
 "sadice": ("cloud_ice_surface_area_density","cm2 cm-3",L2,R),
 "sadsnow": ("snow_surface_area_density","cm2 cm-3",L2,R),
 "prect": ("lwe_large_scale_precipitation_rate_at_surface","m s-1",L1,R),
 "nevapr": ("precipitation_evaporation_rate_wrt_moist_air_and_condensed_water","kg kg-1 s-1",L2,R),
 "prain": ("precipitation_production_rate_wrt_moist_air_and_condensed_water","kg kg-1 s-1",L2,R),
 "lflx": ("cloud_liquid_sedimentation_flux","kg m-2 s-1",LI,R),
 "iflx": ("cloud_ice_sedimentation_flux","kg m-2 s-1",LI,R),
 "gflx": ("graupel_sedimentation_flux","kg m-2 s-1",LI,R),
 "rflx": ("rain_sedimentation_flux","kg m-2 s-1",LI,R),
 "sflx": ("snow_sedimentation_flux","kg m-2 s-1",LI,R),
 "qrout": ("rain"+q+"_of_new_state_in_precipitating_fraction_of_gridcell","kg kg-1",L2,R),
 "qsout": ("snow"+q+"_of_new_state_in_precipitating_fraction_of_gridcell","kg kg-1",L2,R),
 "qgout": ("graupel_water"+q+"_of_new_state_in_precipitating_fraction_of_gridcell","kg kg-1",L2,R),
 "nrout": ("rain_number_concentration_of_new_state_in_precipitating_fraction_of_gridcell","m-3",L2,R),
 "nsout": ("snow_number_concentration_of_new_state_in_precipitating_fraction_of_gridcell","m-3",L2,R),
 "refl": ("analytic_radar_reflectivity_at_94_GHz_in_precipitating_fraction_of_gridcell","dBZ",L2,R),
 "arefl": ("analytic_radar_reflectivity_at_94_GHz","dBZ",L2,R),
 "areflz": ("analytic_radar_reflectivity_z_factor_at_94_GHz","mm6 m-3",L2,R),
 "frefl": ("fraction_of_gridcell_with_nonzero_radar_reflectivity","fraction",L2,R),
 "csrfl": ("analytic_radar_reflectivity_at_94_GHz_with_cloudsat_thresholds_in_precipitating_fraction_of_gridcell","dBZ",L2,R),
 "acsrfl": ("analytic_radar_reflectivity_at_94_GHz_with_cloudsat_thresholds","dBZ",L2,R),
 "fcsrfl": ("fraction_of_gridcell_with_nonzero_radar_reflectivity_with_cloudsat_thresholds","fraction",L2,R),
 "refl10cm": ("analytic_radar_reflectivity_at_10_cm_wavelength","dBZ",L2,R),
 "reflz10cm": ("analytic_radar_reflectivity_z_factor_at_10_cm_wavelength","mm6 m-3",L2,R),
 "rercld": ("effective_radius_of_stratiform_cloud_liquid_plus_rain_drops","um",L2,R),
 "ncai": ("available_ice_nuclei_number_concentration_of_new_state","m-3",L2,R),
 "ncal": ("available_cloud_condensation_nuclei_number_concentration_of_new_state","m-3",L2,R),
 "qrout2": ("rain"+q+"_of_new_state","kg kg-1",L2,R),
 "qsout2": ("snow"+q+"_of_new_state","kg kg-1",L2,R),
 "nrout2": ("rain_number_concentration_of_new_state","m-3",L2,R),
 "nsout2": ("snow_number_concentration_of_new_state","m-3",L2,R),
 "qgout2": ("graupel_water"+q+"_of_new_state","kg kg-1",L2,R),
 "ngout2": ("graupel_number_concentration_of_new_state","m-3",L2,R),
 "freqg": ("fraction_of_gridcell_with_graupel","fraction",L2,R),
 "freqs": ("fraction_of_gridcell_with_snow","fraction",L2,R),
 "freqr": ("fraction_of_gridcell_with_rain","fraction",L2,R),
 "nfice": ("fraction_of_frozen_water_to_total_condensed_water","fraction",L2,R),
 "qcrat": ("fraction_of_cloud_liquid_tendency_applied_to_state","fraction",L2,R),
 "acprecl": ("accumulated_stratiform_precipitation_rate_for_precipitation_efficiency","m s-1",L1,R),
 "acgcme": ("accumulated_cloud_liquid_condensation_rate_for_precipitation_efficiency","m s-1",L1,R),
 "acnum": ("number_of_timesteps_accumulated_for_precipitation_efficiency","count",L1,R),
}
EXTRA = {"rercld": "  # FOR_REVIEWER_FEEDBACK: the values are in meters (calc_rercld in micro_pumas_v1.F90), but\n"
                   "  #    micro_pumas_ccpp.meta and micro_pumas_ccpp_dimensions_post.meta declare 'um'.  The units\n"
                   "  #    here match so that the CCPP framework does not apply a unit conversion.  Should the\n"
                   "  #    units be fixed to 'm' in the PUMAS and dimensions_post metadata (and then here)?\n"}
def vars_(names):
    return [v(x, *V[x][:3], typ=V[x][3], extra=EXTRA.get(x)) for x in names] + err()
pre  = ("ncol pver top_lev gravit ast cldliq cldice pdel "
        "micro_mg_adjust_cpt temp zm tropLev_cpp tropT_cpp tropZ_cpp").split()
mini = ["pver", "top_lev", "micro_mg_warm_rain", "micro_mg_adjust_cpt"]
mrun = ("ncol pver top_lev dtime gravit rair tmelt pi rhoh2o cpairv micro_mg_warm_rain "
        "pmid pdel temp q cldliq cldice numliq numice qr nr qs ns qg ng ast cld concld cmeliq prec_dp prec_sh "
        "naai naai_hom calday lat lon eccen obliqr lambm0 mvelpp micro_mg_adjust_cpt tropLev_cpp tropT_cpp "
        "proc_rates tlat qvlat qcten qiten ncten niten effi rel des sadice sadsnow nevapr prain prect "
        "lflx iflx gflx rflx sflx qrout qsout qgout nrout nsout refl arefl areflz frefl csrfl acsrfl fcsrfl "
        "refl10cm reflz10cm rercld ncai ncal qrout2 qsout2 nrout2 nsout2 qgout2 ngout2 "
        "freqg freqs freqr nfice qcrat acprecl acgcme acnum").split()
out = []
out += [hdr("micro_pumas_pre_diagnostics"), tbl("micro_pumas_pre_diagnostics_run")] + vars_(pre)
out += ["\n" + hdr("micro_pumas_diagnostics"), tbl("micro_pumas_diagnostics_init")] + vars_(mini)
out += ["\n" + tbl("micro_pumas_diagnostics_run")] + vars_(mrun)
if len(sys.argv) > 1:
    outfile = sys.argv[1]
else:
    outfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "micro_pumas_diagnostics.meta")
with open(outfile, "w") as fh:
    fh.write("".join(out))
