module socrates_calc_mod

! Socrates calculation interface modules
! MDH added FMS diagnostics

!----------
!DIAG ExoFMS diagnostics
   use    diag_manager_mod,   only: register_diag_field, send_data

! ExoFMS time
   use    time_manager_mod,   only: time_type, &
                                    operator(+), operator(-), operator(/=)
!----------
implicit none

!----------
!DIAG ExoFMS diagnostic fields

integer :: id_soc_olr, id_soc_olr_spectrum_lw, id_soc_surf_spectrum_sw
integer :: id_soc_heating_sw, id_soc_heating_lw, id_soc_heating_rate

character(len=10), parameter :: soc_mod_name = 'socrates'

real :: missing_value = -999
!----------

contains

! Initialise ExoFMS diagnostics

subroutine socrates_init(is, ie, js, je, num_levels, axes, Time, lat)


!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real, intent(in) , dimension(:,:)   :: lat
!-------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------
    id_soc_olr = &
    register_diag_field ( soc_mod_name, 'soc_olr', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'watts/m2', missing_value=missing_value               )

    id_soc_olr_spectrum_lw = &
    register_diag_field ( soc_mod_name, 'soc_olr_spectrum_lw',(/ axes(1:2), axes(5)/) , Time, &
               'socrates substellar LW OLR spectrum', &
               'watts/m2', missing_value=missing_value               )

    id_soc_surf_spectrum_sw = &
    register_diag_field ( soc_mod_name, 'soc_surf_spectrum_sw',(/ axes(1:2), axes(5)/) , Time, &
               'socrates substellar SW surface spectrum', &
               'watts/m2', missing_value=missing_value               )

    id_soc_heating_lw = &
    register_diag_field ( soc_mod_name, 'soc_heating_lw', axes(1:3), Time, &
               'socrates LW heating rate', &
               'J/s', missing_value=missing_value               )

    id_soc_heating_sw = &
    register_diag_field ( soc_mod_name, 'soc_heating_sw', axes(1:3), Time, &
               'socrates SW heating rate', &
               'J/s', missing_value=missing_value               )

    id_soc_heating_rate = &
    register_diag_field ( soc_mod_name, 'soc_heating_rate', axes(1:3), Time, &
               'socrates total heating rate', &
               'J/s', missing_value=missing_value               )
return
end subroutine socrates_init
! ==================================================================================




! Set up the call to the Socrates radiation scheme
! -----------------------------------------------------------------------------
!DIAG Added Time
subroutine socrates_calc(Time_diag,control, spectrum,                                    &
  n_profile, n_layer, n_cloud_layer, n_aer_mode,                               &
  cld_subcol_gen, cld_subcol_req,                                              &
  p_layer, t_layer, t_layer_boundaries, d_mass, density,                       &
  h2o, o3,                                                                     &
  t_rad_surf, cos_zenith_angle, solar_irrad, orog_corr,                        &
  l_planet_grey_surface, planet_albedo, planet_emissivity,                     &
  layer_heat_capacity,                                                         &
  flux_direct, flux_down, flux_up, heating_rate)

use rad_pcf
use def_control,  only: StrCtrl
use def_spectrum, only: StrSpecData
use def_dimen,    only: StrDim
use def_atm,      only: StrAtm,   deallocate_atm
use def_bound,    only: StrBound, deallocate_bound
use def_cld,      only: StrCld,   deallocate_cld, deallocate_cld_prsc
use def_aer,      only: StrAer,   deallocate_aer, deallocate_aer_prsc
use def_out,      only: StrOut,   deallocate_out

use set_control_mod, only: set_control
use set_dimen_mod,   only: set_dimen
use set_atm_mod,     only: set_atm
use set_bound_mod,   only: set_bound
use set_cld_mod,     only: set_cld
use set_aer_mod,     only: set_aer

use soc_constants_mod,   only: i_def, r_def

implicit none

!DIAG ExoFMS diagnostic time
type(time_type), intent(in)         :: Time_diag

! Spectral data:
type (StrSpecData), intent(in) :: spectrum

! Controlling options:
type (StrCtrl), intent(inout) :: control

integer(i_def), intent(in) :: n_profile
!   Number of columns to operate on
integer(i_def), intent(in) :: n_layer
!   Number of layers for radiation
integer(i_def), intent(in) :: n_cloud_layer
!   Number of potentially cloudy layers
integer(i_def), intent(in) :: n_aer_mode
!   Number of aerosol modes
integer(i_def), intent(in) :: cld_subcol_gen
!   Number of sub-grid cloud columns generated
integer(i_def), intent(in) :: cld_subcol_req
!   Number of sub-grid cloud columns required

real(r_def), intent(in) :: p_layer(n_profile, n_layer)
!   Pressure at layer centres
real(r_def), intent(in) :: t_layer(n_profile, n_layer)
!   Temperature at layer centres
real(r_def), intent(in) :: t_layer_boundaries(n_profile, 0:n_layer)
!   Temperature at layer boundaries
real(r_def), intent(in) :: d_mass(n_profile, n_layer)
!   Mass of layer (kg m-2)
real(r_def), intent(in) :: density(n_profile, n_layer)
!   Density of layer (kg m-3)
real(r_def), intent(in) :: h2o(n_profile, n_layer)
!   Mass mixing ratio of water vapour
real(r_def), intent(in) :: o3(n_profile, n_layer)
!   Mass mixing ratio of ozone

real(r_def), intent(in) :: t_rad_surf(n_profile)
!   Effective radiative temperature over whole grid-box
real(r_def), intent(in) :: cos_zenith_angle(n_profile)
!   Cosine of solar zenith angle
real(r_def), intent(in) :: solar_irrad(n_profile)
!   Solar irradiance at top-of-atmosphere (mean over timestep)
real(r_def), intent(in) :: orog_corr(n_profile)
!   Orographic correction factor

logical, intent(in) :: l_planet_grey_surface
!   Set a single grey albedo / emissivity for the surface
real(r_def), intent(in) :: planet_albedo
!   Surface albedo used for SW calculations
real(r_def), intent(in) :: planet_emissivity
!   Surface emissivity used for LW calculations

real(r_def), intent(in) :: layer_heat_capacity(n_profile, n_layer)
!   Heat capacity of layer

real(r_def), intent(out) :: flux_direct(n_profile, 0:n_layer)
!   Direct (unscattered) downwards flux (Wm-2)
real(r_def), intent(out) :: flux_down(n_profile, 0:n_layer)
!   Downwards flux (Wm-2)
real(r_def), intent(out) :: flux_up(n_profile, 0:n_layer)
!   Upwards flux (Wm-2)
real(r_def), intent(out) :: heating_rate(n_profile, n_layer)
!   Heating rate (Ks-1)


! Dimensions:
TYPE (StrDim) :: dimen

! Atmospheric properties:
TYPE(StrAtm) :: atm

! Boundary conditions:
TYPE(StrBound) :: bound

! Cloud properties:
TYPE(StrCld) :: cld

! Aerosol properties:
TYPE(StrAer) :: aer

! Output fields from core radiation code:
TYPE(StrOut) :: radout

integer(i_def) :: l, i
!   Loop variablesi

!DIAG Diagnostic
logical :: used


call set_control(control)

call set_dimen(control, dimen, spectrum, n_profile, n_layer,                   &
  n_cloud_layer, n_aer_mode, cld_subcol_gen, cld_subcol_req)

call set_atm(control, dimen, spectrum, atm, n_profile, n_layer,                &
  p_layer, t_layer, t_layer_boundaries, d_mass, density, h2o, o3)

call set_bound(control, dimen, spectrum, bound, n_profile,                     &
  t_rad_surf, cos_zenith_angle, solar_irrad, orog_corr,                        &
  l_planet_grey_surface, planet_albedo, planet_emissivity)

call set_cld(control, dimen, spectrum, cld, n_profile)

call set_aer(control, dimen, spectrum, aer, n_profile)


! DEPENDS ON: radiance_calc
call radiance_calc(control, dimen, spectrum, atm, cld, aer, bound, radout)


if (control%isolir == 10) then
PRINT*, 'ook'
do i=1,20
PRINT*, radout%flux_down_band(1,39,i)
PRINT*, ','
end do
end if




! set heating rates and diagnostics
do l=1, n_profile
  do i=1, n_layer
    heating_rate(l, i) = (radout%flux_down(l,i-1,1)-radout%flux_down(l,i,1)    &
                        + radout%flux_up(l,i,1)-radout%flux_up(l,i-1,1))       &
                       / layer_heat_capacity(l, i)
  end do
end do


do l=1, n_profile
  do i=0, n_layer
    flux_direct(l, i) = radout%flux_direct(l, i, 1)
    flux_down(l, i)   = radout%flux_down(l, i, 1)
    flux_up(l, i)     = radout%flux_up(l, i, 1)
  end do
end do

!DIAG

! Send SW diagnostics
if ( control%isolir == 1 ) then
     used = send_data ( id_soc_surf_spectrum_sw, RESHAPE(radout%flux_down_band(:,40,:), (/144,3,20/)), Time_diag)
     used = send_data ( id_soc_heating_sw, RESHAPE(heating_rate, (/144,3,40/)), Time_diag)
endif

! Send LW diagnosticis
if ( control%isolir == 2 ) then
     used = send_data ( id_soc_olr, RESHAPE(flux_up(:,0), (/144,3/)), Time_diag)
     used = send_data ( id_soc_olr_spectrum_lw, RESHAPE(radout%flux_up_band(:,0,:), (/144,3,20/)), Time_diag)
     used = send_data ( id_soc_heating_lw, RESHAPE(heating_rate, (/144,3,40/)), Time_diag)
endif

call deallocate_out(radout)
call deallocate_aer_prsc(aer)
call deallocate_aer(aer)
call deallocate_cld_prsc(cld)
call deallocate_cld(cld)
call deallocate_bound(bound)
call deallocate_atm(atm)

end subroutine socrates_calc
end module socrates_calc_mod
