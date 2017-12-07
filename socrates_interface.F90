module socrates_interface_mod

! Socrates calculation interface module
! MDH 07/12/17

!----------
!DIAG ExoFMS diagnostics
   use    diag_manager_mod,   only: register_diag_field, send_data

! ExoFMS time
   use    time_manager_mod,   only: time_type, &
                                    operator(+), operator(-), operator(/=)
!----------
implicit none

contains

! Set up the call to the Socrates radiation scheme
! -----------------------------------------------------------------------------
!DIAG Added Time
subroutine socrates_interface(Time_diag, spectrum_lw, spectrum_sw,        &
  fms_temp, fms_t_surf, fms_p_full, fms_p_half, n_profile, n_layer        &
  fms_heating_rate, net_surf_sw_down, surf_lw_down )


USE realtype_rd
USE soc_constants_mod
USE read_control_mod
USE socrates_calc_mod
USE compress_spectrum_mod
USE def_spectrum
USE def_dimen,   ONLY: StrDim
USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
USE def_atm,     ONLY: StrAtm,   allocate_atm,       deallocate_atm
USE def_cld,     ONLY: StrCld,   allocate_cld,       deallocate_cld, &
                                 allocate_cld_prsc,  deallocate_cld_prsc, &
                                 allocate_cld_mcica, deallocate_cld_mcica
USE def_aer,     ONLY: StrAer,   allocate_aer,       deallocate_aer, &
                                 allocate_aer_prsc,  deallocate_aer_prsc
USE def_bound,   ONLY: StrBound, allocate_bound,     deallocate_bound
USE def_out,     ONLY: StrOut,                       deallocate_out
!-----------------------------------------------------------------------
implicit none

! Input time
type(time_type), intent(in)         :: Time_diag

INTEGER(i_def), intent(in) :: n_profile
INTEGER(i_def), intent(in) :: n_layer

! Input arrays
real(r_def), intent(in) :: fms_temp(:,:,:)
real(r_def), intent(in) :: fms_p_full(:,:,:)
real(r_def), intent(in) :: fms_p_half(:,:,:)
real(r_def), intent(in) :: fms_t_surf(:,:)

!Output arrays
real(r_def), intent(out) :: fms_heating_rate(:,:,:)
real(r_def), intent(out) :: fms_net_surf_sw_down(:,:)
real(r_def), intent(out) :: fms_surf_lw_down(:,:)

!Input spectra
type (StrSpecData), intent(in) :: spectrum_lw
type (StrSpecData), intent(in) :: spectrum_sw

! Arrays to send to Socrates
real, dimension(n_profile,n_layer) :: input_p, input_t, input_mixing_ratio, &
                                      input_d_mass, input_density, input_layer_heat_capacity, &
                                      soc_heating_rate, output_heating_rate, input_o3_mixing_ratio
real, dimension(n_profile,0:n_layer) :: input_p_level, input_t_level, soc_flux_direct, &
                                        soc_flux_down, soc_flux_up, output_flux_net
real, dimension(n_profile) :: input_t_surf, input_cos_zenith_angle, input_solar_irrad, &
                              input_orog_corr


!SOCRATES options
LOGICAL :: input_l_planet_grey_surface = .FALSE.
REAL(r_def) :: input_planet_albedo = 0.06
REAL(r_def) :: input_planet_emissivity
INTEGER(i_def) :: input_n_cloud_layer
INTEGER(i_def) :: input_n_aer_mode
INTEGER(i_def) :: input_cld_subcol_gen
INTEGER(i_def) :: input_cld_subcol_req

!-------------------------------
! Socrates input options -- to become namelist
real :: soc_stellar_constant = 3500000.0!3000.0
logical :: soc_tide_locked = .TRUE.
!namelist/socrates_nml/ soc_tide_locked, soc_stellar_constant

! Dimensions:
  TYPE(StrDim) :: dimen

! Control options:
  TYPE(StrCtrl) :: control
  TYPE(StrCtrl) :: control_sw

! Spectral information:
  TYPE(StrSpecData) :: spectrum
  TYPE(StrSpecData) :: spectrum_sw

! Atmospheric input:
  TYPE(StrAtm) :: atm_input


!----------------------------

! Set array sizes
input_n_cloud_layer = n_layer
input_n_aer_mode = n_layer
input_cld_subcol_gen = n_layer
input_cld_subcol_req = n_layer

!Set input T, p, p_level, and mixing ratio profiles
input_t = RESHAPE(fms_temp, (/432, 40/))
input_p = RESHAPE(p_full, (/432, 40/))
input_p_level = RESHAPE(p_half, (/432, 41/))
input_mixing_ratio = 1.E-1
input_o3_mixing_ratio = 1.E-1


!-------------
-
!Default parameters
input_cos_zenith_angle = 0.7
!input_solar_irrad = 1370.0
input_orog_corr = 0.0
input_layer_heat_capacity = 29.07

!Set tide-locked flux - should be set by namelist eventually!
fms_stellar_flux = soc_stellar_constant*cos(rlat)*cos(rlon)
WHERE (fms_stellar_flux < 0.0) fms_stellar_flux = 0.0
input_solar_irrad = RESHAPE(fms_stellar_flux, (/432/))


!--------------

!Set input t_level by scaling t - NEEDS TO CHANGE!
DO i = 1, 3
      DO j = 1, 144
            DO k = 0,nlev
                 input_t_level(j + 144*(i-1),k) = 0.5*(input_t(j+144*(i-1),k+1)+input_t(j+144*(i-1),k))
            END DO
            input_t_level(j+144*(i-1),40) = input_t(j+144*(i-1),40) + input_t(j+144*(i-1),40) - input_t_level(j+144*(i-1),39)
            input_t_level(j+144*(i-1),0) = input_t(j+144*(i-1),1) - (input_t_level(j+144*(i-1),1) - input_t(j+144*(i-1),1))
      END DO
END DO

!Set input dry mass, density, and heat capacity profiles
DO i=n_layer, 1, -1
      DO l=1, n_profile
        input_d_mass(l, i) = (input_p_level(l, i)-input_p_level(l, i-1))/23.0
        input_density(l, i) = input_p(l, i)/(8.31*input_t(l, i))!1000.!atm%p(l ,i) / 1000
!KLUDGE
        input_layer_heat_capacity(l,i) = input_d_mass(l,i)*1005.0
      END DO
END DO

!--------------
! LW calculation
control_lw%isolir = 2
CALL read_control(control_lw,spectrum_lw)

CALL socrates_calc(Time, control_lw, spectrum_lw,                                          &
  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
  input_cld_subcol_gen, input_cld_subcol_req,                                  &
  input_p, input_t, input_t_level, input_d_mass, input_density,                &
  input_mixing_ratio, input_o3_mixing_ratio,                                      &
  input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
  input_layer_heat_capacity,                                                   &
  soc_flux_direct, soc_flux_down, soc_flux_up, soc_heating_rate)

! Set output arrays
surf_lw_down = RESHAPE(soc_flux_down(:,40) , (/144,3/))
output_heating_rate =soc_heating_rate(:,:)

!--------------

! SW calculation
control_sw%isolir = 1
CALL read_control(control_sw, spectrum_sw

CALL socrates_calc(Time, control_sw, spectrum_sw,                                          &
  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
  input_cld_subcol_gen, input_cld_subcol_req,                                  &
  input_p, input_t, input_t_level, input_d_mass, input_density,                &
  input_mixing_ratio, input_o3_mixing_ratio,                                      &
  input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
  input_layer_heat_capacity,                                                   &
  soc_flux_direct, soc_flux_down, soc_flux_up, soc_heating_rate))

! Set output arrays
net_surf_sw_down = RESHAPE(soc_flux_down(:,0) , (/144,3/))
output_heating_rate = output_heating_rate + soc_heating_rate
fms_heating_rate = RESHAPE(output_heating_rate, (144, 3, 40/))

!--------------

end subroutine socrates_interface
end module socrates_interface_mod
