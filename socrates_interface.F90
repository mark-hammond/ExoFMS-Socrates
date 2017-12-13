module socrates_interface_mod

! Socrates calculation interface module
! Takes FMS time, spectra, temperature, and pressure
! Outputs FMS heating rate, and downwards surface LW and SW
! MDH 07/12/17

!----------
! ExoFMS diagnostics
   use    diag_manager_mod,   only: register_diag_field, send_data

! ExoFMS time
   use    time_manager_mod,   only: time_type, &
                                    operator(+), operator(-), operator(/=)

USE read_control_mod
USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
USE def_spectrum

implicit none

! Input spectra
type (StrSpecData) :: spectrum_lw
type (StrSpecData) :: spectrum_sw

! Control options:
type(StrCtrl) :: control
type(StrCtrl) :: control_sw
type(StrCtrl) :: control_lw

! Diagnostic IDs, name, and missing value
integer :: id_soc_olr, id_soc_olr_spectrum_lw, id_soc_surf_spectrum_sw
integer :: id_soc_heating_sw, id_soc_heating_lw, id_soc_heating_rate
character(len=10), parameter :: soc_mod_name = 'socrates'
real :: missing_value = -999

! Socrates inputs from namelist
real :: stellar_flux = 1370.0
logical :: tidally_locked = .true.
namelist/socrates_nml/ stellar_flux, tide_locked


contains

subroutine socrates_init(is, ie, js, je, num_levels, axes, Time, lat)
      !! Initialises Socrates spectra, arrays, and constants

! Arguments
integer, intent(in), dimension(4) :: axes
      !! NB axes refers to the handles of the axes defined in fv_diagnostics
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real, intent(in) , dimension(:,:)   :: lat
!-------------------------------------------------------------------------------------

! Read in namelist
unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=radiation_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'socrates_nml')
enddo
10 call close_file (unit)

!-----------------------------------------------------------------------

! Socrates spectral files -- should be set by namelist
control_lw%spectral_file = '~/spec_file_co2_co_lowres'
control_sw%spectral_file = '~/spec_file_co2_co_lowres'

! Read in spectral files
CALL read_spectrum(control_lw%spectral_file,spectrum_lw)
CALL read_spectrum(control_sw%spectral_file,spectrum_sw)

! Set Socrates configuration
CALL read_control(control_lw,spectrum_lw)
CALL read_control(control_sw,spectrum_sw)

! Specify LW and SW setups
control_sw%isolir=1
control_lw%isolir=2


! Register diagostic fields
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
subroutine socrates_interface(Time_diag, rlat, rlon,        &
  fms_temp, fms_t_surf, fms_p_full, fms_p_half, n_profile, n_layer,        &
  output_heating_rate, fms_net_surf_sw_down, fms_surf_lw_down, fms_stellar_flux )

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

INTEGER(i_def), intent(in) :: n_profile, n_layer

! Input arrays
real(r_def), intent(in) :: fms_temp(:,:,:)
real(r_def), intent(in) :: fms_p_full(:,:,:)
real(r_def), intent(in) :: fms_p_half(:,:,:)
real(r_def), intent(in) :: fms_t_surf(:,:)
real(r_def), intent(in) :: fms_stellar_flux(:,:)
real(r_def), intent(in) :: rlon(:,:)
real(r_def), intent(in) :: rlat(:,:)

!Output arrays
real(r_def), intent(out) :: fms_net_surf_sw_down(:,:)
real(r_def), intent(out) :: fms_surf_lw_down(:,:)

!Input spectra
!type (StrSpecData) :: spectrum_lw
!type (StrSpecData) :: spectrum_sw

! Arrays to send to Socrates
real, dimension(n_profile,n_layer) :: input_p, input_t, input_mixing_ratio, &
                                      input_d_mass, input_density, input_layer_heat_capacity, &
                                      soc_heating_rate, output_heating_rate, input_o3_mixing_ratio, &
                                      soc_heating_rate_lw, soc_heating_rate_sw
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
real :: soc_stellar_constant
logical :: soc_tide_locked
namelist/socrates_nml/ soc_tide_locked, soc_stellar_constant

! Dimensions:
  TYPE(StrDim) :: dimen

  TYPE(StrAtm) :: atm_input

! Loop variables
integer(i_def) :: i, j, k, l

!DIAG Diagnostic
logical :: used

!----------------------------

! Set array sizes
input_n_cloud_layer = n_layer
input_n_aer_mode = n_layer
input_cld_subcol_gen = n_layer
input_cld_subcol_req = n_layer

!Set input T, p, p_level, and mixing ratio profiles
input_t = RESHAPE(fms_temp, (/432, 40/))
input_p = RESHAPE(fms_p_full, (/432, 40/))
input_p_level = RESHAPE(fms_p_half, (/432, 41/))
input_mixing_ratio = 1.E-1
input_o3_mixing_ratio = 1.E-1



!-------------

!Default parameters
input_cos_zenith_angle = 0.7
input_orog_corr = 0.0
input_layer_heat_capacity = 29.07

!Set tide-locked flux - should be set by namelist eventually!
input_solar_irrad = RESHAPE(fms_stellar_flux, (/432/))
input_t_surf = RESHAPE(fms_t_surf, (/432/))


!--------------

!Set input t_level by scaling t - NEEDS TO CHANGE!
DO i = 1, 3
      DO j = 1, 144
            DO k = 0,n_layer
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
CALL read_control(control_lw, spectrum_lw)

CALL socrates_calc(Time_diag, control_lw, spectrum_lw,                                          &
  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
  input_cld_subcol_gen, input_cld_subcol_req,                                  &
  input_p, input_t, input_t_level, input_d_mass, input_density,                &
  input_mixing_ratio, input_o3_mixing_ratio,                                      &
  input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
  input_layer_heat_capacity,                                                   &
  soc_flux_direct, soc_flux_down, soc_flux_up, soc_heating_rate_lw)

! Set output arrays
fms_surf_lw_down = RESHAPE(soc_flux_down(:,40) , (/144,3/))




 !   used = send_data ( id_soc_heating_lw, RESHAPE(soc_heating_rate, (/144,3,40/)), Time_diag)
 !   used = send_data ( id_soc_heating_sw, RESHAPE(soc_flux_up, (/144,3,40/)), Time_diag)

!--------------



! SW calculation
control_sw%isolir = 1
CALL read_control(control_sw, spectrum_sw)

CALL socrates_calc(Time_diag, control_sw, spectrum_sw,                                          &
  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
  input_cld_subcol_gen, input_cld_subcol_req,                                  &
  input_p, input_t, input_t_level, input_d_mass, input_density,                &
  input_mixing_ratio, input_o3_mixing_ratio,                                      &
 input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
  input_layer_heat_capacity,                                                   &
  soc_flux_direct, soc_flux_down, soc_flux_up, soc_heating_rate_sw)

! Set output arrays
fms_net_surf_sw_down = RESHAPE(soc_flux_down(:,0) , (/144,3/))
output_heating_rate = output_heating_rate_sw + soc_heating_rate_lw





! Send LW diagnosticis
!     used = send_data ( id_soc_olr, RESHAPE(flux_up(:,0), (/144,3/)), Time_diag)
!     used = send_data ( id_soc_olr_spectrum_lw, RESHAPE(radout%flux_up_band(:,0,:), (/144,3,20/)), Time_diag)
!     used = send_data ( id_soc_heating_sw, RESHAPE(soc_heating_rate, (/144,3,40/)), Time_diag)

endif

!--------------

end subroutine socrates_interface
end module socrates_interface_mod
