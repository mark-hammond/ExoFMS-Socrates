!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Modified by MDH to interface with SOCRATES (UK Met Office) radiative transfer routine


module atmosphere_mod

!-----------------------------------------------------------------------
!
!    interface for FV dynamical core with Held-Suarez forcing
!
!-----------------------------------------------------------------------


use time_manager_mod, only: time_type, get_time, set_time, operator(+)

use fms_mod,          only: file_exist, open_namelist_file,   &
                            error_mesg, FATAL,                &
                            check_nml_error, stdlog,          &
                            write_version_number,             &
                            close_file, set_domain

  use hs_forcing_mod,   only: hs_forcing_init
  use constants_mod, only: omega, cp_air, rdgas, kappa, radius, grav, rvgas, &
                           cp_vapor, cp_water, latent_heat, rho_cp, rho0, r_vir, cp_vir

!------------------
! FV specific codes:
!------------------
!  use fv_pack
!use            fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, endlon, &
!                              u, v, pt, q, ua, va, delp, phis, ps, pe, peln,    &
!                              pk, pkz, omga, u_srf, v_srf, rlonb, rlatb,        &
!                              cold_start, ncnst, pnats, consv_te, ptop,         &
!                              fv_init, fv_domain, fv_end, change_time, map_dt,  &
!                              adiabatic
use fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, endlon, &
     rlonb, rlatb, rlat, rlon,       &
     cold_start, ncnst, pnats, consv_te, ptop,         &
     fv_init, fv_domain, fv_end, change_time, map_dt,  &
     adiabatic, restart_format, &
     get_eta_level, p_var, ak, bk, ng_d, master, &
     nt_phys
use update_fv_phys_mod, only: update_fv_phys
use fv_arrays_mod, only: fv_array_sync

  use fv_diagnostics, only: fv_diag_init, fv_diag, fv_time
  use timingModule,   only: timing_on, timing_off
  use fv_restart_mod, only: fv_restart, write_fv_rst
  use fv_dynamics_mod, only: fv_dynamics
  use fv_phys_mod, only: fv_phys

  use  diag_manager_mod, only: register_diag_field, send_data
!  use  radiation_mod, only: radiation_init, radiation_down, radiation_up, &
!                                    radiation_end
  use  surface_flux_mod, only: surface_flux
  use dry_convection_mod, only: dry_convection
  use  qe_moist_convection_mod, only: moist_convection, compute_k, d622
  use  simple_sat_vapor_pres_mod, only: escomp

!-----------------------------------------------------------------------

!-------
!SOCRATES specific
!-------
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
private

public   atmosphere_init, atmosphere,  atmosphere_end

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 13.0.2.2 2006/05/19 16:44:39 wfc Exp $'
character(len=128) :: tag = '$Name:  $'
character(len=10), parameter :: mod_name='atmosphere'

!-----------------------------------------------------------------------
!---- private data ----

type        (time_type) :: Time_step_atmos
real                    :: dt_atmos
integer :: sec
integer days, seconds

!-----------------------------------------------------------------------
!df 1401010
integer :: i, j, k
integer :: ij, nx, tsiz, isiz
integer :: is, ie
!real :: eff_heat_capacity = rho_cp*1.
real :: mld = 1.
integer :: id_t_surf, id_conv_rain, id_cond_rain, id_pme
integer :: id_conv_rain_profile, id_cond_rain_profile
integer :: id_flux_t, id_flux_q, id_flux_r, id_rh
logical :: used
real, allocatable, dimension(:,:)   :: t_surf, fms_stellar_flux, lats, lons
real, allocatable, dimension(:,:,:) :: p_half, p_full, rh
real, allocatable, dimension(:,:)   :: flux_t, flux_q, flux_r
real, allocatable, dimension(:,:)   :: conv_rain, cond_rain, pme
real, allocatable, dimension(:,:)   :: net_surf_sw_down, surf_lw_down
real, allocatable, dimension(:,:,:) :: conv_rain_profile, cond_rain_profile

!-------------------------------
!SOCRATES specific
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


!  REAL  (RealK), ALLOCATABLE :: flux_net(:,:,:)
!       Net flux
!  REAL  (RealK), ALLOCATABLE :: heating_rate(:,:,:)
!       Heating rates

! Controlling variables:
!  INTEGER :: i, j, l, ic, ll

!----------------------------------

!-----------------------------------------------------------------------
contains

!#######################################################################

 subroutine atmosphere_init ( Time_init, Time, Time_step )

#include <fv_arrays.h>

 type (time_type), intent(in) :: Time_step
 type (time_type), intent(in) :: Time_init
 type (time_type), intent(in) :: Time
 logical cold_start

! local:
 integer axes(4)
 integer ss, ds

#include <fv_point.inc>
!----- write version and namelist to log file -----

    call write_version_number ( version, tag )

!---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    call get_time (Time_step_atmos, sec)
    dt_atmos = real(sec)

!----- initialize FV dynamical core -----

    call fv_init( sec )
    call fv_restart( days, seconds )
   
    if ( .not. cold_start ) then 
!
! Check consistency in Time
!
    fv_time = set_time (seconds, days)
    call get_time (Time, ss,  ds)

    if( seconds /= ss .or. days /= ds ) then
!       write(6,*) 'FMS:', ds, ss
!       write(6,*) 'FV:', days, seconds
        call error_mesg('FV_init:', &
     'Time inconsistent between fv_rst and INPUT/atmos_model.res', &
                         FATAL)
    endif
    else
    fv_time = time
    endif

    call fv_diag_init( axes, Time )

    if( nlev > 1 ) call hs_forcing_init ( axes, Time )

!-----------------------------------------------------------------------
!df 1401010
allocate(t_surf      (1:nlon, beglat:endlat))
allocate(fms_stellar_flux      (1:nlon, beglat:endlat))
allocate(lats      (1:nlon, beglat:endlat))
allocate(lons      (1:nlon, beglat:endlat))
!t_surf(:,:) = tsini 
!----------------------------------------------------------------------
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
if (.not. file_exist('INPUT/atmos_tracers.res.nc')) then
    do ij=1,tsiz
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
       t_surf(is:ie,j) = pt(is:ie,j,nlev) + 1.
    end do
else
    do ij=1,tsiz
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
       t_surf(is:ie,j) = q(is:ie,j,nlev,2) !-20.
       if (master) write (6,*) "restart Ts"
    end do 
end if

    q(:,:,nlev,2) = 0. !1e-4

allocate(p_half      (1:nlon, beglat:endlat, nlev+1))
allocate(p_full      (1:nlon, beglat:endlat, nlev))
allocate(flux_t      (1:nlon, beglat:endlat))
allocate(flux_q      (1:nlon, beglat:endlat))
allocate(flux_r      (1:nlon, beglat:endlat))
allocate(conv_rain   (1:nlon, beglat:endlat))
allocate(cond_rain   (1:nlon, beglat:endlat))
allocate(pme         (1:nlon, beglat:endlat))
allocate(net_surf_sw_down   (1:nlon, beglat:endlat))
allocate(surf_lw_down       (1:nlon, beglat:endlat))
allocate(conv_rain_profile  (1:nlon, beglat:endlat, nlev))
allocate(cond_rain_profile  (1:nlon, beglat:endlat, nlev))
allocate(rh                 (1:nlon, beglat:endlat, nlev))

p_half = 0.; p_full =0.
!     ----- register diagnostic fields -----
id_t_surf = register_diag_field(mod_name, 't_surf',        &
     axes(1:2), Time, 'Surface temperature','Kelvin')

id_flux_t = register_diag_field(mod_name, 'flux_t',        &
     axes(1:2), Time, 'Sensible heat flux','W/m**2')

id_flux_q = register_diag_field(mod_name, 'flux_q',        &
     axes(1:2), Time, 'Evaporation flux','W/m**2')

id_flux_r = register_diag_field(mod_name, 'flux_r',        &
     axes(1:2), Time, 'Upward IR at surface','W/m**2')

id_conv_rain = register_diag_field(mod_name, 'conv_rain',        &
     axes(1:2), Time, 'Convection rain','kg/s/m**2')

id_cond_rain = register_diag_field(mod_name, 'cond_rain',        &
     axes(1:2), Time, 'Condensation rain','kg/s/m**2')

id_pme = register_diag_field(mod_name, 'pme',        &
     axes(1:2), Time, 'P-E','kg/s/m**2')

id_conv_rain_profile = register_diag_field(mod_name, 'conv_profile',  &
     axes(1:3), Time, 'Convection rain profile','kg/s/m**2')

id_cond_rain_profile = register_diag_field(mod_name, 'cond_profile',  &
     axes(1:3), Time, 'Condensation rain profile','kg/s/m**2')

id_rh = register_diag_field(mod_name, 'rh',        &
     axes(1:3), Time, 'Relative humidity','%')

!-----------------------------------------------------------------------
!call radiation_init(1, nlon, beglat, endlat, nlev, axes, Time,rlat(:,:))

!-----------------------------------------------------------------------

!control%spectral_file = '/network/group/aopp/testvol2/plan/fms-scratch-mdh/spec_file_co2_co'
!control_sw%spectral_file = '/network/group/aopp/testvol2/plan/fms-scratch-mdh/spec_file_co2_co'

control%spectral_file = '~/spec_file_co2_co_lowres'
control_sw%spectral_file = '~/spec_file_co2_co_lowres'

!control%spectral_file = '/network/group/aopp/testvol2/plan/fms-scratch-mdh/sp_lw_300_jm2'

!control%spectral_file = '/network/group/aopp/testvol2/plan/fms-scratch-mdh/sp_lw_ga7'

!control%spectral_file = '/network/group/aopp/planetary/RTP001_HAMMOND_55CNCE-C/spec_file_co2_co_lowres'
!control_sw%spectral_file = '/network/group/aopp/planetary/RTP001_HAMMOND_55CNCE-C/spec_file_co2_co_lowres'
CALL read_spectrum(control%spectral_file,spectrum)

!control_sw%spectral_file = '/network/group/aopp/testvol2/plan/fms-scratch-mdh/sp_sw_ga7'
CALL read_spectrum(control_sw%spectral_file,spectrum_sw)

!CALL read_control(control,spectrum)
!CALL compress_spectrum(control,spectrum)

!CALL read_control(control_sw,spectrum_sw)
!CALL compress_spectrum(control_sw,spectrum_sw)

! Read Socrates namelist
!unit = open_file ('input.nml', action='read')
!ierr=1
!do while (ierr /= 0)
!   read  (unit, nml=socrates_nml, iostat=io, end=10)
!   ierr = check_nml_error (io, 'radiation_nml')
!enddo
!10 call close_file (unit)

 end subroutine atmosphere_init


!#######################################################################

 subroutine atmosphere (Time)
#include <fv_arrays.h>
 type(time_type), intent(in) :: Time

! local:
 real    zvir
 logical p_map                      ! Perform only partial remapping
#include <fv_point.inc>

!df 141011 used in my physics schemes---------------------------------
real, dimension(1:nlon, beglat:endlat, nlev) :: tg_tmp, qg_tmp, cp, &
                                        u_tmp, v_tmp, rh_tmp, esat, &
                                        dt_ug, dt_vg, dt_tg, dt_qg, &
                                        diff_u, diff_v, atm_mass, atm_rho
real, dimension(1:nlon, beglat:endlat) :: q1, p1, dp1, dmass, lmass, rl, rt, ct, kappa1, t11, &
                                 zh, Ee, sigma, psg_tmp, Ek, dt_psg

real, dimension(1:nlon, beglat:endlat) ::      &
     gust,                 &   !gustiness constant
     flux_u,               &   ! surface flux of zonal mom.
     flux_v,               &   ! surface flux of meridional mom.
     drag_m,               &   ! momentum drag coefficient
     drag_t,               &   ! heat drag coefficient
     drag_q,               &   ! moisture drag coefficient
     w_atm,                &   ! wind speed
     ustar,                &   ! friction velocity
     bstar,                &   ! buoyancy scale
     qstar,                &   ! moisture scale
     dhdt_surf,            &   ! d(sensible heat flux)/d(surface temp)
     dedt_surf,            &   ! d(latent heat flux)/d(surface temp)???
     dedq_surf,            &   ! d(latent heat flux)/d(surface moisture)???
     drdt_surf,            &   ! d(upward longwave)/d(surface temp)
     dhdt_atm,             &   ! d(sensible heat flux)/d(atmos.temp)
     dedq_atm,             &   ! d(latent heat flux)/d(atmospheric mixing rat.)
     dtaudv_atm                ! d(stress component)/d(atmos wind)
logical, dimension(1:nlon, beglat:endlat) ::    land, avail

real, dimension(nlev) :: t_prev, q_prev, pfull_prev, t_after, q_after, &
                               pfull_after, rain_profile, lnpfull_prev, &
                               lnp_full, water, u_after, v_after, &
                               tmean, qmean
                       
real, dimension(nlev+1) :: phalf_prev, phalf_after, lnp_half

real, dimension(1:nlon, nlev) :: tm, qm, um, vm
real, dimension(1:nlon)       :: psm, tsm
real :: psmean, tsmean

integer :: k,n  
integer :: ngroup, nn, l
!integer :: i, j, k, n
!integer ij, nx, tsiz, isiz
!integer is, ie

!SOCRATES
LOGICAL :: input_l_planet_grey_surface = .FALSE.
REAL(r_def) :: input_planet_albedo = 0.06
REAL(r_def) :: input_planet_emissivity = 0.94
INTEGER(i_def), PARAMETER :: n_profile = 432
INTEGER(i_def), PARAMETER :: n_layer = 40
INTEGER(i_def) :: input_n_cloud_layer = 40
INTEGER(i_def) :: input_n_aer_mode = 40
INTEGER(i_def) :: input_cld_subcol_gen = 40
INTEGER(i_def) :: input_cld_subcol_req = 40
!INTEGER, PARAMETER :: n_profile = 432
!INTEGER, PARAMETER :: n_layer = 40

! Fluxes and radiances calculated
!  REAL  (RealK), ALLOCATABLE :: flux_diffuse_diown(:,:,:)
!       Diffuse downward flux
!  REAL  (RealK), ALLOCATABLE :: flux_net(:,:,:)
!       Net flux
!  REAL  (RealK), ALLOCATABLE :: heating_rate(:,:,:)
!       Heating rates

real, dimension(n_profile,n_layer) :: input_p, input_t, input_mixing_ratio, &
                                      input_d_mass, input_density, input_layer_heat_capacity, &
                                      soc_heating_rate, output_heating_rate, input_o3_mixing_ratio
real, dimension(n_profile,0:n_layer) :: input_p_level, input_t_level, soc_flux_direct, &
                                        soc_flux_down, soc_flux_up, output_flux_net
real, dimension(n_profile) :: input_t_surf, input_cos_zenith_angle, input_solar_irrad, &
                              input_orog_corr


real :: Ep, delta_t_surf, Emass
real :: lh, rain1
real :: k3, k4, Ep_2, water_1, water_2
real :: vcoeff, delta_t
real :: ftop
real gmean

land = .false.; avail = .true.; gust = 1.
dt_ug = 0.; dt_vg = 0.; dt_tg = 0.; dt_qg = 0.; dt_psg = 0.
diff_u = 0.; diff_v = 0.
delta_t = dt_atmos
!-----------------------------------------------------------------------

  fv_time = Time + Time_step_atmos
  call get_time (fv_time, seconds,  days)

!---- call fv dynamics -----
  zvir = 0.         ! no virtual effect if not full physics
!  zvir = d608

  if ( mod(seconds, map_dt) == 0 ) then
       p_map = .false.
  else
       p_map = .true.
  endif
                                call timing_on('fv_dynamics')

#ifndef USE_LIMA
  call fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
                    ncnst,   pnats,  p_map,   consv_te,            &
                    u,       v,      delp,    pt,       q,         &
                    ps,      pe,     pk,      pkz,      phis,      &
                    omga,    peln,   ptop,    omega,    sec,       &  
                    zvir,    cp_air, rdgas,   kappa,  radius, ua, va, fv_time )
#else
  call fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
                    ncnst,   pnats,  p_map,   consv_te,            &
                    u,       v,      delp,    pt,       q,         &
                    ps,      pe,     pk,      pkz,      phis,      &
                    omga,    peln,   ptop,    omega,    sec,       &  
                    zvir,    cp_air, rdgas,   kappa,  radius, ua, va )
#endif

                                call timing_off('fv_dynamics')

      if( nlev /=1 .and. .not. adiabatic ) then
                                call timing_on('FV_PHYS')
!                call fv_phys ( fv_time, sec )
!df physics 141011
!----------------------------------------------------------------------
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    do ij=1,tsiz              
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
!       do k=1,nlev
       do i=is,ie
          call get_eta_level(nlev, ps(i,j), p_full(i,j,:), p_half(i,j,:))
          u_tmp(i,j,:) = ua(i,j,:)
          v_tmp(i,j,:) = va(i,j,:)
          tg_tmp(i,j,:) = pt(i,j,:)
          qg_tmp(i,j,:) = q(i,j,:,1)
          psg_tmp(i,j) = ps(i,j)
       end do
    end do
       n = nlev
       cp = cp_air * (1 - qg_tmp) + cp_vapor * qg_tmp


!-----------
! Socrates interface
! Retrieve output_heating_rate, and downward surface SW and LW fluxes
CALL socrates_interface(Time, spectrum_lw, spectrum_sw,     &
     tg_tmp, t_surf, p_full, p_half, n_profile, n_layer     &
     fms_heating_rate, net_surf_sw_down, surf_lw_down )

tg_tmp

tg_tmp = tg_tmp + fms_heating_rate* delta_t
!-----------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SOCRATES
!Set to LW
control%isolir = 2

!Read in control values

!control%spectral_file = '~/spec_file_co2_co_lowres'
!control_sw%spectral_file = '~/spec_file_co2_co_lowres'

!CALL read_spectrum(control%spectral_file,spectrum)
CALL read_control(control,spectrum)
!CALL compress_spectrum(control,spectrum)

!Set input T, p, p_level, and mixing ratio profiles
input_t = RESHAPE(tg_tmp, (/432, 40/))
input_p = RESHAPE(p_full, (/432, 40/))
input_p_level = RESHAPE(p_half, (/432, 41/))
input_mixing_ratio = 1.E-1
input_o3_mixing_ratio = 1.E-1

!PRINT*, input_p(1,:)
!Optional extra layers for radiative balance
!control%l_extra_top = .TRUE.
!control_sw%l_extra_top = .TRUE.

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

!Default parameters
input_cos_zenith_angle = 0.7 
!input_solar_irrad = 1370.0
input_orog_corr = 0.0
input_layer_heat_capacity = 29.07

!Set tide-locked flux - should be set by namelist eventually!
fms_stellar_flux = soc_stellar_constant*cos(rlat)*cos(rlon)
WHERE (fms_stellar_flux < 0.0) fms_stellar_flux = 0.0
input_solar_irrad = RESHAPE(fms_stellar_flux, (/432/))

!KLUDGE
!!where (input_solar_irrad(:) < 1.e5)
!!   input_solar_irrad(:) = 1.e5
!!endwhere


!Set input surface T
input_t_surf = RESHAPE(t_surf, (/432/))


!Start of a rough T-dependent albedo
!bound%rho_alb(:,:,:) = 0.07
!WHERE (t_surf < 273.0) bound%rho_alb(:,:,:) = 0.7

!Set input dry mass, density, and heat capacity profiles
DO i=n_layer, 1, -1
      DO l=1, n_profile
        input_d_mass(l, i) = (input_p_level(l, i)-input_p_level(l, i-1))/23.0
        input_density(l, i) = input_p(l, i)/(8.31*input_t(l, i))!1000.!atm%p(l ,i) / 1000
!KLUDGE
        input_layer_heat_capacity(l,i) = input_d_mass(l,i)*1005.0
      END DO
END DO


CALL socrates_calc(Time, control, spectrum,                                          &
  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
  input_cld_subcol_gen, input_cld_subcol_req,                                  &
  input_p, input_t, input_t_level, input_d_mass, input_density,                &
  input_mixing_ratio, input_o3_mixing_ratio,                                      &
  input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
  input_layer_heat_capacity,                                                   &
  soc_flux_direct, soc_flux_down, soc_flux_up, soc_heating_rate)


!output_flux_net(:,1) =  soc_flux_down(:,40)!soc_flux_up(1,:) - soc_flux_down(1,:)

!A
surf_lw_down = RESHAPE(soc_flux_down(:,40) , (/144,3/))

!PRINT*, 'ook'
!PRINT*, soc_flux_up(1,:)


!output_flux_net(1,:) = soc_flux_up(1,:) - soc_flux_down(1,:)

!A-PROBLEM HERE!
!soc_heating_rate(:,:10)=0.0
output_heating_rate =soc_heating_rate(:,:)

!PRINT*, 'ook'
!PRINT*, soc_heating_rate(1,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set SW
control_sw%isolir = 1
!CALL read_spectrum(control_sw%spectral_file,spectrum_sw)
CALL read_control(control_sw, spectrum_sw)
!CALL compress_spectrum(control_sw,spectrum_sw)

CALL socrates_calc(Time, control_sw, spectrum_sw,                                          &
  n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,                   &
  input_cld_subcol_gen, input_cld_subcol_req,                                  &
  input_p, input_t, input_t_level, input_d_mass, input_density,                &
  input_mixing_ratio, input_o3_mixing_ratio,                                      &
  input_t_surf, input_cos_zenith_angle, input_solar_irrad, input_orog_corr,    &
  input_l_planet_grey_surface, input_planet_albedo, input_planet_emissivity,   &
  input_layer_heat_capacity,                                                   &
  soc_flux_direct, soc_flux_down, soc_flux_up, soc_heating_rate)

!output_flux_net(:,1) = output_flux_net(:,1) + soc_flux_down(:,40)!soc_flux_up(1,:) - soc_flux_down(1,:)

!A
net_surf_sw_down = RESHAPE(soc_flux_down(:,0) , (/144,3/))
!net_surf_sw_down = fms_stellar_flux(:,:)

!A - 
output_heating_rate = output_heating_rate + soc_heating_rate*0.0
output_heating_rate(:,:6) = 0.0


!PRINT*, 'ook'
!PRINT*, output_heating_rate(:144,30)

!PRINT*, 'ook'
!PRINT*, soc_flux_down(1,:)

!PRINT*, 'ook'
!PRINT*, input_t(1,:)

!PRINT*, 'ook'
!PRINT*, input_t_surf(1)
 
tg_tmp = tg_tmp + RESHAPE(output_heating_rate, (/144, 3, 40/)) * delta_t*0.01



omga=omga*0.0
u=u*0.0
v=v*0.0

!PRINT*, 'ook'
!PRINT*, tg_tmp(:,1,20)

!KLUDGE - FIX!
!tg_tmp(:,:,1) = tg_tmp(:,:,2)
!WHERE (tg_tmp < 1300.0) tg_tmp = 1300.0
!tg_tmp(:,:,:5) = 200.0

!------------------------------------------------------

 
       call surface_flux(       &
        tg_tmp(:,:,nlev),       &
        qg_tmp(:,:,nlev),       &
        u_tmp(:,:,nlev),            &
        v_tmp(:,:,nlev),            &
        p_full(:,:,nlev),       &
        dt_psg(:,:),            & !zatm                          
        p_half(:,:,nlev+1),     &
        t_surf(:,:),            &
        t_surf(:,:),            &
        dt_psg(:,:),            & !q_surf
        dt_psg(:,:),            & !u_surf
        dt_psg(:,:),            & !v_surf
        dt_psg(:,:),            & !rough_mom
        dt_psg(:,:),            & !rough_heat
        dt_psg(:,:),            & !rough_moist
        gust(:,:),              &
        flux_t(:,:),            &
        flux_q(:,:),            &
        flux_r(:,:),            &
                                  flux_u(:,:),                              &
                                  flux_v(:,:),                              &
                                  drag_m(:,:),                              &
                                  drag_t(:,:),                              &
                                  drag_q(:,:),                              &
                                   w_atm(:,:),                              &
                                   ustar(:,:),                              &
                                   bstar(:,:),                              &
                                   qstar(:,:),                              &
                               dhdt_surf(:,:),                              &
                               dedt_surf(:,:),                              &
                               dedq_surf(:,:),                              &
                               drdt_surf(:,:),                              &
                                dhdt_atm(:,:),                              &
                                dedq_atm(:,:),                              &
                              dtaudv_atm(:,:),                              &
                                      delta_t,                              &
                                    land(:,:),                              &
                                   avail(:,:)  )

! boundary layer scheme
tg_tmp(:,:,n) = tg_tmp(:,:,n) + flux_t * delta_t &
           * grav / (p_half(:,:,n+1) - p_half(:,:,n)) &
           / cp(:,:,n)


!dt_ug(:, :, n) = - ug(:,:,n,previous) / 86400.  
!dt_vg(:, :, n) = - vg(:,:,n,previous) / 86400.  !* 2.
!Ek = (p_half(:,:,n+1) - p_half(:,:,n)) / grav/ 86400.  &
!     * (ug(:,:,n,previous)**2 + vg(:,:,n,previous)**2)  !t_surf
!do i = is, ie
!   do j = js, je
!      if (flux_q(i,j) .le. 0) flux_q(i,j) = 0.
!   end do
!end do
!Rayleigh damping
vcoeff = -1./(1. - 0.7) / 86400. 
do k = 1, nlev
   sigma(:,:) = p_full(:,:,k) / psg_tmp(:,:)
   where (sigma(:,:) > 0.7)
      dt_ug(:,:,k) = u_tmp(:,:,k) *vcoeff *(sigma(:,:) - 0.7)
      dt_vg(:,:,k) = v_tmp(:,:,k) *vcoeff *(sigma(:,:) - 0.7)
!   elsewhere (p_full(:,:,k) < 1e2)
!      dt_ug(:,:,k) = -u_tmp(:,:,k) / 86400.
!      dt_vg(:,:,k) = -v_tmp(:,:,k) / 86400.
   elsewhere
      dt_ug(:,:,k) = 0.
      dt_vg(:,:,k) = 0.
   endwhere
end do
      dt_ug(:,:,1:5) = -u_tmp(:,:,1:5) / 86400. /10.
      dt_vg(:,:,1:5) = -v_tmp(:,:,1:5) / 86400. /10.

!!vertical diffusion: constant km, df 14/11/17
!vcoeff =  .5 !m**2/sec
!do k=2,nlev-1
!       diff_u(:,:,k) = (u_tmp(:,:,k+1)-u_tmp(:,:,k-1))/(p_full(:,:,k+1)-p_full(:,:,k-1))* &
!                       (p_full(:,:,k)/rdgas/tg_tmp(:,:,k)/(1.+0.608*qg_tmp(:,:,k)))**2    &
!                       *grav**2 *vcoeff
!       diff_v(:,:,k) = (v_tmp(:,:,k+1)-v_tmp(:,:,k-1))/(p_full(:,:,k+1)-p_full(:,:,k-1))* &
!                       (p_full(:,:,k)/rdgas/tg_tmp(:,:,k)/(1.+0.608*qg_tmp(:,:,k)))**2    &
!                       *grav**2 *vcoeff
!end do
!       k=nlev
!       diff_u(:,:,k) = (0.-u_tmp(:,:,k))/(p_half(:,:,k+1)-p_full(:,:,k))* &
!                       (p_full(:,:,k)/rdgas/tg_tmp(:,:,k)/(1.+0.608*qg_tmp(:,:,k)))**2 &
!                       *grav**2 *vcoeff
!       diff_v(:,:,k) = (0.-v_tmp(:,:,k))/(p_half(:,:,k+1)-p_full(:,:,k))* &
!                       (p_full(:,:,k)/rdgas/tg_tmp(:,:,k)/(1.+0.608*qg_tmp(:,:,k)))**2 &
!                       *grav**2 *vcoeff
!do k=1,nlev-1
!   sigma(:,:) = p_full(:,:,k) / psg_tmp(:,:)
!   where (sigma(:,:) > 0.7)
!!      dt_ug(:,:,k) = u_tmp(:,:,k) *vcoeff *(sigma(:,:) - 0.7)
!!      dt_vg(:,:,k) = v_tmp(:,:,k) *vcoeff *(sigma(:,:) - 0.7)
!       dt_ug(:,:,k)  = (diff_u(:,:,k+1)-diff_u(:,:,k-1))/(p_full(:,:,k+1)-p_full(:,:,k-1))
!       dt_vg(:,:,k)  = (diff_v(:,:,k+1)-diff_v(:,:,k-1))/(p_full(:,:,k+1)-p_full(:,:,k-1))
!    elsewhere (p_full(:,:,k) < 1e2)
!      dt_ug(:,:,k) = -u_tmp(:,:,k) / 86400.
!      dt_vg(:,:,k) = -v_tmp(:,:,k) / 86400.
!   elsewhere
!      dt_ug(:,:,k) = 0.
!      dt_vg(:,:,k) = 0.
!   endwhere
!end do
!       k=nlev
!       dt_ug(:,:,k)  = (diff_u(:,:,k)-diff_u(:,:,k-1))/(p_full(:,:,k)-p_full(:,:,k-1))
!       dt_vg(:,:,k)  = (diff_v(:,:,k)-diff_v(:,:,k-1))/(p_full(:,:,k)-p_full(:,:,k-1))
!


tg_tmp = tg_tmp - delta_t / cp * &
         (dt_ug*(u_tmp(:,:,:)+dt_ug*delta_t/2.)+ &
          dt_vg*(v_tmp(:,:,:)+dt_vg*delta_t/2.) )
u_tmp = u_tmp + dt_ug * delta_t
v_tmp = v_tmp + dt_vg * delta_t
where (flux_q(:,:) < 0.0)
   flux_q(:,:) = 0.
endwhere


!evaporation, convection
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
do ij=1,tsiz
   j  = beglat + (ij-1) / nx
   is = 1 + isiz * mod(ij-1, nx)
   ie = is + isiz - 1

   do i=is,ie
! moist convection------------------------------------------------------
    t_prev = tg_tmp(i,j,:) *1.
    q_prev = qg_tmp(i,j,:) *1.
    phalf_prev = p_half(i,j,:) *1.
    pfull_prev = p_full(i,j,:) *1.
!dry convection
    call dry_convection(t_prev,  pfull_prev, &
            t_after) 

    !large-scale condensation
    !prev is actually after
    ! update surface temperature and pressure

!KLUDGE - flux_t off
    delta_t_surf = (surf_lw_down(i,j) + net_surf_sw_down(i,j)  &
                     - flux_t(i,j) - flux_r(i,j) ) &
                    * delta_t / rho_cp / mld !eff_heat_capacity
!PRINT*, '-----------------------'
!PRINT*, 't_surf'
!PRINT*, t_surf(i,j)
!PRINT*, 'lw'
!    delta_t_surf = (surf_lw_down(i,j) + net_surf_sw_down(i,j))  &
!                    * delta_t / rho_cp / mld !eff_heat_capacity

!    delta_t_surf = (surf_lw_down(i,j) + net_surf_sw_down(i,j)-0.5*5.67e-8*t_surf(i,j)**4) &
!                   * delta_t / 200000.0
    t_surf(i,j) = t_surf(i,j) + delta_t_surf*0.0001
    !correct the energy imbalance due to vertical interpolation

    tg_tmp(i,j,:) = t_after
    !qg_tmp(i,j,:) = q_after
    !u_tmp(i,j,:) = u_after
    !v_tmp(i,j,:) = v_after
!    p_half(i,j,:) = phalf_after
!    p_full(i,j,:) = pfull_after
    call escomp(tg_tmp(i,j,:), esat(i,j,:))
    rh_tmp(i,j,:) = qg_tmp(i,j,:)*pfull_after /&
        (0.622 + 0.378 * qg_tmp(i,j,:)) /esat(i,j,:) * 100
 
   end do

   ps(is:ie,j)    = psg_tmp(is:ie,j)
   pt(is:ie,j,:)  = tg_tmp(is:ie,j,:)
   q(is:ie,j,:,1) = qg_tmp(is:ie,j,:)
   !u(is:ie,j,:)   = u_tmp(is:ie,j,:)
   !v(is:ie,j,:)   = v_tmp(is:ie,j,:)
   rh(is:ie,j,:)  = rh_tmp(is:ie,j,:)
   u_dt(is:ie,j,:) = (u_tmp(is:ie,j,:)-ua(is:ie,j,:))/delta_t
   v_dt(is:ie,j,:) = (v_tmp(is:ie,j,:)-va(is:ie,j,:))/delta_t
   t_dt(is:ie,j,:) = 0.
   q_dt(is:ie,j,:,:) = 0.
end do


      do j=beglat,endlat
!----------------------------------------------------------------------
!high latitude filter
!         if ((j==1) .or. (j==mlat-2)) then ! .or. (j==4) .or. (j==mlat-5)) then
!                 tm(:,:) = sum(pt(:,j:j+2,:), dim=2) / 3.
!                 qm(:,:) = sum( q(:,j:j+2,:, 1), dim=2) / 3.
!                 um(:,:) = sum(ua(:,j:j+2,:), dim=2) / 3.
!                 vm(:,:) = sum(va(:,j:j+2,:), dim=2) / 3.
!                 psm(:)  = sum(ps(:,j:j+2), dim=2) /3. 
!                 tsm(:)  = sum(t_surf(:,j:j+2), dim=2) /3. 
!                 do k=1,3
!                    pt(:,j+k-1,:) = tm(:,:) *1.
!                    q(:,j+k-1,:,1)= qm(:,:)*1.
!                    ua(:,j+k-1,:) = um(:,:) *1.
!                    va(:,j+k-1,:) = vm(:,:) *1.
!                    ps(:,j+k-1)   = psm(:) *1.
!                    t_surf(:,j+k-1)   = tsm(:) *1.
!                 end do
!         end if
!----------------------------------------------------------------------
!high latitude filter
!         if ((j==2) .or. (j==mlat-1) ) then !  .or. (j==5) .or. (j==mlat-4)) then
!         if (((j>1) .and. (j<6)) .or. ((j>mlat-5) .and. (j<mlat))) then
!                    pt(:,j,:) =sum(pt(:,j-1:j+1,:), dim=2) / 3.
!                    q(:,j,:,1)=sum( q(:,j-1:j+1,:, 1), dim=2) / 3.
!                    ua(:,j,:) =sum(ua(:,j-1:j+1,:), dim=2) / 3.
!                    va(:,j,:) =sum(ua(:,j-1:j+1,:), dim=2) / 3. 
!                    ps(:,j)   =sum(ps(:,j-1:j+1), dim=2) /3. 
!                    t_surf(:,j)=sum(t_surf(:,j-1:j+1), dim=2) /3. 
!         end if
!----------------------------------------------------------------------
!high latitude filter
!         if ((j==3) .or. (j==mlat-3) ) then 
!                    pt(:,j,:) =sum(pt(:,j-2:j-1,:), dim=2) / 2.
!                    q(:,j,:,1)=sum( q(:,j-2:j-1,:, 1), dim=2) / 2.
!                    ua(:,j,:) =sum(ua(:,j-2:j-1,:), dim=2) / 2.
!                    va(:,j,:) =sum(ua(:,j-2:j-1,:), dim=2) / 2. 
!                    ps(:,j)   =sum(ps(:,j-2:j-1), dim=2) /2. 
!                    t_surf(:,j)=sum(t_surf(:,j-2:j-1), dim=2) /2. 
!         end if
!         if ((j==4) .or. (j==mlat-2) ) then 
!                    pt(:,j,:) =sum(pt(:,j+1:j+2,:), dim=2) / 2.
!                    q(:,j,:,1)=sum( q(:,j+1:j+2,:, 1), dim=2) / 2.
!                    ua(:,j,:) =sum(ua(:,j+1:j+2,:), dim=2) / 2.
!                    va(:,j,:) =sum(ua(:,j+1:j+2,:), dim=2) / 2. 
!                    ps(:,j)   =sum(ps(:,j+1:j+2), dim=2) /2. 
!                    t_surf(:,j)=sum(t_surf(:,j+1:j+2), dim=2) /2. 
!         end if
!----------------------------------------------------------------------
!          do ngroup=1,3
!          if ((j==ngroup) .or. (j==mlat-ngroup+1)) then
!             nn = 2**(4-ngroup)
!             do i=1,nlon, nn
!                umean(:) = sum(ua(i:i+nn-1,j,:), dim=1) /nn
!                vmean(:) = sum(va(i:i+nn-1,j,:), dim=1) /nn
!                tmean(:) = sum(pt(i:i+nn-1,j,:), dim=1) /nn
!                qmean(:) = sum( q(i:i+nn-1,j,:, 1), dim=1) /nn
!                psmean = sum(ps(i:i+nn-1, j), dim=1) /nn
!                tsmean = sum(t_surf(i:i+nn-1, j), dim=1) /nn
!                do k=1,nn
!                   ua(i+k-1,j,:) = umean(:) *1.
!                   va(i+k-1,j,:) = vmean(:) *1.
!                   pt(i+k-1,j,:) = tmean(:) *1.
!                   q(i+k-1,j,:,1) = qmean(:)*1.
!                   ps(i+k-1,j) = psmean
!                   t_surf(i+k-1,j) = tsmean
!                end do
!             end do
!          end if
!          end do
!----------------------------------------------------------------------
         if (j==1) then
!            psmean   = sum(ps(:,j+1), dim=1)/nlon
            tmean(:) = sum(pt(:,j+1,:), dim=1)/nlon
            qmean(:) = sum(q(:,j+1,:,1),dim=1)/nlon
            do i=1,nlon
 !              ps(i,j)   = psmean
               pt(i,j,:) = tmean(:)*1.
               q(i,j,:,1)= qmean(:)*1.
             end do
          else if (j==mlat) then
 !            psmean   = sum(ps(:,j-1), dim=1)/nlon
             tmean(:) = sum(pt(:,j-1,:), dim=1) / nlon
             qmean(:) = sum(q(:,j-1,:,1),dim=1)/nlon
             do i=1,nlon
!                ps(i,j)   = psmean
                pt(i,j,:) = tmean(:)*1.
                q(i,j,:,1)= qmean(:)*1.
             end do
          end if


         do i=1,nlon
            do k=1,nlev
               delp(i,j,k) = ak(k+1) - ak(k) + ps(i,j) * (bk(k+1) - bk(k))
            enddo
         enddo
      enddo
call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
               pe, peln, pk, pkz, kappa, q, ng_d, ncnst, .false. )
!----------------------------------------------------------------------
    call update_fv_phys ( delta_t, nt_phys, .false., .false., Time )
!--------------------------------------------------------------
!df 141214 change in mixed layer depth
!if(master) write(6,*) 'surface temperature correction (K)=', ftop
!if(master) write(6,*) 'p_map=', p_map
!if(master) write(6,*) r_vir, cp_vir, zvir
!--------------------------------------------------------------
                                call timing_off('FV_PHYS')
      endif

!---- diagnostics for FV dynamics -----

                                call timing_on('FV_DIAG')

  call fv_diag(fv_time, nlon, mlat, nlev, beglat, endlat, ncnst, zvir,   &
               dt_atmos, .true.)

                                call timing_off('FV_DIAG')

!--------------------------------------------------------
!df 141010
if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time)
if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time)
if(id_flux_q > 0) used = send_data(id_flux_q, flux_q, Time)
if(id_flux_r > 0) used = send_data(id_flux_r, flux_r, Time)
if(id_conv_rain > 0) used = send_data(id_conv_rain, conv_rain, Time)
if(id_cond_rain > 0) used = send_data(id_cond_rain, cond_rain, Time)
if(id_pme       > 0) used = send_data(id_pme      , pme      , Time)
if(id_conv_rain_profile > 0) used = send_data(id_conv_rain_profile, conv_rain_profile, Time)
if(id_cond_rain_profile > 0) used = send_data(id_cond_rain_profile, cond_rain_profile, Time)

if(id_rh > 0)     used = send_data(id_rh, rh, Time)

!--------------------------------------------------------
 end subroutine atmosphere


 subroutine atmosphere_end

#include <fv_arrays.h>
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length
! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed
!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    do ij=1,tsiz
       j  = beglat + (ij-1) / nx
       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1
       q(is:ie,j,nlev,2) = t_surf(is:ie,j)
    end do

!----- initialize domains for writing global physics data -----

    call set_domain ( fv_domain )
    call get_time (fv_time, seconds,  days)
    call write_fv_rst( 'RESTART/fv_rst.res', days, seconds, grav, &
         restart_format )


    call fv_end(days, seconds)
   
    deallocate(t_surf)
    deallocate(p_half, p_full)
    deallocate(flux_t, flux_q, flux_r)
    deallocate(conv_rain, cond_rain)
    deallocate(net_surf_sw_down, surf_lw_down)
    deallocate(conv_rain_profile, cond_rain_profile, rh)

!    call radiation_end

 end subroutine atmosphere_end

end module atmosphere_mod
