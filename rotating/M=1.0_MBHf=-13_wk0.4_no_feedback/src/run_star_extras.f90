! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras
      
         use star_lib
         use star_def
         use const_def
         use math_lib
         use utils_lib, only: mesa_error
         
         implicit none
         
         ! these routines are called by the standard run_star check_model
         contains
         
         subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            
            ! this is the place to set any procedure pointers you want to change
            ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
            
            
            ! the extras functions in this file will not be called
            ! unless you set their function pointers as done below.
            ! otherwise we use a null_ version which does nothing (except warn).
            
            s% extras_startup => extras_startup
            s% extras_start_step => extras_start_step
            s% extras_check_model => extras_check_model
            s% extras_finish_step => extras_finish_step
            s% extras_after_evolve => extras_after_evolve
            s% how_many_extra_history_columns => how_many_extra_history_columns
            s% data_for_extra_history_columns => data_for_extra_history_columns
            s% how_many_extra_profile_columns => how_many_extra_profile_columns
            s% data_for_extra_profile_columns => data_for_extra_profile_columns  
            
            s% how_many_extra_history_header_items => how_many_extra_history_header_items
            s% data_for_extra_history_header_items => data_for_extra_history_header_items
            s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
            s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
            
         end subroutine extras_controls
         
         subroutine do1_relax_R_center(s, new_Rcenter, ierr)
            ! adjust all lnR's to keep same density for each cell as 1st guess for next model
            type (star_info), pointer :: s
            real(dp), intent(in) :: new_Rcenter ! cm
            integer, intent(out) :: ierr
            real(dp) :: dm, rho, dr3, rp13
            integer :: k
            ierr = 0
            s% R_center = new_Rcenter
            ! adjust lnR's
            rp13 = s% R_center*s% R_center*s% R_center
            do k = s% nz, 1, -1
               dm = s% dm(k)
               rho = s% rho(k)
               dr3 = dm/(rho*four_thirds_pi) ! dm/rho is cell volume
               s% xh(s% i_lnR,k) = log(rp13 + dr3)*one_third
               rp13 = rp13 + dr3
            end do
         end subroutine do1_relax_R_center
       
         subroutine black_hole_accretion(id, s, startup, ierr)
            integer, intent(in) :: id
            logical, intent(in) :: startup
            type (star_info), pointer :: s
            integer, intent(out) :: ierr
            
            real(dp) :: G, c2, c_s, rho, gamma1, opacity, dt
            real(dp) :: nabla_ad, P_rad, P_gas
            real(dp) :: M_BH, M_BH_new, J_BH, J_BH_new
            real(dp) :: M_dot, dm, R_B, L_BH, new_core_mass, core_avg_rho, core_avg_eps
            real(dp) :: j_B, j_B_div_j_ISCO, j_ISCO_ref
            real(dp) :: rad_eff, con_eff, timestep_factor, max_accretion_fraction
            real(dp) :: transition_width
   
            ierr = 0
            transition_width = 0.5d0 ! Set default transition width
            
            
            ! Get control parameters and physical quantities
            call get_control_parameters(s, rad_eff, con_eff, timestep_factor, max_accretion_fraction)
            call get_physical_quantities(s, dt, G, c_s, rho, opacity, P_rad, P_gas, gamma1, nabla_ad)
            
            ! Get black hole properties
            M_BH = s% xtra(1)
            J_BH = s% xtra(15)
            c2 = clight**2
            
            ! Calculate Bondi radius and reference j_ISCO
            R_B = compute_bondi_radius(G, M_BH, c_s)
            j_ISCO_ref = compute_schwarzschild_j_ISCO(G, M_BH)
            
            ! Find angular momentum at Bondi radius
            call find_j_at_bondi_radius(s, R_B, j_B)

            
            ! Determine accretion regime
            j_B_div_j_ISCO = j_B / j_ISCO_ref



            ! Use smooth blending between the two regimes
            call smooth_accretion_transition(s, M_BH, J_BH, j_B, j_B_div_j_ISCO, &
               G, c_s, rho, opacity, gamma1, dt, &
               rad_eff, con_eff, timestep_factor, &
               max_accretion_fraction, transition_width, &
               M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B, ierr)

            if (ierr /= 0) return
            
            ! if (j_B_div_j_ISCO < 1.0_dp) then
            !    ! Bondi accretion without disk
            !    call bondi_accretion_no_disk(s, M_BH, J_BH, j_B, G, c_s, rho, opacity, dt, &
            !                                 timestep_factor, max_accretion_fraction, &
            !                                 M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B, ierr)
            ! else
            !    ! Disk accretion with feedback
            !    call disk_accretion_with_feedback(s, M_BH, J_BH, j_B, G, c_s, rho, opacity, gamma1, dt, &
            !                                      rad_eff, con_eff, timestep_factor, &
            !                                      M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B)
            ! endif
            
            ! if (ierr /= 0) return
            
            ! Update stellar structure
            call update_stellar_structure(s, id, startup, M_BH_new, J_BH_new, L_BH, R_B, M_dot, &
                                          rad_eff, dt, rho, ierr)
            
            ! Store diagnostic quantities
            call store_diagnostics(s, M_BH_new, J_BH_new, L_BH, R_B, M_dot, dm, dt, &
                                  rad_eff, opacity, P_rad, P_gas, nabla_ad, j_B_div_j_ISCO)
            
            ! Print summary
            call print_bh_summary(s, M_BH_new, L_BH, M_dot, rad_eff, j_B_div_j_ISCO)
            
         end subroutine black_hole_accretion
         
         
         subroutine get_control_parameters(s, rad_eff, con_eff, timestep_factor, max_accretion_fraction)
            type (star_info), pointer :: s
            real(dp), intent(out) :: rad_eff, con_eff, timestep_factor, max_accretion_fraction
            
            rad_eff = s% x_ctrl(1)  ! epsilon - radiative efficiency
            con_eff = s% x_ctrl(2)  ! eta - convection efficiency
            timestep_factor = s% x_ctrl(3)
            max_accretion_fraction = s% x_ctrl(5)
            
            ! Set default if not specified
            if (max_accretion_fraction <= 0d0) max_accretion_fraction = 0.01d0
            
         end subroutine get_control_parameters
         
         
         subroutine get_physical_quantities(s, dt, G, c_s, rho, opacity, P_rad, P_gas, gamma1, nabla_ad)
            type (star_info), pointer :: s
            real(dp), intent(out) :: dt, G, c_s, rho, opacity, P_rad, P_gas, gamma1, nabla_ad
            
            dt      = s% dt
            G       = s% cgrav(s% nz)
            c_s     = s% csound(s% nz)
            rho     = s% rho(s% nz)
            opacity = s% opacity(s% nz)
            P_rad   = s% prad(s% nz)
            P_gas   = s% pgas(s% nz)
            gamma1  = s% gamma1(s% nz)
            nabla_ad = 1d0 - 1d0 / gamma1
            
         end subroutine get_physical_quantities
         
         
         function compute_bondi_radius(G, M_BH, c_s) result(R_B)
            real(dp), intent(in) :: G, M_BH, c_s
            real(dp) :: R_B
            
            R_B = 2d0 * G * M_BH / (c_s * c_s)
            
         end function compute_bondi_radius
         
         
         function compute_schwarzschild_j_ISCO(G, M_BH) result(j_ISCO)
            real(dp), intent(in) :: G, M_BH
            real(dp) :: j_ISCO
            
            j_ISCO = 2d0 * sqrt(3d0) * G * M_BH / clight
            
         end function compute_schwarzschild_j_ISCO
         
         
         subroutine find_j_at_bondi_radius(s, R_B, j_B)
            type (star_info), pointer :: s
            real(dp), intent(in) :: R_B
            real(dp), intent(out) :: j_B
            
            real(dp) :: r_center, r_surface
            integer :: k_B, k_lo, k_hi
            
            r_center = s% r(s% nz)
            r_surface = s% r(1)
            
            if (R_B <= r_center) then
               ! Bondi radius inside innermost cell - extrapolate
               k_B = s% nz
               j_B = (2.0_dp/3.0_dp) * s% omega(k_B) * (R_B**2)
               write(*,*) 'Bondi radius inside inner cell; extrapolated j_B from omega(nz).'
               
            else if (R_B >= r_surface) then
               ! Bondi radius outside star - use surface value
               k_B = 1
               j_B = s% j_rot(k_B)
               write(*,*) 'Bondi radius larger than star; using j at surface.'
               
            else
               ! R_B within stellar model - interpolate
               do k_hi = s% nz - 1, 1, -1
                  if ((s% r(k_hi) >= R_B) .and. (R_B >= s% r(k_hi+1))) then
                     k_lo = k_hi + 1
                     exit
                  end if
               end do
               
               ! Linear interpolation
               j_B = s% j_rot(k_hi) + (s% j_rot(k_lo) - s% j_rot(k_hi)) * &
                     (R_B - s% r(k_hi)) / (s% r(k_lo) - s% r(k_hi))
            end if
            
         end subroutine find_j_at_bondi_radius
         
         
         subroutine compute_kerr_ISCO(M_BH, J_BH, j_B, G, R_ISCO, j_ISCO)
            real(dp), intent(in) :: M_BH, J_BH, j_B, G
            real(dp), intent(out) :: R_ISCO, j_ISCO
            
            real(dp) :: a_star, Z1, Z2, r_tilde, sqrt_r, den, L_tilde
            
            ! Dimensionless spin parameter
            a_star = clight * J_BH / (G * M_BH * M_BH)
            a_star = max(-0.999999d0, min(0.999999d0, a_star))
            
            ! Bardeen-Petterson ISCO calculation
            Z1 = 1d0 + (1d0 - a_star*a_star)**(1d0/3d0) * &
                 ((1d0 + a_star)**(1d0/3d0) + (1d0 - a_star)**(1d0/3d0))
            Z2 = sqrt(3d0 * a_star * a_star + Z1 * Z1)
            
            if (j_B >= 0d0) then
               ! Prograde ISCO
               r_tilde = 3d0 + Z2 - sqrt((3d0 - Z1) * (3d0 + Z1 + 2d0*Z2))
            else
               ! Retrograde ISCO
               r_tilde = 3d0 + Z2 + sqrt((3d0 - Z1) * (3d0 + Z1 + 2d0*Z2))
            end if
            r_tilde = max(r_tilde, 1d0)
            
            ! Specific angular momentum at ISCO
            sqrt_r = sqrt(r_tilde)
            den = sqrt(max(1d-30, r_tilde*r_tilde - 3d0*r_tilde + 2d0*a_star*sign(1d0, j_B)*sqrt_r))
            L_tilde = (r_tilde*r_tilde - 2d0*a_star*sign(1d0, j_B)*sqrt_r + a_star*a_star) / den
            
            R_ISCO = r_tilde * G * M_BH / (clight * clight)
            j_ISCO = L_tilde * G * M_BH / clight
            
         end subroutine compute_kerr_ISCO
         
         
         subroutine bondi_accretion_no_disk(s, M_BH, J_BH, j_B, G, c_s, rho, opacity, dt, &
                                            timestep_factor, max_accretion_fraction, &
                                            M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B, ierr)
            type (star_info), pointer :: s
            real(dp), intent(in) :: M_BH, J_BH, j_B, G, c_s, rho, opacity, dt
            real(dp), intent(in) :: timestep_factor, max_accretion_fraction
            real(dp), intent(out) :: M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B
            integer, intent(out) :: ierr
            
            real(dp) :: M_dot_BH, available_mass, max_dm, j_acc, dj
            
            ierr = 0
            !write(*,*) 'Bondi Accretion without a disk'
            
            ! Bondi accretion rate
            M_dot_BH = 4d0 * pi / sqrt(2d0) * rho * (G * M_BH)**2 / (c_s**3)
            M_dot = M_dot_BH
            L_BH = 0d0  ! No feedback
            
            ! Check available mass
            available_mass = s% xmstar
            if (available_mass <= 0d0) then
               write(*,*) 'ERROR: No available envelope mass for accretion'
               write(*,*) 's% xmstar (envelope mass in g):', available_mass
               write(*,*) 's% mstar (total mass in g):', s% mstar
               write(*,*) 's% M_center (core mass in Msun):', s% M_center / Msun
               ierr = -1
               return
            endif
            
            ! Limit accretion by available mass
            max_dm = max_accretion_fraction * available_mass
            dm = min(M_dot * dt, max_dm)
            
            if (M_dot * dt > max_dm) then
               write(*,*) 'WARNING: Accretion limited by available mass'
               write(*,*) 'Requested dm/Msun:', M_dot * dt / Msun
               write(*,*) 'Available envelope mass/Msun:', available_mass / Msun
               write(*,*) 'Actual dm/Msun:', dm / Msun
               M_dot = dm / dt
               s% max_timestep = timestep_factor * max_dm / M_dot
            else
               s% max_timestep = timestep_factor * M_BH / M_dot
            endif
            
            ! Update black hole properties (spherical inflow uses local j_B)
            j_acc = j_B
            M_BH_new = M_BH + dm
            dj = dm * j_acc
            J_BH_new = J_BH + dj
            
            ! Update Bondi radius with new BH mass
            R_B = compute_bondi_radius(G, M_BH_new, c_s)
            
         end subroutine bondi_accretion_no_disk
         
         
         subroutine disk_accretion_with_feedback(s, M_BH, J_BH, j_B, G, c_s, rho, opacity, gamma1, dt, &
                                                 rad_eff, con_eff, timestep_factor, &
                                                 M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B)
            type (star_info), pointer :: s
            real(dp), intent(in) :: M_BH, J_BH, j_B, G, c_s, rho, opacity, gamma1, dt
            real(dp), intent(in) :: rad_eff, con_eff, timestep_factor
            real(dp), intent(out) :: M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B
            
            real(dp) :: c2, M_dot_BH, L_Bondi, L_Edd, R_ISCO, j_ISCO, j_acc, dj
            
            ! write(*,*) 'Disk Accretion'
            c2 = clight**2
            
            ! Accretion rate limited by convection
            M_dot_BH = 16d0 * pi / (rad_eff / (1d0 - rad_eff)) * con_eff / gamma1 / c_s * &
                       rho * (G * M_BH)**2 / c2
            L_Bondi = (rad_eff / (1d0 - rad_eff)) * M_dot_BH * c2
            
            ! Eddington limit (optional)
            L_Edd = 4d0 * pi * clight * G * M_BH / opacity
            if (s% x_logical_ctrl(2)) then
               L_BH = min(L_Bondi, L_Edd)
            else
               L_BH = L_Bondi
            endif
            
            ! Actual mass accretion onto BH
            M_dot = L_BH / (rad_eff * c2)
            dm = (1d0 - rad_eff) * M_dot * dt
            s% max_timestep = timestep_factor * M_BH / ((1d0 - rad_eff) * M_dot)
            
            ! Compute Kerr ISCO properties
            call compute_kerr_ISCO(M_BH, J_BH, j_B, G, R_ISCO, j_ISCO)
            
            ! Accretion occurs at ISCO
            j_acc = sign(1d0, j_B) * abs(j_ISCO)
            
            ! Update black hole properties
            M_BH_new = M_BH + dm
            dj = dm * j_acc
            J_BH_new = J_BH + dj
            
            ! Update Bondi radius with new BH mass
            R_B = compute_bondi_radius(G, M_BH_new, c_s)
            
         end subroutine disk_accretion_with_feedback
         

         function smoothstep(x, x_min, x_max) result(f)
            ! Smooth interpolation function (Hermite cubic)
            ! Returns 0 for x <= x_min, 1 for x >= x_max, smooth transition in between
            real(dp), intent(in) :: x, x_min, x_max
            real(dp) :: f, t
            
            if (x <= x_min) then
               f = 0d0
            else if (x >= x_max) then
               f = 1d0
            else
               t = (x - x_min) / (x_max - x_min)
               f = t * t * (3d0 - 2d0 * t)  ! 3t^2 - 2t^3
            endif
            
         end function smoothstep


         subroutine smooth_accretion_transition(s, M_BH, J_BH, j_B, j_B_div_j_ISCO, &
                                                G, c_s, rho, opacity, gamma1, dt, &
                                                rad_eff, con_eff, timestep_factor, &
                                                max_accretion_fraction, transition_width, &
                                                M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B, ierr)
            type (star_info), pointer :: s
            real(dp), intent(in) :: M_BH, J_BH, j_B, j_B_div_j_ISCO
            real(dp), intent(in) :: G, c_s, rho, opacity, gamma1, dt
            real(dp), intent(in) :: rad_eff, con_eff, timestep_factor, max_accretion_fraction
            real(dp), intent(inout) :: transition_width
            real(dp), intent(out) :: M_BH_new, J_BH_new, M_dot, dm, L_BH, R_B
            integer, intent(out) :: ierr
            
            ! Variables for Bondi regime (no disk)
            real(dp) :: M_dot_bondi, dm_bondi, L_BH_bondi, M_BH_new_bondi, J_BH_new_bondi, R_B_bondi
            
            ! Variables for disk regime (with feedback)
            real(dp) :: M_dot_disk, dm_disk, L_BH_disk, M_BH_new_disk, J_BH_new_disk, R_B_disk
            
            ! Blending variables
            real(dp) :: alpha, j_min, j_max
            integer :: ierr_bondi
            
            ierr = 0
            
            ! Define transition region
            ! Default: transition from 0.5 to 1.5 in j_B/j_ISCO
            if (transition_width <= 0d0) transition_width = 0.5d0
            j_min = 1d0 - transition_width
            j_max = 1d0 + transition_width
            
            ! Compute blending factor (0 = pure Bondi, 1 = pure disk)
            alpha = smoothstep(j_B_div_j_ISCO, j_min, j_max)
            
            ! Compute both accretion regimes
            call bondi_accretion_no_disk(s, M_BH, J_BH, j_B, G, c_s, rho, opacity, dt, &
                                       timestep_factor, max_accretion_fraction, &
                                       M_BH_new_bondi, J_BH_new_bondi, M_dot_bondi, &
                                       dm_bondi, L_BH_bondi, R_B_bondi, ierr_bondi)
            
            if (ierr_bondi /= 0) then
               ierr = ierr_bondi
               return
            endif
            
            call disk_accretion_with_feedback(s, M_BH, J_BH, j_B, G, c_s, rho, opacity, gamma1, dt, &
                                             rad_eff, con_eff, timestep_factor, &
                                             M_BH_new_disk, J_BH_new_disk, M_dot_disk, &
                                             dm_disk, L_BH_disk, R_B_disk)
            
            ! Blend the results
            M_dot   = (1d0 - alpha) * M_dot_bondi   + alpha * M_dot_disk
            dm      = (1d0 - alpha) * dm_bondi      + alpha * dm_disk
            L_BH    = (1d0 - alpha) * L_BH_bondi    + alpha * L_BH_disk
            R_B     = (1d0 - alpha) * R_B_bondi     + alpha * R_B_disk
            
            ! For BH properties, use dm-weighted average to conserve mass/angular momentum
            M_BH_new = M_BH + dm
            
            ! Blend angular momentum update
            J_BH_new = (1d0 - alpha) * J_BH_new_bondi + alpha * J_BH_new_disk
            
            ! Set timestep appropriately
            if (alpha < 0.5d0) then
               ! More Bondi-like
               s% max_timestep = timestep_factor * M_BH / max(M_dot_bondi, 1d-99)
            else
               ! More disk-like
               s% max_timestep = timestep_factor * M_BH / ((1d0 - rad_eff) * max(M_dot_disk, 1d-99))
            endif
            
            ! Print transition information
            if (alpha > 0d0 .and. alpha < 1d0) then
               write(*,*) 'SMOOTH TRANSITION MODE'
               write(*,*) '  j_B/j_ISCO = ', j_B_div_j_ISCO
               write(*,*) '  Blending factor (alpha) = ', alpha
               write(*,*) '  Bondi weight = ', 1d0 - alpha
               write(*,*) '  Disk weight = ', alpha
            else if (alpha == 0d0) then
               write(*,*) 'Pure Bondi Accretion (no disk)'
            else
               write(*,*) 'Pure Disk Accretion'
            endif
            
         end subroutine smooth_accretion_transition
         
         subroutine update_stellar_structure(s, id, startup, M_BH_new, J_BH_new, L_BH, R_B, M_dot, &
                                             rad_eff, dt, rho, ierr)
            type (star_info), pointer :: s
            integer, intent(in) :: id
            logical, intent(in) :: startup
            real(dp), intent(in) :: M_BH_new, J_BH_new, L_BH, R_B, M_dot, rad_eff, dt, rho
            integer, intent(out) :: ierr
            
            real(dp) :: M_cav, new_core_mass, core_avg_rho, core_avg_eps
            
            ! Calculate cavity mass and new core mass
            M_cav = 8d0 * pi / 3d0 * rho * (R_B**3)
            new_core_mass = (M_BH_new + M_cav) / Msun
            
            ! Average core properties
            core_avg_eps = L_BH / (new_core_mass * Msun)
            core_avg_rho = 1d0 / (4d0 / 3d0 * pi) * (new_core_mass * Msun) / (R_B**3)
            
            if (startup) then
               ! Relax core during startup
               call star_relax_core( &
                    id, new_core_mass, s% job% dlg_core_mass_per_step, &
                    s% job% relax_core_years_for_dt, &
                    core_avg_rho, core_avg_eps, ierr)
            else
               ! Update stellar properties
               s% M_center = new_core_mass * Msun
               s% mstar = s% mstar - rad_eff * M_dot * dt
               
               if (s% mstar - s% M_center < 0d0) then
                  write(*,*) 'M_Center > M_Star: Stopping Calculation'
                  return
               endif
               
               s% xmstar = s% mstar - s% M_center
               s% L_center = L_BH
               call do1_relax_R_center(s, R_B, ierr)
            endif
            
         end subroutine update_stellar_structure
         
         
         subroutine store_diagnostics(s, M_BH_new, J_BH_new, L_BH, R_B, M_dot, dm, dt, &
                                      rad_eff, opacity, P_rad, P_gas, nabla_ad, j_B_div_j_ISCO)
            type (star_info), pointer :: s
            real(dp), intent(in) :: M_BH_new, J_BH_new, L_BH, R_B, M_dot, dm, dt
            real(dp), intent(in) :: rad_eff, opacity, P_rad, P_gas, nabla_ad, j_B_div_j_ISCO
            
            real(dp) :: c2, L_Bondi, L_Edd, M_cav, rho
            
            rho = s% rho(s% nz)
            c2 = clight**2
            
            ! Recompute some diagnostic quantities
            L_Bondi = (rad_eff / (1d0 - rad_eff)) * M_dot * c2
            L_Edd = 4d0 * pi * clight * s% cgrav(s% nz) * M_BH_new / opacity
            M_cav = 8d0 * pi / 3d0 * rho * (R_B**3)
            
            s% xtra(1)  = M_BH_new
            s% xtra(2)  = L_BH
            s% xtra(3)  = R_B
            s% xtra(4)  = M_dot
            s% xtra(5)  = safe_log10(dm) - safe_log10(dt)
            s% xtra(6)  = rad_eff
            s% xtra(7)  = opacity
            s% xtra(8)  = L_Bondi
            s% xtra(9)  = L_Edd
            s% xtra(10) = M_cav
            s% xtra(11) = P_rad
            s% xtra(12) = P_gas
            s% xtra(13) = nabla_ad
            s% xtra(14) = safe_log10(j_B_div_j_ISCO)
            s% xtra(15) = J_BH_new
            
         end subroutine store_diagnostics
         
         
         subroutine print_bh_summary(s, M_BH_new, L_BH, M_dot, rad_eff, j_B_div_j_ISCO)
            type (star_info), pointer :: s
            real(dp), intent(in) :: M_BH_new, L_BH, M_dot, rad_eff, j_B_div_j_ISCO
            
            real(dp) :: new_core_mass
            
            new_core_mass = s% M_center / Msun
            
            print*, '--- Black Hole Properties ---'
            print*, 'M/M_sun: ',              s% mstar / Msun
            print*, 'M_BH/M_sun: ',           M_BH_new / Msun
            print*, 'L_BH/L_sun: ',           L_BH / Lsun
            print*, 'new_core_mass/M_sun: ',  new_core_mass
            print*, 'M_dot (M_sun/yr):',      M_dot / Msun * secyer
            print*, 'radiative efficiency: ', rad_eff
            print*, 'log (j_Bondi / j_ISCO)', safe_log10(j_B_div_j_ISCO)
            print*, '-----------------------------'
            
         end subroutine print_bh_summary
         
         
         subroutine extras_startup(id, restart, ierr)
             integer, intent(in) :: id
             logical, intent(in) :: restart
             integer, intent(out) :: ierr
             type (star_info), pointer :: s
   
             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return
   
             s%lxtra(16:22) = .false.
   
             if (s% x_logical_ctrl(1)) then
                 s% xtra(1) = s% job% new_core_mass * Msun
                 call black_hole_accretion(id, s, .true., ierr) 
             end if
         end subroutine extras_startup
         
         
         integer function extras_start_step(id)
             integer, intent(in) :: id
             integer :: ierr
             type (star_info), pointer :: s
   
             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return
             extras_start_step = 0
   
             if (s% x_logical_ctrl(1)) then
                 call black_hole_accretion(id, s, .false., ierr) 
             end if
         end function extras_start_step
         
         
         ! returns either keep_going, retry, or terminate.
         integer function extras_check_model(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_check_model = keep_going
            if (.false. .and. s% star_mass_h1 < 0.35d0) then
               ! stop when star hydrogen mass drops to specified level
               extras_check_model = terminate
               write(*, *) 'have reached desired hydrogen mass'
               return
            end if
            
            ! terminate if the BH mass exceeds x_ctrl(4)
           ! if (s% x_ctrl(4) > 0 .and. s% xtra(1) / Msun > s% x_ctrl(4)) then
           !    extras_check_model = terminate
           !    termination_code_str(t_xtra1) = 'black hole'
           ! end if
            
            ! terminate if the BH + Cavity mass exceeds the stellar mass 
            if (s% xtra(1) + s% xtra(10) >= s% mstar) then
               extras_check_model = terminate
               write(*, *) 'M_BH >= Mstar'
               termination_code_str(t_xtra1) = 'black hole'
            end if
            
            ! if you want to check multiple conditions, it can be useful
            ! to set a different termination code depending on which
            ! condition was triggered.  MESA provides 9 customizeable
            ! termination codes, named t_xtra1 .. t_xtra9.  You can
            ! customize the messages that will be printed upon exit by
            ! setting the corresponding termination_code_str value.
            ! termination_code_str(t_xtra1) = 'my termination condition'
   
            ! by default, indicate where (in the code) MESA terminated
            if (extras_check_model == terminate) s% termination_code = t_extras_check_model
         end function extras_check_model
   
   
         integer function how_many_extra_history_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_columns = 16
         end function how_many_extra_history_columns
         
         
         subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            
            integer :: i
            real(dp) :: X0, mX0, rX0
            
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            
            names(1:16) = 'empty'
            vals(1:16) = -1d99
            
   
            if (s% x_logical_ctrl(1)) then 
                names(1) = "M_BH"
                vals(1) = s% xtra(1) / Msun   ! M_BH / Msun
                names(2) = "L_BH"
                vals(2)  = s% xtra(2) / Lsun  ! L_BH / Lsun
                names(3) = "R_B"
                vals(3)  = s% xtra(3) / Rsun  ! R_B  / Rsun
                names(4) = "M_dot"
                vals(4)  = s% xtra(4) / Msun  ! M_dot / Msun
                names(5) = "log10(dm/dt)"
                vals(5)  = s% xtra(5)         ! g/s
                names(6) = "rad_eff"
                vals(6)  = s% xtra(6)         ! epsilon  
                names(7) = "kap_center"
                vals(7)  = s% xtra(7)         ! cm^2/g
                names(8) = "L_B"
                vals(8) = s% xtra(8) / Lsun   ! Bondi luminosity 
                names(9) = "L_E"
                vals(9) = s% xtra(9) / Lsun ! Eddington luminosity 
                names(10) = "M_cav"
                vals(10) = s% xtra(10) / Msun ! M_cav / Msun
                names(11) = "prad_center"
                vals(11)  = s% xtra(11)
                names(12) = "pgas_center"
                vals(12)  = s% xtra(12)
                names(13) = "nabla_ad_center"
                vals(13)  = s% xtra(13)
                names(14) = "log10_j_B_div_j_ISCO"
                vals(14) =  s% xtra(14)   
                names(15) = "cs_center"
                vals(15) = s% csound(s% nz)   ! cm/s
                names(16) = "J_BH"
                vals(16) = s% xtra(15)   ! cm^2/s  
            end if
            
            ! note: do NOT add the extras names to history_columns.list
            ! the history_columns.list is only for the built-in history column options.
            ! it must not include the new column names you are adding here.
            
   
         end subroutine data_for_extra_history_columns
   
   
         integer function how_many_extra_profile_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_columns = 0
         end function how_many_extra_profile_columns
         
         
         subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
            use star_def, only: star_info, maxlen_profile_column_name
            use const_def, only: dp
            integer, intent(in) :: id, n, nz
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(nz,n)
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            integer :: k, op_err, net_lwork
            logical :: okay
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
       end subroutine data_for_extra_profile_columns
   
   
         integer function how_many_extra_history_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_header_items = 0
         end function how_many_extra_history_header_items
   
   
         subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return
   
            ! here is an example for adding an extra history header item
            ! also set how_many_extra_history_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha
   
         end subroutine data_for_extra_history_header_items
   
   
         integer function how_many_extra_profile_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_header_items = 0
         end function how_many_extra_profile_header_items
   
   
         subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return
   
            ! here is an example for adding an extra profile header item
            ! also set how_many_extra_profile_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha
   
         end subroutine data_for_extra_profile_header_items
   
   
         ! returns either keep_going or terminate.
         ! note: cannot request retry; extras_check_model can do that.
         integer function extras_finish_step(id)
            integer, intent(in) :: id
            integer :: ierr
      
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_finish_step = keep_going
   
   
           
            !s% xtra(1) / Msun  ! Mbh
            !s% xtra(2) ! Lbh
   
            !s% xtra(8)  = L_Bondi
            !s% xtra(9)  = L_Edd
   
            s% xtra(16) = s% L(1)  ! L_star  
            ! s% xtra(2) ! L_BH
            
            s% profile_data_suffix = '.data'
   
   
   
            
            if (s%xtra(2) / s%xtra(16) >= 0.001 .and. s%xtra(2) / s%xtra(16) < 0.01 .and. .not. s%lxtra(16)) then
                ! L_BH = 1/1000 L_star
                s% profile_data_suffix = '_Lbh_div_L_0.001.data'
                s%need_to_save_profiles_now = .true.
                s%need_to_update_history_now = .true.
                write(*,*) 'Save profile at L_BH = 1/1000 L_star'
                s%lxtra(16) = .true.
            else if (s%xtra(2) / s%xtra(16) >= 0.01 .and. s%xtra(2) / s%xtra(16) < 0.1 .and. .not. s%lxtra(17)) then
                ! L_BH = 1/100 L_star
                s% profile_data_suffix = '_Lbh_div_L_0.01.data'
                s%need_to_save_profiles_now = .true.
                s%need_to_update_history_now = .true.
                write(*,*) 'Save profile at L_BH = 1/100 L_star'
                s%lxtra(17) = .true.
            else if (s%xtra(2) / s%xtra(16) >= 0.1 .and. s%xtra(2) / s%xtra(16) < 1 .and. .not. s%lxtra(18)) then
                ! L_BH = 1/10 L_star
                s% profile_data_suffix = '_Lbh_div_L_0.1.data'
                s%need_to_save_profiles_now = .true.
                s%need_to_update_history_now = .true.
                write(*,*) 'Save profile at L_BH = 1/10 L_star'
                s%lxtra(18) = .true.
            else if (s%xtra(2) / s%xtra(16) >= 1 .and. s%xtra(2) / s%xtra(16) < 10 .and. .not. s%lxtra(19)) then
                ! L_BH = L_star
                s% profile_data_suffix = '_Lbh_div_L_1.data'
                s%need_to_save_profiles_now = .true.
                s%need_to_update_history_now = .true.
                write(*,*) 'Save profile at L_BH = L_star'
                s%lxtra(19) = .true.
            else if (s%xtra(2) / s%xtra(16) >= 10 .and. s%xtra(2) / s%xtra(16) < 100 .and. .not. s%lxtra(20)) then
                ! L_BH = 10 L_star
                s% profile_data_suffix = '_Lbh_div_L_10.data'
                s%need_to_save_profiles_now = .true.
                s%need_to_update_history_now = .true.
                write(*,*) 'Save profile at L_BH = 10 L_star'
                s%lxtra(20) = .true.
            else if (s%xtra(2) / s%xtra(16) >= 100 .and. s%xtra(2) / s%xtra(16) < 1000 .and. .not. s%lxtra(21)) then
                ! L_BH = 100 L_star
                s% profile_data_suffix = '_Lbh_div_L_100.data'
                s%need_to_save_profiles_now = .true.
                s%need_to_update_history_now = .true.
                write(*,*) 'Save profile at L_BH = 100 L_star'
                s%lxtra(21) = .true.
            else if (s%xtra(2) / s%xtra(16) >= 1000 .and. .not. s%lxtra(22)) then
                ! L_BH = 1000 L_star
                s% profile_data_suffix = '_Lbh_div_L_1000.data'
                s%need_to_save_profiles_now = .true.
                s%need_to_update_history_now = .true.
                write(*,*) 'Save profile at L_BH = 1000 L_star'
                s%lxtra(22) = .true.
            endif
            ! Save a profile at first episode of disk accretion
   
            if (s% xtra(14) >= 0.0 .and. .not. s%lxtra(23)) then 
               s% profile_data_suffix = '_Disk_Formation.data'
               s%need_to_save_profiles_now = .true.
               s%need_to_update_history_now = .true.
               write(*,*) 'Save profile at First Disk Formation. log (J/Jisco):',  s% xtra(14)
               s%lxtra(23) = .true.
            endif 

   
            ! see extras_check_model for information about custom termination codes
            ! by default, indicate where (in the code) MESA terminated
            if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
         end function extras_finish_step
         
         
         subroutine extras_after_evolve(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
         end subroutine extras_after_evolve
   
         end module run_star_extras
   