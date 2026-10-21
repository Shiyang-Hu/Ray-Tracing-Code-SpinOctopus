module metric_parameters
implicit none
real*8::spin=0.94d0
end module metric_parameters
    
module task_parameters
implicit none
integer::task_model=11
integer::cpu=28
end module task_parameters    
    
module physics_parameters
implicit none
real*8::mass_ratio=4.14d0*1d6
real*8::distance_obs=8.127d0*1d-3
real*8::light_speed=299792458d0
real*8::gravitational_constant=6.69d0*1d-11
real*8::sun_mass=1.989d0*1d30
real*8::math_pi=dacos(-1d0)
end module physics_parameters    
    
module observer_parameters
implicit none
real*8::t_obs=0d0
real*8::r_obs=1000.0d0
real*8::theta_obs=0.29670597283903605d0
real*8::phi_obs=0.0d0
integer::resolution_x=200
integer::resolution_y=200
real*8::x_ini=-15.0d0,x_end=15.0d0
real*8::y_ini=-15.0d0,y_end=15.0d0
real*8::delta_x=0d0
real*8::delta_y=0d0
real*8::vt_obs=1d0
real*8::vr_obs=0d0
real*8::vtheta_obs=0d0
real*8::vphi_obs=0d0
real*8::x_start=4.0d0
real*8::y_start=4.0d0
real*8::z_start=0d0
end module observer_parameters
    
module integration_parameters
implicit none
real*8::solution_error=1d-8
end module integration_parameters
     
module raytracing_parameters
implicit none
real*8::hit_error=0.001d0
real*8::observer=1200d0
real*8::ini_step=0.0025d0
real*8::step_ratio=1.8d0
real*8::theta_correction=0.01d0
end module raytracing_parameters    
    
module disk_range_parameters
implicit none
real*8::r_inner=4.0d0
real*8::r_outer=500.0d0
integer::disk_model=1
end module disk_range_parameters    
    
module emission_scale_parameters
implicit none
real*8::j0=1.0d0
real*8::a0=0.0d0
end module emission_scale_parameters
  
module thick_equatorial_disk_parameters
implicit none
real*8::disk_intensity=1.0d0
real*8::thick_emit_decay1=-1.5d0
real*8::thick_emit_decay2=-0.5d0
real*8::thickness=0.2d0
real*8::disk_extend=0.1d0
real*8::platform_intensity=1.0d0
real*8::platform_location=5.0d0
real*8::platform_width=3.0d0
real*8::bump_intensity=1.0d0
real*8::bump_location=9.0d0
real*8::sigma_bump_r=1.0d0
real*8::sigma_bump_theta=0.2d0
real*8::blob_intensity=0.0d0
real*8::blob_r=11.0d0
real*8::blob_theta=1.5707963267948966d0
real*8::blob_phi=2.792526803190927d0
real*8::sigma_blob_r=2.0d0
real*8::sigma_blob_theta=0.17453292519943295d0
real*8::sigma_blob_phi=0.3490658503988659d0
real*8::inspect_azimuth=4.71238898038469d0
end module thick_equatorial_disk_parameters
   
module thin_equatorial_disk_parameters
implicit none
real*8::fudge_factor=0.6666666666666666d0
real*8::thin_emit_decay1=-2.0d0
real*8::thin_emit_decay2=-0.5d0
integer::max_hit=3
end module thin_equatorial_disk_parameters

module equatorial_accretion_dynamics_parameters
implicit none
real*8::velocity_max_r=0.05d0
real*8::velocity_max_theta=0.2d0
real*8::velocity_max_phi=0.4d0
real*8::velocity_decay_r=0.5d0
real*8::velocity_decay_theta=0.5d0
real*8::velocity_decay_phi=0.3d0
end module equatorial_accretion_dynamics_parameters

module point_source_parameters
implicit none
real*8::point_source_x=-6.0d0
real*8::point_source_y=0.0d0
real*8::point_source_z=0.0d0
real*8::point_source_radius=0.5d0
integer::show_line=0
integer::delta_latitude=2
integer::delta_longitude=2
real*8::point_source_sigma=0.2d0
real*8::point_source_vt=1.0d0
real*8::point_source_vr=0.2d0
real*8::point_source_vtheta=0.0d0
real*8::point_source_vphi=0.5d0
end module point_source_parameters

module jet_parameters
implicit none
real*8::jet_intensity=1.0d0
real*8::jet_base=2.0d0
real*8::jet_decay=2.0d0
real*8::jet_angle_ini=0.05d0
real*8::jet_angle_inf=0.26d0
real*8::jet_angle_expansion=20.0d0
real*8::jet_velocity_ini=0.1d0
real*8::jet_velocity_inf=0.8d0
real*8::jet_range=200.0d0
real*8::jet_omega=0.0d0
real*8::jet_omega_decay=0.5d0
real*8::gyration_radius=0.01d0
end module jet_parameters

module massive_particle_orbit_parameters
implicit none
real*8::t_mp=0.0d0
real*8::r_mp=15.0d0
real*8::theta_mp=1.5707963267948966d0
real*8::phi_mp=0.0d0
real*8::pt_mp=-0.96d0
real*8::pr_mp=0.05d0
real*8::ptheta_mp=0.1d0
real*8::pphi_mp=3.55d0
real*8::step_mp=-0.1d0
real*8::particle_radius=0.25d0
real*8::motion_time=1000d0
integer::max_point=500
integer::sample=50
integer::condition_type=3
end module massive_particle_orbit_parameters

module spherically_accretion_flow_parameters
implicit none
real*8::spherically_flow_intensity=1.0d0
real*8::spherically_flow_outer=10.0d0
real*8::spherically_flow_inner=6.0d0
real*8::spherically_flow_peak=8.0d0
real*8::infalling_flow_max=0.2d0
real*8::infalling_accelerating=0.5d0
real*8::emission_brake=0.5d0
real*8::spherically_flow_sigma=1.0d0
end module spherically_accretion_flow_parameters

module tilted_thin_disk_parameters
implicit none
real*8::disk_inclination=0.5235987755982988d0
end module tilted_thin_disk_parameters

subroutine covariant_metric(r,theta,g_ab) 
use metric_parameters
implicit none
real*8,dimension(4,4)::g_ab
real*8 r,theta
g_ab(1,1) = 2.0d0*r/(r**2 + spin**2*dcos(theta)**2) - 1.0d0
g_ab(1,2) = 0.0d0
g_ab(1,3) = 0.0d0
g_ab(1,4) = -2.0d0*r*spin*dsin(theta)**2/(r**2 + spin**2*dcos(theta)**2)
g_ab(2,1) = 0.0d0
g_ab(2,2) = (r**2 + spin**2*dcos(theta)**2)/(r**2 - 2.0d0*r + spin**2)
g_ab(2,3) = 0.0d0
g_ab(2,4) = 0.0d0
g_ab(3,1) = 0.0d0
g_ab(3,2) = 0.0d0
g_ab(3,3) = r**2 + spin**2*dcos(theta)**2
g_ab(3,4) = 0.0d0
g_ab(4,1) = -2.0d0*r*spin*dsin(theta)**2/(r**2 + spin**2*dcos(theta)**2)
g_ab(4,2) = 0.0d0
g_ab(4,3) = 0.0d0
g_ab(4,4) = (r**2 + 2.0d0*r*spin**2*dsin(theta)**2/(r**2 + spin**2*dcos(theta)**2) + &
      spin**2)*dsin(theta)**2
return
end subroutine
    
subroutine contravariant_metric(r,theta,gab) 
use metric_parameters
implicit none
real*8,dimension(4,4)::gab
real*8 r,theta
gab(1,1) = (-1.0d0*r**4 - 1.0d0*r**2*spin**2*dcos(theta)**2 - 1.0d0*r**2*spin**2 - &
      2.0d0*r*spin**2*dsin(theta)**2 - 1.0d0*spin**4*dcos(theta)**2)/( &
      1.0d0*r**4 - 2.0d0*r**3 + 1.0d0*r**2*spin**2*dcos(theta)**2 + &
      1.0d0*r**2*spin**2 + 2.0d0*r*spin**2*dsin(theta)**2 - 2.0d0*r* &
      spin**2 + 1.0d0*spin**4*dcos(theta)**2)
gab(1,2) = 0
gab(1,3) = 0
gab(1,4) = -2.0d0*r*spin/(1.0d0*r**4 - 2.0d0*r**3 + 1.0d0*r**2*spin**2*dcos(theta) &
      **2 + 1.0d0*r**2*spin**2 + 2.0d0*r*spin**2*dsin(theta)**2 - 2.0d0 &
      *r*spin**2 + 1.0d0*spin**4*dcos(theta)**2)
gab(2,1) = 0
gab(2,2) = (1.0d0*r**2 - 2.0d0*r + 1.0d0*spin**2)/(1.0d0*r**2 + 1.0d0*spin**2*dcos( &
      theta)**2)
gab(2,3) = 0
gab(2,4) = 0
gab(3,1) = 0
gab(3,2) = 0
gab(3,3) = 1d0/(r**2 + spin**2*dcos(theta)**2)
gab(3,4) = 0
gab(4,1) = -2.0d0*r*spin/(1.0d0*r**4 - 2.0d0*r**3 + 1.0d0*r**2*spin**2*dcos(theta) &
      **2 + 1.0d0*r**2*spin**2 + 2.0d0*r*spin**2*dsin(theta)**2 - 2.0d0 &
      *r*spin**2 + 1.0d0*spin**4*dcos(theta)**2)
gab(4,2) = 0
gab(4,3) = 0
gab(4,4) = (1.0d0*r**2 - 2.0d0*r + 1.0d0*spin**2*dcos(theta)**2)/(1.0d0*r**4*dsin( &
      theta)**2 - 2.0d0*r**3*dsin(theta)**2 + 1.0d0*r**2*spin**2*dsin( &
      theta)**2*dcos(theta)**2 + 1.0d0*r**2*spin**2*dsin(theta)**2 + &
      2.0d0*r*spin**2*dsin(theta)**4 - 2.0d0*r*spin**2*dsin(theta)**2 + &
      1.0d0*spin**4*dsin(theta)**2*dcos(theta)**2)
return
end subroutine    
    
subroutine contravariant_metric_diff(r,theta,gab_r,gab_theta)
use metric_parameters
implicit none
real*8,dimension(4,4)::gab_r,gab_theta
real*8 r,theta
gab_r(1,1) = (-4.0d0*r**3 - 2.0d0*r*spin**2*dcos(theta)**2 - 2.0d0*r*spin**2 - 2.0d0* &
      spin**2*dsin(theta)**2)/(1.0d0*r**4 - 2.0d0*r**3 + 1.0d0*r**2* &
      spin**2*dcos(theta)**2 + 1.0d0*r**2*spin**2 + 2.0d0*r*spin**2* &
      dsin(theta)**2 - 2.0d0*r*spin**2 + 1.0d0*spin**4*dcos(theta)**2) &
      + 0.25d0*(-1.0d0*r**4 - 1.0d0*r**2*spin**2*dcos(theta)**2 - 1.0d0 &
      *r**2*spin**2 - 2.0d0*r*spin**2*dsin(theta)**2 - 1.0d0*spin**4* &
      dcos(theta)**2)*(-4.0d0*r**3 + 6.0d0*r**2 - 2.0d0*r*spin**2*dcos( &
      theta)**2 - 2.0d0*r*spin**2 - 2.0d0*spin**2*dsin(theta)**2 + &
      2.0d0*spin**2)/(0.5d0*r**4 - r**3 + 0.5d0*r**2*spin**2*dcos(theta &
      )**2 + 0.5d0*r**2*spin**2 + r*spin**2*dsin(theta)**2 - r*spin**2 &
      + 0.5d0*spin**4*dcos(theta)**2)**2
gab_r(1,2) = 0
gab_r(1,3) = 0
gab_r(1,4) = -0.5d0*r*spin*(-4.0d0*r**3 + 6.0d0*r**2 - 2.0d0*r*spin**2*dcos(theta)**2 &
      - 2.0d0*r*spin**2 - 2.0d0*spin**2*dsin(theta)**2 + 2.0d0*spin**2) &
      /(0.5d0*r**4 - r**3 + 0.5d0*r**2*spin**2*dcos(theta)**2 + 0.5d0*r &
      **2*spin**2 + r*spin**2*dsin(theta)**2 - r*spin**2 + 0.5d0*spin** &
      4*dcos(theta)**2)**2 - 2.0d0*spin/(1.0d0*r**4 - 2.0d0*r**3 + &
      1.0d0*r**2*spin**2*dcos(theta)**2 + 1.0d0*r**2*spin**2 + 2.0d0*r* &
      spin**2*dsin(theta)**2 - 2.0d0*r*spin**2 + 1.0d0*spin**4*dcos( &
      theta)**2)
gab_r(2,1) = 0
gab_r(2,2) = -2.0d0*r*(1.0d0*r**2 - 2.0d0*r + 1.0d0*spin**2)/(1.0d0*r**2 + 1.0d0*spin &
      **2*dcos(theta)**2)**2 + (2.0d0*r - 2.0d0)/(1.0d0*r**2 + 1.0d0* &
      spin**2*dcos(theta)**2)
gab_r(2,3) = 0
gab_r(2,4) = 0
gab_r(3,1) = 0
gab_r(3,2) = 0
gab_r(3,3) = -2*r/(r**2 + spin**2*dcos(theta)**2)**2
gab_r(3,4) = 0
gab_r(4,1) = -0.5d0*r*spin*(-4.0d0*r**3 + 6.0d0*r**2 - 2.0d0*r*spin**2*dcos(theta)**2 &
      - 2.0d0*r*spin**2 - 2.0d0*spin**2*dsin(theta)**2 + 2.0d0*spin**2) &
      /(0.5d0*r**4 - r**3 + 0.5d0*r**2*spin**2*dcos(theta)**2 + 0.5d0*r &
      **2*spin**2 + r*spin**2*dsin(theta)**2 - r*spin**2 + 0.5d0*spin** &
      4*dcos(theta)**2)**2 - 2.0d0*spin/(1.0d0*r**4 - 2.0d0*r**3 + &
      1.0d0*r**2*spin**2*dcos(theta)**2 + 1.0d0*r**2*spin**2 + 2.0d0*r* &
      spin**2*dsin(theta)**2 - 2.0d0*r*spin**2 + 1.0d0*spin**4*dcos( &
      theta)**2)
gab_r(4,2) = 0
gab_r(4,3) = 0
gab_r(4,4) = (2.0d0*r - 2.0d0)/(1.0d0*r**4*dsin(theta)**2 - 2.0d0*r**3*dsin(theta)**2 &
      + 1.0d0*r**2*spin**2*dsin(theta)**2*dcos(theta)**2 + 1.0d0*r**2* &
      spin**2*dsin(theta)**2 + 2.0d0*r*spin**2*dsin(theta)**4 - 2.0d0*r &
      *spin**2*dsin(theta)**2 + 1.0d0*spin**4*dsin(theta)**2*dcos(theta &
      )**2) + 0.25d0*(1.0d0*r**2 - 2.0d0*r + 1.0d0*spin**2*dcos(theta) &
      **2)*(-4.0d0*r**3*dsin(theta)**2 + 6.0d0*r**2*dsin(theta)**2 - &
      2.0d0*r*spin**2*dsin(theta)**2*dcos(theta)**2 - 2.0d0*r*spin**2* &
      dsin(theta)**2 - 2.0d0*spin**2*dsin(theta)**4 + 2.0d0*spin**2* &
      dsin(theta)**2)/(0.5d0*r**4*dsin(theta)**2 - r**3*dsin(theta)**2 &
      + 0.5d0*r**2*spin**2*dsin(theta)**2*dcos(theta)**2 + 0.5d0*r**2* &
      spin**2*dsin(theta)**2 + r*spin**2*dsin(theta)**4 - r*spin**2* &
      dsin(theta)**2 + 0.5d0*spin**4*dsin(theta)**2*dcos(theta)**2)**2
gab_theta(1,1) = 0.25d0*(2.0d0*r**2*spin**2*dsin(theta)*dcos(theta) - 4.0d0*r*spin**2* &
      dsin(theta)*dcos(theta) + 2.0d0*spin**4*dsin(theta)*dcos(theta))* &
      (-1.0d0*r**4 - 1.0d0*r**2*spin**2*dcos(theta)**2 - 1.0d0*r**2* &
      spin**2 - 2.0d0*r*spin**2*dsin(theta)**2 - 1.0d0*spin**4*dcos( &
      theta)**2)/(0.5d0*r**4 - r**3 + 0.5d0*r**2*spin**2*dcos(theta)**2 &
      + 0.5d0*r**2*spin**2 + r*spin**2*dsin(theta)**2 - r*spin**2 + &
      0.5d0*spin**4*dcos(theta)**2)**2 + (2.0d0*r**2*spin**2*dsin(theta &
      )*dcos(theta) - 4.0d0*r*spin**2*dsin(theta)*dcos(theta) + 2.0d0* &
      spin**4*dsin(theta)*dcos(theta))/(1.0d0*r**4 - 2.0d0*r**3 + 1.0d0 &
      *r**2*spin**2*dcos(theta)**2 + 1.0d0*r**2*spin**2 + 2.0d0*r*spin &
      **2*dsin(theta)**2 - 2.0d0*r*spin**2 + 1.0d0*spin**4*dcos(theta) &
      **2)
gab_theta(1,2) = 0
gab_theta(1,3) = 0
gab_theta(1,4) = -0.5d0*r*spin*(2.0d0*r**2*spin**2*dsin(theta)*dcos(theta) - 4.0d0*r*spin &
      **2*dsin(theta)*dcos(theta) + 2.0d0*spin**4*dsin(theta)*dcos( &
      theta))/(0.5d0*r**4 - r**3 + 0.5d0*r**2*spin**2*dcos(theta)**2 + &
      0.5d0*r**2*spin**2 + r*spin**2*dsin(theta)**2 - r*spin**2 + 0.5d0 &
      *spin**4*dcos(theta)**2)**2
gab_theta(2,1) = 0
gab_theta(2,2) = 2.0d0*spin**2*(1.0d0*r**2 - 2.0d0*r + 1.0d0*spin**2)*dsin(theta)*dcos( &
      theta)/(1.0d0*r**2 + 1.0d0*spin**2*dcos(theta)**2)**2
gab_theta(2,3) = 0
gab_theta(2,4) = 0
gab_theta(3,1) = 0
gab_theta(3,2) = 0
gab_theta(3,3) = 2*spin**2*dsin(theta)*dcos(theta)/(r**2 + spin**2*dcos(theta)**2)**2
gab_theta(3,4) = 0
gab_theta(4,1) = -0.5d0*r*spin*(2.0d0*r**2*spin**2*dsin(theta)*dcos(theta) - 4.0d0*r*spin &
      **2*dsin(theta)*dcos(theta) + 2.0d0*spin**4*dsin(theta)*dcos( &
      theta))/(0.5d0*r**4 - r**3 + 0.5d0*r**2*spin**2*dcos(theta)**2 + &
      0.5d0*r**2*spin**2 + r*spin**2*dsin(theta)**2 - r*spin**2 + 0.5d0 &
      *spin**4*dcos(theta)**2)**2
gab_theta(4,2) = 0
gab_theta(4,3) = 0
gab_theta(4,4) = -2.0d0*spin**2*dsin(theta)*dcos(theta)/(1.0d0*r**4*dsin(theta)**2 - &
      2.0d0*r**3*dsin(theta)**2 + 1.0d0*r**2*spin**2*dsin(theta)**2* &
      dcos(theta)**2 + 1.0d0*r**2*spin**2*dsin(theta)**2 + 2.0d0*r*spin &
      **2*dsin(theta)**4 - 2.0d0*r*spin**2*dsin(theta)**2 + 1.0d0*spin &
      **4*dsin(theta)**2*dcos(theta)**2) + 0.25d0*(1.0d0*r**2 - 2.0d0*r &
      + 1.0d0*spin**2*dcos(theta)**2)*(-2.0d0*r**4*dsin(theta)*dcos( &
      theta) + 4.0d0*r**3*dsin(theta)*dcos(theta) + 2.0d0*r**2*spin**2* &
      dsin(theta)**3*dcos(theta) - 2.0d0*r**2*spin**2*dsin(theta)*dcos( &
      theta)**3 - 2.0d0*r**2*spin**2*dsin(theta)*dcos(theta) - 8.0d0*r* &
      spin**2*dsin(theta)**3*dcos(theta) + 4.0d0*r*spin**2*dsin(theta)* &
      dcos(theta) + 2.0d0*spin**4*dsin(theta)**3*dcos(theta) - 2.0d0* &
      spin**4*dsin(theta)*dcos(theta)**3)/(0.5d0*r**4*dsin(theta)**2 - &
      r**3*dsin(theta)**2 + 0.5d0*r**2*spin**2*dsin(theta)**2*dcos( &
      theta)**2 + 0.5d0*r**2*spin**2*dsin(theta)**2 + r*spin**2*dsin( &
      theta)**4 - r*spin**2*dsin(theta)**2 + 0.5d0*spin**4*dsin(theta) &
      **2*dcos(theta)**2)**2
return
end subroutine
    
subroutine delta_function(r,theta,delta)
use metric_parameters
implicit none
real*8 r,theta
real*8 delta
delta = r**2 - 2.0d0*r + spin**2
return
end subroutine 
    
subroutine delta_function_diff(r,theta,delta_diff)
use metric_parameters
implicit none
real*8 r,theta
real*8 delta_diff
delta_diff = 2*r - 2.0d0
return
end subroutine
