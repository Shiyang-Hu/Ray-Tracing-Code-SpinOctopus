#-------------------------------------------------------------------------------------------------------------------------------#
import sympy as sp
import numpy as np
#-------------------------------------------------------------------------------------------------------------------------------#
# Define the symbolic variables used in the metric here.
# "t", "r", "theta", and "phi" are the generalized coordinates, and "spin" represents the spin parameter of the central object. 
# DO NOT modify these 5 core variables.
# If the metric involves additional parameters, append them sequentially after "spin".
# For instance, if the metric includes a deformation parameter denoted by "deform", this line should be written as:
# t, r, theta, phi, spin, deform = sp.symbols('t r theta phi spin deform')
# Important compatibility note:
# If the target spacetime is spherically symmetric, DO NOT delete the "spin" variable.
# Instead, simply set its numerical value to 0.0 in the assignment section below.
#-------------------------------------------------------------------------------------------------------------------------------#
t, r, theta, phi, spin = sp.symbols('t r theta phi spin')
#-------------------------------------------------------------------------------------------------------------------------------#
# Define the covariant metric components of the target spacetime below.
# Indices 1 through 4 correspond to the coordinates "t", "r", "theta", and "phi", respectively.
# Note: Please input all numerical constants as floating-point values (e.g., write 1.0 instead of 1, 15.0 instead of 15).
# This ensures seamless conversion into double-precision Fortran code.
#-------------------------------------------------------------------------------------------------------------------------------#
g_11=-(1.0-2.0*r/(r*r+spin*spin*sp.cos(theta)*sp.cos(theta)))
g_12=0.0
g_13=0.0
g_14=-(2.0*spin*r*sp.sin(theta)*sp.sin(theta)/(r*r+spin*spin*sp.cos(theta)*sp.cos(theta)))
#-------------------------------------------------------------------------------------------------------------------------------#
g_21=0.0
g_22=(r*r+spin*spin*sp.cos(theta)*sp.cos(theta))/(r*r-2.0*r+spin*spin)
g_23=0.0
g_24=0.0
#-------------------------------------------------------------------------------------------------------------------------------#
g_31=0.0
g_32=0.0
g_33=r*r+spin*spin*sp.cos(theta)*sp.cos(theta)
g_34=0.0
#-------------------------------------------------------------------------------------------------------------------------------#
g_41=g_14
g_42=0.0
g_43=0.0
g_44=sp.sin(theta)*sp.sin(theta)*(r*r+spin*spin+2.0*r*spin*spin*sp.sin(theta)*sp.sin(theta)/(r*r+spin*spin*sp.cos(theta)*sp.cos(theta)))
#-------------------------------------------------------------------------------------------------------------------------------#
# Please specify the denominator of the g_rr (g_22) metric component below.
# The program utilizes this expression to calculate the event horizon radius.
# For instance, in the Kerr spacetime, this denominator is delta = r*r + spin*spin - 2.0*r.
# DO NOT change the variable name "delta".
#-------------------------------------------------------------------------------------------------------------------------------#
delta=r*r+spin*spin-2.0*r
#-------------------------------------------------------------------------------------------------------------------------------#
# Please specify the numerical value for the spin of the central object, as well as the values for any additional metric parameters.
# In this framework, numerical variables are defined by appending the suffix "_value" to their corresponding symbolic names.
# For example, "spin_value" and "deform_value" represent the numerical values of the spin and deformation parameters, respectively.
# Important: Please ensure that "real*8::symbolic={symbolic_value}d0" is added 
# to the "module metric_parameters" block within the Fortran code generation section.
# For example: "real*8::spin={spin_value}d0" and "real*8::deform={deform_value}d0".
# To summarize the mapping logic: for every metric parameter whose symbolic name is defined above as "symbolic-name",
# and whose numerical value "symbolic-name_value" is assigned below, a double-precision declaration must appear 
# in the Fortran "module metric_parameters" as:
# real*8::symbolic-name={symbolic-name_value}d0
#-------------------------------------------------------------------------------------------------------------------------------#
spin_value=0.94
#{symbolic-name}_value=
#-------------------------------------------------------------------------------------------------------------------------------#
# Please specify the execution mode flag and the number of OpenMP threads.
# Note: Both inputs must be strictly integers. 
# The mode flag determines the computational task.
# task_model = 1 : Black hole critical curve
# task_model = 2 : Emission profiles of the accretion environment
# task_model = 3 : Image of a geometrically thick accretion disk
# task_model = 4 : Image of a geometrically thin accretion disk
# task_model = 5 : Four-color screen against a celestial background
# task_model = 6 : Jet emission profiles
# task_model = 7 : Jet image
# task_model = 8 : Null geodesics (photon trajectories)
# task_model = 9 : Timelike geodesics (massive particle orbits)
# task_model = 10: Point source image (gravitational lensing of a point source)
# task_model = 11: Timelike geodesics and corresponding light curves
# task_model = 12: Spherical shell accretion image
# task_model = 13: Equatorial grid image in curved spacetime (Under optimization)
# task_model = 14: Image of a tilted geometrically thin accretion disk
#-------------------------------------------------------------------------------------------------------------------------------#
task_model=3
cpu=28
#-------------------------------------------------------------------------------------------------------------------------------#
# Please specify the observer's coordinates.
# Note: All angular coordinates (such as inclination and azimuth) MUST be provided in radians.
#-------------------------------------------------------------------------------------------------------------------------------#
r_obs=1000.0
theta_obs=17.0*np.pi/180.0
phi_obs=0.0*np.pi/180.0
#-------------------------------------------------------------------------------------------------------------------------------#
# Define the image resolution (in pixels) and the size of the field of view (FOV).
# Important: The FOV size is expressed in geometric units, rather than physical units.
# |---------y_end---------|
# |           |           |
# |           |           |
# |x_ini------|------x_end|
# |           |           |
# |           |           |
# |---------y_ini---------|
#-------------------------------------------------------------------------------------------------------------------------------#
resolution_x=200
resolution_y=200
x_ini=-15.0
x_end=15.0
y_ini=-15.0
y_end=15.0
#-------------------------------------------------------------------------------------------------------------------------------#
# Please define the parameters for the adaptive step-size strategy in the ray-tracing.
# Note: These parameters significantly impact both tracing accuracy and computational efficiency.
# Troubleshooting: If abnormal pixels (artifacts) appear in the rendered image, 
# please decrease "ini_step" or "step_ratio", or increase "theta_correction".
#-------------------------------------------------------------------------------------------------------------------------------#
ini_step=0.0025
step_ratio=1.8
theta_correction=0.01
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the scaling parameters for the emission (j0) and absorption (a0) of the accreting material.
# Note: Setting a0 = 0.0 corresponds to an optically thin medium.
# In general cases, it is recommended to keep j0 = 1.0 and a0 = 0.0.
#-------------------------------------------------------------------------------------------------------------------------------#
j0=1.0
a0=0.0
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the inner and outer boundaries of the accretion environment (Thin disk or Thick disk).
# Note: When disk_model = 1, any provided value for "r_inner" is ignored, and the inner boundary is rigidly fixed at the event horizon.
# To enforce the user-defined "r_inner" value, set disk_model = 2.
#-------------------------------------------------------------------------------------------------------------------------------#
r_inner=4.0
r_outer=500.0
disk_model=1
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the geometric and emission parameters for the geometrically thick accretion environment.
# Note: The parameters in this section are strictly required for simulation tasks task_model = 2 and task_model = 3.
#------------------------------------------------------------------------------------------#
# 1. Global Disk Emission Parameters
#------------------------------------------------------------------------------------------#
# disk_intensity     : Weighting factor for the emission intensity of the thick disk.
# thick_emit_decay1  : Primary radial decay factor of the thick disk emission (must be strictly negative).
# thick_emit_decay2  : Secondary radial decay factor of the thick disk emission (must be strictly negative).
# thickness          : Initial geometrical thickness of the accretion disk.
# disk_extend        : Expansion factor dictating the disk thickness profile.
#------------------------------------------------------------------------------------------#
# 2. Emission Platform (Radiation Compensation) Parameters
#------------------------------------------------------------------------------------------#
# platform_intensity : Weighting factor for the emission platform (radiation compensation).
# platform_location  : Radial position of the emission platform.
# platform_width     : Radial width of the compensated emission region.
#------------------------------------------------------------------------------------------#
# 3. Gaussian Bump Parameters (Axisymmetric Enhancements)
#------------------------------------------------------------------------------------------#
# bump_intensity     : Weighting factor for the Gaussian bump emission.
# bump_location      : Radial position of the Gaussian bump.
# sigma_bump_r       : Radial standard deviation (spatial extent) of the Gaussian bump.
# sigma_bump_theta   : Polar standard deviation (angular extent) of the Gaussian bump.
#------------------------------------------------------------------------------------------#
# 4. Localized Hotspot Parameters (Non-axisymmetric Spot)
#------------------------------------------------------------------------------------------#
# blob_intensity     : Emission weighting factor for the localized hotspot.
# blob_r             : Radial coordinate of the hotspot.
# blob_theta         : Polar coordinate of the hotspot.
# blob_phi           : Azimuthal coordinate of the hotspot.
# sigma_blob_r       : Radial extent of the hotspot.
# sigma_blob_theta   : Polar angular extent of the hotspot.
# sigma_blob_phi     : Azimuthal angular extent of the hotspot.
# inspect_azimuth    : The line-of-sight azimuthal angle adopted when inspecting the emission profile.
#-------------------------------------------------------------------------------------------------------------------------------#
disk_intensity=1.0
thick_emit_decay1=-1.5
thick_emit_decay2=-0.5
thickness=0.2
disk_extend=0.1
#------------------------------------------------------------------------------------------#
platform_intensity=1.0
platform_location=5.0
platform_width=3.0
#------------------------------------------------------------------------------------------#
bump_intensity=1.0
bump_location=9.0
sigma_bump_r=1.0
sigma_bump_theta=0.2
#------------------------------------------------------------------------------------------#
blob_intensity=0.0
blob_r=11.0
blob_theta=90.0*np.pi/180.0
blob_phi=160.0*np.pi/180.0
sigma_blob_r=2.0
sigma_blob_theta=10.0*np.pi/180.0
sigma_blob_phi=20.0*np.pi/180.0
inspect_azimuth=270.0*np.pi/180.0
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the kinematic parameters for the accretion environment.
# Note: This section is strictly required when rendering images of both geometrically thick and geometrically thin accretion disks.
# velocity_max_r       : Maximum value of the radial velocity.
# velocity_max_theta   : Maximum value of the polar (theta) velocity.
# velocity_max_phi     : Maximum value of the azimuthal (phi) velocity.
# velocity_decay_r     : Radial acceleration profile factor (must be strictly positive).
# velocity_decay_theta : Polar acceleration profile factor (must be strictly positive).
# velocity_decay_phi   : Azimuthal acceleration profile factor (must be strictly positive).
#-------------------------------------------------------------------------------------------------------------------------------#
velocity_max_r=0.05
velocity_max_theta=0.2
velocity_max_phi=0.4
velocity_decay_r=0.5
velocity_decay_theta=0.5
velocity_decay_phi=0.3
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the emission and ray-tracing parameters for the geometrically thin accretion disk.
# fudge_factor     : Fudge factor for the emission model (typically fixed at 2.0/3.0).
# thin_emit_decay1 : Primary radial decay factor of the thin disk emission (must be strictly negative).
# thin_emit_decay2 : Secondary radial decay factor of the thin disk emission (must be strictly negative).
# max_hit          : Maximum number of times a ray is permitted to intersect the disk.
#                    (Integer type. Generally set to 3; MUST NOT exceed 4.)
#-------------------------------------------------------------------------------------------------------------------------------#
fudge_factor=2.0/3.0
thin_emit_decay1=-2.0
thin_emit_decay2=-0.5
max_hit=3   
#-------------------------------------------------------------------------------------------------------------------------------#
# Please specify the spatial coordinates, finite radius, and the standard deviation of the point source.
#-------------------------------------------------------------------------------------------------------------------------------#
point_source_x=-6.0
point_source_y=0.0
point_source_z=0.0
point_source_radius=0.5
point_source_sigma=0.2
#-------------------------------------------------------------------------------------------------------------------------------#
# Please specify the instantaneous physical velocity vector of the point source.
#-------------------------------------------------------------------------------------------------------------------------------#
point_source_vt=1.0
point_source_vr=0.2
point_source_vtheta=0.0
point_source_vphi=0.5
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the intervals for the longitude and latitude grid lines on the celestial background.
# Note: These intervals must be provided in degrees.
# Additionally, set the "show_line" flag to toggle the display of these grid lines on the four-color image (0 = hide grid lines, 1 = display grid lines).
#-------------------------------------------------------------------------------------------------------------------------------#
show_line=0
delta_latitude=2
delta_longitude=2
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the geometric, emission, and kinematic parameters for the jet source.
# jet_intensity       : Weighting factor for the jet emission intensity.
# jet_base            : Radial coordinate of the jet base.
# jet_range           : Artificially imposed truncation radius of the jet.
# jet_decay           : Radial decay factor of the jet emission (must be strictly positive).
# jet_angle_ini       : Initial opening angle of the jet (at the base).
# jet_angle_inf       : Asymptotic opening angle of the jet (after expansion stabilizes).
# jet_angle_expansion : Characteristic scale for the expansion of the jet opening angle.
# jet_velocity_ini    : Initial radial velocity of the jet.
# jet_velocity_inf    : Asymptotic (saturation) radial velocity of the jet.
#                       (Typically reaches saturation concurrently with the opening angle).
# jet_omega           : Azimuthal velocity (angular frequency) of the particles within the jet.
# jet_omega_decay     : Decay factor for the azimuthal velocity of the jet particles (must be strictly positive).
# gyration_radius     : Minimum gyration radius of the particles within the jet region.
#-------------------------------------------------------------------------------------------------------------------------------#
jet_intensity=1.0
jet_base=2.0
jet_range=200.0
jet_decay=2.0
jet_angle_ini=0.05
jet_angle_inf=0.26
jet_angle_expansion=20.0
#------------------------------------------------------------------------------------------#
jet_velocity_ini=0.1
jet_velocity_inf=0.8
#------------------------------------------------------------------------------------------#
jet_omega=0.0
jet_omega_decay=0.5
gyration_radius=0.01
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the simulation parameters for timelike particles (massive particles).
# Note: The first eight parameters define the initial conditions of the particle.
# step_mp         : Integration step-size for the RK6 solver. 
#                   (Note: This value MUST be strictly negative.)
# motion_time     : Total proper time for the particle's motion.
# max_point       : Maximum number of data points to be recorded for the trajectory (Applicable exclusively to light curve simulations).
# sample          : Sampling frequency for data export (i.e., record trajectory data once every "sample" integration steps).
#                   Applicable exclusively to light curve simulations.
# particle_radius : Radius of the discretized spheres used to render the particle trajectory (Applicable exclusively to light curve simulations).
# condition_type  : Method used to compute the initial conditions of the particle.
#-------------------------------------------------------------------------------------------------------------------------------#
t_mp=0.0
r_mp=15.0
theta_mp=90.0*np.pi/180.0
phi_mp=0.0
pt_mp=-0.96
pr_mp=0.05
ptheta_mp=0.1
pphi_mp=3.55
#------------------------------------------------------------------------------------------#
step_mp=-0.1
motion_time=1000
max_point=500
sample=50
particle_radius=0.25
condition_type=3
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the geometric, emission, and kinematic parameters for the spherically symmetric accretion flow (spherical shell).
# spherically_flow_intensity : Weighting factor for the spherical shell emission.
# spherically_flow_outer     : Outer boundary radius of the spherical shell.
# spherically_flow_inner     : Inner boundary radius of the spherical shell.
# spherically_flow_peak      : Radial position of the peak emission within the shell.
# infalling_flow_max         : Maximum radial velocity of the infalling accretion flow.
# infalling_accelerating     : Radial acceleration profile factor for the infalling flow.
# emission_brake             : Braking coefficient for the infalling flow.
# spherically_flow_sigma     : Standard deviation (radial extent) of the shell emission.
#-------------------------------------------------------------------------------------------------------------------------------#
spherically_flow_intensity=1.0
spherically_flow_outer=10.0
spherically_flow_inner=6.0
spherically_flow_peak=8.0
#------------------------------------------------------------------------------------------#
infalling_flow_max=0.2
infalling_accelerating=0.5
emission_brake=0.5
spherically_flow_sigma=1.0
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the tilt angle of the geometrically thin tilted accretion disk.
# Note: This angular value MUST be provided in radians.
#-------------------------------------------------------------------------------------------------------------------------------#
disk_inclination=30.0*np.pi/180.0
#-------------------------------------------------------------------------------------------------------------------------------#
# Specify the pixel coordinates for simulating a single ray trajectory.
#-------------------------------------------------------------------------------------------------------------------------------#
x_start=4.0
y_start=4.0
#-------------------------------------------------------------------------------------------------------------------------------#
# The remainder of this script is strictly dedicated to evaluating contravariant metric components, metric derivatives, and generating the Fortran model code.
# DO NOT modify any code below this point, EXCEPT for injecting additional spacetime metric parameters.
# 
# Instructions for adding custom metric parameters (e.g., deformation parameters):
# Approximately 110 lines below, locate the following Fortran code block generation:
#
#     fortran_code = f"""module metric_parameters
#     implicit none
#     real*8::spin={spin_value}d0
#     end module metric_parameters
#
# Insert any additional target metric parameters directly beneath the "spin" declaration line. 
# For instance, if the metric incorporates a deformation parameter named 'deform', append the following line:
#     real*8::deform={deform_value}d0
#
# Important: Ensure that the corresponding Python variable (e.g., 'deform_value') has been properly defined and assigned a numerical value earlier in this script.
#-------------------------------------------------------------------------------------------------------------------------------#
covariantmetric = sp.Matrix([
    [g_11, g_12, g_13, g_14],
    [g_21, g_22, g_23, g_24],
    [g_31, g_32, g_33, g_34],
    [g_41, g_42, g_43, g_44]
])
print("------------------------------------------------------------------")
print("Covariant Metric:")
sp.pprint(covariantmetric)
print("------------------------------------------------------------------")
contravariantmetric = covariantmetric.inv()
sp.simplify(contravariantmetric)
print("Contravariant Metric:")
sp.pprint(contravariantmetric)
#-------------------------------------------------------------------------#
g11=contravariantmetric[0, 0]
g12=contravariantmetric[0, 1]
g13=contravariantmetric[0, 2]
g14=contravariantmetric[0, 3]
#-------------------------------------------------------------------------#
g21=contravariantmetric[1, 0]
g22=contravariantmetric[1, 1]
g23=contravariantmetric[1, 2]
g24=contravariantmetric[1, 3]
#-------------------------------------------------------------------------#
g31=contravariantmetric[2, 0]
g32=contravariantmetric[2, 1]
g33=contravariantmetric[2, 2]
g34=contravariantmetric[2, 3]
#-------------------------------------------------------------------------#
g41=contravariantmetric[3, 0]
g42=contravariantmetric[3, 1]
g43=contravariantmetric[3, 2]
g44=contravariantmetric[3, 3]
#-------------------------------------------------------------------------#
g11_dr=sp.diff(g11,r)
g12_dr=sp.diff(g12,r)
g13_dr=sp.diff(g13,r)
g14_dr=sp.diff(g14,r)
#-------------------------------------------------------------------------#
g21_dr=sp.diff(g21,r)
g22_dr=sp.diff(g22,r)
g23_dr=sp.diff(g23,r)
g24_dr=sp.diff(g24,r)
#-------------------------------------------------------------------------#
g31_dr=sp.diff(g31,r)
g32_dr=sp.diff(g32,r)
g33_dr=sp.diff(g33,r)
g34_dr=sp.diff(g34,r)
#-------------------------------------------------------------------------#
g41_dr=sp.diff(g41,r)
g42_dr=sp.diff(g42,r)
g43_dr=sp.diff(g43,r)
g44_dr=sp.diff(g44,r)
#-------------------------------------------------------------------------#
g11_dtheta=sp.diff(g11,theta)
g12_dtheta=sp.diff(g12,theta)
g13_dtheta=sp.diff(g13,theta)
g14_dtheta=sp.diff(g14,theta)
#-------------------------------------------------------------------------#
g21_dtheta=sp.diff(g21,theta)
g22_dtheta=sp.diff(g22,theta)
g23_dtheta=sp.diff(g23,theta)
g24_dtheta=sp.diff(g24,theta)
#-------------------------------------------------------------------------#
g31_dtheta=sp.diff(g31,theta)
g32_dtheta=sp.diff(g32,theta)
g33_dtheta=sp.diff(g33,theta)
g34_dtheta=sp.diff(g34,theta)
#-------------------------------------------------------------------------#
g41_dtheta=sp.diff(g41,theta)
g42_dtheta=sp.diff(g42,theta)
g43_dtheta=sp.diff(g43,theta)
g44_dtheta=sp.diff(g44,theta)
#-------------------------------------------------------------------------#
delta_dr=sp.diff(delta,r)
#-------------------------------------------------------------------------#
sp.pprint(g11_dr)
sp.pprint(g12_dr)
sp.pprint(g13_dr)
sp.pprint(g14_dr)
sp.pprint(g21_dr)
sp.pprint(g22_dr)
sp.pprint(g23_dr)
sp.pprint(g24_dr)
sp.pprint(g31_dr)
sp.pprint(g32_dr)
sp.pprint(g33_dr)
sp.pprint(g34_dr)
sp.pprint(g41_dr)
sp.pprint(g42_dr)
sp.pprint(g43_dr)
sp.pprint(g44_dr)
#-------------------------------------------------------------------------#
sp.pprint(g11_dtheta)
sp.pprint(g12_dtheta)
sp.pprint(g13_dtheta)
sp.pprint(g14_dtheta)
sp.pprint(g21_dtheta)
sp.pprint(g22_dtheta)
sp.pprint(g23_dtheta)
sp.pprint(g24_dtheta)
sp.pprint(g31_dtheta)
sp.pprint(g32_dtheta)
sp.pprint(g33_dtheta)
sp.pprint(g34_dtheta)
sp.pprint(g41_dtheta)
sp.pprint(g42_dtheta)
sp.pprint(g43_dtheta)
sp.pprint(g44_dtheta)
#-------------------------------------------------------------------------#
custom_funcs = {
        'sin': 'dsin',
        'cos': 'dcos', 
        'log': 'dlog'
}
#-------------------------------------------------------------------------#
# The Fortran module for appending additional target metric parameters is located below.
#-------------------------------------------------------------------------#
fortran_code = f"""module metric_parameters
implicit none
real*8::spin={spin_value}d0
end module metric_parameters
    
module task_parameters
implicit none
integer::task_model={task_model}
integer::cpu={cpu}
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
real*8::r_obs={r_obs}d0
real*8::theta_obs={theta_obs}d0
real*8::phi_obs={phi_obs}d0
integer::resolution_x={resolution_x}
integer::resolution_y={resolution_y}
real*8::x_ini={x_ini}d0,x_end={x_end}d0
real*8::y_ini={y_ini}d0,y_end={y_end}d0
real*8::delta_x=0d0
real*8::delta_y=0d0
real*8::vt_obs=1d0
real*8::vr_obs=0d0
real*8::vtheta_obs=0d0
real*8::vphi_obs=0d0
real*8::x_start={x_start}d0
real*8::y_start={y_start}d0
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
real*8::ini_step={ini_step}d0
real*8::step_ratio={step_ratio}d0
real*8::theta_correction={theta_correction}d0
end module raytracing_parameters    
    
module disk_range_parameters
implicit none
real*8::r_inner={r_inner}d0
real*8::r_outer={r_outer}d0
integer::disk_model={disk_model}
end module disk_range_parameters    
    
module emission_scale_parameters
implicit none
real*8::j0={j0}d0
real*8::a0={a0}d0
end module emission_scale_parameters
  
module thick_equatorial_disk_parameters
implicit none
real*8::disk_intensity={disk_intensity}d0
real*8::thick_emit_decay1={thick_emit_decay1}d0
real*8::thick_emit_decay2={thick_emit_decay2}d0
real*8::thickness={thickness}d0
real*8::disk_extend={disk_extend}d0
real*8::platform_intensity={platform_intensity}d0
real*8::platform_location={platform_location}d0
real*8::platform_width={platform_width}d0
real*8::bump_intensity={bump_intensity}d0
real*8::bump_location={bump_location}d0
real*8::sigma_bump_r={sigma_bump_r}d0
real*8::sigma_bump_theta={sigma_bump_theta}d0
real*8::blob_intensity={blob_intensity}d0
real*8::blob_r={blob_r}d0
real*8::blob_theta={blob_theta}d0
real*8::blob_phi={blob_phi}d0
real*8::sigma_blob_r={sigma_blob_r}d0
real*8::sigma_blob_theta={sigma_blob_theta}d0
real*8::sigma_blob_phi={sigma_blob_phi}d0
real*8::inspect_azimuth={inspect_azimuth}d0
end module thick_equatorial_disk_parameters
   
module thin_equatorial_disk_parameters
implicit none
real*8::fudge_factor={fudge_factor}d0
real*8::thin_emit_decay1={thin_emit_decay1}d0
real*8::thin_emit_decay2={thin_emit_decay2}d0
integer::max_hit={max_hit}
end module thin_equatorial_disk_parameters

module equatorial_accretion_dynamics_parameters
implicit none
real*8::velocity_max_r={velocity_max_r}d0
real*8::velocity_max_theta={velocity_max_theta}d0
real*8::velocity_max_phi={velocity_max_phi}d0
real*8::velocity_decay_r={velocity_decay_r}d0
real*8::velocity_decay_theta={velocity_decay_theta}d0
real*8::velocity_decay_phi={velocity_decay_phi}d0
end module equatorial_accretion_dynamics_parameters

module point_source_parameters
implicit none
real*8::point_source_x={point_source_x}d0
real*8::point_source_y={point_source_y}d0
real*8::point_source_z={point_source_z}d0
real*8::point_source_radius={point_source_radius}d0
integer::show_line={show_line}
integer::delta_latitude={delta_latitude}
integer::delta_longitude={delta_longitude}
real*8::point_source_sigma={point_source_sigma}d0
real*8::point_source_vt={point_source_vt}d0
real*8::point_source_vr={point_source_vr}d0
real*8::point_source_vtheta={point_source_vtheta}d0
real*8::point_source_vphi={point_source_vphi}d0
end module point_source_parameters

module jet_parameters
implicit none
real*8::jet_intensity={jet_intensity}d0
real*8::jet_base={jet_base}d0
real*8::jet_decay={jet_decay}d0
real*8::jet_angle_ini={jet_angle_ini}d0
real*8::jet_angle_inf={jet_angle_inf}d0
real*8::jet_angle_expansion={jet_angle_expansion}d0
real*8::jet_velocity_ini={jet_velocity_ini}d0
real*8::jet_velocity_inf={jet_velocity_inf}d0
real*8::jet_range={jet_range}d0
real*8::jet_omega={jet_omega}d0
real*8::jet_omega_decay={jet_omega_decay}d0
real*8::gyration_radius={gyration_radius}d0
end module jet_parameters

module massive_particle_orbit_parameters
implicit none
real*8::t_mp={t_mp}d0
real*8::r_mp={r_mp}d0
real*8::theta_mp={theta_mp}d0
real*8::phi_mp={phi_mp}d0
real*8::pt_mp={pt_mp}d0
real*8::pr_mp={pr_mp}d0
real*8::ptheta_mp={ptheta_mp}d0
real*8::pphi_mp={pphi_mp}d0
real*8::step_mp={step_mp}d0
real*8::particle_radius={particle_radius}d0
real*8::motion_time={motion_time}d0
integer::max_point={max_point}
integer::sample={sample}
integer::condition_type={condition_type}
end module massive_particle_orbit_parameters

module spherically_accretion_flow_parameters
implicit none
real*8::spherically_flow_intensity={spherically_flow_intensity}d0
real*8::spherically_flow_outer={spherically_flow_outer}d0
real*8::spherically_flow_inner={spherically_flow_inner}d0
real*8::spherically_flow_peak={spherically_flow_peak}d0
real*8::infalling_flow_max={infalling_flow_max}d0
real*8::infalling_accelerating={infalling_accelerating}d0
real*8::emission_brake={emission_brake}d0
real*8::spherically_flow_sigma={spherically_flow_sigma}d0
end module spherically_accretion_flow_parameters

module tilted_thin_disk_parameters
implicit none
real*8::disk_inclination={disk_inclination}d0
end module tilted_thin_disk_parameters

subroutine covariant_metric(r,theta,g_ab) 
use metric_parameters
implicit none
real*8,dimension(4,4)::g_ab
real*8 r,theta
g_ab(1,1) = {sp.fcode(g_11, source_format='free', user_functions=custom_funcs)}
g_ab(1,2) = {sp.fcode(g_12, source_format='free', user_functions=custom_funcs)}
g_ab(1,3) = {sp.fcode(g_13, source_format='free', user_functions=custom_funcs)}
g_ab(1,4) = {sp.fcode(g_14, source_format='free', user_functions=custom_funcs)}
g_ab(2,1) = {sp.fcode(g_21, source_format='free', user_functions=custom_funcs)}
g_ab(2,2) = {sp.fcode(g_22, source_format='free', user_functions=custom_funcs)}
g_ab(2,3) = {sp.fcode(g_23, source_format='free', user_functions=custom_funcs)}
g_ab(2,4) = {sp.fcode(g_24, source_format='free', user_functions=custom_funcs)}
g_ab(3,1) = {sp.fcode(g_31, source_format='free', user_functions=custom_funcs)}
g_ab(3,2) = {sp.fcode(g_32, source_format='free', user_functions=custom_funcs)}
g_ab(3,3) = {sp.fcode(g_33, source_format='free', user_functions=custom_funcs)}
g_ab(3,4) = {sp.fcode(g_34, source_format='free', user_functions=custom_funcs)}
g_ab(4,1) = {sp.fcode(g_41, source_format='free', user_functions=custom_funcs)}
g_ab(4,2) = {sp.fcode(g_42, source_format='free', user_functions=custom_funcs)}
g_ab(4,3) = {sp.fcode(g_43, source_format='free', user_functions=custom_funcs)}
g_ab(4,4) = {sp.fcode(g_44, source_format='free', user_functions=custom_funcs)}
return
end subroutine
    
subroutine contravariant_metric(r,theta,gab) 
use metric_parameters
implicit none
real*8,dimension(4,4)::gab
real*8 r,theta
gab(1,1) = {sp.fcode(g11, source_format='free', user_functions=custom_funcs)}
gab(1,2) = {sp.fcode(g12, source_format='free', user_functions=custom_funcs)}
gab(1,3) = {sp.fcode(g13, source_format='free', user_functions=custom_funcs)}
gab(1,4) = {sp.fcode(g14, source_format='free', user_functions=custom_funcs)}
gab(2,1) = {sp.fcode(g21, source_format='free', user_functions=custom_funcs)}
gab(2,2) = {sp.fcode(g22, source_format='free', user_functions=custom_funcs)}
gab(2,3) = {sp.fcode(g23, source_format='free', user_functions=custom_funcs)}
gab(2,4) = {sp.fcode(g24, source_format='free', user_functions=custom_funcs)}
gab(3,1) = {sp.fcode(g31, source_format='free', user_functions=custom_funcs)}
gab(3,2) = {sp.fcode(g32, source_format='free', user_functions=custom_funcs)}
gab(3,3) = {sp.fcode(g33, source_format='free', user_functions=custom_funcs)}
gab(3,4) = {sp.fcode(g34, source_format='free', user_functions=custom_funcs)}
gab(4,1) = {sp.fcode(g41, source_format='free', user_functions=custom_funcs)}
gab(4,2) = {sp.fcode(g42, source_format='free', user_functions=custom_funcs)}
gab(4,3) = {sp.fcode(g43, source_format='free', user_functions=custom_funcs)}
gab(4,4) = {sp.fcode(g44, source_format='free', user_functions=custom_funcs)}
return
end subroutine    
    
subroutine contravariant_metric_diff(r,theta,gab_r,gab_theta)
use metric_parameters
implicit none
real*8,dimension(4,4)::gab_r,gab_theta
real*8 r,theta
gab_r(1,1) = {sp.fcode(g11_dr, source_format='free', user_functions=custom_funcs)}
gab_r(1,2) = {sp.fcode(g12_dr, source_format='free', user_functions=custom_funcs)}
gab_r(1,3) = {sp.fcode(g13_dr, source_format='free', user_functions=custom_funcs)}
gab_r(1,4) = {sp.fcode(g14_dr, source_format='free', user_functions=custom_funcs)}
gab_r(2,1) = {sp.fcode(g21_dr, source_format='free', user_functions=custom_funcs)}
gab_r(2,2) = {sp.fcode(g22_dr, source_format='free', user_functions=custom_funcs)}
gab_r(2,3) = {sp.fcode(g23_dr, source_format='free', user_functions=custom_funcs)}
gab_r(2,4) = {sp.fcode(g24_dr, source_format='free', user_functions=custom_funcs)}
gab_r(3,1) = {sp.fcode(g31_dr, source_format='free', user_functions=custom_funcs)}
gab_r(3,2) = {sp.fcode(g32_dr, source_format='free', user_functions=custom_funcs)}
gab_r(3,3) = {sp.fcode(g33_dr, source_format='free', user_functions=custom_funcs)}
gab_r(3,4) = {sp.fcode(g34_dr, source_format='free', user_functions=custom_funcs)}
gab_r(4,1) = {sp.fcode(g41_dr, source_format='free', user_functions=custom_funcs)}
gab_r(4,2) = {sp.fcode(g42_dr, source_format='free', user_functions=custom_funcs)}
gab_r(4,3) = {sp.fcode(g43_dr, source_format='free', user_functions=custom_funcs)}
gab_r(4,4) = {sp.fcode(g44_dr, source_format='free', user_functions=custom_funcs)}
gab_theta(1,1) = {sp.fcode(g11_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(1,2) = {sp.fcode(g12_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(1,3) = {sp.fcode(g13_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(1,4) = {sp.fcode(g14_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(2,1) = {sp.fcode(g21_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(2,2) = {sp.fcode(g22_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(2,3) = {sp.fcode(g23_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(2,4) = {sp.fcode(g24_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(3,1) = {sp.fcode(g31_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(3,2) = {sp.fcode(g32_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(3,3) = {sp.fcode(g33_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(3,4) = {sp.fcode(g34_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(4,1) = {sp.fcode(g41_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(4,2) = {sp.fcode(g42_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(4,3) = {sp.fcode(g43_dtheta, source_format='free', user_functions=custom_funcs)}
gab_theta(4,4) = {sp.fcode(g44_dtheta, source_format='free', user_functions=custom_funcs)}
return
end subroutine
    
subroutine delta_function(r,theta,delta)
use metric_parameters
implicit none
real*8 r,theta
real*8 delta
delta = {sp.fcode(delta, source_format='free', user_functions=custom_funcs)}
return
end subroutine 
    
subroutine delta_function_diff(r,theta,delta_diff)
use metric_parameters
implicit none
real*8 r,theta
real*8 delta_diff
delta_diff = {sp.fcode(delta_dr, source_format='free', user_functions=custom_funcs)}
return
end subroutine
"""
output_filename = 'model.f90'
with open(output_filename, 'w', encoding='utf-8') as f:
    f.write(fortran_code)
print(fortran_code)