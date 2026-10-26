subroutine thick_equatorial_disk_emission(r,theta,phi,j_local,a_local)
    use disk_range_parameters
    use thick_equatorial_disk_parameters
    use physics_parameters
    use emission_scale_parameters
    implicit none
    real*8 r,theta,phi
    real*8 j_local,a_local
    real*8 j_base,j_angle,j_bump,j_blob
    real*8 r_eff
    r_eff=r-platform_intensity*platform_width*(dtanh((r-platform_location)/platform_width)-dtanh((r_inner-platform_location)/platform_width))
    j_angle=dexp(-(theta-math_pi/2d0)**2d0/(2d0*(thickness*dexp(disk_extend*r))**2d0))
    j_base=disk_intensity*j_angle*dexp(thick_emit_decay1*dlog(r_eff/r_inner)+thick_emit_decay2*(dlog(r_eff/r_inner))**2d0)
    j_bump=bump_intensity*dexp(-(r-bump_location)**2d0/2d0/sigma_bump_r/sigma_bump_r-(theta-math_pi/2d0)**2d0/2d0/sigma_bump_theta/sigma_bump_theta)
    j_blob=blob_intensity*dexp(-(r-blob_r)**2d0/2d0/sigma_blob_r/sigma_blob_r-(theta-blob_theta)**2d0/2d0/sigma_blob_theta/sigma_blob_theta-(1d0-dcos(phi-blob_phi))/2d0/sigma_blob_phi/sigma_blob_phi)
    j_local=j0*(j_base+j_bump+j_blob)
    a_local=a0*(j_base+j_bump+j_blob)
    return
    end subroutine
    
subroutine thick_equatorial_disk_velocity(r,theta,phi,flow_velocity,r_local_velocity,theta_local_velocity,phi_local_velocity)
    use disk_range_parameters
    use equatorial_accretion_dynamics_parameters
    use physics_parameters
    use metric_parameters
    implicit none
    real*8 r,theta,phi
    real*8,dimension(4,4)::g_ab
    real*8,dimension(4)::flow_velocity
    real*8 r_local_velocity,theta_local_velocity,phi_local_velocity
    real*8 lapse_alpha,lorentz_gamma,drag_omega
    call covariant_metric(r,theta,g_ab)
    r_local_velocity=-velocity_max_r*(r_inner/r)**velocity_decay_r
    phi_local_velocity=velocity_max_phi*(r_inner/r)**velocity_decay_phi
    theta_local_velocity=-(theta-math_pi/2d0)*velocity_max_theta*(r_inner/r)**velocity_decay_theta
    lapse_alpha=dsqrt(-g_ab(1,1)+g_ab(1,4)*g_ab(1,4)/g_ab(4,4))
    drag_omega=-g_ab(1,4)/g_ab(4,4)
    lorentz_gamma=1d0/dsqrt(1d0-r_local_velocity*r_local_velocity-phi_local_velocity*phi_local_velocity-theta_local_velocity*theta_local_velocity)
    flow_velocity(1)=lorentz_gamma/lapse_alpha
    flow_velocity(2)=lorentz_gamma*r_local_velocity/dsqrt(g_ab(2,2))
    flow_velocity(3)=lorentz_gamma*theta_local_velocity/dsqrt(g_ab(3,3))
    flow_velocity(4)=lorentz_gamma*(drag_omega/lapse_alpha+phi_local_velocity/dsqrt(g_ab(4,4)))
    return
    end subroutine 
    
subroutine thin_equatorial_disk_emission(r,j_local,r_eh)
    use emission_scale_parameters
    use thin_equatorial_disk_parameters
    implicit none
    real*8 z
    real*8 r,j_local,r_eh
    z=dlog(r/r_eh)
    j_local=j0*dexp(thin_emit_decay1*z+thin_emit_decay2*z*z)
    return 
    end subroutine
    
subroutine jet_emission(r,theta,phi,j_local,a_local)
    use jet_parameters
    use physics_parameters
    use emission_scale_parameters
    implicit none
    real*8 r,theta,phi,theta_pole
    real*8 j_local,a_local
    real*8 jet_radius_emission,jet_angle_emission,jet_core_emission,jet_core_intensity
    real*8 sigma_function
    if(r.lt.jet_base)then
        j_local=0d0
        a_local=0d0
    else
        theta_pole=min(theta,math_pi-theta)
        jet_radius_emission=(jet_base/r)**(jet_decay)
        sigma_function=jet_angle_inf-(jet_angle_inf-jet_angle_ini)*dexp(-r/jet_angle_expansion)
        jet_angle_emission=dexp(-theta_pole*theta_pole/2d0/sigma_function/sigma_function)
        !jet_core_emission=dexp(-(r-jet_core_location)**2d0/2d0/sigma_core_r/sigma_core_r-theta*theta/2d0/sigma_core_theta/sigma_core_theta)
        jet_core_emission=0d0
        jet_core_intensity=0d0
        j_local=j0*(jet_intensity*jet_radius_emission*jet_angle_emission+jet_core_intensity*jet_core_emission)
        a_local=a0*(jet_intensity*jet_radius_emission*jet_angle_emission+jet_core_intensity*jet_core_emission)
    end if
    return
    end subroutine    
    
subroutine jet_dynamic(r,theta,phi,jet_velocity,r_local_velocity)   
    use jet_parameters
    use physics_parameters
    implicit none
    real*8 r,theta,phi,theta_pole
    real*8 r_local_velocity,theta_local_velocity,phi_local_velocity
    real*8 lapse_alpha,lorentz_gamma,drag_omega,sigma_function,jet_omega_function,gyration_radius_function
    real*8,dimension(4)::jet_velocity
    real*8,dimension(4,4)::g_ab
    call covariant_metric(r,theta,g_ab)
    theta_pole=min(theta,math_pi-theta)
    sigma_function=jet_angle_inf-(jet_angle_inf-jet_angle_ini)*dexp(-r/jet_angle_expansion)
    r_local_velocity=(jet_velocity_inf-(jet_velocity_inf-jet_velocity_ini)*dexp(-(r-jet_base)/jet_angle_expansion))*dexp(-theta_pole*theta_pole/2d0/sigma_function/sigma_function)
    jet_omega_function=jet_omega*(jet_base/r)**jet_omega_decay
    gyration_radius_function=gyration_radius*sigma_function/jet_angle_ini
    phi_local_velocity=sign(1d0,dcos(theta))*gyration_radius_function*jet_omega_function*dexp(-theta_pole*theta_pole/2d0/sigma_function/sigma_function)
    theta_local_velocity=0d0
    lapse_alpha=dsqrt(-g_ab(1,1)+g_ab(1,4)*g_ab(1,4)/g_ab(4,4))
    drag_omega=-g_ab(1,4)/g_ab(4,4)
    lorentz_gamma=1d0/dsqrt(1d0-r_local_velocity*r_local_velocity-phi_local_velocity*phi_local_velocity-theta_local_velocity*theta_local_velocity)
    jet_velocity(1)=lorentz_gamma/lapse_alpha
    jet_velocity(2)=lorentz_gamma*r_local_velocity/dsqrt(g_ab(2,2))
    jet_velocity(3)=lorentz_gamma*theta_local_velocity/dsqrt(g_ab(3,3))
    jet_velocity(4)=lorentz_gamma*(drag_omega/lapse_alpha+phi_local_velocity/dsqrt(g_ab(4,4)))
    !jet_velocity(1)=dsqrt((-1d0-g_ab(2,2)*jet_velocity(2)*jet_velocity(2))/g_ab(1,1))
    !if(isnan(jet_velocity(1)))then
    !    write(*,*)(-1d0-g_ab(2,2)*jet_velocity(2)*jet_velocity(2))/g_ab(1,1)
    !    write(*,*)g_ab(1,1)
    !    write(*,*)g_ab(2,2)
    !    pause
    !else
    !end if
    return
end subroutine
    
subroutine point_source_emission(distance,j_local,a_local)
    use point_source_parameters
    use emission_scale_parameters
    implicit none
    real*8 distance,j_local,a_local
    j_local=j0*dexp(-distance*distance/2d0/point_source_sigma/point_source_sigma)
    a_local=j_local*a0/j0
    !a_local=a0*(1d0-distance/point_source_radius)**2d0
    return
end subroutine

subroutine point_source_dynamic(r,theta,phi,flow_velocity,r_local_velocity,theta_local_velocity,phi_local_velocity)
    use point_source_parameters
    implicit none
    real*8 r,theta,phi,r_local_velocity,theta_local_velocity,phi_local_velocity
    real*8,dimension(4)::flow_velocity
    real*8,dimension(4,4)::g_ab
    real*8 lapse_alpha,lorentz_gamma,drag_omega
    call covariant_metric(r,theta,g_ab)
    r_local_velocity=point_source_vr
    theta_local_velocity=point_source_vtheta
    phi_local_velocity=point_source_vphi
    lapse_alpha=dsqrt(-g_ab(1,1)+g_ab(1,4)*g_ab(1,4)/g_ab(4,4))
    drag_omega=-g_ab(1,4)/g_ab(4,4)
    lorentz_gamma=1d0/dsqrt(1d0-r_local_velocity*r_local_velocity-phi_local_velocity*phi_local_velocity-theta_local_velocity*theta_local_velocity)
    flow_velocity(1)=lorentz_gamma/lapse_alpha
    flow_velocity(2)=lorentz_gamma*r_local_velocity/dsqrt(g_ab(2,2))
    flow_velocity(3)=lorentz_gamma*theta_local_velocity/dsqrt(g_ab(3,3))
    flow_velocity(4)=lorentz_gamma*(drag_omega/lapse_alpha+phi_local_velocity/dsqrt(g_ab(4,4)))
    return
end subroutine


subroutine spherically_accretion_flow_velocity(r,theta,phi,flow_velocity,r_local_velocity,r_eh)
    use spherically_accretion_flow_parameters
    implicit none
    real*8 r,theta,phi,r_local_velocity,theta_local_velocity,phi_local_velocity,r_eh
    real*8 lapse_alpha,lorentz_gamma,drag_omega
    real*8,dimension(4)::flow_velocity
    real*8,dimension(4,4)::g_ab
    call covariant_metric(r,theta,g_ab)
    r_local_velocity=-infalling_flow_max*(r_eh/r)**infalling_accelerating*dexp(1d0-emission_brake*dexp(-(r-spherically_flow_peak)**2d0/2d0/spherically_flow_sigma/spherically_flow_sigma))
    phi_local_velocity=0d0
    theta_local_velocity=0d0
    lapse_alpha=dsqrt(-g_ab(1,1)+g_ab(1,4)*g_ab(1,4)/g_ab(4,4))
    drag_omega=-g_ab(1,4)/g_ab(4,4)
    lorentz_gamma=1d0/dsqrt(1d0-r_local_velocity*r_local_velocity-phi_local_velocity*phi_local_velocity-theta_local_velocity*theta_local_velocity)
    flow_velocity(1)=lorentz_gamma/lapse_alpha
    flow_velocity(2)=lorentz_gamma*r_local_velocity/dsqrt(g_ab(2,2))
    flow_velocity(3)=lorentz_gamma*theta_local_velocity/dsqrt(g_ab(3,3))
    flow_velocity(4)=lorentz_gamma*(drag_omega/lapse_alpha+phi_local_velocity/dsqrt(g_ab(4,4)))
    !flow_velocity(1)=dsqrt((-1d0-g_ab(2,2)*flow_velocity(2)*flow_velocity(2))/g_ab(1,1))
    return
end subroutine

subroutine spherically_accretion_flow_emission(r,j_local,a_local)  
    use spherically_accretion_flow_parameters
    use emission_scale_parameters
    implicit none
    real*8 r,j_local,a_local
    j_local=j0*spherically_flow_intensity*dexp(-(r-spherically_flow_peak)**2d0/2d0/spherically_flow_sigma/spherically_flow_sigma)
    a_local=j_local*a0/j0
    return
end subroutine




     
    

    
