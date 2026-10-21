subroutine photon_initial_values(x_obs,y_obs,vector)
    use metric_parameters
    use observer_parameters
    implicit none
    real*8,dimension(8)::vector
    real*8,dimension(4,4)::g_ab
    real*8 x_obs,y_obs,z_obs
    real*8 coefficient_d,coefficient_w
    real*8 coefficient_a,coefficient_b,coefficient_c
    real*8 x_bh,y_bh,z_bh
    real*8 t,r,theta,phi,p_t,p_r,p_theta,p_phi
    real*8 t_dot,r_dot,theta_dot,phi_dot
    real*8 radius,psi
    real*8 energy
    z_obs=z_start
    !************************************************************************************************************************!
    coefficient_d=(dsqrt(r_obs*r_obs+spin*spin)-z_obs)*dsin(theta_obs)-y_obs*dcos(theta_obs)
    x_bh=coefficient_d*dcos(phi_obs)-x_obs*dsin(phi_obs)
    y_bh=coefficient_d*dsin(phi_obs)+x_obs*dcos(phi_obs)
    z_bh=(r_obs-z_obs)*dcos(theta_obs)+y_obs*dsin(theta_obs)
    !************************************************************************************************************************!
    coefficient_w=x_bh*x_bh+y_bh*y_bh+z_bh*z_bh-spin*spin
    r=dsqrt((coefficient_w+dsqrt(coefficient_w*coefficient_w+4d0*spin*spin*z_bh*z_bh))/2d0)
    theta=dacos(z_bh/r)
    phi=datan2(y_bh,x_bh)
    !************************************************************************************************************************!
    call covariant_metric(r,theta,g_ab)
    radius=dsqrt(r*r+spin*spin)
    psi=phi-phi_obs
    r_dot=(r*radius*dsin(theta)*dsin(theta_obs)*dcos(psi)+radius*radius*dcos(theta)*dcos(theta_obs))/g_ab(3,3)
    theta_dot=(radius*dcos(theta)*dsin(theta_obs)*dcos(psi)-r*dsin(theta)*dcos(theta_obs))/g_ab(3,3)
    phi_dot=-dsin(theta_obs)*dsin(psi)/radius/dsin(theta)
    !************************************************************************************************************************!
    p_r=r_dot*g_ab(2,2)
    p_theta=theta_dot*g_ab(3,3)
    !************************************************************************************************************************!
    coefficient_a=g_ab(1,1)
    coefficient_b=2d0*g_ab(1,4)*phi_dot
    coefficient_c=g_ab(2,2)*r_dot*r_dot+g_ab(3,3)*theta_dot*theta_dot+g_ab(4,4)*phi_dot*phi_dot
    t_dot=(-coefficient_b+dsqrt(coefficient_b*coefficient_b-4d0*coefficient_a*coefficient_c))/2d0/coefficient_a
    !************************************************************************************************************************!
    p_t=g_ab(1,1)*t_dot+g_ab(1,4)*phi_dot
    p_phi=g_ab(4,4)*phi_dot+g_ab(1,4)*t_dot
    !************************************************************************************************************************!
    energy=dabs(p_t)
    p_r=p_r/energy
    p_theta=p_theta/energy
    p_phi=p_phi/energy
    p_t=-1d0
    !************************************************************************************************************************!
    call dimensional(vector,t,r,theta,phi,p_t,p_r,p_theta,p_phi)
    return
    end subroutine
    
subroutine massive_particle_initialization(vector)
    use massive_particle_orbit_parameters
    implicit none
    real*8,dimension(8)::vector
    real*8,dimension(4)::flow_velocity
    real*8,dimension(4,4)::gab,g_ab
    real*8 coefficient_a,coefficient_b,coefficient_c
    real*8 hami
    real*8 r_local_v,theta_local_v,phi_local_v
    call contravariant_metric(r_mp,theta_mp,gab)
    call covariant_metric(r_mp,theta_mp,g_ab)
    select case(condition_type)
    case(1)
        call dimensional(vector,t_mp,r_mp,theta_mp,phi_mp,pt_mp,pr_mp,ptheta_mp,pphi_mp)
    case(2)
        ptheta_mp=0d0
        pr_mp=dsqrt((-1d0-2d0*gab(1,4)*pt_mp*pphi_mp-gab(4,4)*pphi_mp*pphi_mp-gab(1,1)*pt_mp*pt_mp)/gab(2,2))
        call dimensional(vector,t_mp,r_mp,theta_mp,phi_mp,pt_mp,pr_mp,ptheta_mp,pphi_mp)
    case(3)
        pr_mp=0d0
        ptheta_mp=dsqrt((-1d0-2d0*gab(1,4)*pt_mp*pphi_mp-gab(4,4)*pphi_mp*pphi_mp-gab(1,1)*pt_mp*pt_mp)/gab(3,3))
        call dimensional(vector,t_mp,r_mp,theta_mp,phi_mp,pt_mp,pr_mp,ptheta_mp,pphi_mp)
    case(4)
        pr_mp=0d0
        ptheta_mp=0d0
        coefficient_a=gab(1,1)
        coefficient_b=2d0*pphi_mp*gab(1,4)
        coefficient_c=1d0+pphi_mp*pphi_mp*gab(4,4)
        pt_mp=-dabs((-coefficient_b+dsqrt(coefficient_b*coefficient_b-4d0*coefficient_a*coefficient_c))/2d0/coefficient_a)
        call dimensional(vector,t_mp,r_mp,theta_mp,phi_mp,pt_mp,pr_mp,ptheta_mp,pphi_mp)
    case(5)
        call thick_equatorial_disk_velocity(r_mp,theta_mp,phi_mp,flow_velocity,r_local_v,theta_local_v,phi_local_v)
        pt_mp=g_ab(1,1)*flow_velocity(1)+g_ab(1,4)*flow_velocity(4)
        pr_mp=g_ab(2,2)*flow_velocity(2)
        ptheta_mp=g_ab(3,3)*flow_velocity(3)
        pphi_mp=g_ab(4,4)*flow_velocity(4)+g_ab(1,4)*flow_velocity(1)
        call dimensional(vector,t_mp,r_mp,theta_mp,phi_mp,pt_mp,pr_mp,ptheta_mp,pphi_mp)
    end select
    call hamiltonian(vector,hami)
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    write(6,*)"Particle Hamiltonian:",hami
    write(6,*)"Warning: Large deviations from -1/2 violate the Hamiltonian constraint, leading to non-geodesic motion."
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    return
end subroutine

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    