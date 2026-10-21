subroutine anti_dimensional(u,c1,c2,c3,c4,p1,p2,p3,p4)
    implicit none
    real*8,dimension(8)::u
    real*8 c1,c2,c3,c4,p1,p2,p3,p4
    c1=u(1)
    c2=u(2)
    c3=u(3)
    c4=u(4)
    p1=u(5)
    p2=u(6)
    p3=u(7)
    p4=u(8)
    return
    end subroutine
    
subroutine dimensional(u,c1,c2,c3,c4,p1,p2,p3,p4)
    implicit none
    real*8,dimension(8)::u
    real*8 c1,c2,c3,c4,p1,p2,p3,p4
    u(1)=c1
    u(2)=c2
    u(3)=c3
    u(4)=c4
    u(5)=p1
    u(6)=p2
    u(7)=p3
    u(8)=p4
    return
    end subroutine
    
subroutine timestep(t,h)
    implicit none
    real*8 t,h
    t=t+h
    return
    end subroutine
    
subroutine record_values(u,u_0)
    implicit none
    real*8,dimension(8)::u,u_0
    u_0=u
    return
    end subroutine    
    
!subroutine adjust_angle(theta,phi)
!    implicit none
!    real*8 theta,phi
!    real*8 math_pi
!    math_pi=dacos(-1d0)
!    do while(theta.gt.math_pi)
!        theta=theta-2d0*math_pi
!    end do
!    do while(theta.lt.(-1d0)*math_pi)
!        theta=theta+2d0*math_pi
!    end do
!    if(theta.lt.0d0)then
!        theta=-theta
!        phi=phi+math_pi
!    else
!    end if
!    do while(phi.ge.2d0*math_pi)
!        phi=phi-2d0*math_pi
!    end do
!    do while(phi.le.0d0)
!        phi=phi+2d0*math_pi
!    end do
!    return
!    end subroutine 

subroutine adjust_angle(theta,phi)
    implicit none
    real*8 theta,phi
    real*8 math_pi
    math_pi=dacos(-1d0)
    theta=modulo(theta,2d0*math_pi)
    phi=modulo(phi,2d0*math_pi)
    if(theta.gt.math_pi)then
        theta=2d0*math_pi-theta
        phi=phi+math_pi
    end if
    phi=modulo(phi,2d0*math_pi)
    return
    end subroutine     


subroutine change_step(r,theta,step,r_eh)
    use raytracing_parameters
    implicit none
    real*8 step,r,r_eh,theta
    step=ini_step*(r/r_eh)**step_ratio*dabs(dsin(theta))/dsqrt(dsin(theta)*dsin(theta)+theta_correction*theta_correction)
    !step=ini_step*(r/r_eh)**step_ratio*dabs(dsin(theta))
    !step=ini_step*(r/r_eh)**step_ratio
    return
    end subroutine    
    
subroutine spherical_transfer(u,x,y,z)
    use metric_parameters
    implicit none
    real*8,dimension(8)::u
    real*8 x,y,z
    !x=u(2)*dsin(u(3))*dcos(u(4))
    !y=u(2)*dsin(u(3))*dsin(u(4))
    !z=u(2)*dcos(u(3))
    x=dsqrt(u(2)*u(2)+spin*spin)*dsin(u(3))*dcos(u(4))
    y=dsqrt(u(2)*u(2)+spin*spin)*dsin(u(3))*dsin(u(4))
    z=u(2)*dcos(u(3))
    return
    end subroutine    
    
    
    
    
    
    
    
    