subroutine rkfourth(u,h)
implicit none
real*8,dimension(8)::u,u0,diff_u
real*8,dimension(8)::k1,k2,k3,k4
real*8 h 
u0=u
call difffunction(u,diff_u)
k1=h*diff_u
u=u0+k1/2d0
call difffunction(u,diff_u)
k2=h*diff_u
u=u0+k2/2d0
call difffunction(u,diff_u)
k3=h*diff_u
u=u0+k3
call difffunction(u,diff_u)
k4=h*diff_u
u=u0+(1d0/6d0)*(k1+2d0*k2+2d0*k3+k4)
return
end subroutine 

subroutine rkfifth(u,h)
implicit none
real*8,dimension(8)::u,u0,diff_u
real*8,dimension(8)::k1,k2,k3,k4,k5,k6
real*8 h 
u0=u
call difffunction(u,diff_u)
k1=h*diff_u
u=u0+k1/4d0
call difffunction(u,diff_u)
k2=h*diff_u
u=u0+3d0*k1/32d0+9d0*k2/32d0
call difffunction(u,diff_u)
k3=h*diff_u
u=u0+1932d0*k1/2197d0-7200d0*k2/2197d0+7296d0*k3/2197d0
call difffunction(u,diff_u)
k4=h*diff_u
u=u0+439d0*k1/216d0-8d0*k2+3680d0*k3/513d0-845d0*k4/4104d0
call difffunction(u,diff_u)
k5=h*diff_u
u=u0-8d0*k1/27d0+2d0*k2-3544d0*k3/2565d0+1859d0*k4/4104d0-11d0*k5/40d0
call difffunction(u,diff_u)
k6=h*diff_u
u=u0+16d0*k1/135d0+6656d0*k3/12825d0+28561d0*k4/56430d0-9d0*k5/50d0+2d0*k6/55d0
return
end subroutine    
    
subroutine rksixth(u,h)
implicit none
real*8,dimension(8)::u,u0,diff_u
real*8,dimension(8)::k1,k2,k3,k4,k5,k6,k7
real*8 h 
real*8 c
c=dsqrt(21d0)
u0=u
call difffunction(u,diff_u)
k1=h*diff_u
u=u0+k1
call difffunction(u,diff_u)
k2=h*diff_u
u=u0+(1d0/8d0)*(3d0*k1+k2)
call difffunction(u,diff_u)
k3=h*diff_u
u=u0+(1d0/27d0)*(8d0*k1+2d0*k2+8d0*k3)
call difffunction(u,diff_u)
k4=h*diff_u
u=u0+(1d0/392d0)*(3d0*(3d0*c-7d0)*k1-8d0*(7d0-c)*k2+48d0*(7d0-c)*k3-3d0*(21d0-c)*k4)
call difffunction(u,diff_u)
k5=h*diff_u
u=u0+(1d0/1960d0)*(-5d0*(231d0+51d0*c)*k1-40d0*(7d0+c)*k2-320d0*c*k3+3d0*(21d0+121d0*c)*k4+392d0*(6d0+c)*k5)
call difffunction(u,diff_u)
k6=h*diff_u
u=u0+(1d0/180d0)*(15d0*(22d0+7d0*c)*k1+120d0*k2+40d0*(7d0*c-5d0)*k3-63d0*(3d0*c-2d0)*k4-14d0*(49d0+9d0*c)*k5+70d0*(7d0-c)*k6)
call difffunction(u,diff_u)
k7=h*diff_u
u=u0+(1d0/180d0)*(9d0*k1+64d0*k3+49d0*k5+49d0*k6+9d0*k7)
return
end subroutine
    
subroutine rkf(u,step)
use integration_parameters
implicit none
real*8,dimension(8)::u,u_1,u_2,eps
real*8 max_eps
real*8 step
integer i
u_1=0d0
u_2=0d0
eps=0d0
do while(.true.)
    u_1=u
    u_2=u
    !call rkfourth(u_1,step)
    !call rkfifth(u_2,step)
    call rkfifth(u_1,step)
    call rksixth(u_2,step)
    do i=1,8,1
        eps(i)=dabs(u_1(i)-u_2(i))/u_2(i)
    end do
    max_eps=maxval(eps)
    if(max_eps.gt.solution_error)then
        step=0.5d0*step
    else
        u=u_1
        exit
    end if
end do
return
end subroutine