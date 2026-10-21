subroutine code_check(r_eh)
    use observer_parameters
    use disk_range_parameters
    use task_parameters
    use metric_parameters
    implicit none
    real*8 r_eh
    delta_x=dabs(x_end-x_ini)/(float(resolution_x)-1d0)
    delta_y=dabs(y_end-y_ini)/(float(resolution_y)-1d0)
    call event_horizon(r_eh)
    select case(disk_model)
    case(1)
        r_inner=r_eh
    case(2)
        r_inner=r_inner
    end select
    if(task_model==3.or.task_model==4)then
        call disk_velocity_check(r_eh)
        call thick_equatorial_disk_radiation(r_eh)
    else
    end if
    if(task_model==7)then
        call jet_velocity_check(r_eh)
        call jet_radiation(r_eh)
    else
    end if
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    write(6,*)"Black Hole Parameters"
    write(6,*)"Black Hole Spin:",spin
    write(6,*)"Radius of the Event Horizon:",r_eh
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    return
    end subroutine
    
    subroutine event_horizon(r_e)
    implicit none
    real*8 r_e
    real*8 r_0,r
    real*8 delta_r,delta_1,delta_0,delta,delta_diff
    real*8::error=1d-15
    real*8::theta=dacos(-1d0)
    integer i,max_n
    i=0
    max_n=1000
    r=15d0
    delta_r=0.1d0
    call delta_function(r,theta,delta_0)
    do while(.true.)
        r=r-delta_r
        call delta_function(r,theta,delta_1)
        if(r.le.0d0.or.isnan(delta_1))then
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            write(6,*)"Warning: Newton iteration failed to converge." 
            write(6,*)"The specified metric parameters might not yield a valid event horizon radius."
            write(6,*)"Please modify the spacetime parameters and retry."
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            stop 
            else
            end if
        if(delta_1*delta_0.le.0d0)exit
        delta_0=delta_1
    end do
    !write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    !write(6,*)"The guess value of the event horizon radius calculated by Newton iteration is",r
    !write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    do while(.true.)
        r_0=r
        call delta_function(r,theta,delta)
        call delta_function_diff(r,theta,delta_diff)
        r=r-delta/delta_diff
        i=i+1
        if(i.gt.max_n)exit
        if(dabs(r-r_0).le.error)exit
    end do
    r_e=dabs(r)
    return
    end subroutine

subroutine shadow_profile(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call shadow_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "
    return
    end subroutine    
 
!subroutine thick_equatorial_disk_image(r_eh)
!    use observer_parameters
!    implicit none
!    real*8 x_screen,y_screen
!    real*8 r_eh
!    integer i,j
!    !$omp parallel do schedule(dynamic,1) default(none) firstprivate(x_ini,y_ini,delta_x,delta_y,resolution_x,resolution_y,x_screen,y_screen,r_eh)
!    do i=1,resolution_y,1
!        y_screen=y_ini+(float(i)-1d0)*delta_y
!        do j=1,resolution_x,1
!            x_screen=x_ini+(float(j)-1d0)*delta_x
!            call imaging_equatorial_thick_task(x_screen,y_screen,r_eh)
!        end do
!    end do
!    !$omp end parallel do
!    return
!    end subroutine  
    
subroutine thick_equatorial_disk_image(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call imaging_equatorial_thick_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "
        return
    end subroutine   
    
subroutine thin_equatorial_disk_image(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call imaging_equatorial_thin_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "
        return
    end subroutine

subroutine four_color_image(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call four_color_image_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "
        return
    end subroutine
    
subroutine jet_image(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call jet_image_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "
        return
    end subroutine
    
subroutine thick_equatorial_disk_radiation(r_eh)    
    use disk_range_parameters
    use physics_parameters
    use thick_equatorial_disk_parameters
    implicit none
    real*8 r_eh
    real*8 ini_r,end_r,delta_r,r
    real*8 ini_theta,delta_theta,theta
    real*8 ini_phi,delta_phi,phi
    real*8 j_v,a_v
    integer i,j
    integer point_r,point_theta,point_phi
    phi=inspect_azimuth
    point_theta=360
    point_phi=720
    point_r=1000
    end_r=50d0
    ini_theta=(10d0/180d0)*math_pi
    ini_phi=0d0
    delta_r=dabs(end_r-r_inner)/(float(point_r)-1d0)
    delta_theta=dabs(math_pi-2d0*ini_theta)/(float(point_theta)-1d0)
    delta_phi=(2d0*math_pi)/(float(point_phi)-1d0)
    ini_r=r_inner
    do i=1,point_r,1
        r=ini_r+delta_r*(float(i)-1d0)
        do j=1,point_theta,1
            theta=ini_theta+delta_theta*(float(j)-1d0)
            call thick_equatorial_disk_emission(r,theta,phi,j_v,a_v)
            write(4,*)r*dsin(theta),r*dcos(theta),j_v
        end do
    end do
    theta=math_pi/2d0
    phi=inspect_azimuth
    do i=1,point_r,1
        r=ini_r+delta_r*(float(i)-1d0)
        call thick_equatorial_disk_emission(r,theta,phi,j_v,a_v)
        write(7,*)r,j_v
    end do
    phi=0d0
    do i=1,point_r,1
        r=ini_r+delta_r*(float(i)-1d0)
        do j=1,point_phi-1,1
            phi=ini_phi+delta_phi*(float(j)-1d0)
            call thick_equatorial_disk_emission(r,theta,phi,j_v,a_v)
            write(5,*)r*dsin(theta)*dcos(phi),r*dsin(theta)*dsin(phi),j_v
        end do
    end do
    return
    end subroutine
   
subroutine jet_radiation(r_eh)
    use physics_parameters
    use jet_parameters
    implicit none
    real*8 r_eh,phi
    real*8 ini_r,end_r,delta_r,r
    real*8 ini_theta,delta_theta,theta
    real*8 j_v,a_v
    integer i,j
    integer point_r,point_theta
    point_theta=100
    point_r=500
    end_r=100d0
    ini_theta=0d0
    delta_r=dabs(end_r-jet_base)/(float(point_r)-1d0)
    delta_theta=dabs(math_pi-ini_theta)/(float(point_theta)-1d0)
    ini_r=jet_base
    do i=1,point_r,1
        r=ini_r+delta_r*(float(i)-1d0)
        do j=1,point_theta,1
            theta=ini_theta+delta_theta*(float(j)-1d0)
            call jet_emission(r,theta,phi,j_v,a_v)
            write(4,*)r*dsin(theta),r*dcos(theta),j_v
            write(4,*)-r*dsin(theta),r*dcos(theta),j_v
        end do
    end do
    return
    end subroutine

subroutine hamiltonian(vector,hami)
    implicit none
    real*8,dimension(8)::vector
    real*8,dimension(4,4)::gab
    real*8 t,r,theta,phi,p_t,p_r,p_theta,p_phi
    real*8 hami
    call anti_dimensional(vector,t,r,theta,phi,p_t,p_r,p_theta,p_phi)
    call contravariant_metric(r,theta,gab)
    hami=(1d0/2d0)*(gab(1,1)*p_t*p_t+gab(2,2)*p_r*p_r+gab(3,3)*p_theta*p_theta+gab(4,4)*p_phi*p_phi+2d0*gab(1,4)*p_t*p_phi)
    return
    end subroutine
    
subroutine difffunction(u,diff_u)
    implicit none
    real*8,dimension(8)::u,diff_u
    real*8,dimension(4,4)::gab,gab_r,gab_theta
    real*8 t,r,theta,phi,p_t,p_r,p_theta,p_phi
    call anti_dimensional(u,t,r,theta,phi,p_t,p_r,p_theta,p_phi)
    call contravariant_metric(r,theta,gab)
    call contravariant_metric_diff(r,theta,gab_r,gab_theta)
    diff_u(1)=-(gab(1,1)*p_t+gab(1,4)*p_phi)
    diff_u(2)=-(gab(2,2)*p_r)
    diff_u(3)=-(gab(3,3)*p_theta)
    diff_u(4)=-(gab(4,4)*p_phi+gab(1,4)*p_t)
    diff_u(5)=0d0
    diff_u(6)=(1d0/2d0)*(gab_r(1,1)*p_t*p_t+gab_r(2,2)*p_r*p_r+gab_r(3,3)*p_theta*p_theta+gab_r(4,4)*p_phi*p_phi+2d0*gab_r(1,4)*p_t*p_phi)
    diff_u(7)=(1d0/2d0)*(gab_theta(1,1)*p_t*p_t+gab_theta(2,2)*p_r*p_r+gab_theta(3,3)*p_theta*p_theta+gab_theta(4,4)*p_phi*p_phi+2d0*gab_theta(1,4)*p_t*p_phi)
    diff_u(8)=0d0
    return
    end subroutine       
    
subroutine photon_path(r_eh)
    use observer_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(8)::vector
    real*8 r_eh
    real*8 step
    real*8 x_c,y_c,z_c
    real*8 temp_hami
    integer integration_step
    integration_step=0
    call photon_initial_values(x_start,y_start,vector)
    do while(.true.)
        call change_step(vector(2),vector(3),step,r_eh)
        integration_step=integration_step+1
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))then
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            write(6,*)"Photon captured by the event horizon."
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            exit
        else
        end if
        if(vector(2).ge.observer)then
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            write(6,*)"Photon escaped to infinity."
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            exit
        else
        end if
        call hamiltonian(vector,temp_hami)
        call spherical_transfer(vector,x_c,y_c,z_c)
        write(14,*)x_c,y_c,z_c
        write(15,*)integration_step,dlog10(dabs(temp_hami)+1d-16)
    end do
    return
end subroutine

subroutine massive_particle_path(r_eh)
    use raytracing_parameters
    use massive_particle_orbit_parameters
    real*8,dimension(8)::vector
    real*8 r_eh
    real*8 x_c,y_c,z_c
    real*8 temp_hami
    real*8 t
    integer integration_step
    integration_step=0
    t=0d0
    call massive_particle_initialization(vector)
    do while(.true.)
        integration_step=integration_step+1
        call timestep(t,-step_mp)
        call rksixth(vector,step_mp)
        if(vector(2).le.(r_eh+0.1d0).or.isnan(vector(2)))then
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            write(6,*)"Particle captured by the event horizon."
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            exit
        else
        end if
        if(vector(2).ge.500d0)then
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            write(6,*)"Particle escaped to infinity."
            write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
            exit
        else
        end if
        !if(integration_step.gt.max_point)exit
        if(t.gt.motion_time)exit
        call hamiltonian(vector,temp_hami)
        call spherical_transfer(vector,x_c,y_c,z_c)
        write(14,*)x_c,y_c,z_c
        write(15,*)integration_step,dlog10(dabs(temp_hami+0.5d0)+1d-16)
    end do
    return
end subroutine

subroutine point_source_image(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call point_source_image_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "
        return
    end subroutine
    
 subroutine spherically_accretion_flow_image(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call spherically_accretion_image_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "
        return
    end subroutine
    
subroutine disk_velocity_check(r_eh)
    use physics_parameters
    use raytracing_parameters
    use disk_range_parameters
    implicit none   
    real*8,dimension(4)::flow_velocity
    real*8 r_local_v,theta_local_v,phi_local_v
    real*8 r_eh
    real*8 ini_r,end_r,delta_r,r
    real*8 ini_theta,end_theta,delta_theta,theta
    real*8 phi
    integer i,j
    integer point_r,point_theta
    ini_r=r_eh+hit_error
    end_r=r_outer
    point_r=1d4
    ini_theta=0.005d0
    end_theta=math_pi-0.005d0
    point_theta=1000
    delta_r=dabs(end_r-ini_r)/(float(point_r)-1d0)
    delta_theta=dabs(end_theta-ini_theta)/(float(point_theta)-1d0)
    do i=1,point_r,1
        r=ini_r+delta_r*(float(i)-1d0)
        do j=1,point_theta,1
            theta=ini_theta+delta_theta*(float(j)-1d0)
            call thick_equatorial_disk_velocity(r,theta,phi,flow_velocity,r_local_v,theta_local_v,phi_local_v)
            if(isnan(flow_velocity(1)).or.isnan(flow_velocity(2)).or.isnan(flow_velocity(3)).or.isnan(flow_velocity(4)))then
                write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
                write(6,*)"Warning: Unphysical dynamical parameters detected for the accretion matter."
                write(6,*)"Current settings could lead to superluminal motion."
                write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
                stop
            else
            end if
        end do
    end do
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    write(6,*)"Sanity check passed: The dynamical parameters are physically valid."
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    return
    end subroutine

    subroutine jet_velocity_check(r_eh)
        use physics_parameters
        use raytracing_parameters
        use jet_parameters
        implicit none   
        real*8,dimension(4)::jet_velocity
        real*8 r_local_v,theta_local_v,phi_local_v
        real*8 r_eh
        real*8 ini_r,end_r,delta_r,r
        real*8 ini_theta,end_theta,delta_theta,theta
        real*8 phi
        integer i,j
        integer point_r,point_theta
        phi=0d0
        if(jet_base.le.r_eh)then
            ini_r=r_eh+hit_error
        else
            ini_r=jet_base+hit_error
        end if
        end_r=jet_range
        point_r=1d4
        ini_theta=0.005d0
        end_theta=math_pi/2d0-0.005d0
        point_theta=1000
        delta_r=dabs(end_r-ini_r)/(float(point_r)-1d0)
        delta_theta=dabs(end_theta-ini_theta)/(float(point_theta)-1d0)
        do i=1,point_r,1
            r=ini_r+delta_r*(float(i)-1d0)
            do j=1,point_theta,1
                theta=ini_theta+delta_theta*(float(j)-1d0)
                call jet_dynamic(r,theta,phi,jet_velocity,r_local_v)
                if(isnan(jet_velocity(1)).or.isnan(jet_velocity(2)).or.isnan(jet_velocity(3)).or.isnan(jet_velocity(4)))then
                    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
                    write(6,*)"Warning: Unphysical dynamical parameters detected for the jets."
                    write(6,*)"Current settings could lead to superluminal motion."
                    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
                    stop
                else
                end if
            end do
        end do
        write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
        write(6,*)"Sanity check passed: The dynamical parameters are physically valid."
        write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
        return
        end subroutine

subroutine grid_image(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call grid_image_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "
        return
    end subroutine    
    
subroutine tilted_disk_four_color_image(r_eh)
    use observer_parameters
    use omp_lib
    implicit none
    real*8 x_screen,y_screen
    real*8 r_eh
    real*8 frac
    integer i,j
    integer task_total,task_done,task_local_done,task_temp
    task_total=resolution_x*resolution_y
    task_done=0
    !$omp parallel default(none) &
    !$omp& private(i,j,x_screen,y_screen,task_local_done,task_temp,frac) &
    !$omp& shared(task_done) &
    !$omp& firstprivate(x_ini,y_ini,delta_x,delta_y, &
    !$omp&              resolution_x,resolution_y,r_eh,task_total)
    task_local_done=0
    !$omp do schedule(dynamic, 1)
        do i=1,resolution_y,1
            y_screen=y_ini+(float(i)-1d0)*delta_y
            do j=1,resolution_x,1
                x_screen=x_ini+(float(j)-1d0)*delta_x
                call tilted_disk_four_color_image_task(x_screen,y_screen,r_eh)
                task_local_done=task_local_done+1
                if(task_local_done.ge.task_total/200) then
                    !$omp atomic
                    task_done=task_done+task_local_done
                    !$omp end atomic
                    task_local_done=0
                    !$omp critical (progress)
                    !$omp atomic read
                    task_temp=task_done
                    !$omp end atomic
                    frac=100d0*task_temp/float(task_total)
                    write(6,'(A,F6.2,A)',advance='no') &
                        char(13)//"Progress: ", frac, " %   "
                    flush(6)
                    !$omp end critical (progress)
                else
                end if
            end do
        end do
        !$omp end do
        if(task_local_done.gt.0)then
            !$omp atomic
            task_done=task_done+task_local_done
            !$omp end atomic
        else
        end if
        !$omp end parallel
        write(6,'(A)')char(13)//"Progress: 100.00 %   Done.     "    
    return
end subroutine