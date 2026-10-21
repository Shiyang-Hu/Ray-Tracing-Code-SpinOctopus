subroutine imaging_equatorial_thick_task(x_screen,y_screen,r_eh)
    use disk_range_parameters
    use observer_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(8)::vector
    real*8,dimension(4)::flow_velocity,vector_obs
    real*8 r_local_v,theta_local_v,phi_local_v
    real*8 x_screen,y_screen,r_eh
    real*8 redshift_factor
    real*8 deepth,i_obs
    real*8 step
    real*8 j_local,a_local
    deepth=0d0
    i_obs=0d0
    call photon_initial_values(x_screen,y_screen,vector)
    vector_obs(1)=vector(5)
    vector_obs(2)=vector(6)
    vector_obs(3)=vector(7)
    vector_obs(4)=vector(8)
    do while(.true.)
        call change_step(vector(2),vector(3),step,r_eh)
        call rkf(vector,step)
        !call rkfifth(vector,step)
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))exit
        if(vector(2).ge.observer)exit
        if(vector(2).le.r_outer)then
            call thick_equatorial_disk_velocity(vector(2),vector(3),vector(4),flow_velocity,r_local_v,theta_local_v,phi_local_v)
            redshift_factor=(vector_obs(1)*vt_obs+vector_obs(2)*vr_obs+vector_obs(3)*vtheta_obs+vector_obs(4)*vphi_obs)/(vector(5)*flow_velocity(1)+vector(6)*flow_velocity(2)+vector(7)*flow_velocity(3)+vector(8)*flow_velocity(4))
            call thick_equatorial_disk_emission(vector(2),vector(3),vector(4),j_local,a_local)
            deepth=deepth+a_local*step/redshift_factor
            i_obs=i_obs+(redshift_factor*redshift_factor*j_local*dexp(-deepth))*step
        else
        end if
    end do
    write(2,*)x_screen,y_screen,i_obs
    return
    end subroutine
    
subroutine shadow_task(x_screen,y_screen,r_eh)
    use raytracing_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh
    real*8,dimension(8)::vector
    real*8 step 
    call photon_initial_values(x_screen,y_screen,vector)
        do while(.true.)
            call change_step(vector(2),vector(3),step,r_eh)
            call rkf(vector,step)
            if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))then
                write(1,*)x_screen,y_screen
                exit
            else
            end if
            if(vector(2).ge.observer)exit
        end do
    return
    end subroutine      
    
subroutine disk_velocity_field(r_eh)
    use disk_range_parameters
    implicit none
    real*8,dimension(4)::flow_velocity
    real*8 r_local_v,theta_local_v,phi_local_v
    real*8 r_eh
    real*8 ini_r,end_r,delta_r,r
    real*8 phi,delta_phi
    real*8 x_velocity,y_velocity,velocity_total
    real*8 angle_direction
    integer point_r,point_phi
    integer i,j
    point_r=10
    point_phi=36
    ini_r=r_inner+0.1d0
    end_r=20d0
    delta_r=dabs(ini_r-end_r)/(float(point_r)-1d0)
    delta_phi=dabs(dacos(-1d0)*2d0)/(float(point_phi)-1d0)
    do i=2,point_r,1
        r=ini_r+(float(i)-1d0)*delta_r
        do j=1,point_phi-1,1
            phi=(float(j)-1d0)*delta_phi
            call thick_equatorial_disk_velocity(r,dacos(-1d0)/2d0,phi,flow_velocity,r_local_v,theta_local_v,phi_local_v)
            x_velocity=r_local_v*dcos(phi)-phi_local_v*dsin(phi)
            y_velocity=r_local_v*dsin(phi)+phi_local_v*dcos(phi)
            velocity_total=dsqrt(x_velocity*x_velocity+y_velocity*y_velocity)
            angle_direction=datan2(y_velocity,x_velocity)
            write(3,'(4(F12.6,1X))')r*dcos(phi),r*dsin(phi),angle_direction,velocity_total
        end do
    end do
    return
    end subroutine
            
subroutine imaging_equatorial_thin_task(x_screen,y_screen,r_eh)
    use disk_range_parameters
    use observer_parameters
    use raytracing_parameters
    use thin_equatorial_disk_parameters
    implicit none
    real*8,dimension(8)::vector,vector0
    real*8,dimension(4)::flow_velocity,redshift_factor,intensity,vector_obs
    real*8 r_local_v,theta_local_v,phi_local_v
    real*8 x_screen,y_screen,r_eh
    real*8 intensity_obs,rsource
    real*8 step
    real*8 j_local
    integer hit,file_name,i
    hit=0
    call photon_initial_values(x_screen,y_screen,vector)
    vector_obs(1)=vector(5)
    vector_obs(2)=vector(6)
    vector_obs(3)=vector(7)
    vector_obs(4)=vector(8)
    do while(.true.)
        call record_values(vector,vector0)
        call change_step(vector(2),vector(3),step,r_eh)
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))exit
        if(vector(2).ge.observer)exit
        if(dcos(vector(3))*dcos(vector0(3)).le.0d0.and.vector(2).le.r_outer.and.vector(2).gt.r_inner)then
            hit=hit+1
            rsource=(vector(2)+vector0(2))/2d0
            call thick_equatorial_disk_velocity(rsource,dacos(-1d0)/2d0,vector(4),flow_velocity,r_local_v,theta_local_v,phi_local_v)
            redshift_factor(hit)=(vector_obs(1)*vt_obs+vector_obs(2)*vr_obs+vector_obs(3)*vtheta_obs+vector_obs(4)*vphi_obs)/(vector(5)*flow_velocity(1)+vector(6)*flow_velocity(2)+vector(7)*flow_velocity(3)+vector(8)*flow_velocity(4))
            call thin_equatorial_disk_emission(rsource,j_local,r_eh) 
            intensity(hit)=j_local*redshift_factor(hit)**3d0
        else
        end if
        if(hit.ge.max_hit)exit
    end do
    if(hit==0)then
        write(1,*)x_screen,y_screen
    else
        do i=1,hit,1
            file_name=i+7
            write(file_name,*)x_screen,y_screen,redshift_factor(i)
        end do
    end if
    select case(hit)
    case(0)
        intensity_obs=0d0
    case(1)
        intensity_obs=intensity(1)
    case(2)
        intensity_obs=intensity(1)+fudge_factor*intensity(2)
    case(3)
        intensity_obs=intensity(1)+fudge_factor*intensity(2)+fudge_factor*intensity(3)
    case(4)
        intensity_obs=intensity(1)+fudge_factor*intensity(2)+fudge_factor*intensity(3)+fudge_factor*intensity(4)
    end select 
    write(2,*)x_screen,y_screen,intensity_obs
    write(12,*)x_screen,y_screen,hit
    return
    end subroutine

subroutine four_color_image_task(x_screen,y_screen,r_eh)
    use point_source_parameters
    use observer_parameters
    use raytracing_parameters
    use physics_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh,step
    real*8,dimension(8)::vector,vector0
    real*8 x,y,z
    real*8 x0,y0,z0
    real*8 radius,radius0
    real*8 latitude,longitude
    integer i,j
    integer point_latitude,point_longitude
    point_longitude=(360-delta_longitude)/delta_longitude+1
    point_latitude=(180-delta_latitude-delta_latitude)/delta_latitude+1
    call photon_initial_values(x_screen,y_screen,vector)
    do while(.true.)
        call change_step(vector(2),vector(3),step,r_eh)
        call record_values(vector,vector0)
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))then
            write(13,*)x_screen,y_screen,0
            exit
        else
        end if
        call spherical_transfer(vector,x,y,z)
        call spherical_transfer(vector0,x0,y0,z0)
        radius=dsqrt((x-point_source_x)**2d0+(y-point_source_y)**2d0+(z-point_source_z)**2d0)   
        radius0=dsqrt((x0-point_source_x)**2d0+(y0-point_source_y)**2d0+(z0-point_source_z)**2d0)
        if(radius.le.point_source_radius.and.radius0.ge.point_source_radius)then
            write(13,*)x_screen,y_screen,5
            exit
        else
        end if
        if(vector(2).ge.observer)then
            call adjust_angle(vector(3),vector(4))
            if(show_line==1)then
                do j=1,point_longitude,1
                    longitude=float(delta_longitude)*math_pi/180d0*(float(j)-1d0)
                    if(dabs(vector(4)-longitude)*dsin(vector(3)).le.0.002)then
                        write(13,*)x_screen,y_screen,0
                        go to 1
                    else
                    end if
                end do
                do i=1,point_latitude,1
                    latitude=delta_latitude*math_pi/180d0+delta_latitude*math_pi/180d0*(float(i)-1d0)
                    if(dabs((vector(3)-latitude)/dsin(vector(3))).le.0.002)then
                        write(13,*)x_screen,y_screen,0
                        go to 1
                    else
                    end if
                end do
            else
            end if
            if(z.gt.0d0.and.y.gt.0d0)write(13,*)x_screen,y_screen,1
            if(z.gt.0d0.and.y.lt.0d0)write(13,*)x_screen,y_screen,2
            if(z.lt.0d0.and.y.gt.0d0)write(13,*)x_screen,y_screen,3
            if(z.lt.0d0.and.y.lt.0d0)write(13,*)x_screen,y_screen,4
            exit
        else
        end if
    end do
1   return
    end subroutine
    
subroutine jet_image_task(x_screen,y_screen,r_eh)   
    use jet_parameters
    use observer_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(8)::vector
    real*8,dimension(4)::jet_velocity,vector_obs
    real*8,dimension(4,4)::g_ab
    real*8 x_screen,y_screen,r_eh
    real*8 redshift_factor
    real*8 deepth,i_obs
    real*8 step,r_local_velocity
    real*8 j_local,a_local 
    deepth=0d0
    i_obs=0d0
    call photon_initial_values(x_screen,y_screen,vector)
    vector_obs(1)=vector(5)
    vector_obs(2)=vector(6)
    vector_obs(3)=vector(7)
    vector_obs(4)=vector(8)
    do while(.true.)
        call change_step(vector(2),vector(3),step,r_eh)
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))exit
        if(vector(2).ge.observer)exit
        call covariant_metric(vector(2),vector(3),g_ab) 
        if(vector(2).le.jet_range.and.vector(2).ge.jet_base)then
            call jet_dynamic(vector(2),vector(3),vector(4),jet_velocity,r_local_velocity)
            redshift_factor=(vector_obs(1)*vt_obs+vector_obs(2)*vr_obs+vector_obs(3)*vtheta_obs+vector_obs(4)*vphi_obs)/(vector(5)*jet_velocity(1)+vector(6)*jet_velocity(2)+vector(7)*jet_velocity(3)+vector(8)*jet_velocity(4))
            call jet_emission(vector(2),vector(3),vector(4),j_local,a_local)
            deepth=deepth+a_local*step/redshift_factor
            i_obs=i_obs+(redshift_factor*redshift_factor*j_local*dexp(-deepth))*step
        else
        end if
    end do
    write(2,*)x_screen,y_screen,i_obs
    return
    end subroutine
    
 subroutine point_source_image_task(x_screen,y_screen,r_eh)
    use point_source_parameters
    use observer_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(8)::vector,vector0
    real*8,dimension(4)::vector_obs,flow_velocity
    real*8 r_local_v,theta_local_v,phi_local_v
    real*8 x_screen,y_screen,r_eh
    real*8 deepth,i_obs,step
    real*8 radius1,radius0
    real*8 x1,y1,z1,x0,y0,z0
    real*8 redshift_factor
    real*8 j_local,a_local
    deepth=0d0
    i_obs=0d0
    call photon_initial_values(x_screen,y_screen,vector)
    vector_obs(1)=vector(5)
    vector_obs(2)=vector(6)
    vector_obs(3)=vector(7)
    vector_obs(4)=vector(8)
    do while(.true.)
        call change_step(vector(2),vector(3),step,r_eh)
        call record_values(vector,vector0)
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))exit
        if(vector(2).ge.observer)exit 
        call spherical_transfer(vector,x1,y1,z1)
        call spherical_transfer(vector0,x0,y0,z0)
        radius1=dsqrt((x1-point_source_x)**2d0+(y1-point_source_y)**2d0+(z1-point_source_z)**2d0)
        radius0=dsqrt((x0-point_source_x)**2d0+(y0-point_source_y)**2d0+(z0-point_source_z)**2d0)
        if(radius1.le.point_source_radius)then
            call point_source_dynamic(vector0(2),vector0(3),vector0(4),flow_velocity,r_local_v,theta_local_v,phi_local_v)
            redshift_factor=(vector_obs(1)*vt_obs+vector_obs(2)*vr_obs+vector_obs(3)*vtheta_obs+vector_obs(4)*vphi_obs)/(vector0(5)*flow_velocity(1)+vector0(6)*flow_velocity(2)+vector0(7)*flow_velocity(3)+vector0(8)*flow_velocity(4))
            call point_source_emission(radius0,j_local,a_local)
            deepth=deepth+a_local*step/redshift_factor
            i_obs=i_obs+(redshift_factor*redshift_factor*j_local*dexp(-deepth))*step
        else
        end if
    end do
    write(2,*)x_screen,y_screen,i_obs
    return
end subroutine

subroutine spherically_accretion_image_task(x_screen,y_screen,r_eh)
    use spherically_accretion_flow_parameters
    use observer_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(8)::vector,vector0
    real*8,dimension(4)::flow_velocity,vector_obs
    real*8 x_screen,y_screen,r_eh
    real*8 deepth,i_obs,step
    real*8 redshift_factor
    real*8 j_local,a_local
    real*8 r_local_velocity
    deepth=0d0
    i_obs=0d0
    call photon_initial_values(x_screen,y_screen,vector)
    vector_obs(1)=vector(5)
    vector_obs(2)=vector(6)
    vector_obs(3)=vector(7)
    vector_obs(4)=vector(8)
    do while(.true.)
        call change_step(vector(2),vector(3),step,r_eh)
        call record_values(vector,vector0)
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))exit
        if(vector(2).ge.observer)exit 
        if(vector(2).ge.spherically_flow_inner.and.vector(2).le.spherically_flow_outer)then
            call spherically_accretion_flow_velocity(vector0(2),vector0(3),vector0(4),flow_velocity,r_local_velocity,r_eh)
            redshift_factor=(vector_obs(1)*vt_obs+vector_obs(2)*vr_obs+vector_obs(3)*vtheta_obs+vector_obs(4)*vphi_obs)/(vector0(5)*flow_velocity(1)+vector0(6)*flow_velocity(2)+vector0(7)*flow_velocity(3)+vector0(8)*flow_velocity(4))
            call spherically_accretion_flow_emission(vector0(2),j_local,a_local)
            deepth=deepth+a_local*step/redshift_factor
            i_obs=i_obs+(redshift_factor*redshift_factor*j_local*dexp(-deepth))*step
        else
        end if
    end do
    write(2,*)x_screen,y_screen,i_obs
    return
end subroutine

subroutine light_curve(r_eh)
    use observer_parameters
    use massive_particle_orbit_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(8)::vector,vector0
    real*8,dimension(max_point,8)::path_data
    real*8,dimension(max_point)::spot_radius
    real*8 r_eh
    real*8 x_mp,y_mp,z_mp
    real*8 max_radius,min_radius
    real*8 integration_time
    real*8 temp_hami
    real*8 x_screen,y_screen
    integer i,j,k
    integer loop,code_stop
    x_screen=0d0
    y_screen=0d0
    integration_time=0d0
    loop=0
    code_stop=0
    path_data=0d0
    call massive_particle_initialization(vector)
    do i=1,max_point,1
        do j=1,sample,1
            call timestep(integration_time,-step_mp)
            call record_values(vector,vector0)
            call rksixth(vector,step_mp)
            if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))then
                write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
                write(6,*)"Particle captured by the event horizon."
                write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
                code_stop=1
                exit
            else
            end if
            if(vector(2).ge.500d0)then
                write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
                write(6,*)"Particle escaped to infinity."
                write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
                code_stop=1
                exit
            else
            end if
        end do
        if(code_stop==1)exit
        loop=loop+1
        do k=1,8,1
            path_data(loop,k)=vector(k)
        end do
        spot_radius(loop)=vector(2)
        call hamiltonian(vector,temp_hami)
        call spherical_transfer(vector,x_mp,y_mp,z_mp)
        if(vector(2).le.200d0)write(14,*)x_mp,y_mp,z_mp
        if(vector(2).le.200d0)write(15,*)loop,dlog10(dabs(temp_hami+0.5d0)+1d-16)
    end do
    max_radius=maxval(spot_radius)+2d0
    min_radius=minval(spot_radius)-2d0
    !$omp parallel do schedule(dynamic,1) default(none) firstprivate(x_ini,y_ini,delta_x,delta_y,resolution_x,resolution_y,x_screen,y_screen,r_eh,max_radius,min_radius,path_data,loop)
    do i=1,resolution_y,1
        y_screen=y_ini+(float(i)-1d0)*delta_y
        do j=1,resolution_x,1
            x_screen=x_ini+(float(j)-1d0)*delta_x
            call orbit_tracing_task(x_screen,y_screen,r_eh,max_radius,min_radius,path_data,loop)
        end do
    end do
    !$omp end parallel do
    return
    end subroutine

subroutine orbit_tracing_task(x_screen,y_screen,r_eh,max_radius,min_radius,path_data,loop)
    use massive_particle_orbit_parameters
    use raytracing_parameters
    use metric_parameters
    implicit none
    real*8,dimension(8)::vector,vector0
    real*8,dimension(max_point,8)::path_data
    real*8 x_screen,y_screen,r_eh,max_radius,min_radius
    real*8 step
    real*8 x1,y1,z1,x0,y0,z0,hotspot_x,hotspot_y,hotspot_z
    real*8 radius0,radius1
    integer i,loop
    call photon_initial_values(x_screen,y_screen,vector)
    do while(.true.)
        call change_step(vector(2),vector(3),step,r_eh)
        call record_values(vector,vector0)
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))exit
        if(vector(2).ge.observer)exit 
        if(vector0(2).le.max_radius.and.vector0(2).ge.min_radius)then
            do i=1,loop,1
                call spherical_transfer(vector,x1,y1,z1)
                call spherical_transfer(vector0,x0,y0,z0)
                hotspot_x=dsqrt(path_data(i,2)**2d0+spin*spin)*dsin(path_data(i,3))*dcos(path_data(i,4))
                hotspot_y=dsqrt(path_data(i,2)**2d0+spin*spin)*dsin(path_data(i,3))*dsin(path_data(i,4))
                hotspot_z=path_data(i,2)*dcos(path_data(i,3))
                radius1=dsqrt((x1-hotspot_x)**2d0+(y1-hotspot_y)**2d0+(z1-hotspot_z)**2d0)
                radius0=dsqrt((x0-hotspot_x)**2d0+(y0-hotspot_y)**2d0+(z0-hotspot_z)**2d0)
                if(radius0.le.particle_radius.and.radius1.ge.particle_radius)then
                    write(16,*)dabs(path_data(i,1))+dabs(vector0(1))
                else
                end if
            end do
        else
        end if
    end do
    return
end subroutine

subroutine grid_image_task(x_screen,y_screen,r_eh)
    use raytracing_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh
    real*8,dimension(8)::vector,vector0
    real*8 step
    real*8 x,y,z,x0,y0,z0,x_hit,y_hit,z_hit
    real*8 loop_x,loop_y,loop_z
    real*8 hit_plane
    integer i,j,k
    hit_plane=0d0
    call photon_initial_values(x_screen,y_screen,vector)
        do while(.true.)
            call change_step(vector(2),vector(3),step,r_eh)
            call record_values(vector,vector0)
            call rkf(vector,step)
            if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))then
                write(1,*)x_screen,y_screen
                exit
            else
            end if
            if(vector(2).ge.observer)exit
            call spherical_transfer(vector,x,y,z)
            call spherical_transfer(vector0,x0,y0,z0)
            if(dabs(x).le.10d0.and.dabs(y).le.10d0)then
                if((z-hit_plane)*(z0-hit_plane).le.0d0)then
                    x_hit=x0+(x-x0)*(hit_plane-z0)/(z-z0)
                    y_hit=y0+(y-y0)*(hit_plane-z0)/(z-z0)
                    do i=1,21,1
                        loop_x=-10d0+1d0*(float(i)-1d0)
                        if((x_hit-loop_x)**2d0.le.0.005d0)then
                            write(17,*)x_screen,y_screen,hit_plane
                            goto 1
                        else
                        end if
                    end do
                    do j=1,21,1
                        loop_y=-10d0+1d0*(float(j)-1d0)
                        if((y_hit-loop_y)**2d0.le.0.005d0)then
                            write(17,*)x_screen,y_screen,hit_plane
                            goto 1
                        else
                        end if
                    end do
                    goto 1
                else
                end if
            else
            end if
        end do
1   return
    end subroutine

subroutine tilted_disk_four_color_image_task(x_screen,y_screen,r_eh)
    use tilted_thin_disk_parameters
    use raytracing_parameters
    use disk_range_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh
    real*8,dimension(8)::vector,vector0
    real*8 step
    real*8 x,y,z,x0,y0,z0
    integer color
    color=-1
    call photon_initial_values(x_screen,y_screen,vector)
    do while(.true.)
        call change_step(vector(2),vector(3),step,r_eh)
        call record_values(vector,vector0)
        call rkf(vector,step)
        if(vector(2).ge.observer)then
            color=5
            write(13,*)x_screen,y_screen,color
            exit
        else
        end if
        if(vector(2).le.(r_eh+hit_error).or.isnan(vector(2)))then
            color=0
            write(1,*)x_screen,y_screen
            write(13,*)x_screen,y_screen,color
            exit
        else
        end if
        call spherical_transfer(vector,x,y,z)
        call spherical_transfer(vector0,x0,y0,z0)
        if(((x*dtan(disk_inclination)-z)*(x0*dtan(disk_inclination)-z0)).le.0d0.and.vector(2).le.r_outer)then
            if(x.ge.0d0.and.y.gt.0d0)color=4
            if(x.gt.0d0.and.y.le.0d0)color=2
            if(x.lt.0d0.and.y.ge.0d0)color=1
            if(x.le.0d0.and.y.lt.0d0)color=3
            write(13,*)x_screen,y_screen,color
            exit
        else
        end if
    end do
    return
    end subroutine