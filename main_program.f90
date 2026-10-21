program octopus_spin
    use omp_lib 
    use task_parameters
    real*8 ini_t,over_t
    real*8 r_eh
    open(1,file='data/shadow.dat')
    open(2,file='data/image.dat')
    open(3,file='data/velocity-field.dat')
    open(4,file='data/j-profile-roz.dat')
    open(5,file='data/j-profile-xoy.dat')
    open(7,file='data/j-profile-r.dat')
    open(8,file='data/redshift-1.dat')
    open(9,file='data/redshift-2.dat')
    open(10,file='data/redshift-3.dat')
    open(11,file='data/redshift-4.dat')
    open(12,file='data/thin-disk-hit.dat')
    open(13,file='data/four-color-image.dat')
    open(14,file='data/particle-path.dat')
    open(15,file='data/energy-error.dat')
    open(16,file='data/time-series.dat')
    open(17,file='data/grid.dat')
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    write(6,*)"Initializing SpinOctopus..."
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    call code_check(r_eh)
    ini_t=omp_get_wtime()
    call omp_set_num_threads(cpu)
    select case(task_model)
    case(1) 
        call shadow_profile(r_eh)
    case(2)
        call thick_equatorial_disk_radiation(r_eh)
    case(3)
        call thick_equatorial_disk_image(r_eh)
    case(4)
        call thin_equatorial_disk_image(r_eh)
    case(5)
        call four_color_image(r_eh)
    case(6)
        call jet_radiation(r_eh)
    case(7)
        call jet_image(r_eh)
    case(8)
        call photon_path(r_eh)
    case(9)
        call massive_particle_path(r_eh)
    case(10)
        call point_source_image(r_eh)
    case(11)
        call light_curve(r_eh)
    case(12)
        call spherically_accretion_flow_image(r_eh)
    case(13)
        call grid_image(r_eh)
    case(14)
        call tilted_disk_four_color_image(r_eh)
    end select
    over_t=omp_get_wtime()
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    write(6,*)"Cpu Cost:",dabs(over_t-ini_t),"seconds"
    write(6,*)"SpinOctopus execution finished."
    write(6,*)"!-----------------------------------------------------------------------------------------------------------!"
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(7)
    close(8)
    close(9)
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    stop
    end 
        