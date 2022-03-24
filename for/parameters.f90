module parameters
  use fft
  use domain
  implicit none
  save


  real(rkind)  time
  integer time_step, rk_step
  real(rkind) save_flow_time, save_stats_time, save_movie_time

  integer previous_time_step

  ! Input.dat
  character(len=35)   flavor
  logical             use_LES
  real(rkind)         nu
  real(rkind)         nu_v_scale
  real(rkind)         beta
  real(rkind)         delta_t, dt, delta_t_next_event, kick, ubulk0, px0

  logical             create_new_flow
  real(rkind)         wall_time_limit, time_limit
  real(rkind)         start_wall_time, previous_wall_time, end_wall_time

  logical             variable_dt, first_time
  logical             reset_time
  real(rkind)         CFL
  integer             update_dt

  real(rkind)         save_flow_dt, save_stats_dt
  real(rkind)         save_movie_dt

  logical             create_new_th(1:N_th)
  real(rkind)         Ri(1:N_th), Pr(1:N_th)
  logical             filter_th(1:N_th)
  integer             filter_int(1:N_th)

  integer     num_read_th
  integer     read_th_index(1:N_th)
  real(rkind) dTHdX(1:N_th), dTHdZ(1:N_th)
  real(rkind) dWdX ! Background vorticity


  integer     IC_type, f_type
  logical     physical_noise
  logical     homogeneousX


  ! Periodic
  real(rkind) ek0, ek, epsilon_target
  logical     background_grad(1:N_th)

  ! Rotating Flows
  real(rkind) Ro_inv, grav_x, grav_y, grav_z, delta

  ! Parameters for oscillatory forcing
  real(rkind) omega0, amp_omega0, force_start
  real(rkind) w_BC_Ymax_c1_transient


  ! Numerical parameters
  real(rkind)  h_bar(3), beta_bar(3), zeta_bar(3) ! For RK
  integer  time_ad_meth
  integer les_model_type






contains

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine init_parameters
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real version, current_version

    integer i, j, k, n
    real(rkind) Re
    logical start_file_exists

    open (11, file='input.dat', form='formatted', status='old')

    ! Read input file.
    !   (Note - if you change the following section of code, update the
    !    CURRENT_VERSION number to make obsolete previous input files !)

    current_version = 3.4
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *) flavor, version
    if (version /= current_version) stop 'Wrong input data format.'
    read (11, *)
    read (11, *) use_mpi, use_LES
    if (use_mpi .eqv. .false.) stop 'Serial processing has been deprecated in diablo3.'
    read (11, *)
    read (11, *) Re, beta, Lx, Lz, Ly
    nu = 1.d0 / Re
    read (11, *)
    read (11, *) nu_v_scale
    read (11, *)
    read (11, *) num_per_dir, create_new_flow
    read (11, *)
    read (11, *) wall_time_limit, time_limit, delta_t, reset_time, &
      variable_dt, CFL, update_dt
    read (11, *)
    read (11, *) verbosity, save_flow_dt, save_stats_dt, save_movie_dt, XcMovie, ZcMovie, YcMovie
    read (11, *)
    ! Read in the parameters for the N_th scalars
    do n = 1, N_th
      read (11, *)
      read (11, *) create_new_th(n)
      read (11, *)
      read (11, *) filter_th(n), filter_int(n)
      read (11, *)
      read (11, *) Ri(n), Pr(n)
    end do


    ! Initialize MPI Variables
    call init_mpi


    if (rank == 0) then
      write (*, *)
      write (*, *) '             ****** WELCOME TO DIABLO ******'
      write (*, *)
    end if

    inquire (file="start.h5", exist=start_file_exists)
    if (start_file_exists) then
      create_new_flow = .false.
      create_new_th(1) = .false.
    end if


    ! Initialize case-specific packages
    if (num_per_dir == 3) then
      stop 'Error: Triply-Periodic Box has been deprecated!'
    elseif (num_per_dir == 2) then
      call input_chan
      call create_grid_chan
      call init_chan_mpi
      if (save_movie_dt /= 0) then
        call init_chan_movie
      end if
    elseif (num_per_dir == 1) then
      stop 'Error: Duct not implemented!'
    elseif (num_per_dir == 0) then
      stop 'Error: Cavity not implemented!'
    end if

    if (rank == 0) &
      write (*, '("Use LES: " L1)') use_LES

    ! Initialize grid
    if (rank == 0) then
      write (*, '("Flavor:  ", A35)') flavor
      write (*, '("Nx    =  ", I10)') Nx
      write (*, '("Ny    =  ", I10)') Nz
      write (*, '("Nz(p) =  ", I10)') Nyp
      do n = 1, N_th
        write (*, '("Scalar Number: ", I2)') n
        write (*, '("  Richardson number = ", ES12.5)') Ri(n)
        write (*, '("  Prandtl number    = ", ES12.5)') Pr(n)
      end do
      write (*, '("Nu   = ", ES12.5)') nu
      write (*, '("Beta = ", ES12.5)') beta
    end if

    ! Initialize RKW3 parameters.
    h_bar(1) = delta_t * (8.d0 / 15.d0)
    h_bar(2) = delta_t * (2.d0 / 15.d0)
    h_bar(3) = delta_t * (5.d0 / 15.d0)
    beta_bar(1) = 1.d0
    beta_bar(2) = 25.d0 / 8.d0
    beta_bar(3) = 9.d0 / 4.d0
    zeta_bar(1) = 0.d0
    zeta_bar(2) = -17.d0 / 8.d0
    zeta_bar(3) = -5.d0 / 4.d0


    return
  end




  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine input_chan
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real version, current_version
    integer i, j, k, n
    real(rkind) ro

    ! Read in input parameters specific for channel flow case
    open (11, file='input_chan.dat', form='formatted', status='old')
    ! Read input file.

    current_version = 3.4
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *) version
    if (version /= current_version) &
      stop 'Wrong input data format input_chan'
    read (11, *)
    read (11, *) time_ad_meth
    read (11, *)
    read (11, *) les_model_type
    read (11, *)
    read (11, *) IC_type, kick, physical_noise
    read (11, *)
    read (11, *) ro
    Ro_inv = 1.d0 / ro
    read (11, *)
    read (11, *) delta
    read (11, *)
    !read (11, *) dWdX ! Background vorticity
    !read (11, *)
    read (11, *) grav_x, grav_z, grav_y
    read (11, *)
    read (11, *) f_type, ubulk0, px0, omega0, amp_omega0, force_start
    read (11, *)
    read (11, *)
    read (11, *) u_BC_Ymin, u_BC_Ymin_c1
    read (11, *)
    read (11, *) w_BC_Ymin, w_BC_Ymin_c1
    read (11, *)
    read (11, *) v_BC_Ymin, v_BC_Ymin_c1
    read (11, *)
    read (11, *) u_BC_Ymax, u_BC_Ymax_c1
    read (11, *)
    read (11, *) w_BC_Ymax, w_BC_Ymax_c1
    read (11, *)
    read (11, *) v_BC_Ymax, v_BC_Ymax_c1
    read (11, *)
    ! Read in boundary conditions and background gradients for the N_th scalars
    do n = 1, N_th
      read (11, *)
      read (11, *) th_BC_Ymin(n), th_BC_Ymin_c1(n)
      read (11, *)
      read (11, *) th_BC_Ymax(n), th_BC_Ymax_c1(n)
    end do

    if (rank == 0) write (*, '("Ro Inverse = " ES26.18)') Ro_inv


    ! Compensate no-slip BC in the GS flow direction due to dTHdx
    !   AND also define dTHdx & dTHdz
    if (IC_Type == 4 .or. IC_Type == 5) then ! Infinite Front
      if (w_BC_Ymin == 1) then
        w_BC_Ymin_c1 = w_BC_Ymin_c1 - 1.d0
      end if
      if (w_BC_Ymax == 1) then
        w_BC_Ymax_c1 = w_BC_Ymax_c1 - 1.d0
      end if
      dTHdX(1) = Ro_inv / delta
      dTHdZ(1) = 0.d0

    else if (IC_Type == 6 .or. IC_Type == 7 .or. IC_Type == 8) then ! Finite Front
      if (w_BC_Ymin == 1) then
        w_BC_Ymin_c1 = w_BC_Ymin_c1 - 2.d0 * delta / Lx
      end if
      if (w_BC_Ymax == 1) then
        w_BC_Ymax_c1 = w_BC_Ymax_c1 - 2.d0 * delta / Lx
      end if
      dTHdX(1) = 2.d0 / Lx * Ro_inv
      dTHdZ(1) = 0.d0

    end if

    w_BC_Ymax_c1_transient = w_BC_Ymax_c1 ! Mean of surface forcing (compensating for GS flow)


    ! Set the valid averaging directions depending on the IC
    if (IC_Type == 5 .or. IC_Type == 6 .or. IC_Type == 7 .or. IC_Type == 8) then
      homogeneousX = .false.
    else ! Infinite, homogeneous front (or other IC...)
      homogeneousX = .true. ! Assume the x-direction is a valid averaging dimension
    endif


    return
  end




end module parameters
