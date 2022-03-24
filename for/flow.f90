module flow
  use fft
  use domain
  use parameters
  use tools
  implicit none
  save

  ! Flow !

  ! 3D
  real(rkind), pointer, contiguous, dimension(:,:,:) :: u1,u2,u3,p,r1,r2,r3,f1,f2,f3,s1
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: cu1,cu2,cu3,cp,cr1,cr2,cr3, &
                                                           cf1,cf2,cf3,cs1

  ! 4D
  real(rkind), pointer, contiguous, dimension(:,:,:,:) :: th,fth,rth
  complex(rkind), pointer, contiguous, dimension(:,:,:,:) :: cth,cfth,crth



  ! LES !
  real(rkind), pointer, contiguous, dimension(:,:,:) :: Sij1,Sij2,Sij3,Sij4,Sij5,Sij6,temp
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: CSij1,CSij2,CSij3,CSij4,CSij5,CSij6,ctemp
  ! AMD !
  real(rkind), pointer, contiguous, dimension(:,:,:) :: temp_th, s1_th
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: Ctemp_th, cs1_th
  real(rkind), pointer, contiguous, dimension(:,:,:) :: Oij4,Oij5,Oij6
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: COij4,COij5,COij6
  real(rkind), pointer, contiguous, dimension(:,:,:) :: SIij1,SIij2,SIij3,SIij4,SIij5,SIij6
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: CSIij1,CSIij2,CSIij3,CSIij4,CSIij5,CSIij6
  real(rkind), pointer, contiguous, dimension(:,:,:) :: du1dx,du2dy,du3dz, &
                                                           du2dx,du3dx, &
                                                           du1dz,du2dz, &
                                                           du1dy,du3dy
  complex(rkind), pointer, contiguous, dimension(:,:,:) :: Cdu1dx,Cdu2dy,Cdu3dz, &
                                                           Cdu2dx,Cdu3dx, &
                                                           Cdu1dz,Cdu2dz, &
                                                           Cdu1dy,Cdu3dy
  real(rkind), pointer, contiguous, dimension(:,:,:,:) :: dthetadx,dthetady,dthetadz
  complex(rkind), pointer, contiguous, dimension(:,:,:,:) :: Cdthetadx,Cdthetady,Cdthetadz


  real(rkind) nu_u1(0:Nyp+1) ! For plane-averaged momentum budget
  real(rkind) nu_u3(0:Nyp+1)

  real(rkind) cross
  logical flag_save_LES
  real(rkind) :: nu_t(0:Nx+1, 0:Nzp+1, 0:Nyp+1) = 0.d0
  real(rkind) :: kappa_t(0:Nx+1, 0:Nzp+1, 0:Nyp+1, 1:N_th) = 0.d0
  integer j1, j2



  ! Diagnostics !
  real(rkind), dimension(0:Nyp+1) :: urms, vrms, wrms, &
                                     ume,  vme,  wme,  &
                                     uv,   uw,   wv,   &
                                     dudy, dwdy, shear,&
                                     omega_x, omega_y, omega_z, &
                                     u1y_left, mke, &
                                     uu_dudx, vv_dvdy, wu_dwdx, &
                                     vu_dvdx, uv_dudy, wv_dwdy

  real(rkind) vme_xy(0:Nxm1, 1:Nyp + 1)
  real(rkind), dimension(0:Nxm1, 1:Nyp) :: uvar_xy, vvar_xy, wvar_xy, thvar_xy, &
                                           uu_xy, vv_xy, ww_xy, &
                                           ume_xy, var_xy, wme_xy, &
                                           dudx_m, dwdx_m, dvdx_m ! For sharing mean gradients with Production from Dissipation
  real(rkind), dimension(0:Nxm1, 1:Nyp, 1:N_th) :: thme_xy

  real(rkind) urms_b, vrms_b, wrms_b, tke_b, u1y_left_b

  real(rkind), dimension(0:Nyp+1,1:N_th) :: thrms, thme, &
                                            thv, thv_m, dthdy, pe_diss
  real(rkind) :: thrms_b(1:N_th)

  complex(rkind)  cuu1_yx(0:Nyp+1,0:Nxp-1)

  real(rkind) epsilon(0:Nyp+1), epsilon_m(0:Nyp+1)









contains

  include 'solver.f90'
  include 'les.f90'
  include 'diagnostics.f90'
  include 'phdf5.f90'



  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine init_flow
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    integer j, n

    call allocate_all
    if (use_LES) call allocate_all_les

    ! Initialize FFT package (includes defining the wavenumber vectors).
    call init_fft(u1,cu1)

    ! Initialize values for reading of scalars
    num_read_th = 0
    do n = 1, N_th
      if (create_new_th(n)) then
        num_read_th = num_read_th
      else
        num_read_th = num_read_th + 1
        read_th_index(num_read_th) = n
      end if
    end do
    call create_th_chan

    ! Create flow.
    if (create_new_flow) then
      if (num_per_dir == 3) then
        stop 'Error: Triply-Periodic Box has been deprecated!'
      elseif (num_per_dir == 2) then
        call create_flow_chan
      elseif (num_per_dir == 1) then
        stop 'Error: Duct not implemented yet!'
      elseif (num_per_dir == 0) then
        stop 'Error: Cavity not implemented yet!'
      end if
      if (rank == 0) &
        write (*, '("A new flowfield has been created.")')
    else
      if (rank == 0) &
        write (*, '("Reading flow...")')
      call read_flow
      if (rank == 0) &
        write (*, '("Starting at time step ", I10)') time_step


      call ghost_chan_mpi

    end if

    ! Initialize flow
    if (reset_time .or. create_new_flow) then
      save_flow_time = save_flow_dt
      save_stats_time = save_stats_dt
      save_movie_time = save_movie_dt
      time_step = 0
      time = 0

      call save_stats(save_movie_dt/=0,.false.)
      if (use_LES) call save_stats_LES_OOL(.true.)
    end if

    previous_time_step = time_step
    call wall_time(previous_wall_time)

    if (rank == 0) then
      write (*, *)
      write (*, *) '             ****** Done Initialising ******'
      write (*, *)
    end if

    return
  end




  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine allocate_all
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|


    call alloc_array3D(u1,cu1)
    call alloc_array3D(u2,cu2)
    call alloc_array3D(u3,cu3)
    call alloc_array3D(p ,cp)
    call alloc_array3D(r1,cr1)
    call alloc_array3D(r2,cr2)
    call alloc_array3D(r3,cr3)
    call alloc_array3D(f1,cf1)
    call alloc_array3D(f2,cf2)
    call alloc_array3D(f3,cf3)
    call alloc_array3D(s1,cs1)

    call alloc_array4D(th, cth)   ! Not using the same memory!
    call alloc_array4D(fth,cfth)
    call alloc_array4D(rth,crth)

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine allocate_all_les
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    call alloc_array3D(Sij1,CSij1)
    call alloc_array3D(Sij2,CSij2)
    call alloc_array3D(Sij3,CSij3)
    call alloc_array3D(Sij4,CSij4)
    call alloc_array3D(Sij5,CSij5)
    call alloc_array3D(Sij6,CSij6)
    call alloc_array3D(temp,ctemp)


    ! Need Extra AMD Arrays
    call alloc_array3D(temp_th,ctemp_th)
    call alloc_array3D(s1_th,cs1_th)

    call alloc_array3D(Oij4,COij4)
    call alloc_array3D(Oij5,COij5)
    call alloc_array3D(Oij6,COij6)

    call alloc_array3D(SIij1,CSIij1)
    call alloc_array3D(SIij2,CSIij2)
    call alloc_array3D(SIij3,CSIij3)
    call alloc_array3D(SIij4,CSIij4)
    call alloc_array3D(SIij5,CSIij5)
    call alloc_array3D(SIij6,CSIij6)

    ! Can we actually point some of these to the SIij ?
    call alloc_array3D(du1dx,Cdu1dx)
    call alloc_array3D(du2dy,Cdu2dy)
    call alloc_array3D(du3dz,Cdu3dz)

    call alloc_array3D(du2dx,Cdu2dx)
    call alloc_array3D(du3dx,Cdu3dx)

    call alloc_array3D(du1dz,Cdu1dz)
    call alloc_array3D(du2dz,Cdu2dz)

    call alloc_array3D(du1dy,Cdu1dy)
    call alloc_array3D(du3dy,Cdu3dy)

    call alloc_array4D(dthetadx,Cdthetadx)   ! Not using the same memory!
    call alloc_array4D(dthetady,Cdthetady)
    call alloc_array4D(dthetadz,Cdthetadz)

  end



  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine deallocate_all
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    ! call fftw_free(p) ! at pointer p...

  end




  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine create_flow_chan
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|


    integer i, j, k, n
    real(rkind) rnum1, rnum2, rnum3
    integer, dimension(:), allocatable :: seed

    ! Initialize the random number generator
    call random_seed(size=k)
    allocate (seed(1:k))
    do i = 1, k
      seed(i) = rank*k + i + 999
    end do
    call random_seed(put=seed)

    ! UBULK0 and kick should be set in input.dat

    ! IC_type is set in input_chan.dat and can be used to easily
    ! control which initial condition is used.  A few examples
    ! are given here. These can be modified, or new types can be
    ! added

    if (IC_type == 0) then
      ! Parabolic profile for laminar closed channel flow
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = (3./2.) * ubulk0 * (1.d0 - gyf(j)**2.)
            u2(i, k, j) = 0.
            u3(i, k, j) = 0.
          end do
        end do
      end do
    else if (IC_type == 1) then
      ! Laminar profile for open channel flow
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          do j = 1, Nyp
            u1(i, k, j) = -(3./2.) * ubulk0 * gyf(j)**2.+3.*ubulk0 * gyf(j)
            u2(i, k, j) = 0.
            u3(i, k, j) = 0.
          end do
          u1(i, k, 0) = 0.
          u3(i, k, 0) = 0.
          u1(i, k, Nyp + 1) = 0.
          u3(i, k, Nyp + 1) = 0.
        end do
      end do
    else if (IC_type == 2) then
      ! Linear profile for laminar Couette flow
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = gyf(j)
            u2(i, k, j) = 0.
            u3(i, k, j) = 0.
          end do
        end do
      end do
    else if (IC_type == 3) then
      ! Tanh shear layer
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = tanh(gyf(j))
            u2(i, k, j) = 0.d0
            u3(i, k, j) = 0.d0
          end do
        end do
      end do
    else if (IC_type == 4) then
      ! Infinite Front
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            u3(i, k, j) = 0.d0 ! TWS exactly balanced in solver
          end do
        end do
      end do
    else if (IC_type == 5) then
      ! Balanced Infinite Front + Barotropic Background vorticity (W)
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            u3(i, k, j) = (tanh((gx(i) - 0.5 * Lx) / delta) - &
                              (tanh(0.5*Lx/delta) - tanh(-0.5*Lx/delta)) / Lx * (gx(i) - 0.5 * Lx)) & ! Make periodic
                            / (delta*(1.d0/delta - (tanh(0.5*Lx/delta)-tanh(-0.5*Lx/delta))/Lx)) ! Scale so middle gradient is exactly 1/delta
          end do
        end do
      end do
    else if (IC_type == 6) then
      ! Finite Sloping Front
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            if (Ri(1) == 0) then
              u3(i, k, j) = (gyf(j) - 0.5 * Ly) / &
                            cosh((gx(i) - 0.5 * Lx) / delta)**2.d0 &
                            / tanh(0.5d0 * Lx / delta) &
                            - 2.d0 * delta / Lx * (gyf(j) - 0.5 * Ly)
            else
              u3(i, k, j) = delta * sqrt(Ri(1)**2.d0 + Ro_inv**2.d0) / Ri(1) * &
                            tanh(((gx(i) - 0.5 * Lx) * Ro_inv + &
                                  (gyf(j) - 0.5 * Ly) * Ri(1)) / &
                                 (delta * sqrt(Ri(1)**2.d0 + Ro_inv**2.d0))) &
                            - 2.d0 * delta / Lx * (gyf(j) - 0.5 * Ly)
            end if
          end do
        end do
      end do
    else if (IC_type == 7) then
      ! Finite Front, with u3 = 0 at bottom (for no-slip BC)
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            u3(i, k, j) = (gyf(j)) / &
                            cosh((gx(i) - 0.5 * Lx) / delta)**2.d0 &
                            / tanh(0.5d0 * Lx / delta) &
                            - 2.d0 * delta / Lx * (gyf(j))
          end do
        end do
      end do
    else if (IC_type == 8) then
      ! Finite Front Unbalanced
      do j = 0, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            u1(i, k, j) = 0.d0
            u2(i, k, j) = 0.d0
            u3(i, k, j) = -2.d0 * delta / Lx * (gyf(j) - 0.5 * Ly) ! Remove the TWS from solver
          end do
        end do
      end do
    else
      write (*, '("Warning, unsupported IC_type in create_flow.")')
    end if


    !call add_SI_Mode


    if (physical_noise) then
      ! Add random noise in physical space
      call random_number(rnum1)
      call random_number(rnum1)
      call random_number(rnum1)
      do j = 0, Nyp + 1
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            call random_number(rnum1)
            u1(i, k, j) = u1(i, k, j) + kick * (rnum1 - 0.5d0)
            call random_number(rnum1)
            u2(i, k, j) = u2(i, k, j) + kick * (rnum1 - 0.5d0)
            call random_number(rnum1)
            u3(i, k, j) = u3(i, k, j) + kick * (rnum1 - 0.5d0)
          end do
        end do
      end do


      ! Convert to Fourier space
      call fft_xz_to_fourier(u1, cu1)
      call fft_xz_to_fourier(u2, cu2)
      call fft_xz_to_fourier(u3, cu3)
      call fft_xz_to_fourier(p, cp)

    else ! Optionally, add random noise in Fourier space instead

      ! Convert to Fourier space
      call fft_xz_to_fourier(u1, cu1)
      call fft_xz_to_fourier(u2, cu2)
      call fft_xz_to_fourier(u3, cu3)
      call fft_xz_to_fourier(p, cp)


      do i = 0, Nxp - 1
        do j = 1, Nyp
          do k = 0, twoNkz
            ! Now, give the velocity field a random perturbation
            call random_number(rnum1)
            call random_number(rnum2)
            cu1(i,k,j) = cu1(i,k,j) + cmplx((rnum1 - 0.5d0), (rnum2 - 0.5d0)) * kick
            call random_number(rnum1)
            call random_number(rnum2)
            cu2(i,k,j) = cu2(i,k,j) + cmplx((rnum1 - 0.5d0), (rnum2 - 0.5d0)) * kick
            call random_number(rnum1)
            call random_number(rnum2)
            cu3(i,k,j) = cu3(i,k,j) + cmplx((rnum1 - 0.5d0), (rnum2 - 0.5d0)) * kick
          end do
          if (twoNkz == 0) then
            ! Here, In the 2d case we want to add a kick to the mean in z
            k = 0
            call random_number(rnum1)
            call random_number(rnum2)
            call random_number(rnum3)

            if (IC_type == 3) then
              cu1(i,k,j) = cu1(i,k,j) + (rnum1 - 0.5) * kick * exp(-(gyf(j) * 20.d0)**2.d0)
              cu2(i,k,j) = cu2(i,k,j) + (rnum1 - 0.5) * kick * exp(-(gyf(j) * 20.d0)**2.d0)
              cu3(i,k,j) = cu3(i,k,j) + (rnum1 - 0.5) * kick * exp(-(gyf(j) * 20.d0)**2.d0)
            else
              cu1(i,k,j) = cu1(i,k,j) + (rnum1 - 0.5) * kick
              cu2(i,k,j) = cu2(i,k,j) + (rnum2 - 0.5) * kick
              cu3(i,k,j) = cu3(i,k,j) + (rnum3 - 0.5) * kick
            end if
          end if
        end do
      end do

    end if


    ! Apply Boundary conditions to velocity field
    call apply_BC_vel_mpi_post
    call ghost_chan_mpi


    ! Remove the divergence of the velocity field
    call rem_div_chan
    call ghost_chan_mpi

    ! Get the pressure from the poisson equation
    ! call poisson_p_chan
    ! Fix for the pressure
    ! call ghost_chan_mpi

    ! Save various statistics to keep track of the initial condition
    ! call save_stats_chan(.false.)

    return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine Add_SI_Mode
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Add exact SI Mode Perturbations


    integer i, j, k, n
    real(rkind) lambda1, lambda2, k_max, l_max, sigma
    complex(rkind) amplitude_scaled

    if (Ro_inv == 1/1.d0) then
      ! Ro = 1
      lambda1 = -19.943848405400693
      lambda2 = -13.660663098433014
      k_max = 26.922065852762007
      l_max = 17.093431565855951
      sigma = 0.765902833138420
    else if (Ro_inv == 1/10.d0) then
      ! Ro = 10
      lambda1 = -17.235267629576086
      lambda2 = -10.952082328361845
      k_max = 14.461814598675044
      l_max = 14.439573333787491
      sigma = 0.300145893642542
    else if (Ro_inv == 1/0.1d0) then
      ! Ro = 0.1
      lambda1 = -11.590816435595643
      lambda2 = -5.307631128443722
      k_max = 85.214180782654623
      l_max = 9.014376678159723
      sigma = 0.850936873842555
    else
      stop 'Exact Mode not defined for this Ro!'
    end if


    do j = 0, Nyp + 1
      do k = 0, Nzp - 1
        do i = 0, Nxm1

          amplitude_scaled = (0, 1) * (lambda1 * exp((0, 1) * lambda1 * gyf(j)) - &
                           lambda2 * exp((0, 1) * lambda2 * gyf(j))) * &
                 exp((0, -1) * k_max * gx(i))
          u1(i, k, j) = u1(i, k, j) + real(amplitude_scaled) * kick * 100

          amplitude_scaled = (0, -1) * k_max * &
                 (exp((0, 1) * lambda1 * gyf(j)) - exp((0, 1) * lambda2 * gyf(j))) * &
                 exp((0, -1) * k_max * gx(i))
          u2(i, k, j) = u2(i, k, j) + real(amplitude_scaled) * kick * 100

          amplitude_scaled = -(u1(i, k, j) * Ro_inv + u2(i, k, j)) / &
                 (sigma + nu * (k_max * k_max + l_max * l_max))
          u3(i, k, j) = u3(i, k, j) + real(amplitude_scaled) * kick * 100
        end do
      end do
    end do

    Return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine create_th_chan
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Initialize the scalar fields
    ! In this subroutine, you should initialize each scalar field for the
    ! particular problem of interest


    integer i, j, k, n

    do n = 1, N_th
      if (create_new_th(n)) then

        if (IC_type == 0) then
          ! As an example, initialize TH1 with a sine in x
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              do j = 1, Nyp
                th(i, k, j, n) = Ri(n) * sin(2.d0 * pi * gx(i) / Lx) / (4.d0 * pi**2.d0)
              end do
            end do
          end do
        else if ((IC_type == 1) .or. (IC_type == 2)) then
          ! Initialize with a linear profile using the BCs
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              if ((th_BC_Ymin(n) == 0) .and. (th_BC_Ymax(n) == 0)) then
                do j = 1, Nyp
                  if (gyf(j) <= 2.0) then
                    th(i, k, j, n) = (th_BC_Ymax_c1(n) - th_BC_Ymin_c1(n)) &
                                     * (gyf(j) + 1.) / 2.0 + th_BC_Ymin_c1(n)
                  else
                    th(i, k, j, n) = th_BC_Ymax_c1(n)
                  end if
                end do
              else if ((th_BC_Ymin(n) == 1) &
                       .and. (th_BC_Ymax(n) == 1)) then
                do j = 1, Nyp
                  ! Linear profile with slope corresponding to lower value
                  th(i, k, j, n) = th_BC_Ymin_c1(n) * gyf(j)
                end do
              else
                if (rank == 0) then
                  write (*, '("Warning, THETA INITIALIZED TO ZERO ...")')
                  write (*, '("Create an initial value in create_flow_chan")')
                end if
              end if
            end do
          end do
        else if (IC_type == 3) then
          ! Tanh vertical profile (for shear layer)
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              do j = 1, Nyp
                th(i, k, j, n) = Ri(n) * tanh(gyf(j))
              end do
            end do
          end do
        else if (IC_type == 4 .or. IC_type == 5) then
          ! Infinite Front (5 includes background sinusoidal lateral W shear)
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              th(i, k, 0, n) = 0.d0
              do j = 1, Nyp
                th(i, k, j, n) = Ri(n) * (gyf(j) - 0.5 * Ly)
              end do
            end do
          end do
        else if (IC_type == 6 .or. IC_type == 7 .or. IC_type == 8) then
          ! Finite Sloped Front
          do k = 0, Nzp - 1
            do i = 0, Nxm1
              do j = 1, Nyp
                if (Ri(n) == 0) then
                  th(i, k, j, n) = Ro_inv &
                                   * tanh((gx(i) - 0.5 * Lx) / delta) &
                                   / tanh(0.5d0 * Lx / delta) &
                                   + Ro_inv * (1.d0 - 2.d0 / Lx * gx(i))
                else
                  th(i, k, j, n) = sqrt(Ri(n)**2.d0 + Ro_inv**2.d0) &
                                   * tanh(((gx(i) - 0.5 * Lx) * Ro_inv + &
                                           (gyf(j) - 0.5 * Ly) * Ri(n)) / &
                                          (delta * sqrt(Ri(n)**2.d0 + Ro_inv**2.d0))) &
                                   + Ro_inv * (1.d0 - 2.d0 / Lx * gx(i))
                end if
              end do
            end do
          end do
        else
          write (*, '("Warning, unsupported IC_type in create_flow.")')
        end if

        call fft_xz_to_fourier(th(:, :, :, n), cth(:, :, :, n))

      end if
    end do

    return
  end






  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_flow
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    character(len=55) fname
    integer i, j, k, n

    fname = 'start.h5'
    if (rank == 0) &
      write (*, '("Reading flow from " A10)') fname

    call mpi_barrier(mpi_comm_world, ierror)
    call ReadHDF5(fname)

    ! Apply initial boundary conditions, set ghost cells
    call apply_BC_vel_mpi_post

    return
  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine save_flow(final)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    character(len=55) fname
    integer i, j, k, n
    logical final, save_pressure
    real(rkind) wall_begin

    call wall_time(wall_begin)

    if (rank == 0) &
      write (*, '("Saving Output File")')

    save_pressure = .false.
    if (final) then
      fname = 'end.h5'
      save_pressure = .true.
    else
      write (fname,'(A4, I0.6, A3)') 'out.', time_step, '.h5'
    end if

    if (num_per_dir /= 2) then
      if (rank == 0) &
        write (*, '("Saving to HDF5 only implemented for Channel.")')
      stop
    end if
    call mpi_barrier(mpi_comm_world, ierror)
    call WriteHDF5(fname, save_pressure)

    call wall_time(end_wall_time)
    if (rank == 0) &
      write (*,'("Elapsed Wall Time to Save Flow: ", ES12.5)') (end_wall_time - wall_begin)

    return
  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine save_stats(movie,final)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    logical movie, final
    real(rkind) wall_begin

    call wall_time(wall_begin)
    call save_stats_chan(movie,final)
    call wall_time(end_wall_time)
    if (rank == 0) &
      write (*,'("Elapsed Wall Time to Save Stats: ", ES13.3)') (end_wall_time - wall_begin)



    return
  end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine courant
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine sets the timestep based on the specified CFL number
  ! The subroutine should be called with the velocity in physical space

  real(rkind) dt_x, dt_y, dt_z
  real(rkind) Nmax
  integer i, j, k, n
  integer imin, jmin, kmin
  real(rkind) r_max ! Maximum fractional change in dt

  ! Set the initial dt to some arbitrary large number
  dt = 1.d0

  ! Set the timestep based on viscosity and diffusivity
  dt = min(dt, 0.5d0 * min(dx(1), dz(1))**(2.d0 * beta) / nu)
  do n = 1, N_th
    dt = min(dt, dt * nu / (nu / Pr(n)))
  end do
  ! Make sure that we capture the inertial period (for rotating flows)
  if (Ro_inv /= 0.d0) then
    dt = min(dt, 2.d0 * pi / abs((Ro_inv / delta)) / 20.d0)
  end if
  ! Make sure that we capture the buoyancy period (for stratified flows)
  do n = 1, N_th
    Nmax = sqrt(abs(th_BC_Ymin_c1(n)))
    dt = min(dt, 0.1 * 2.d0 * pi / Nmax)
  end do

  ! Make sure we capture the Geostrophic Flow-through time
  !   because the split advection component may still be at most Vg !
  !dt=min(dt,min(dx(1),dy(1))/(LY*delta/Ro_inv*dTHdX))

  ! Use the model velocity to calculate the CFL number

  if (flavor == 'Front') then
    ! Add thermal wind to velocity when calculating the CFL number
    do n = 1, N_th
      do j = jstart, jend
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            dt_x = CFL * dx(i) / abs(u1(i, k, j) - dTHdZ(n) &
                                     * gyf(j) / (Ro_inv / delta))
            dt_y = CFL * dy(j) / abs(u2(i, k, j))
            dt_z = CFL * dz(k) / abs(u3(i, k, j) + (1.d0 / (Ro_inv / delta)) &
                                     * dTHdX(n) * (gyf(j) - 0.5d0*Ly))
            dt = min(dt, dt_x, dt_y, dt_z)
          end do
        end do
      end do
    end do
  else
    do j = 1, Nyp
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          dt_x = CFL * dx(i) / abs(u1(i, k, j))
          dt_y = CFL * dy(j) / abs(u2(i, k, j))
          dt_z = CFL * dz(k) / abs(u3(i, k, j))
          dt = min(dt, dt_x, dt_y, dt_z)
        end do
      end do
    end do
  end if

  call get_minimum_mpi(dt)

  if (dt <= 0) then
    if (rank == 0) &
      write (*, '("Error: dt <= 0 in courant")')
    ! Set DELTA_T to some small default value
    delta_t = 0.00001d0
  else if (dt >= 1000.) then
    write (*, '("Warning: dt > 1000, value capped at 1000")')
    delta_t = 1000.d0
  else
    delta_t = dt
  end if

  ! Find the next event we need to meet at
  delta_t_next_event = min(time_limit, &
                           save_stats_time - time, &
                           save_flow_time - time,  &
                           time_limit - time)

  ! Try to schedule the next 2 time-steps to not be too small (keep CFL ~ 0.5)
  r_max = (1. - 0.025)
  if (delta_t_next_event < dt) then ! Directly go there (collect $200)
    delta_t = delta_t_next_event + 1.d-14
  else if (delta_t_next_event < (1.+r_max)*dt) then ! Two steps this time
    delta_t = 1./(1.+r_max)*delta_t_next_event
  else
    delta_t = dt
  end if

  ! if (time + delta_t > save_stats_time) then
  !   delta_t = save_stats_time - time + 1.d-14
  ! else if (time + delta_t > save_flow_time) then
  !   delta_t = save_flow_time - time + 1.d-14
  ! else if (time + delta_t > time_limit) then
  !   delta_t = time_limit - time
  ! end if

  if (rank == 0 .and. delta_t < 1.d-10) &
    write (*, '("Warning: dt < 1e-10, so something may be going wrong!")')


  h_bar(1) = delta_t * (8.d0 / 15.d0)
  h_bar(2) = delta_t * (2.d0 / 15.d0)
  h_bar(3) = delta_t * (5.d0 / 15.d0)

  return
end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine ghost_chan_mpi
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine is part of the MPI package for the channel flow
  ! Diablo package.
  ! Here, we define a set of ghost cells on each process
  ! the ghost cells contain information from the neighboring nodes
  ! and allow us to compute finite differences over the local gridpoints.
  ! We need to update the contents of the ghost cells at the start of
  ! each Runge-Kutta substep



  integer i, j, k, n

  ! Define the arrays that will be used for data packing.  This makes the
  ! communication between processes more efficient by only requiring one
  ! send and recieve.
  ! The communication will be done in Fourier space, so these arrays should
  ! be complex arrays to match the velocity
  ! The size of the buffer array is 0:Nkx,0:twoNkz,# of variables
  complex(rkind) ocpack(0:Nxp - 1, 0:twoNkz, 4 + N_th)
  complex(rkind) icpack(0:Nxp - 1, 0:twoNkz, 4 + N_th)

  ! If we are using more than one processor, then we need to pass data

  if (NprocY > 1) then

    ! First, Pass data up the chain to higher ranked processes

    if (rankY == 0) then
      ! If we are the lowest ranked process, then we don't need to recieve
      ! data at the lower ghost cells, these will be filled with boundary
      ! condition information
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ocpack(i, k, 1) = cu1(i, k, Nyp)
          ocpack(i, k, 2) = cu2(i, k, Nyp)
          ocpack(i, k, 3) = cu3(i, k, Nyp)
          ocpack(i, k, 4) = cp(i, k, Nyp)
          do n = 1, N_th
            ocpack(i, k, 4 + n) = cth(i, k, Nyp, n)
          end do
        end do
      end do
      ! Now, we have packed the data into a compact array, pass the data up
      call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY + 1, 1, mpi_comm_y, ierror)

      ! End if RANK=0
    else if (rankY < NprocY - 1) then
      ! Here, we are one of the middle processes and we need to pass data
      ! up and recieve data from below
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ocpack(i, k, 1) = cu1(i, k, Nyp)
          ocpack(i, k, 2) = cu2(i, k, Nyp)
          ocpack(i, k, 3) = cu3(i, k, Nyp)
          ocpack(i, k, 4) = cp(i, k, Nyp)
          do n = 1, N_th
            ocpack(i, k, 4 + n) = cth(i, k, Nyp, n)
          end do
        end do
      end do
      ! Use MPI_SENDRECV since we need to recieve and send data
      call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY + 1, 1, mpi_comm_y, ierror)

      call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY - 1, 1, mpi_comm_y, status, ierror)
      ! Now, unpack the data that we have recieved
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, 1) = icpack(i, k, 1)
          cu2(i, k, 1) = icpack(i, k, 2)
          cu3(i, k, 1) = icpack(i, k, 3)
          cp(i, k, 1) = icpack(i, k, 4)
          do n = 1, N_th
            cth(i, k, 1, n) = icpack(i, k, 4 + n)
          end do
        end do
      end do

    else
      ! Otherwise, we must be the uppermost process with RANK=Nprocs-1
      ! Here, we need to recieve data from below, but don't need to send data up
      call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY - 1, 1, mpi_comm_y, status, ierror)
      ! Unpack the data that we have recieved
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, 1) = icpack(i, k, 1)
          cu2(i, k, 1) = icpack(i, k, 2)
          cu3(i, k, 1) = icpack(i, k, 3)
          cp(i, k, 1) = icpack(i, k, 4)
          do n = 1, N_th
            cth(i, k, 1, n) = icpack(i, k, 4 + n)
          end do
        end do
      end do
    end if

    ! AT this point we have passed data up the chain
    if (rankY == NprocY - 1) then
      ! If we are the higest ranked process, then we don't need to recieve
      ! data at the upper ghost cells, these will be filled with boundary
      ! condition information
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ocpack(i, k, 1) = cu1(i, k, 2)
          ocpack(i, k, 2) = cu2(i, k, 2)
          ocpack(i, k, 3) = cu3(i, k, 2)
          ocpack(i, k, 4) = cp(i, k, 2)
          do n = 1, N_th
            ocpack(i, k, 4 + n) = cth(i, k, 2, n)
          end do
        end do
      end do
      ! Now, we have packed the data into a compact array, pass the data up
      call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY - 1, 3, mpi_comm_y, ierror)
    else if (rankY > 0) then
      ! Here, we are one of the middle processes and we need to pass data
      ! down and recieve data from above us
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ocpack(i, k, 1) = cu1(i, k, 2)
          ocpack(i, k, 2) = cu2(i, k, 2)
          ocpack(i, k, 3) = cu3(i, k, 2)
          ocpack(i, k, 4) = cp(i, k, 2)
          do n = 1, N_th
            ocpack(i, k, 4 + n) = cth(i, k, 2, n)
          end do
        end do
      end do

      call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY - 1, 3, mpi_comm_y, ierror)

      call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY + 1, 3, mpi_comm_y, status, ierror)
      ! Now, unpack the data that we have recieved
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, Nyp + 1) = icpack(i, k, 1)
          cu2(i, k, Nyp + 1) = icpack(i, k, 2)
          cu3(i, k, Nyp + 1) = icpack(i, k, 3)
          cp(i, k, Nyp + 1) = icpack(i, k, 4)
          do n = 1, N_th
            cth(i, k, Nyp + 1, n) = icpack(i, k, 4 + n)
          end do
        end do
      end do
    else
      ! Here, we must be the lowest process (RANK=0) and we need to recieve
      ! data from above
      call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY + 1, 3, mpi_comm_y, status, ierror)
      ! Unpack the data that we have recieved
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, Nyp + 1) = icpack(i, k, 1)
          cu2(i, k, Nyp + 1) = icpack(i, k, 2)
          cu3(i, k, Nyp + 1) = icpack(i, k, 3)
          cp(i, k, Nyp + 1) = icpack(i, k, 4)
          do n = 1, N_th
            cth(i, k, Nyp + 1, n) = icpack(i, k, 4 + n)
          end do
        end do
      end do
    end if

  end if

  return
end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine ghost_chan_mpi_j0
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! ghost_chan_mpi subroutine to update the j = 0 cell
  ! Needed for les output and epsilon calculation



  integer i, j, k, n

  ! Define the arrays that will be used for data packing.  This makes the
  ! communication between processes more efficient by only requiring one
  ! send and recieve.
  ! The communication will be done in Fourier space, so these arrays should
  ! be complex arrays to match the velocity
  ! The size of the buffer array is 0:Nkx,0:twoNkz,# of variables
  complex(rkind) ocpack(0:Nxp - 1, 0:twoNkz, 4 + N_th)
  complex(rkind) icpack(0:Nxp - 1, 0:twoNkz, 4 + N_th)

  ! If we are using more than one processor, then we need to pass data

  if (NprocY > 1) then

    if (rankY == 0) then
      ! If we are the lowest ranked process, then we don't need to recieve
      ! data at the lower ghost cells, these will be filled with boundary
      ! condition information
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ocpack(i, k, 1) = cu1(i, k, Nyp - 1)
          ocpack(i, k, 2) = cu2(i, k, Nyp - 1)
          ocpack(i, k, 3) = cu3(i, k, Nyp - 1)
          ocpack(i, k, 4) = cp(i, k, Nyp - 1)
          do n = 1, N_th
            ocpack(i, k, 4 + n) = cth(i, k, Nyp - 1, n)
          end do
        end do
      end do
      ! Now, we have packed the data into a compact array, pass the data up
      call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY + 1, 1, mpi_comm_y, ierror)

      ! End if RANK=0
    else if (rankY < NprocY - 1) then
      ! Here, we are one of the middle processes and we need to pass data
      ! up and recieve data from below
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ocpack(i, k, 1) = cu1(i, k, Nyp - 1)
          ocpack(i, k, 2) = cu2(i, k, Nyp - 1)
          ocpack(i, k, 3) = cu3(i, k, Nyp - 1)
          ocpack(i, k, 4) = cp(i, k, Nyp - 1)
          do n = 1, N_th
            ocpack(i, k, 4 + n) = cth(i, k, Nyp - 1, n)
          end do
        end do
      end do
      ! Use MPI_SENDRECV since we need to recieve and send data
      call mpi_send(ocpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY + 1, 1, mpi_comm_y, ierror)

      call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY - 1, 1, mpi_comm_y, status, ierror)
      ! Now, unpack the data that we have recieved
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, 0) = icpack(i, k, 1)
          cu2(i, k, 0) = icpack(i, k, 2)
          cu3(i, k, 0) = icpack(i, k, 3)
          cp(i, k, 0) = icpack(i, k, 4)
          do n = 1, N_th
            cth(i, k, 0, n) = icpack(i, k, 4 + n)
          end do
        end do
      end do

    else
      ! Otherwise, we must be the uppermost process with RANK=Nprocs-1
      ! Here, we need to recieve data from below, but don't need to send data up
      call mpi_recv(icpack, (4 + N_th) * (Nxp) * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY - 1, 1, mpi_comm_y, status, ierror)
      ! Unpack the data that we have recieved
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, 0) = icpack(i, k, 1)
          cu2(i, k, 0) = icpack(i, k, 2)
          cu3(i, k, 0) = icpack(i, k, 3)
          cp(i, k, 0) = icpack(i, k, 4)
          do n = 1, N_th
            cth(i, k, 0, n) = icpack(i, k, 4 + n)
          end do
        end do
      end do

    end if

  end if

  return
end






!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_vel_lower_post
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  ! This subroutine is called after initializing the flow
  ! It sets the appropriate boundary conditions including ghost cell values
  !  on the velocity field in Fourier space

  integer i, k

  ! Now, apply the boundary conditions depending on the type specified
  if (u_BC_Ymin == 0) then
    ! Dirichlet
    ! Start with zero
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu1(i, k, 0) = 0.d0
        cu1(i, k, 1) = 0.d0
      end do
    end do
    ! Now, set only the mean
    if (rankZ == 0) then
      cu1(0, 0, 1) = u_BC_Ymin_c1
      cu1(0, 0, 0) = u_BC_Ymin_c1
    end if
  else if (u_BC_Ymin == 1) then
    ! Neumann
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu1(i, k, 1) = cu1(i, k, 2)
        cu1(i, k, 0) = cu1(i, k, 1)
      end do
    end do
    ! Now, Apply BC to mean
    if (rankZ == 0) then
      cu1(0, 0, 1) = cu1(0, 0, 2) - dy(2) * u_BC_Ymin_c1
      cu1(0, 0, 0) = cu1(0, 0, 1) - dy(1) * u_BC_Ymin_c1
    end if
  else
    stop 'Error: u_BC_Zmin must be 0, or 1'
  end if

  if (w_BC_Ymin == 0) then
    ! Dirichlet
    ! Start with zero
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu3(i, k, 0) = 0.d0
        cu3(i, k, 1) = 0.d0
      end do
    end do
    ! Now, set only the mean
    if (rankZ == 0) then
      cu3(0, 0, 1) = w_BC_Ymin_c1
      cu3(0, 0, 0) = w_BC_Ymin_c1
    end if
  else if (w_BC_Ymin == 1) then
    ! Neumann
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu3(i, k, 1) = cu3(i, k, 2)
        cu3(i, k, 0) = cu3(i, k, 1)
      end do
    end do
    ! Now, Apply BC to mean
    if (rankZ == 0) then
      cu3(0, 0, 1) = cu3(0, 0, 2) - dy(2) * w_BC_Ymin_c1
      cu3(0, 0, 0) = cu3(0, 0, 1) - dy(1) * w_BC_Ymin_c1
    end if
  else
    stop 'Error: v_BC_Zmin must be 0, 1, or 2'
  end if

  if (v_BC_Ymin == 0) then
    ! Dirichlet
    ! Set the vertical velocity at GYF(1) (halfway between GY(2) and GY(1))
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu2(i, k, 1) = 2.d0 * v_BC_Ymin_c1 - cu2(i, k, 2)
        cu2(i, k, 0) = cu2(i, k, 1)
      end do
    end do
  else if (v_BC_Ymin == 1) then
    ! Neumann
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu2(i, k, 1) = cu2(i, k, 2)
        cu2(i, k, 0) = cu2(i, k, 1)
      end do
    end do
    if (rankZ == 0) then
      cu2(0, 0, 1) = cu2(0, 0, 2) - dyf(1) * v_BC_Ymin_c1
      cu2(0, 0, 0) = cu2(0, 0, 1) - dyf(1) * v_BC_Ymin_c1
    end if
  else
    stop 'Error: w_BC_Zmin must be 0 or 1'
  end if

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_vel_upper_post
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  ! This subroutine is called after initializing the flow
  ! It sets the appropriate boundary conditions including ghost cell values
  !  on the velocity field in Fourier space

  integer i, k

  ! Now, apply boundary conditions to the top of the domain
  if (u_BC_Ymax == 0) then
    ! Dirichlet
    ! Start with zero
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu1(i, k, Nyp) = 0.d0
        cu1(i, k, Nyp + 1) = 0.d0
      end do
    end do
    ! Now, set only the mean
    if (rankZ == 0) then
      cu1(0, 0, Nyp) = u_BC_Ymax_c1
      cu1(0, 0, Nyp + 1) = u_BC_Ymax_c1
    end if
  else if (u_BC_Ymax == 1) then
    ! Neumann
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu1(i, k, Nyp) = cu1(i, k, Nyp - 1)
        cu1(i, k, Nyp + 1) = cu1(i, k, Nyp)
      end do
    end do
    ! Now, Apply BC to mean
    if (rankZ == 0) then
      cu1(0, 0, Nyp) = cu1(0, 0, Nyp - 1) + dy(Nyp) * u_BC_Ymax_c1
      cu1(0, 0, Nyp + 1) = cu1(0, 0, Nyp) + dy(Nyp) * u_BC_Ymax_c1
    end if
  else
    stop 'Error: u_BC_Zmax must be 0 or 1'
  end if

  if (w_BC_Ymax == 0) then
    ! Dirichlet
    ! Start with zero
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu3(i, k, Nyp) = 0.d0
        cu3(i, k, Nyp + 1) = 0.d0
      end do
    end do
    ! Now, set only the mean
    if (rankZ == 0) then
      cu3(0, 0, Nyp) = w_BC_Ymax_c1
      cu3(0, 0, Nyp + 1) = w_BC_Ymax_c1
    end if
    ! Ghost cell not used
    cu3(0, 0, Nyp + 1) = 0.d0
  else if (w_BC_Ymax == 1) then
    ! Neumann
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu3(i, k, Nyp) = cu3(i, k, Nyp - 1)
        cu3(i, k, Nyp + 1) = cu3(i, k, Nyp)
      end do
    end do
    ! Now, Apply BC to mean
    if (rankZ == 0) then
      if (f_type == 4)  w_BC_Ymax_c1_transient = w_BC_Ymax_c1 + amp_omega0 * sin(omega0 * (Ro_inv/delta * time - force_start) ) ! force_start is the phase to start _down-front_ forcing
      cu3(0, 0, Nyp) = cu3(0, 0, Nyp - 1) + dy(Nyp) * w_BC_Ymax_c1_transient
      cu3(0, 0, Nyp + 1) = cu3(0, 0, Nyp) + dy(Nyp) * w_BC_Ymax_c1_transient
    end if
  else
    stop 'Error: v_BC_Zmax must be 0 or 1'
  end if

  if (v_BC_Ymax == 0) then
    ! Dirichlet
    ! Set the vertical velocity at GYF(Nyp) (halfway between GY(Nyp) and GY(Nyp+1))
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu2(0, 0, Nyp + 1) = 2.d0 * v_BC_Ymax_c1 - cu2(0, 0, Nyp)
      end do
    end do
  else if (v_BC_Ymax == 1) then
    ! Neumann
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu2(i, k, Nyp) = cu2(i, k, Nyp - 1)
        cu2(i, k, Nyp + 1) = cu2(i, k, Nyp)
      end do
    end do
    ! Now, Apply BC to mean
    if (rankZ == 0) then
      cu2(0, 0, Nyp + 1) = cu2(0, 0, Nyp) + dy(Nyp) * v_BC_Ymax_c1
    end if
  else
    stop 'Error: w_BC_Zmax must be 0 or 1'
  end if

  return
end






!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_th_phys_lower
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  ! This subroutine is called after initializing the flow
  ! It sets the appropriate boundary conditions including ghost cell values
  !  on the velocity field in Physical space

  integer i, k, n

  do n = 1, N_th

    ! Now, apply the boundary conditions depending on the type specified
    if (th_BC_Ymin(n) == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          th(i, k, 0, n) = th_BC_Ymin_c1(n)
          th(i, k, 1, n) = th_BC_Ymin_c1(n)
        end do
      end do
    else if (th_BC_Ymin(n) == 1) then
      ! Neumann
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          th(i, k, 1, n) = th(i, k, 2, n) - dy(2) * th_BC_Ymin_c1(n)
          th(i, k, 0, n) = th(i, k, 1, n) - dy(1) * th_BC_Ymin_c1(n)
        end do
      end do
    else
      stop 'Error: TH_BC_Zmin must be 0 or 1'
    end if

  end do

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_th_phys_upper
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  ! This subroutine is called after initializing the flow
  ! It sets the appropriate boundary conditions including ghost cell values
  !  on the velocity field in Fourier space

  integer i, k, n

  do n = 1, N_th

    ! Now, apply boundary conditions to the top of the domain
    if (th_BC_Ymax(n) == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          th(i, k, Nyp, n) = th_BC_Ymax_c1(n)
          th(i, k, Nyp + 1, n) = th_BC_Ymax_c1(n)
        end do
      end do
    else if (th_BC_Ymax(n) == 1) then
      ! Neumann
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          th(i, k, Nyp, n) = th(i, k, Nyp - 1, n) + dy(Nyp) * th_BC_Ymax_c1(n)
          th(i, k, Nyp + 1, n) = th(i, k, Nyp, n) + dy(Nyp) * th_BC_Ymax_c1(n)
        end do
      end do
    else
      stop 'Error: TH_BC_Zmax must be 0 or 1'
    end if

  end do

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_vel_phys_lower
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  ! This subroutine is called after initializing the flow
  ! It sets the appropriate boundary conditions including ghost cell values
  !  on the velocity field in Physical space

  integer i, k

  ! Now, apply the boundary conditions depending on the type specified
  if (u_BC_Ymin == 0) then
    ! Dirichlet
    ! Start with zero
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u1(i, k, 0) = u_BC_Ymin_c1
        u1(i, k, 1) = u_BC_Ymin_c1
      end do
    end do
  else if (u_BC_Ymin == 1) then
    ! Neumann
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u1(i, k, 1) = u1(i, k, 2) - dy(2) * u_BC_Ymin_c1
        u1(i, k, 0) = u1(i, k, 1) - dy(1) * u_BC_Ymin_c1
      end do
    end do
  else
    stop 'Error: u_BC_Zmin must be 0 or 1'
  end if

  if (w_BC_Ymin == 0) then
    ! Dirichlet
    ! Start with zero
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u3(i, k, 0) = w_BC_Ymin_c1
        u3(i, k, 1) = w_BC_Ymin_c1
      end do
    end do
  else if (w_BC_Ymin == 1) then
    ! Neumann
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u3(i, k, 1) = u3(i, k, 2) - dy(2) * w_BC_Ymin_c1
        u3(i, k, 0) = u3(i, k, 1) - dy(1) * w_BC_Ymin_c1
      end do
    end do
  else
    stop 'Error: v_BC_Zmin must be 0 or 1'
  end if

  if (v_BC_Ymin == 0) then
    ! Dirichlet
    ! Set the vertical velocity at GYF(1) (halfway between GY(2) and GY(1))
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u2(i, k, 1) = 2.d0 * v_BC_Ymin_c1 - u2(i, k, 2)
        u2(i, k, 0) = u2(i, k, 1)
      end do
    end do
  else if (v_BC_Ymin == 1) then
    ! Neumann
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u2(i, k, 1) = u2(i, k, 2) - dyf(1) * v_BC_Ymin_c1
        u2(i, k, 0) = u2(i, k, 1) - dyf(1) * v_BC_Ymin_c1
      end do
    end do
  else
    stop 'Error: w_BC_Zmin must be 0 1'
  end if

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_vel_phys_upper
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  ! This subroutine is called after initializing the flow
  ! It sets the appropriate boundary conditions including ghost cell values
  !  on the velocity field in Fourier space

  integer i, k

  ! Now, apply boundary conditions to the top of the domain
  if (u_BC_Ymax == 0) then
    ! Dirichlet
    ! Start with zero
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u1(i, k, Nyp) = u_BC_Ymax_c1
        u1(i, k, Nyp + 1) = u_BC_Ymax_c1
      end do
    end do
  else if (u_BC_Ymax == 1) then
    ! Neumann
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u1(i, k, Nyp) = u1(i, k, Nyp - 1) + dy(Nyp) * u_BC_Ymax_c1
        u1(i, k, Nyp + 1) = u1(i, k, Nyp) + dy(Nyp) * u_BC_Ymax_c1
      end do
    end do
  else
    stop 'Error: u_BC_Zmax must be 0 or 1'
  end if

  if (w_BC_Ymax == 0) then
    ! Dirichlet
    ! Start with zero
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u3(i, k, Nyp) = w_BC_Ymax_c1
        u3(i, k, Nyp + 1) = w_BC_Ymax_c1
      end do
    end do
  else if (w_BC_Ymax == 1) then
    ! Neumann
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        if (f_type == 4)  w_BC_Ymax_c1_transient = w_BC_Ymax_c1 + amp_omega0 * sin(omega0 * (Ro_inv/delta * time - force_start) ) ! force_start is the phase to start _down-front_ forcing
        u3(i, k, Nyp) = u3(i, k, Nyp - 1) + dy(Nyp) * w_BC_Ymax_c1_transient
        u3(i, k, Nyp + 1) = u3(i, k, Nyp) + dy(Nyp) * w_BC_Ymax_c1_transient
      end do
    end do
  else
    stop 'Error: v_BC_Zmax must be 0 or 1'
  end if

  if (v_BC_Ymax == 0) then
    ! Dirichlet
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u2(i, k, Nyp + 1) = v_BC_Ymax_c1
        u2(i, k, Nyp) = v_BC_Ymax_c1
      end do
    end do
  else if (v_BC_Ymax == 1) then
    ! Neumann
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        u2(i, k, Nyp + 1) = u2(i, k, Nyp) + dy(Nyp) * v_BC_Ymax_c1
      end do
    end do
  else
    stop 'Error: w_BC_Zmax must be 0 1'
  end if

  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_th_lower_post
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  ! This subroutine is called after initializing the flow
  ! It sets the appropriate boundary conditions including ghost cell values
  !  on the scalar fields in Fourier space

  integer i, k, n

  do n = 1, N_th

    ! Now, apply the boundary conditions depending on the type specified
    if (th_BC_Ymin(n) == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cth(i, k, 0, n) = 0.d0
          cth(i, k, 1, n) = 0.d0
        end do
      end do
      ! Now, set only the mean
      if (rankZ == 0) then
        cth(0, 0, 1, n) = th_BC_Ymin_c1(n)
        cth(0, 0, 0, n) = th_BC_Ymin_c1(n)
      end if
    else if (th_BC_Ymin(n) == 1) then
      ! Neumann
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cth(i, k, 1, n) = cth(i, k, 2, n)
          cth(i, k, 0, n) = cth(i, k, 1, n)
        end do
      end do
      ! Now, Apply BC to mean
      if (rankZ == 0) then
        cth(0, 0, 1, n) = cth(0, 0, 2, n) - dy(2) * th_BC_Ymin_c1(n)
        cth(0, 0, 0, n) = cth(0, 0, 1, n) - dy(1) * th_BC_Ymin_c1(n)
      end if
    else
      stop 'Error: TH_BC_Zmin must be 0, or 1'
    end if

  end do

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_th_upper_post
  !----*|--.---------.---------.---------.---------.---------.---------.-|--
  ! This subroutine is called after initializing the flow
  ! It sets the appropriate boundary conditions including ghost cell values
  !  on the scalar fields in Fourier space

  integer i, k, n

  do n = 1, N_th

    ! Now, apply boundary conditions to the top of the domain
    if (th_BC_Ymax(n) == 0) then
      ! Dirichlet
      ! Start with zero
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cth(i, k, Nyp, n) = 0.d0
          cth(i, k, Nyp + 1, n) = 0.d0
        end do
      end do
      ! Now, set only the mean
      if (rankZ == 0) then
        cth(0, 0, Nyp, n) = th_BC_Ymax_c1(n)
        cth(0, 0, Nyp + 1, n) = th_BC_Ymax_c1(n)
      end if
    else if (th_BC_Ymax(n) == 1) then
      ! Neumann
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cth(i, k, Nyp, n) = cth(i, k, Nyp - 1, n)
          cth(i, k, Nyp + 1, n) = cth(i, k, Nyp, n)
        end do
      end do
      ! Now, Apply B! to mean
      if (rankZ == 0) then
        cth(0, 0, Nyp, n) = cth(0, 0, Nyp - 1, n) + dy(Nyp) * th_BC_Ymax_c1(n)
        cth(0, 0, Nyp + 1, n) = cth(0, 0, Nyp, n) + dy(Nyp) * th_BC_Ymax_c1(n)
      end if
    else
      stop 'Error: TH_BC_Zmax must be 0 or 1'
    end if

  end do

  return
end







!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_vel_mpi_post
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions for the Poisson Eq.
  ! Note, MATL, MATD, etc. are dimensioned in header


  ! Apply Boundary conditions to velocity field
  if (rankY == 0) then
    call apply_BC_vel_lower_post
  end if
  if (rankY == NprocY - 1) then
    call apply_BC_vel_upper_post
  end if

  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_th_mpi_post
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions for the Poisson Eq.
  ! Note, MATL, MATD, etc. are dimensioned in header


  ! Apply Boundary conditions to scalar field
  if (rankY == 0) then
    call apply_BC_th_lower_post
  end if
  if (rankY == NprocY - 1) then
    call apply_BC_th_upper_post
  end if

  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_vel_phys_mpi
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions for the Poisson Eq.
  ! Note, MATL, MATD, etc. are dimensioned in header


  ! Apply Boundary conditions to velocity field
  if (rankY == 0) then
    call apply_BC_vel_phys_lower
  end if
  if (rankY == NprocY - 1) then
    call apply_BC_vel_phys_upper
  end if

  return
end



end module flow
