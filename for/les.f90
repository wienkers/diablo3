
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine les_chan
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine models the terms owing to the subgrid scale stress
  ! if the computation is to be treated as an LES not a DNS
  ! This subroutine should be called when the velocity is in fourier space
  ! in the periodic directions, on output, the velocity will be
  ! in physical space.
  ! It is assumed that the test filter and the LES filter are performed
  ! by the same operation
  ! On output S1 should contain |S| which may be used again in les_chan_th
  ! if for the subgrid scalar dissipation


  integer i, j, k, l, m, ij

  real(rkind) S1_mean(0:Nyp + 1)
  real(rkind) nu_t_mean(0:Nyp + 1)
  real(rkind) eps_sgs1_mean(0:Nyp + 1)
  real(rkind) U3_bar(0:Nyp + 1)
  real(rkind) U1_bar(0:Nyp + 1)

  real(rkind), parameter :: c_smag = 0.17d0
  real(rkind), parameter :: c_amd = 0.2887d0
  real(rkind) delta_y(0:Nyp + 1), delta_yf(0:Nyp + 1)
  real(rkind) deltax, deltay, deltaz
  real(rkind) denominator_sum

  ! Array for writing HDF5 files
  real(rkind) Diag(1:Nyp)
  character(len=20) gname

  character(len=35) fname

  ! Array to store the velocity index for each component of the strain rate tensor
  integer U_index1(6)
  integer U_index2(6)

  ! Here, alpha is the test/LES filter width ratio (Not actually used...)
  real(rkind), parameter :: alpha_sgs = 2.449d0
  ! beta is the LES/grid filter width ratio
  real(rkind), parameter :: beta_sgs = 3.d0

  ! Set the indices that are used when adding the off-diagnoal SGS stress terms
  if (rankY == NprocY - 1) then
    ! We are at the upper wall
    j1 = jstart
    j2 = Nyp - 1
  else if (rankY == 0) then
    ! We are at the lower wall
    j1 = 2
    j2 = jend
  else
    ! We are on a middle process
    j1 = 2
    j2 = Nyp
  end if

  ! First, for all models, apply boundary conditions to the velocity field
  ! (fill ghost cells) to ensure accurate calculation of gradients
  ! Apply Boundary conditions to velocity field
  call apply_BC_vel_mpi_post
  call apply_BC_th_mpi_post ! Pre-emptively for les_chan_th
  call ghost_chan_mpi
  if (flag_save_LES .and. rk_step == 1) then
    call ghost_chan_mpi_j0 ! Need u1/u3 at j = 0 for saving stats
  end if


  ! If we are using Neuman boundary conditions, over-write the values of the
  ! velocity at the ghost cells so that the LES model doesn't use the large
  ! velocity gradient
  call apply_BC_les
  call apply_BC_th_les ! Pre-emptively for les_chan_th

  if (les_model_type == 1) then
    ! Constant Smagorinsky model

    ! First, compute the rate of strain tensor S_ij
    call compute_strain_chan

    ! To remove mean shear: first get the mean horizontal velocity as in save_stats
    ! Save the mean velocity (maybe JSTART to JEND)
    if (rankZ == 0) then
      do j = 0, Nyp + 1
        U1_bar(j) = dble(cu1(0, 0, j))
        U3_bar(j) = dble(cu3(0, 0, j))
      end do
    end if
    call mpi_bcast(U1_bar, Nyp + 2, mpi_double_precision, 0, &
                   mpi_comm_z, ierror)
    call mpi_bcast(U3_bar, Nyp + 2, mpi_double_precision, 0, &
                   mpi_comm_z, ierror)

    ! Convert the velocity to physical space
    call fft_xz_to_physical(cu1, u1)
    call fft_xz_to_physical(cu2, u2)
    call fft_xz_to_physical(cu3, u3)

    ! Compute |S| at GYF points, store in S1
    ! Interpolation to GYF points is easy since by definition
    ! GYF points are exactly midway between neighboring GY points
    ! Sij4 and Sij6 have the mean shear
    ! remove the zero value U1,3bar as in compute strain
    do j = 1, Nyp ! Need j = 1 for saving stats
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          s1(i, k, j) = sqrt( &
                        2.d0 * Sij1(i, k, j)**2.d0 &
                        + 4.d0 * (0.5d0 * (Sij4(i, k, j + 1) + Sij4(i, k, j) &
                                           ! Optionally remove mean shear
                                           !- 0.5d0 * (U1_bar(j) - U1_bar(j - 1)) / dy(j) &
                                           !- 0.5d0 * (U1_bar(j + 1) - U1_bar(j)) / dy(j + 1) &
                                           ))**2.d0 &
                        + 4.d0 * Sij5(i, k, j)**2.d0 &
                        + 2.d0 * Sij2(i, k, j)**2.d0 &
                        + 4.d0 * (0.5d0 * (Sij6(i, k, j + 1) + Sij6(i, k, j) &
                                           ! Optionally remove mean shear
                                           !- 0.5d0 * (U3_bar(j) - U3_bar(j - 1)) / dy(j) &
                                           !- 0.5d0 * (U3_bar(j + 1) - U3_bar(j)) / dy(j + 1) &
                                           ))**2.d0 &
                        + 2.d0 * Sij3(i, k, j)**2.d0)
        end do
      end do
    end do

    ! Compute |S| at GY points and store in temp
    do j = 1, Nyp + 1  ! Need j = 1 for saving stats
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          temp(i, k, j) = sqrt( &
                          2.d0 * ((Sij1(i, k, j) * dyf(j - 1) + Sij1(i, k, j - 1) * dyf(j)) &
                                  / (2.d0 * dy(j)))**2.d0 &
                          + 4.d0 * (Sij4(i, k, j) &
                                    ! Optionally remove mean shear
                                    !- 0.5d0 * (U1_bar(j) - U1_bar(j - 1)) / dy(j) &
                                    )**2.d0 &
                          + 4.d0 * ((Sij5(i, k, j) * dyf(j - 1) + Sij5(i, k, j - 1) * dyf(j)) &
                                    / (2.d0 * dy(j)))**2.d0 &
                          + 2.d0 * ((Sij2(i, k, j) * dyf(j - 1) + Sij2(i, k, j - 1) * dyf(j)) &
                                    / (2.d0 * dy(j)))**2.d0 &
                          + 4.d0 * (Sij6(i, k, j) &
                                    ! Optionally remove mean shear
                                    !- 0.5d0 * (U3_bar(j) - U3_bar(j - 1)) / dy(j) &
                                    )**2.d0 &
                          + 2.d0 * ((Sij3(i, k, j) * dyf(j - 1) + Sij3(i, k, j - 1) * dyf(j)) &
                                    / (2.d0 * dy(j)))**2.d0 )
        end do
      end do
    end do

    ! Now, compute |S|*S_ij, storing in Sij
    ! First compute at GYF points
    do j = 1, Nyp  ! Need j = 1 for saving stats
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          Sij1(i, k, j) = s1(i, k, j) * Sij1(i, k, j)
          Sij5(i, k, j) = s1(i, k, j) * Sij5(i, k, j)
          ! Sij2 is added through an implicit eddy viscosity
          Sij2(i, k, j) = 0.d0
          Sij3(i, k, j) = s1(i, k, j) * Sij3(i, k, j)
        end do
      end do
    end do

    ! Now, compute at |S|*S_ij at GY points
    do j = 1, Nyp + 1  ! Need j = 1 for saving stats
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! The terms dU1/dy and dU3/dy in CSij4(:,:,:) and CSij6(:,:,:) respectively
          ! are subtracted out from Sij here since they are treated implicitly
          ! in eddy viscosity terms
          Sij4(i, k, j) = temp(i, k, j) &
                          * (Sij4(i, k, j) - 0.5 * (u1(i, k, j) - u1(i, k, j - 1)) / dy(j))
          Sij6(i, k, j) = temp(i, k, j) &
                          * (Sij6(i, k, j) - 0.5 * (u3(i, k, j) - u3(i, k, j - 1)) / dy(j))
        end do
      end do
    end do

    ! We now have |S|*S_ij stored in Sij in Physical space

    ! Compute the filter lengthscale
    ! Absorb -2.d0*C_SMAG**2.d0 here for efficiency
    do j = 1, Nyp  ! Need j = 1 for saving stats
      ! At GYF points:
      ! Constant Smagorinsky
      delta_yf(j) = -2.d0 * c_smag**2.d0 &
                    * (dx(1) * beta_sgs * dyf(j) * 1.d0 * dz(1) * beta_sgs)**(2.d0 / 3.d0)
      ! Wall Damping
      ! delta_yf(j) = -2.d0 * (0.1d0 * (1.d0 - exp((-gyf(j)/(nu*25.d0))**3.d0)))**2.d0 &
      !                 * (dx(1)*beta_sgs*dyf(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)

    end do

    do j = 1, Nyp + 1  ! Need j = 1 for saving stats
      ! At GY points:
      ! Constant Smagorinsky
      delta_y(j) = -2.d0 * c_smag**2.d0 &
                   * (dx(1) * beta_sgs * dy(j) * 2.d0 * dz(1) * beta_sgs)**(2.d0 / 3.d0)
      ! Wall Damping
      ! delta_y(j) = -2.d0 * (0.1d0 * (1.d0 - exp((-gy(j)/(nu*25.d0))**3.d0)))**2.d0 &
      !                * (dx(1)*beta_sgs*dy(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
    end do

    ! Get the eddy viscosity at GY points
    ! NU_T = (C_S^2 * DELTA^2)*|S|
    do j = 2, Nyp
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          nu_t(i, k, j) = -0.5d0 * delta_y(j) * temp(i, k, j)
        end do
      end do
    end do

    ! Now that we have calculated NU_T, set the value at the ghost cells
    ! by sharing with neighboring processes.  This subroutine also sets
    ! the value of NU_T to zero at both walls
    call ghost_les_mpi

    ! Convert the stress tensor to Fourier space

    call fft_xz_to_fourier(Sij1, CSij1)
    ! Sij2 is added through an implicit eddy viscosity
    ! call fft_xz_to_fourier(Sij2, CSij2)
    call fft_xz_to_fourier(Sij3, CSij3)
    call fft_xz_to_fourier(Sij4, CSij4)
    call fft_xz_to_fourier(Sij5, CSij5)
    call fft_xz_to_fourier(Sij6, CSij6)

    ! Now, compute TAU, store in the corresponging Sij
    do j = 1, Nyp  ! Need j = 1 for saving stats
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          CSij1(i, k, j) = delta_yf(j) * CSij1(i, k, j)
          CSij5(i, k, j) = delta_yf(j) * CSij5(i, k, j)
          ! CSij2 is added through an implicit eddy viscosity
          ! CSij2(i, k, j) = delta_yf(j) * CSij2(i, k, j)
          CSij3(i, k, j) = delta_yf(j) * CSij3(i, k, j)
        end do
      end do
    end do
    do j = 1, Nyp + 1  ! Need j = 1 for saving stats
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          CSij4(i, k, j) = delta_y(j) * CSij4(i, k, j)
          CSij6(i, k, j) = delta_y(j) * CSij6(i, k, j)
        end do
      end do
    end do



  else if ((les_model_type == 2) .or. (les_model_type == 3)) then
    ! Here, use a dynamic smagorinsky model with or without scale similar part

    stop 'Error: LES_MODEL_TYPE=2, 3 not supported yet in MPI'

  else if (les_model_type == 4) then
    ! Anisotropic Minimum Dissipation (AMD) Model (Rozema)

    ! First, compute the rate of strain tensor S_ij
    call compute_strain_chan

    ! Then compute the rate of strain tensor omega_ij
    call compute_rotation_chan

    ! Convert the velocity to physical space
    call fft_xz_to_physical(cu1, u1)
    call fft_xz_to_physical(cu2, u2)
    call fft_xz_to_physical(cu3, u3)

    ! Set filter length (based on grid size) in x,z directions
    ! Based on constant Smag code
    deltax = (dx(1) * beta_sgs)**2.d0
    deltaz = (dz(1) * beta_sgs)**2.d0


    ! Compute max{0,-delta_k*(I3-I4)}/(I1-I2) at GYF points, store in S1
    ! Interpolation to GYF points is easy since by definition
    ! GYF points are exactly midway between neighboring GY points

    do j = 1, Nyp ! Need j = 1 for saving stats
      ! Set filter length (based on grid size) in y direction
      ! Based on constant Smag code
      deltay = (dy(j) * 2.d0)**2.d0
      do k = 0, Nzp - 1
        do i = 0, Nxm1

          ! First calculate delta_k*I3 and store it in S1.
          s1(i, k, j) = &
                        deltax * Sij1(i, k, j)**3.d0 &
                      + deltay * Sij2(i, k, j)**3.d0 &
                      + deltaz * Sij3(i, k, j)**3.d0 &
              + (2.d0*deltax + deltay) * Sij1(i, k, j) &
                     * (0.5d0 * (Sij4(i, k, j + 1) + Sij4(i, k, j)))**2.d0 &
              + (2.d0*deltax + deltaz) * Sij1(i, k, j) * Sij5(i, k, j)**2.d0 &
              + (2.d0*deltay + deltax) * Sij2(i, k, j) &
                     * (0.5d0 * (Sij4(i, k, j + 1) + Sij4(i, k, j)))**2.d0 &
              + (2.d0*deltay + deltaz) * Sij2(i, k, j) &
                     * (0.5d0 * (Sij6(i, k, j + 1) + Sij6(i, k, j)))**2.d0 &
              + (2.d0*deltaz + deltax) * Sij3(i, k, j) * Sij5(i, k, j)**2.d0 &
              + (2.d0*deltaz + deltay)*Sij3(i, k, j) &
                     * (0.5d0 * (Sij6(i, k, j + 1)+Sij6(i, k, j)))**2.d0 &
              + 2.d0*(deltax + deltay + deltaz) &
                     * (0.5d0 * (Sij4(i, k, j + 1) + Sij4(i, k, j))) &
                     * Sij5(i, k, j) * (0.5d0 * (Sij6(i, k, j + 1) + Sij6(i, k, j)))

          ! Then calculate -delta_k*(I3-I4) = -s1+delta_k*I4 and store in s1.

          s1(i, k, j) = (-s1(i, k, j) &
               - deltay * Sij1(i, k, j) &
                     * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j)))**2.d0 &
               - deltaz * Sij1(i, k, j) * Oij5(i, k, j)**2.d0 &
               - deltax * Sij2(i, k, j) &
                     * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j)))**2.d0 &
               - deltaz * Sij2(i, k, j) &
                     * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j)))**2.d0 &
               - deltax * Sij3(i, k, j) * Oij5(i, k, j)**2.d0 &
               - deltay * Sij3(i, k, j) &
                     * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j)))**2.d0 &
             - 2.d0*deltaz * (0.5d0*(Sij4(i, k, j + 1) + Sij4(i, k, j))) &
                     * Oij5(i, k, j) * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j))) &
             + 2.d0*deltay * Sij5(i, k, j) &
                     * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j))) &
                     * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j))) &
             - 2.d0*deltax * (0.5d0*(Sij6(i, k, j + 1) + Sij6(i, k, j))) &
                     * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j))) * Oij5(i, k, j)  )

          if (s1(i, k, j) <= 0.0d0) then
            s1(i, k, j) = 0.0d0
          else
            s1(i, k, j) = s1(i, k, j) / (Sij1(i, k, j)**2.d0 &
                  + Sij2(i, k, j)**2.d0 &
                  + Sij3(i, k, j)**2.d0 &
                  + 2.d0 * (0.5d0*(Sij4(i, k, j + 1) + Sij4(i, k, j)))**2.d0 &
                  + 2.d0 * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j)))**2.d0 &
                  + 2.d0 * Sij5(i, k, j)**2.d0 &
                  + 2.d0 * Oij5(i, k, j)**2.d0 &
                  + 2.d0 * (0.5d0*(Sij6(i, k, j + 1) + Sij6(i, k, j)))**2.d0 &
                  + 2.d0 * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j)))**2.d0)
          end if

        end do
      end do
    end do



    ! Compute max{0,-delta_k*(I3-I4)}/(I1-I2) at GY points and store in TEMP
    do j = 2, Nyp + 1
      ! Set filter length (based on grid size) in y direction
      ! Based on constant Smag code
      deltay = (dy(j) * 2.d0)**2.d0

      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! First calculate delta_k*I3 and store it in TEMP.
          temp(i, k, j) = &
              deltax * ((Sij1(i, k, j) * dyf(j-1) + Sij1(i, k, j - 1) * dyf(j)) &
                                    / (2.d0*dy(j)))**3.d0 &
            + deltay * ((Sij2(i, k, j) * dyf(j-1) + Sij2(i, k, j - 1) * dyf(j)) &
                                    / (2.d0*dy(j)))**3.d0 &
            + deltaz * ((Sij3(i, k, j) * dyf(j-1) + Sij3(i, k, j - 1) * dyf(j)) &
                                    / (2.d0*dy(j)))**3.d0 &
            + (2.d0*deltax + deltay) &
              * ((Sij1(i, k, j) * dyf(j-1) + Sij1(i, k, j - 1) * dyf(j)) / (2.d0*dy(j))) &
              * Sij4(i, k, j)**2.d0 &
            + (2.d0*deltax + deltaz) &
              * ((Sij1(i, k, j) * dyf(j-1) + Sij1(i, k, j - 1) * dyf(j)) / (2.d0*dy(j))) &
              * ((Sij5(i, k, j) * dyf(j-1) + Sij5(i, k, j - 1) * dyf(j)) / (2.d0*dy(j)))**2.d0 &
            + (2.d0*deltay+deltax) &
              * ((Sij2(i, k, j) * dyf(j-1) + Sij2(i, k, j - 1) * dyf(j)) / (2.d0*dy(j))) &
              * Sij4(i, k, j)**2.d0 &
            + (2.d0*deltay + deltaz) &
              * ((Sij2(i, k, j) * dyf(j-1) + Sij2(i, k, j - 1) * dyf(j)) / (2.d0*dy(j))) &
              * Sij6(i, k, j)**2.d0 &
            + (2.d0*deltaz + deltax) &
              * ((Sij3(i, k, j) * dyf(j-1) + Sij3(i, k, j - 1) * dyf(j)) / (2.d0*dy(j))) &
              * ((Sij5(i, k, j) * dyf(j-1) + Sij5(i, k, j - 1) * dyf(j)) / (2.d0*dy(j)))**2.d0 &
            + (2.d0*deltaz + deltay) &
              * ((Sij3(i, k, j) * dyf(j-1) + Sij3(i, k, j - 1) * dyf(j)) / (2.d0*dy(j))) &
              * Sij6(i, k, j)**2.d0 &
            + 2.d0*(deltax + deltay + deltaz) * Sij4(i, k, j) * Sij6(i, k, j) &
              * ((Sij5(i, k, j) * dyf(j-1) + Sij5(i, k, j - 1) * dyf(j)) / (2.d0*dy(j))  )

          ! Then calculate -(I3-I4) = -TEMP+delta_k*I4 and store in TEMP.
          temp(i, k, j) = -temp(i, k, j) &
            - deltay * ((Sij1(i, k, j) * dyf(j-1) + Sij1(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j))) &
              * Oij4(i, k, j)**2.d0 &
            - deltaz * ((Sij1(i, k, j) * dyf(j-1) + Sij1(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j))) &
              * ((Oij5(i, k, j) * dyf(j-1) + Oij5(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j)))**2.d0 &
            - deltax * ((Sij2(i, k, j) * dyf(j-1) + Sij2(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j))) &
              * Oij4(i, k, j)**2.d0 &
            - deltaz * ((Sij2(i, k, j) * dyf(j-1) + Sij2(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j))) &
              * Oij6(i, k, j)**2.d0 &
            - deltax * ((Sij3(i, k, j) * dyf(j-1) + Sij3(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j))) &
              * ((Oij5(i, k, j) * dyf(j-1) + Oij5(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j)))**2.d0 &
            - deltay * ((Sij3(i, k, j) * dyf(j-1) + Sij3(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j))) &
              * Oij6(i, k, j)**2.d0 &
            - 2.d0*deltaz * Sij4(i, k, j) &
              * ((Oij5(i, k, j) * dyf(j-1) + Oij5(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j))) &
              * Oij6(i, k, j) &
            + 2.d0*deltay * ((Sij5(i, k, j)*dyf(j-1)+Sij5(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j))) &
              * Oij4(i, k, j) * Oij6(i, k, j) &
            - 2.d0*deltax * Sij6(i, k, j) * Oij4(i, k, j) &
              * ((Oij5(i, k, j) * dyf(j-1) + Oij5(i, k, j - 1) * dyf(j)) &
                                     / (2.d0*dy(j)))


          if (temp(i, k, j) <= 0.0d0) then
            temp(i, k, j) = 0.0d0
          else
            temp(i, k, j) = temp(i, k, j) /  &
                (((Sij1(i, k, j) * dyf(j-1) + Sij1(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j)))**2.d0 &
              + ((Sij2(i, k, j) * dyf(j-1) + Sij2(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j)))**2.d0 &
              + ((Sij3(i, k, j) * dyf(j-1) + Sij3(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j)))**2.d0 &
              + 2.d0*Sij4(i, k, j)**2.d0 &
              + 2.d0*Oij4(i, k, j)**2.d0 &
              + 2.d0 * ((Sij5(i, k, j) * dyf(j-1) + Sij5(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j)))**2.d0 &
              + 2.d0 * ((Oij5(i, k, j) * dyf(j-1) + Oij5(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j)))**2.d0 &
              + 2.d0*Sij6(i, k, j)**2.d0 &
              + 2.d0*Oij6(i, k, j)**2.d0 )

          end if

        end do
      end do
    end do

    ! Now, compute max{0,-delta_k*(I3-I4)}/(I1-I2)*S_ij, storing in Sij
    ! First compute at GYF points
    do j = 1, Nyp + 1
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          Sij1(i, k, j) = s1(i, k, j) * Sij1(i, k, j)
          Sij5(i, k, j) = s1(i, k, j) * Sij5(i, k, j)
          ! Sij2 is added through an implicit eddy viscosity
          Sij2(i, k, j) = 0.d0
          Sij3(i, k, j) = s1(i, k, j) * Sij3(i, k, j)
        end do
      end do
    end do


    ! Now, compute at max{0,-delta_k*(I3-I4)}/(I1-I2)*S_ij at GY points
    do j = 2, Nyp + 1
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! The terms dU1/dy and dU3/dy in CSij4(:,:,:) and CSij6(:,:,:) respectively
          ! are subtracted out from Sij here since they are treated implicitly
          ! in eddy viscosity terms
          Sij4(i, k, j) = temp(i, k, j) &
              * (Sij4(i, k, j) - 0.5*(u1(i, k, j) - u1(i, k, j - 1)) / dy(j))
          Sij6(i, k, j) = temp(i, k, j) &
              * (Sij6(i, k, j) - 0.5*(u3(i, k, j) - u3(i, k, j - 1)) / dy(j))
        end do
      end do
    end do

    ! We now have max{0,-delta_k*(I3-I4)}/(I1-I2)*S_ij stored in Sij in Physical space

    ! Compute  -2.d0*c_amd**2.d0 here for efficiency
    do j = 1, Nyp + 1
      ! At GYF points:
      ! AMD (based off Constant Smagorinsky)
      delta_yf(j) = -2.d0 * c_amd**2.d0
      !       * (dx(1)*beta_sgs*dyf(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
      ! Wall Damping
      ! delta_yf(j) = -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
      !       * (dx(1)*beta_sgs*dyf(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)

    end do

    do j = 1, Nyp + 1
      ! At GY points:
      ! AMD (based off Constant Smagorinsky)
      delta_y(j) = -2.d0 * c_amd**2.d0
      !        * (dx(1)*beta_sgs*dy(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
      ! Wall Damping
      ! delta_y(j) = -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
      !        * (dx(1)*beta_sgs*dy(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
    end do

    ! Get the eddy viscosity at GY points
    ! NU_T = (C_S^2 * DELTA^2)*|S|
    do j = 2, Nyp
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          nu_t(i, k, j) = -0.5d0 * delta_y(j) * temp(i, k, j)
        end do
      end do
    end do

    ! Now that we have calculated NU_T, set the value at the ghost cells
    ! by sharing with neighboring processes.  This subroutine also sets
    ! the value of NU_T to zero at both walls
    call ghost_les_mpi

    ! Convert the stress tensor to Fourier space

    call fft_xz_to_fourier(Sij1, CSij1)
    ! Sij2 is added through an implicit eddy viscosity
    ! call fft_xz_to_fourier(Sij2, CSij2)
    call fft_xz_to_fourier(Sij3, CSij3)
    call fft_xz_to_fourier(Sij4, CSij4)
    call fft_xz_to_fourier(Sij5, CSij5)
    call fft_xz_to_fourier(Sij6, CSij6)

    ! Now, compute TAU, store in the corresponging Sij
    do j = 1, Nyp  ! Need j = 1 for saving stats
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          CSij1(i, k, j) = delta_yf(j) * CSij1(i, k, j)
          CSij5(i, k, j) = delta_yf(j) * CSij5(i, k, j)
          ! CSij2 is added through an implicit eddy viscosity
          ! CSij2(i, k, j) = delta_yf(j) * CSij2(i, k, j)
          CSij3(i, k, j) = delta_yf(j) * CSij3(i, k, j)
        end do
      end do
    end do
    do j = 1, Nyp + 1  ! Need j = 1 for saving stats
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          CSij4(i, k, j) = delta_y(j) * CSij4(i, k, j)
          CSij6(i, k, j) = delta_y(j) * CSij6(i, k, j)
        end do
      end do
    end do





  else if (les_model_type == 5) then
    ! Anisotropic Minimum Dissipation (AMD) Model (Verstappen)

    ! First, compute the rate of strain tensor S_ij
    call compute_strain_chan

    ! Then compute the rate of strain tensor invariant SI_ij
    call compute_strain_chan_invariant(beta_sgs)

    ! Then compute the rate of rotation tensor invariant omega_ij
    call compute_rotation_chan_invariant(beta_sgs)

    ! (Pre-emptively, while in FF) Compute all the velocity and scalar  gradients
    call compute_all_gradients_chan

    ! Convert the velocity to physical space
    call fft_xz_to_physical(cu1, u1)
    call fft_xz_to_physical(cu2, u2)
    call fft_xz_to_physical(cu3, u3)

    ! Compute max{0,-(I3-I4)}/(I1-I2) at GYF points, store in S1
    ! Interpolation to GYF points is easy since by definition
    ! GYF points are exactly midway between neighboring GY points
    do j = jstart, jend
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! First calculate I3 and store it in S1.
          s1(i, k, j) = &
                   SIij1(i, k, j)**3.d0 &
                 + SIij2(i, k, j)**3.d0 &
                 + SIij3(i, k, j)**3.d0 &
               + 3.d0 * SIij1(i, k, j) * (0.5d0*(SIij4(i, k, j + 1) + SIij4(i, k, j)))**2.d0 &
               + 3.d0 * SIij1(i, k, j) * SIij5(i, k, j)**2.d0 &
               + 3.d0 * SIij2(i, k, j) * (0.5d0*(SIij4(i, k, j + 1) + SIij4(i, k, j)))**2.d0 &
               + 3.d0 * SIij2(i, k, j) * (0.5d0*(SIij6(i, k, j + 1) + SIij6(i, k, j)))**2.d0 &
               + 3.d0 * SIij3(i, k, j) * SIij5(i, k, j)**2.d0 &
               + 3.d0 * SIij3(i, k, j) * (0.5d0*(SIij6(i, k, j + 1) + SIij6(i, k, j)))**2.d0 &
               + 6.d0 * (0.5d0*(SIij4(i, k, j + 1) + SIij4(i, k, j))) * SIij5(i, k, j) &
                        * (0.5d0*(SIij6(i, k, j + 1) + SIij6(i, k, j)))

          ! Then calculate -(I3-I4) = -S1+I4 and store in S1.

          s1(i, k, j) = ( -s1(i, k, j) &
               - SIij1(i, k, j) * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j)))**2.d0 &
               - SIij1(i, k, j) * Oij5(i, k, j)**2.d0 &
               - SIij2(i, k, j) * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j)))**2.d0 &
               - SIij2(i, k, j) * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j)))**2.d0 &
               - SIij3(i, k, j) * Oij5(i, k, j)**2.d0 &
               - SIij3(i, k, j) * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j)))**2.d0 &
               - 2.d0 * (0.5d0*(SIij4(i, k, j + 1) + SIij4(i, k, j))) * Oij5(i, k, j) &
                    * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j)))  &
               + 2.d0 * SIij5(i, k, j) * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j))) &
                    * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j))) &
               - 2.d0 * (0.5d0*(SIij6(i, k, j + 1) + SIij6(i, k, j))) &
                    * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j))) * Oij5(i, k, j) )

          if (s1(i, k, j) <= 0.0d0) then
            s1(i, k, j) = 0.0d0
          else
            s1(i, k, j) = s1(i, k, j) / (SIij1(i, k, j)**2.d0 &
                        + SIij2(i, k, j)**2.d0 &
                        + SIij3(i, k, j)**2.d0 &
               + 2.d0 * (0.5d0*(SIij4(i, k, j + 1) + SIij4(i, k, j)))**2.d0 &
               + 2.d0 * (0.5d0*(Oij4(i, k, j + 1) + Oij4(i, k, j)))**2.d0 &
               + 2.d0 * SIij5(i, k, j)**2.d0 &
               + 2.d0*Oij5(i, k, j)**2.d0 &
               + 2.d0 * (0.5d0*(SIij6(i, k, j + 1) + SIij6(i, k, j)))**2.d0 &
               + 2.d0 * (0.5d0*(Oij6(i, k, j + 1) + Oij6(i, k, j)))**2.d0)
          end if

        end do
      end do
    end do



    ! Compute max{0,-(I3-I4)}/(I1-I2) at GY points and store in TEMP
    do j = 2, Nyp
      do k = 0, Nzp - 1
        do i = 0, Nxm1

          ! First calculate I3 and store it in TEMP.
          temp(i, k, j) = &
                  ((SIij1(i, k, j)*dyf(j-1) + SIij1(i, k, j - 1)*dyf(j)) &
                                 /(2.d0*dy(j)))**3.d0 &
                 + ((SIij2(i, k, j)*dyf(j-1) + SIij2(i, k, j - 1)*dyf(j)) &
                                 /(2.d0*dy(j)))**3.d0 &
                 + ((SIij3(i, k, j)*dyf(j-1) + SIij3(i, k, j - 1)*dyf(j)) &
                                 /(2.d0*dy(j)))**3.d0 &
               + 3.d0*((SIij1(i, k, j)*dyf(j-1) + SIij1(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j))) &
                    * SIij4(i, k, j)**2.d0 &
               + 3.d0*((SIij1(i, k, j)*dyf(j-1) + SIij1(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j))) &
                   * ((SIij5(i, k, j)*dyf(j-1) + SIij5(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j)))**2.d0 &
               + 3.d0*((SIij2(i, k, j)*dyf(j-1) + SIij2(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j))) &
                    * SIij4(i, k, j)**2.d0 &
               + 3.d0*((SIij2(i, k, j)*dyf(j-1) + SIij2(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j))) &
                    * SIij6(i, k, j)**2.d0 &
               + 3.d0*((SIij3(i, k, j)*dyf(j-1) + SIij3(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j))) &
                   * ((SIij5(i, k, j)*dyf(j-1) + SIij5(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j)))**2.d0 &
               + 3.d0*((SIij3(i, k, j)*dyf(j-1) + SIij3(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j))) &
                    * SIij6(i, k, j)**2.d0 &
               + 6.d0 * SIij4(i, k, j) * SIij6(i, k, j) &
                   * ((SIij5(i, k, j)*dyf(j-1) + SIij5(i, k, j - 1)*dyf(j)) &
                                   /(2.d0*dy(j)))

          ! Then calculate -(I3-I4) = -TEMP+I4 and store in TEMP.
          temp(i, k, j) = -temp(i, k, j) &
               - ((SIij1(i, k, j)*dyf(j-1) + SIij1(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j))) &
               * Oij4(i, k, j)**2.d0 &
               - ((SIij1(i, k, j)*dyf(j-1) + SIij1(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j))) &
               * ((Oij5(i, k, j)*dyf(j-1) + Oij5(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j)))**2.d0 &
               - ((SIij2(i, k, j)*dyf(j-1) + SIij2(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j))) &
               * Oij4(i, k, j)**2.d0 &
               - ((SIij2(i, k, j)*dyf(j-1) + SIij2(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j))) &
               * Oij6(i, k, j)**2.d0 &
               - ((SIij3(i, k, j)*dyf(j-1) + SIij3(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j))) &
               * ((Oij5(i, k, j)*dyf(j-1) + Oij5(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j)))**2.d0 &
               - ((SIij3(i, k, j)*dyf(j-1) + SIij3(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j))) &
               * Oij6(i, k, j)**2.d0 &
               - 2.d0 * SIij4(i, k, j) &
               * ((Oij5(i, k, j)*dyf(j-1) + Oij5(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j))) &
               * Oij6(i, k, j) &
               + 2.d0*((SIij5(i, k, j)*dyf(j-1) + SIij5(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j))) &
               * Oij4(i, k, j)*Oij6(i, k, j) &
               - 2.d0 * SIij6(i, k, j)*Oij4(i, k, j) &
               * ((Oij5(i, k, j)*dyf(j-1) + Oij5(i, k, j - 1)*dyf(j)) &
                                     /(2.d0*dy(j)))

          if (temp(i, k, j) <= 0.0d0) then
            temp(i, k, j) = 0.0d0
          else
            temp(i, k, j) = temp(i, k, j)/ &
                (((SIij1(i, k, j)*dyf(j-1) + SIij1(i, k, j - 1)*dyf(j)) &
                                 /(2.d0*dy(j)))**2.d0 &
                + ((SIij2(i, k, j)*dyf(j-1) + SIij2(i, k, j - 1)*dyf(j)) &
                                 /(2.d0*dy(j)))**2.d0 &
                + ((SIij3(i, k, j)*dyf(j-1) + SIij3(i, k, j - 1)*dyf(j)) &
                                 /(2.d0*dy(j)))**2.d0 &
               + 2.d0 * SIij4(i, k, j)**2.d0 &
               + 2.d0*Oij4(i, k, j)**2.d0 &
               + 2.d0*((SIij5(i, k, j)*dyf(j-1) + SIij5(i, k, j - 1)*dyf(j)) &
                                 /(2.d0*dy(j)))**2.d0 &
               + 2.d0*((Oij5(i, k, j)*dyf(j-1) + Oij5(i, k, j - 1)*dyf(j)) &
                                 /(2.d0*dy(j)))**2.d0 &
               + 2.d0 * SIij6(i, k, j)**2.d0 &
               + 2.d0*Oij6(i, k, j)**2.d0 )
          end if

      end do
    end do
  end do


  ! Now, compute max{0,-(I3-I4)}/(I1-I2)*S_ij, storing in Sij
  ! First compute at GYF points
  do j = jstart, jend
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        Sij1(i, k, j) = s1(i, k, j) * Sij1(i, k, j)
        Sij5(i, k, j) = s1(i, k, j) * Sij5(i, k, j)
        ! Sij2 is added through an implicit eddy viscosity
        Sij2(i, k, j) = 0.d0
        Sij3(i, k, j) = s1(i, k, j) * Sij3(i, k, j)
      end do
    end do
  end do


  ! Now, compute at max{0,-(I3-I4)}/(I1-I2)*S_ij at GY points
  do j = 2, Nyp + 1  ! Need j = 1 for saving stats
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        ! The terms dU1/dy and dU3/dy in CSij4(:,:,:) and CSij6(:,:,:) respectively
        ! are subtracted out from Sij here since they are treated implicitly
        ! in eddy viscosity terms
        Sij4(i, k, j) = temp(i, k, j) &
            * (Sij4(i, k, j) - 0.5*(u1(i, k, j) - u1(i, k, j - 1))/dy(j))
        Sij6(i, k, j) = temp(i, k, j) &
            * (Sij6(i, k, j) - 0.5*(u3(i, k, j) - u3(i, k, j - 1))/dy(j))
      end do
    end do
  end do

  ! We now have max{0,-(I3-I4)}/(I1-I2)*S_ij stored in Sij in Physical space

  ! Compute the filter lengthscale
  ! Absorb -2.d0*C_AMD**2.d0 here for efficiency
  do j = 1, Nyp + 1
    ! At GYF points:
    ! AMD (based off constant Smag)
    delta_yf(j) = -2.d0*C_AMD**2.d0 &
       * 3.d0 / (1.d0 / (dx(1)*beta_sgs)**2.d0 + 1.d0 / (dyf(j)*2.d0)**2.d0 &
          + 1.d0 / (dz(1)*beta_sgs)**2.d0)
    !           *(dx(1)*beta_sgs*dyf(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
    ! Wall Damping
    !        delta_yf(j) =
    !          -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
    !                  *(dx(1)*beta_sgs*dyf(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
  end do

  do j = 1, Nyp + 1
    ! At GY points:
    ! AMD (based off Constant Smagorinsky)
    delta_y(j) = -2.d0*C_AMD**2.d0 &
       * 3.d0 / (1.d0 / (dx(1)*beta_sgs)**2.d0 + 1.d0 / (dy(j)*2.d0)**2.d0 &
          + 1.d0 / (dz(1)*beta_sgs)**2.d0)
    !              *(dx(1)*beta_sgs*dy(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
    ! Wall Damping
    !        delta_y(j) =
    !          -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
    !                  *(dx(1)*beta_sgs*dy(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
  end do

  ! Get the eddy viscosity at GY points
  ! NU_T = (C_S^2 * DELTA^2)*|S|
  do j = 2, Nyp + 1
    do k = 0, Nzp - 1
      do i = 0, Nxm1
       nu_t(i, k, j) = -0.5d0*delta_y(j) * temp(i, k, j)
      end do
    end do
  end do

  ! Now that we have calculated NU_T, set the value at the ghost cells
  ! by sharing with neighboring processes.  This subroutine also sets
  ! the value of NU_T to zero at both walls
  call ghost_les_mpi

  ! Convert the stress tensor to Fourier space


  call fft_xz_to_fourier(Sij1, CSij1)
  ! Sij2 is added through an implicit eddy viscosity
  ! call fft_xz_to_fourier(Sij2, CSij2)
  call fft_xz_to_fourier(Sij3, CSij3)
  call fft_xz_to_fourier(Sij4, CSij4)
  call fft_xz_to_fourier(Sij5, CSij5)
  call fft_xz_to_fourier(Sij6, CSij6)

  ! Now, compute TAU, store in the corresponging Sij
  do j = 1, Nyp  ! Need j = 1 for saving stats
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        CSij1(i, k, j) = delta_yf(j) * CSij1(i, k, j)
        CSij5(i, k, j) = delta_yf(j) * CSij5(i, k, j)
        ! CSij2 is added through an implicit eddy viscosity
        ! CSij2(i, k, j) = delta_yf(j) * CSij2(i, k, j)
        CSij3(i, k, j) = delta_yf(j) * CSij3(i, k, j)
      end do
    end do
  end do
  do j = 1, Nyp + 1  ! Need j = 1 for saving stats
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        CSij4(i, k, j) = delta_y(j) * CSij4(i, k, j)
        CSij6(i, k, j) = delta_y(j) * CSij6(i, k, j)
      end do
    end do
  end do



  end if

  ! tau_ij is now contained in CSij in Fourier space




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Now, add the subgrid scale forcing to CFi
  ! (This includes the subgrid scale stress as an explicit R-K term)

  do j = j1, j2
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf1(i, k, j) = cf1(i, k, j) &
                       - cikx(i) * CSij1(i, k, j) &
                       - (CSij4(i, k, j + 1) - CSij4(i, k, j)) / dyf(j) &
                       - cikz(k) * CSij5(i, k, j)
        cf3(i, k, j) = cf3(i, k, j) &
                       - cikx(i) * CSij5(i, k, j) &
                       - (CSij6(i, k, j + 1) - CSij6(i, k, j)) / dyf(j) &
                       - cikz(k) * CSij3(i, k, j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf2(i, k, j) = cf2(i, k, j) &
                      - cikx(i) * CSij4(i, k, j) &
        ! Sij2 is added through an implict eddy viscosity
        !             - (CSij2(i, k, j) - CSij2(i, k, j - 1)) / dy(j) &
                      - cikz(k) * CSij6(i, k, j)
      end do
    end do
  end do

  ! Add back the Cross-terms using an explicit treatment
  !   These are the ones subtracted out above (in Sij4 and Sij6)
  !     (but which won't be treated in the C-N solve in solver.f90)
  do k = 0, Nzp - 1
    do i = 0, Nxm1
      do j = 2, Nyp
        temp(i, k, j) = nu_t(i, k, j) * (u1(i, k, j) - u1(i, k, j - 1)) / dy(j)
      end do
    end do
  end do
  call fft_xz_to_fourier(temp, ctemp)
  do k = 0, twoNkz
    do i = 0, Nxp - 1
      do j = 2, Nyp
        cf2(i, k, j) = cf2(i, k, j) + cikx(i) * ctemp(i, k, j)
      end do
    end do
  end do
  do k = 0, Nzp - 1
    do i = 0, Nxm1
      do j = 2, Nyp
        temp(i, k, j) = nu_t(i, k, j) * (u3(i, k, j) - u3(i, k, j - 1)) / dy(j)
      end do
    end do
  end do
  call fft_xz_to_fourier(temp, ctemp)
  do k = 0, twoNkz
    do i = 0, Nxp - 1
      do j = 2, Nyp
        cf2(i, k, j) = cf2(i, k, j) + cikz(k) * ctemp(i, k, j)
      end do
    end do
  end do



  ! Periodically, output mean quantities
  if (flag_save_LES .and. rk_step == 1) then
    !flag_save_LES = .false. ! Set this false in les_chan_th !
    ! NOTE eps_sgs1 now is consistently evaluated at GY points

    ! Get plane averages
    do j = 0, Nyp + 1
      S1_mean(j) = 0.d0
      nu_t_mean(j) = 0.d0
      do i = 0, Nxm1
        do k = 0, Nzp - 1
          S1_mean(j) = S1_mean(j) + s1(i, k, j)
          nu_t_mean(j) = nu_t_mean(j) + nu_t(i, k, j)
        end do
      end do
    end do

    eps_sgs1_mean = 0.d0 ! NOTE: eps_sgs1 is actually negative of the convention
    ! compute part of sgs term
    ! U1*F1+U2*F2+U3*F3
    ! compute part of CF1, then interpolate onto GY
    do j = 1, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ctemp(i, k, j) = &
            - cikx(i) * CSij1(i, k, j) &
            - (CSij4(i, k, j + 1) - CSij4(i, k, j)) / dyf(j) &
            - cikz(k) * CSij5(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_physical(ctemp, temp)
    do j = 2, Nyp
      do i = 0, Nxm1
        do k = 0, Nzp - 1
          eps_sgs1_mean(j) = eps_sgs1_mean(j) + &
                              0.5d0 * (u1(i, k, j) * temp(i, k, j) &
                                     + u1(i, k, j - 1) * temp(i, k, j - 1))
        end do
      end do
    end do
    ! compute part of CF3, then interpolate onto GY
    do j = 1, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ctemp(i, k, j) = &
            - cikx(i) * CSij5(i, k, j) &
            - (CSij6(i, k, j + 1) - CSij6(i, k, j)) / dyf(j) &
            - cikz(k) * CSij3(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_physical(ctemp, temp)
    do j = 2, Nyp
      do i = 0, Nxm1
        do k = 0, Nzp - 1
          eps_sgs1_mean(j) = eps_sgs1_mean(j) + &
                              0.5d0 * (u3(i, k, j) * temp(i, k, j) &
                                     + u3(i, k, j - 1) * temp(i, k, j - 1))
        end do
      end do
    end do
    ! compute part of CF2, on GY
    do j = 2, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ctemp(i, k, j) = &
            - cikx(i) * CSij4(i, k, j) &
            ! Sij2 is added through an implict eddy viscosity
            ! - (CSij2(i, k, j) - CSij2(i, k, j - 1)) / dy(j)
            - cikz(k) * CSij6(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_physical(ctemp, temp)
    do j = 2, Nyp
      do i = 0, Nxm1
        do k = 0, Nzp - 1
          eps_sgs1_mean(j) = eps_sgs1_mean(j) + temp(i, k, j) * u2(i, k, j)
        end do
      end do
    end do
    ! Add back the Cross-terms (previously subtracted from Sij, thus tau_ij)
    ! using an explicit treatment
    ! compute second part of CF2, on GY
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        do j = 2, Nyp
          temp(i, k, j) = nu_t(i, k, j) * (u1(i, k, j) - u1(i, k, j - 1)) / dy(j)
        end do
      end do
    end do
    call fft_xz_to_fourier(temp, ctemp)
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        do j = 2, Nyp
          ctemp(i, k, j) = cikx(i) * ctemp(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_physical(ctemp, temp)
    do j = 2, Nyp
      do i = 0, Nxm1
        do k = 0, Nzp - 1
          eps_sgs1_mean(j) = eps_sgs1_mean(j) + temp(i, k, j) * u2(i, k, j)
        end do
      end do
    end do
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        do j = 2, Nyp
          temp(i, k, j) = nu_t(i, k, j) * (u3(i, k, j) - u3(i, k, j - 1)) / dy(j)
        end do
      end do
    end do
    call fft_xz_to_fourier(temp, ctemp)
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        do j = 2, Nyp
          ctemp(i, k, j) = cikz(k) * ctemp(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_physical(ctemp, temp)
    do j = 2, Nyp
      do i = 0, Nxm1
        do k = 0, Nzp - 1
          eps_sgs1_mean(j) = eps_sgs1_mean(j) + temp(i, k, j) * u2(i, k, j)
        end do
      end do
    end do

    call mpi_allreduce(mpi_in_place, S1_mean, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, nu_t_mean, Nyp + 2 &
                       , mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, eps_sgs1_mean, Nyp + 2 &
                       , mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)

    do j = 0, Nyp + 1
      S1_mean(j) = S1_mean(j) / float(Nx * Nz)
      nu_t_mean(j) = nu_t_mean(j) / float(Nx * Nz)
      eps_sgs1_mean(j) = eps_sgs1_mean(j) / float(Nx * Nz)
    end do

    fname = 'mean.h5'

    if (rankZ == 0) then

      gname = 'nu_sgs'
      Diag = nu_t_mean(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'eps_sgs1'
      Diag = eps_sgs1_mean(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

    end if

  end if

  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine les_chan_th
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine models the subgrid scale diffusive terms
  ! Velocity should already be in physical space (from les_chan)

  integer i, j, k, l, m, ij, n

  real(rkind) kappa_t_mean(0:Nyp + 1)

  real(rkind), parameter :: c_amd = 0.2887d0
  real(rkind) delta_y(0:Nyp + 1), delta_yf(0:Nyp + 1)
  real(rkind) deltax, deltay, deltaz
  real(rkind) denominator_sum

  ! Array for writing HDF5 files
  real(rkind) Diag(1:Nyp)
  character(len=20) gname

  character(len=35) fname

  ! Array to store the velocity index for each component of the strain rate tensor
  integer U_index1(6)
  integer U_index2(6)

  ! Here, alpha is the test/LES filter width ratio
  real(rkind), parameter :: alpha_sgs = 2.449d0
  ! beta is the LES/grid filter width ratio
  real(rkind), parameter :: beta_sgs = 3.d0

  ! First, for all models, apply boundary conditions to the velocity field
  ! (fill ghost cells) to ensure accurate calculation of gradients
  ! Apply Boundary conditions to velocity field
  !call apply_BC_vel_mpi_post
  !call apply_BC_th_mpi_post ! Already done when called les_chan
  !call ghost_chan_mpi


  ! If we are using Neuman boundary conditions, over-write the values of the
  ! velocity at the ghost cells so that the LES model doesn't use the large
  ! velocity gradient
  ! Does not yet include Neuman boundary condition for scalar
  !call apply_BC_les ! Already done when called les_chan
  !call apply_BC_th_les ! Already done when called les_chan



  if (les_model_type == 1 .or. les_model_type == 4) then
    ! Constant Smagorinsky model
    !  or  Anisotropic Minimum Dissipation (AMD) Model (Rozema)

    ! APPLY constant SGS Prandlt number
    do n = 1, N_th
      do j = 1, Nyp + 1
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            kappa_t(i, k, j, n) = 1.d0 * nu_t(i, k, j)
          end do
        end do
      end do
    end do


    ! Add the horizontal diffusive terms using explicit timestepping

    do n = 1, N_th

    do j = 1, Nyp + 1
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cs1(i, k, j) = cikx(i) * cth(i, k, j, n)
        end do
      end do
    end do
    call fft_xz_to_physical(cs1, s1)
    do j = 1, Nyp + 1
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          s1(i, k, j) = kappa_t(i, k, j, n) * s1(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_fourier(s1, cs1)
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cfth(i, k, j, n) = cfth(i, k, j, n) + cikx(i)*cs1(i, k, j)
        end do
      end do
    end do

    do j = 1, Nyp + 1
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cs1(i, k, j) = cikz(k) * cth(i, k, j, n)
        end do
      end do
    end do
    call fft_xz_to_physical(cs1, s1)
    do j = 1, Nyp + 1
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          s1(i, k, j) = kappa_t(i, k, j, n) * s1(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_fourier(s1, cs1)
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cfth(i, k, j, n) = cfth(i, k, j, n) + cikz(k) * cs1(i, k, j)
        end do
      end do
    end do
    end do ! end do n

    ! Now, convert TH to physical space for calculation of nonlinear terms
    do n = 1, N_th
      call fft_xz_to_physical(cth(:, :, :, n), th(:, :, :, n))
    end do




  else if ((les_model_type == 2) .or. (les_model_type == 3)) then
    ! Here, use a dynamic smagorinsky model with or without scale similar part

    stop 'Error: LES_MODEL_TYPE=2, 3 not supported yet in MPI'


  else if (les_model_type == 5) then
    ! Anisotropic Minimum Dissipation (AMD) Model (Verstappen)


    ! Compute all the velocity and scalar  gradients (Now done in les_chan)
    !call compute_all_gradients_chan

    ! Convert the scalar to physical space
    do n = 1, N_th
      call fft_xz_to_physical(cth(:, :, :, n), th(:, :, :, n))
    end do

    deltax = (dx(1)*beta_sgs)
    deltaz = (dz(1)*beta_sgs)

    do n = 1, N_th
      ! First compute at GYF points
      do j = jstart_th(n), jend_th(n)
        ! Set filter length (based on grid size) in y direction
        ! Based on Constant Smagorinsky code
        deltay = (dyf(j) * 2.d0)
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            ! First calculate numerator and store it in s1_th.
            s1_th(i, k, j) = -( &
                 deltax**2 * du1dx(i, k, j) * dthetadx(i, k, j, n)**2 &
               + deltay**2 * du2dy(i, k, j) &
                 * (0.5d0 * (dthetady(i, k, j + 1, n) + dthetady(i, k, j, n)))**2 &
               + deltaz**2 * du3dz(i, k, j) * dthetadz(i, k, j, n)**2 &
               + (deltay**2 * 0.5d0 * (du1dy(i, k, j + 1) + du1dy(i, k, j)) &
                 + deltax**2 * 0.5d0 * (du2dx(i, k, j + 1) + du2dx(i, k, j))) &
                 * dthetadx(i, k, j, n) &
                 * 0.5d0 * (dthetady(i, k, j + 1, n) + dthetady(i, k, j, n)) &
               + (deltaz**2 * du1dz(i, k, j) + deltax**2 * du3dx(i, k, j)) &
                 * dthetadx(i, k, j, n) &
                 * dthetadz(i, k, j, n) &
               + (deltaz**2 * 0.5d0 * (du2dz(i, k, j + 1) + du2dz(i, k, j)) &
                 + deltay**2 * 0.5d0 * (du3dy(i, k, j + 1) + du3dy(i, k, j))) &
                 * 0.5d0 * (dthetady(i, k, j + 1, n) + dthetady(i, k, j, n)) &
                 * dthetadz(i, k, j, n) &
                )

            if (s1_th(i, k, j) <= 0.0d0) then
              s1_th(i, k, j) = 0.0d0
            else
              s1_th(i, k, j) = s1_th(i, k, j) /   &
                ((deltax * dthetadx(i, k, j, n))**2 &
              + (deltay * 0.5d0 * (dthetady(i, k, j + 1, n) + dthetady(i, k, j, n)))**2 &
              + (deltaz * dthetadz(i, k, j, n))**2)

            end if

          end do
        end do
      end do

      ! Compute kappa_e at GY points and store in temp_th
      do j = 2, Nyp + 1
        ! Set filter length (based on grid size) in y direction
        ! Based on constant Smag code
        deltay = (dy(j)*2.d0)
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            ! First calculate numerator and store it in temp_th.
            temp_th(i, k, j) = -(  &
                 deltax**2 * (du1dx(i, k, j) * dyf(j-1) + du1dx(i, k, j - 1) * dyf(j)) &
                                      / (2.d0*dy(j)) &
               * ((dthetadx(i, k, j, n) * dyf(j-1) + dthetadx(i, k, j - 1, n) * dyf(j)) &
                                      / (2.d0*dy(j)))**2 &
               + deltay**2*(du2dy(i, k, j) * dyf(j-1) + du2dy(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j)) &
                 * (dthetady(i, k, j, n))**2 &
               + deltaz**2*(du3dz(i, k, j)*dyf(j-1)+du3dz(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j)) &
                 * ((dthetadz(i, k, j, n) * dyf(j-1) + dthetadz(i, k, j - 1, n) * dyf(j)) &
                                       / (2.d0*dy(j)))**2 &
               + (deltay**2 * du1dy(i, k, j) &
                 + deltax**2 * du2dx(i, k, j)) &
                 * (dthetadx(i, k, j, n) * dyf(j-1) + dthetadx(i, k, j - 1, n) * dyf(j)) &
                                       / (2.d0*dy(j)) &
                 * (dthetady(i, k, j, n)) &
               + (deltaz**2*(du1dz(i, k, j) * dyf(j-1) + du1dz(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j)) &
                 + deltax**2*(du3dx(i, k, j) * dyf(j-1) + du3dx(i, k, j - 1) * dyf(j)) &
                                       / (2.d0*dy(j))) &
                 * (dthetadx(i, k, j, n) * dyf(j-1) + dthetadx(i, k, j - 1, n) * dyf(j)) &
                                       / (2.d0*dy(j)) &
                 * (dthetadz(i, k, j, n) * dyf(j-1) + dthetadz(i, k, j - 1, n) * dyf(j)) &
                                       / (2.d0*dy(j)) &
               + (deltaz**2 * du2dz(i, k, j) &
                 + deltay**2 * du3dy(i, k, j)) &
                 * (dthetady(i, k, j, n)) &
                 * (dthetadz(i, k, j, n)*dyf(j-1) + dthetadz(i, k, j - 1, n) * dyf(j)) &
                                       / (2.d0*dy(j)) &
                 )

            if (temp_th(i, k, j) <= 0.0d0) then
              temp_th(i, k, j) = 0.0d0
            else
              temp_th(i, k, j) = temp_th(i, k, j) /  &
                ((deltax * (dthetadx(i, k, j, n) * dyf(j-1) + dthetadx(i, k, j - 1, n) * dyf(j)) &
                                           / (2.d0*dy(j)))**2 &
              + (deltay * dthetady(i, k, j, n))**2 &
              + (deltaz * (dthetadz(i, k, j, n) * dyf(j-1) + dthetadz(i, k, j - 1, n) * dyf(j)) &
                                           / (2.d0*dy(j)))**2)
            end if

          end do
        end do
      end do

      ! Now, compute s1_th*dthetadx_i, storing in dthetadx_i
      ! Need only to compute at GYF points
      do j = jstart_th(n), jend_th(n)
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            dthetadx(i, k, j, n) = s1_th(i, k, j) * dthetadx(i, k, j, n)
            ! dthetady is added through an implicit eddy diffusivity
            dthetadz(i, k, j, n) = s1_th(i, k, j) * dthetadz(i, k, j, n)
          end do
        end do
      end do

      ! Compute the filter lengthscale
      ! Absorb -2.d0*C_AMD**2.d0 here for efficiency
      do j = 1, Nyp + 1
        ! At GYF points:
        ! AMD (based off constant Smag)
        delta_yf(j) = -2.d0 * C_AMD**2.d0 &
            * 3.d0 / (1.d0 / (dx(1) * beta_sgs)**2.d0 + 1.d0 / (dyf(j)*2.d0)**2.d0 &
              + 1.d0 / (dz(1) * beta_sgs)**2.d0)
        !           * (dx(1)*beta_sgs*dyf(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
        ! Wall Damping
        !        delta_yf(J) =
        !          -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
        !                  * (dx(1)*beta_sgs*dyf(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
      end do

      do j = 1, Nyp + 1
        ! At GY points:
        ! AMD (based off Constant Smagorinsky)
        delta_y(j) = -2.d0 * C_AMD**2.d0 &
           * 3.d0 / (1.d0 / (dx(1) * beta_sgs)**2.d0 + 1.d0 / (dy(j)*2.d0)**2.d0 &
              + 1.d0 / (dz(1)*beta_sgs)**2.d0)
        !              * (dx(1)*beta_sgs*dy(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
        ! Wall Damping
        !        delta_y(J) =
        !          -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
        !                  * (dx(1)*beta_sgs*dy(j)*2.d0*dz(1)*beta_sgs)**(2.d0/3.d0)
      end do

      ! Get the eddy diffusivity at GY points
      do j = 2, Nyp
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            kappa_t(i, k, j, n) = -0.5d0 * delta_y(j) * temp_th(i, k, j)
          end do
        end do
      end do

      ! Now that we have calculated kappa_t, set the value at the ghost cells
      ! by sharing with neighboring processes.  This subroutine also sets
      ! the value of kappa_t to zero at both walls
      call ghost_les_th_mpi

      ! Convert the scalar flux tensor to Fourier space
      call fft_xz_to_fourier(dthetadx(:, :, :, n), Cdthetadx(:, :, :, n))
      ! dthetady is added through an implicit eddy diffusivity
      call fft_xz_to_fourier(dthetadz(:, :, :, n), Cdthetadz(:, :, :, n))

      ! Now, compute TAU, store in the corresponging dthetadx_i
      do j = 1, Nyp + 1
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            Cdthetadx(i, k, j, n) = 0.5d0 * delta_yf(j) * Cdthetadx(i, k, j, n)
            ! cthetady(:,:,:,N) is added through an implicit eddy diffusivity
            ! Cdthetady(i, k, j, n) = delta_yf(J)*Cdthetady(i, k, j, n)
            Cdthetadz(i, k, j, n) = 0.5d0 * delta_yf(j) * Cdthetadz(i, k, j, n)
          end do
        end do
      end do

      ! Add to forcing on RHS of th equation (RK)
      do j = jstart_th(n), jend_th(n)
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cfth(i, k, j, n) = cfth(i, k, j, n) &
                      - cikx(i) * Cdthetadx(i, k, j, n) &
                ! dthetady is added through implicit eddy diffusivity
                      - cikz(k) * Cdthetadz(i, k, j, n)
          end do
        end do
      end do

    end do

  end if


  ! Periodically, output mean quantities
  if (flag_save_LES .and. rk_step == 1) then
    flag_save_LES = .false.

    ! Get plane averages
    do j = 0, Nyp + 1
      kappa_t_mean(j) = 0.d0
      do i = 0, Nxm1
        do k = 0, Nzp - 1
          kappa_t_mean(j) = kappa_t_mean(j) + kappa_t(i, k, j, 1)
        end do
      end do
    end do


    call mpi_allreduce(mpi_in_place, kappa_t_mean, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)

    do j = 0, Nyp + 1
      kappa_t_mean(j) = kappa_t_mean(j) / float(Nx * Nz)
    end do

    fname = 'mean.h5'

    if (rankZ == 0) then

      gname = 'kappa_sgs'
      Diag = kappa_t_mean(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

    end if

  end if

  return
end






!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine compute_strain_chan
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine computes S_ij for the filtered velocity field
  ! The input velocity field should be in fourier space in the periodic
  ! directions.
  ! For use in the LES model in channel flow (2 periodic directions)

  integer i, j, k, ij

  do j = 0, Nyp + 1 ! Need j = 0 for saving stats
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        CSij1(i, k, j) = cikx(i) * cu1(i, k, j)
        CSij3(i, k, j) = cikz(k) * cu3(i, k, j)
        CSij5(i, k, j) = 0.5d0 * (cikz(k) * cu1(i, k, j) &
                                  + cikx(i) * cu3(i, k, j))
        if (j /= Nyp + 1) then
          CSij2(i, k, j) = (cu2(i, k, j + 1) - cu2(i, k, j)) / dyf(j)
        end if
      end do
    end do
  end do
  do j = 1, Nyp + 1 ! Need j = 1 for saving stats
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        CSij4(i, k, j) = 0.5d0 * ((cu1(i, k, j) - cu1(i, k, j - 1)) / dy(j) &
                                  + cikx(i) * cu2(i, k, j))
        CSij6(i, k, j) = 0.5d0 * (cikz(k) * cu2(i, k, j) &
                                  + (cu3(i, k, j) - cu3(i, k, j - 1)) / dy(j))
      end do
    end do
  end do

  ! Need to communicate down CSij2 to each j = Nyp + 1
  call ghost_CS2_mpi(CSij2)


  call fft_xz_to_physical(CSij1, Sij1)
  call fft_xz_to_physical(CSij2, Sij2)
  call fft_xz_to_physical(CSij3, Sij3)
  call fft_xz_to_physical(CSij4, Sij4)
  call fft_xz_to_physical(CSij5, Sij5)
  call fft_xz_to_physical(CSij6, Sij6)

  ! We now have S_ij in Physical space

  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine compute_strain_chan_invariant(beta_sgs)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine computes S_ij for the filtered velocity field
  ! The input velocity field should be in fourier space in the periodic
  ! directions.
  ! For use in the LES model in channel flow (2 periodic directions)

  integer i, j, k, ij
  real(rkind) deltax, deltay, deltaz
  real(rkind) beta_sgs

  ! Set filter length (based on grid size) in x,z directions
  ! Based on constant Smag code
  deltax = (dx(1) * beta_sgs)
  deltaz = (dz(1) * beta_sgs)

  do j = 0, Nyp + 1 ! Need j = 0 for saving stats
    ! Set filter length (based on grid size) in y direction
    ! Based on constant Smag code
    deltay = (dyf(j)*2.d0)
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        CSIij1(i, k, j) = cikx(i) * cu1(i, k, j)
        CSIij3(i, k, j) = cikz(k)*cu3(i, k, j)
        CSIij5(i, k, j) = 0.5d0*(cikz(k) * cu1(i, k, j) &
                           * (deltaz / deltax) &
                           + cikx(i) * cu3(i, k, j) &
                             * (deltax / deltaz) )
        if (j /= Nyp + 1) then
          CSIij2(i, k, j) = (cu2(i, k, j + 1) - cu2(i, k, j)) / dyf(j)
        end if
      end do
    end do
  end do
  do j = 1, Nyp + 1
    deltay = (dy(j)*2.d0)
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        CSIij4(i, k, j) = 0.5d0*((cu1(i, k, j) - cu1(i, k, j - 1))/dy(j) &
                            * (deltay / deltax) &
                            + cikx(i) * cu2(i, k, j) &
                             * (deltax / deltay) )
        CSIij6(i, k, j) = 0.5d0*(cikz(k) * cu2(i, k, j) &
                            * (deltaz/deltay) &
                            + (cu3(i, k, j) - cu3(i, k, j - 1))/dy(j) &
                             * (deltay / deltaz) )
      end do
    end do
  end do

  ! Need to communicate down CSIij2 to each j = Nyp + 1
  call ghost_CS2_mpi(CSIij2)

  call fft_xz_to_physical(CSIij1, SIij1)
  call fft_xz_to_physical(CSIij2, SIij2)
  call fft_xz_to_physical(CSIij3, SIij3)
  call fft_xz_to_physical(CSIij4, SIij4)
  call fft_xz_to_physical(CSIij5, SIij5)
  call fft_xz_to_physical(CSIij6, SIij6)

  ! We now have S_ij in Physical space

  return
end





!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine compute_rotation_chan
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine computes Oij (omega_ij) for the filtered velocity field
  ! The input velocity field should be in fourier space in the periodic directions.
  ! For use in the LES model in channel flow (2 periodic directions)

  integer i, j, k, ij

  do j = 0, Nyp + 1 ! Need j = 0 for saving stats
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        ! omega_ik, this is anti-cyclic permutation
        COij5(i, k, j) = 0.5d0*(cikz(k) * cu1(i, k, j) - cikx(i) * cu3(i, k, j))
      end do
    end do
  end do
  do j = 1, Nyp + 1 ! Need j = 1 for saving stats
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        ! omega_ij, cyclic permutation
        COij4(i, k, j) = 0.5d0*((cu1(i, k, j) - cu1(i, k, j - 1)) / dy(j) - cikx(i) * cu2(i, k, j) )
        ! omega_jk, cyclic permutation
        COij6(i, k, j) = 0.5d0*(cikz(k) * cu2(i, k, j) - (cu3(i, k, j) - cu3(i, k, j - 1)) / dy(j) )
      end do
    end do
  end do


  call fft_xz_to_physical(COij4, Oij4)
  call fft_xz_to_physical(COij5, Oij5)
  call fft_xz_to_physical(COij6, Oij6)

  ! We now have omega_ij in Physical space
  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine compute_rotation_chan_invariant(beta_sgs)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine computes Oij (omega_ij) for the filtered velocity field
  ! The input velocity field should be in fourier space in the periodic directions.
  ! For use in the LES model in channel flow (2 periodic directions)

  integer i, j, k, ij
  real(rkind) deltax, deltay, deltaz
  real(rkind) beta_sgs

  ! Set filter length (based on grid size) in x,z directions
  ! Based on constant Smag code
  deltax = (dx(1) * beta_sgs)
  deltaz = (dz(1) * beta_sgs)

  do j = 0, Nyp + 1 ! Need j = 0 for saving stats
    ! Set filter length (based on grid size) in y direction
    ! Based on constant Smag code
    deltay = (dyf(j)*2.d0)
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        COij5(i, k, j) = 0.5d0*(cikz(k) * cu1(i, k, j) &
                          * (deltaz / deltax) &
                          - cikx(i) * cu3(i, k, j) &
                          * (deltax / deltaz) )
                          ! omega_ik, this is anti-cyclic permutation
      end do
    end do
  end do
  do j = 1, Nyp + 1
    deltay = (dy(j)*2.d0)
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        COij4(i, k, j) = 0.5d0*((cu1(i, k, j) - cu1(i, k, j - 1))/dy(j) &
                          * (deltay / deltax) &
                          - cikx(i) * cu2(i, k, j) &
                          * (deltax / deltay) )
                          ! omega_ij, cyclic permutation
        COij6(i, k, j) = 0.5d0*(cikz(k) * cu2(i, k, j) &
                          * (deltaz / deltay) &
                          - (cu3(i, k, j) - cu3(i, k, j - 1))/dy(j) &
                          * (deltay / deltaz) )
                          ! omega_jk, cyclic permutation
      end do
    end do
  end do

  call fft_xz_to_physical(COij4, Oij4)
  call fft_xz_to_physical(COij5, Oij5)
  call fft_xz_to_physical(COij6, Oij6)

  ! We now have omega_ij in Physical space

  return
end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine compute_all_gradients_chan
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine computes all gradients for the filtered velocity field
  ! and scalar field (testing with one scalar field only).
  ! The input velocity + scalar field should be in fourier space in the periodic
  ! directions.
  ! For use in the LES model in channel flow (2 periodic directions)

  integer i, j, k, ij, n

  do j = 0, Nyp + 1 ! Need j = 0 for saving stats
    do k = 0, twoNkz
      do i = 0, Nxp - 1

        Cdu1dx(i, k, j) = cikx(i) * cu1(i, k, j)
        Cdu1dz(i, k, j) = cikz(k) * cu1(i, k, j)

        Cdu2dx(i, k, j) = cikx(i) * cu2(i, k, j)
        Cdu2dy(i, k, j) = (cu2(i, k, j + 1) - cu2(i, k, j))/dyf(j)
        Cdu2dz(i, k, j) = cikz(k) * cu2(i, k, j)

        Cdu3dx(i, k, j) = cikx(i) * cu3(i, k, j)
        Cdu3dz(i, k, j) = cikz(k) * cu3(i, k, j)

      end do
    end do
  end do

  do j = 1, Nyp + 1 ! Need j = 1 for saving stats
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        Cdu1dy(i, k, j) = (cu1(i, k, j) - cu1(i, k, j - 1)) / dy(j)
        Cdu3dy(i, k, j) = (cu3(i, k, j) - cu3(i, k, j - 1)) / dy(j)
      end do
    end do
  end do


  do n = 1, N_th
    do j = 1, Nyp + 1 ! Need j = 1 for saving stats
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          Cdthetadx(i, k, j, n) = cikx(i) * cth(i, k, j, n)
          Cdthetadz(i, k, j, n) = cikz(k) * cth(i, k, j, n)
        end do
      end do
    end do
    do j = 1, Nyp + 1 ! Need j = 1 for saving stats
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          Cdthetady(i, k, j, n) = (cth(i, k, j, n) - cth(i, k, j - 1, n)) / dy(j)
        end do
      end do
    end do

  end do


  call fft_xz_to_physical(Cdu1dx,du1dx)
  call fft_xz_to_physical(Cdu2dy,du2dy)
  call fft_xz_to_physical(Cdu3dz,du3dz)

  call fft_xz_to_physical(Cdu2dx,du2dx)
  call fft_xz_to_physical(Cdu3dx,du3dx)

  call fft_xz_to_physical(Cdu1dz,du1dz)
  call fft_xz_to_physical(Cdu2dz,du2dz)

  call fft_xz_to_physical(Cdu1dy,du1dy)
  call fft_xz_to_physical(Cdu3dy,du3dy)


  do n = 1, N_th
    call fft_xz_to_physical(Cdthetadx(:, :, :, n), dthetadx(:, :, :, n))
  end do
  do n = 1, N_th
    call fft_xz_to_physical(Cdthetady(:, :, :, n), dthetady(:, :, :, n))
  end do
  do n = 1, N_th
    call fft_xz_to_physical(Cdthetadz(:, :, :, n), dthetadz(:, :, :, n))
  end do

  ! We now have all the gradients in physical space

  return
end




! !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! subroutine les_filter_chan(a, jstart, jend)
!   !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!   ! This subroutine applies the LES filter to the input field
!   ! The indices to the start and end of the array in the y-direction
!   ! are also inputted to make the routine cablable of filtering fields
!   ! at either GYF or GY points.
!   ! The array that is passed should be in physical space
!
!
!   integer i, j, k, n
!
!   real(rkind) a(0:Nx + 1, 0:Nz + 1, 0:Nyp + 1)
!   real(rkind) b(0:Nx - 1, 0:Nz - 1, 0:Nyp + 1)
!
!   integer im2(0:Nx - 1), im1(0:Nx - 1), ip1(0:Nx + 1), ip2(0:Nx + 2)
!   integer km2(0:Nz - 1), km1(0:Nz - 1), kp1(0:Nz + 1), kp2(0:Nz + 2)
!
!   ! These are the weights for the filtering operation used
!   real(rkind) w0, w1, w2, Wm1, Wm2, Wm1_j, W0_j, W1_j
!
!   ! The following is for the 3-point trapezoidal rule, alpha*beta=sqrt(6)
!   !      Wm2=0.d0
!   !      Wm1=1.d0/4.d0
!   !      W0=1.d0/2.d0
!   !      W1=1.d0/4.d0
!   !      W2=0.d0
!   Wm1_j = 1.d0 / 4.d0
!   W0_j = 1.d0 / 2.d0
!   W1_j = 1.d0 / 4.d0
!   ! The following is for the 5-point trapezoidal rule, alpha*beta=9
!   Wm2 = 1.d0 / 8.d0
!   Wm1 = 1.d0 / 4.d0
!   w0 = 1.d0 / 4.d0
!   w1 = 1.d0 / 4.d0
!   w2 = 1.d0 / 8.d0
!
!   Nxm1 = Nx - 1
!   Nzm1 = Nz - 1
!
!   !      do j=0,Nyp+1
!   !        do k=0,Nzm1
!   !          do i=0,Nxm1
!   !            B(i, k, j) = A(i, k, j)
!   !          end do
!   !        end do
!   !      end do
!
!   ! Filter in the periodic directions using cshift
!   ! Note, cshift if not used since it appears to be slower
!   ! Apply filter in the x-direction
!   !      B=Wm2*CSHIFT(B,-2,1)+Wm1*CSHIFT(B,-1,1)+W0*B+W1*CSHIFT(B,1,1)
!   !             + W2*CSHIFT(B,2,1)
!
!   ! Filter using more efficient F77 syntax:
!   ! Set up array to loop around periodic directions
!   do i = 2, Nxm1
!     im2(i) = i - 2
!   end do
!   im2(1) = Nxm1
!   im2(0) = Nx - 2
!   do i = 1, Nxm1
!     im1(i) = i - 1
!   end do
!   im1(0) = Nxm1
!   do i = 0, Nx - 2
!     ip1(i) = i + 1
!   end do
!   ip1(Nxm1) = 0
!   do i = 0, Nx - 3
!     ip2(i) = i + 2
!   end do
!   ip2(Nx - 2) = 0
!   ip2(Nxm1) = 1
!
!   do j = jstart, jend
!     do k = 0, Nzm1
!       do i = 0, Nxm1
!         b(i, k, j) = Wm2 * a(im2(i), k, j) + Wm1 * a(im1(i), k, j) + w0 * a(i, k, j) &
!                      + w1 * a(ip1(i), k, j) + w2 * a(ip2(i), k, j)
!       end do
!     end do
!   end do
!
!   ! Apply filter in the z-diretion
!   !      B=Wm2*CSHIFT(B,-2,2)+Wm1*CSHIFT(B,-1,2)+W0*B+W1*CSHIFT(B,1,2)
!   !             + W2*CSHIFT(B,2,2)
!   ! Filter using more efficient F77 syntax:
!   ! Set up array to loop around periodic directions
!   do k = 2, Nzm1
!     km2(k) = k - 2
!   end do
!   km2(1) = Nzm1
!   km2(0) = Nz - 2
!   do k = 1, Nzm1
!     km1(k) = k - 1
!   end do
!   km1(0) = Nzm1
!   do k = 0, Nz - 2
!     kp1(k) = k + 1
!   end do
!   kp1(Nzm1) = 0
!   do k = 0, Nz - 3
!     kp2(k) = k + 2
!   end do
!   kp2(Nz - 2) = 0
!   kp2(Nzm1) = 1
!
!   do j = jstart, jend
!     do k = 0, Nzm1
!       do i = 0, Nxm1
!         a(i, k, j) = Wm2 * b(i, km2(k), j) + Wm1 * b(i, km1(k), j) + w0 * b(i, k, j) &
!                      + w1 * b(i, kp1(k), j) + w2 * b(i, kp2(k), j)
!       end do
!     end do
!   end do
!
!   ! Apply filter in the vertical direction at all physical cells
!   ! (filter is not applied to ghost cells, but the values of the ghost cells
!   ! is used for averaging)
!   !      B(:,:,jstart+1:jend-1) = W0*B(:,:,jstart:jend-2)
!   !                            + W1*B(:,:,jstart+1:jend-1)
!   !                            + W2*B(:,:,jstart+2:jend)
!   ! Use more efficient F77 syntax:
!   !       do j=jstart+1,jend-1
!   !         do k=0,Nzm1
!   !           do i=0,Nxm1
!   !             B(i, k, j) = Wm1_j*B(i, k, j - 1)+W0_j*B(i, k, j)+W1_j*B(i, k, j + 1)
!   !           end do
!   !         end do
!   !       end do
!
!   !      do j=jstart,jend
!   !        do k=0,Nzm1
!   !          do i=0,Nxm1
!   !            A(i, k, j) = B(i, k, j)
!   !          end do
!   !        end do
!   !      end do
!
!   return
! end

! !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! subroutine les_filter_chan_fourier(a, jstart, jend)
!   !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!   ! This subroutine applies the LES filter to the input field
!   ! The filter is a spectral cutoff filter
!   ! The indices to the start and end of the array in the y-direction
!   ! are also inputted to make the routine cablable of filtering fields
!   ! at either GYF or GY points.
!   ! The array that is passed should be in physical space
!
!
!   integer i, j, k
!
!   real(rkind) pi, Lx, Lz
!   integer Nkx, Nkz, twoNkz
!
!   real(rkind) kx(0:Nx / 3), kz(0:2 * (Nz / 3))
!
!   real(rkind) a(0:Nx + 1, 0:Nz + 1, 0:Nyp + 1)
!
!   real(rkind), parameter :: alpha = 2.d0 ! Ratio of filter scales
!
!   real(rkind) b(0:Nx + 1, 0:Nz + 1, 0:Nyp + 1)
!
!   complex(rkind) cb(0:Nx / 2, 0:Nz + 1, 0:Nyp + 1)
!
!   equivalence(b, cb)
!
!   pi = 4.*atan(1.0)
!
!   Lx = pi
!   Lz = 2.d0 * pi
!
!   ! Get the wavenumber vectors:
!   Nkx = Nx / 3
!   do i = 0, Nkx
!     kx(i) = i * (2.*pi) / Lx
!   end do
!
!   Nkz = Nz / 3
!   twoNkz = Nkz * 2
!   do k = 0, Nkz
!     kz(k) = k * (2.*pi) / Lz
!   end do
!   do k = 1, Nkz
!     kz(twoNkz + 1 - k) = -k * (2.*pi) / Lz
!   end do
!
!   do j = 0, Nyp + 1
!     do k = 0, Nzm1
!       do i = 0, Nxm1
!         b(i, k, j) = a(i, k, j)
!       end do
!     enddo
!   end do
!
!   ! Convert to fourier space
!   call fft_xz_to_fourier(b, cb)
!
!   ! Perform the filtering
!   do j = jstart, jend
!     do k = 0, twoNkz
!       do i = 0, Nkx
!         if (sqrt(kx(i)**2.d0 + kz(k)**2.d0) &
!             > sqrt(kx(Nkx)**2.d0 + kz(Nkz)**2.d0) / alpha) then
!           cb(i, k, j) = 0.d0
!         end if
!       end do
!     end do
!   end do
!
!   ! Now, convert back to physical space
!   call fft_xz_to_physical(cb, b)
!
!   do j = jstart, jend
!     do k = 0, Nzm1
!       do i = 0, Nxm1
!         a(i, k, j) = b(i, k, j)
!       end do
!     end do
!   end do
!
!   return
! end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine apply_BC_les
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

  integer i, j, k

  ! If we are using Neuman boundary conditions, over-write the values of the
  ! velocity at the ghost cells so that the LES model doesn't use the large
  ! velocity gradient
  if (u_BC_Ymax == 1) then
    if (rankY == NprocY - 1) then
      ! We are on process at the upper wall
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, Nyp) = cu1(i, k, Nyp - 1)
          cu1(i, k, Nyp + 1) = cu1(i, k, Nyp)
        end do
      end do
    end if
  end if

  if (w_BC_Ymax == 1) then
    if (rankY == NprocY - 1) then
      ! We are on process at the upper wall
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu3(i, k, Nyp) = cu3(i, k, Nyp - 1)
          cu3(i, k, Nyp + 1) = cu3(i, k, Nyp)
        end do
      end do
    end if
  end if

  if (u_BC_Ymin == 1) then
    if (rankY == 0) then
      ! We are on a process at the bottom wall
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu1(i, k, 1) = cu1(i, k, 2)
          cu1(i, k, 0) = cu1(i, k, 1)
        end do
      end do
    end if
  end if

  if (w_BC_Ymin == 1) then
    if (rankY == 0) then
      ! We are on a process at the bottom wall
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cu3(i, k, 1) = cu3(i, k, 2)
          cu3(i, k, 0) = cu3(i, k, 1)
        end do
      end do
    end if
  end if

  return
end





!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine apply_BC_th_les
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

  integer i, j, k, n

  ! If we are using Neuman boundary conditions, over-write the values of the
  ! scalar at the ghost cells so that the LES model doesn't use the large
  ! scalar gradient
  do n = 1, N_th
    if (th_BC_Ymax(n) == 1) then
      if (rankY == NprocY - 1) then
        ! We are on process at the upper wall
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cth(i, k, Nyp, n) = cth(i, k, Nyp - 1, n)
            cth(i, k, Nyp + 1, n) = cth(i, k, Nyp, n)
          end do
        end do
      end if
    end if
    if (th_BC_Ymin(n) == 1) then
      if (rankY == 0) then
        ! We are on process at the upper wall
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            cth(i, k, 1, n) = cth(i, k, 2, n)
            cth(i, k, 0, n) = cth(i, k, 1, n)
          end do
        end do
      end if
    end if

  end do

  return
end





!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine ghost_les_mpi
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine is part of the MPI package for the LES subroutine
  ! Here, after calculating the SGS viscosity, NU_T on each core,
  ! We need to share the ghost cells between neighboring processors



  integer i, j, k, n

  ! Define the arrays that will be used for data packing.  This makes the
  ! communication between processes more efficient by only requiring one
  ! send and recieve.
  ! The communication will be done in Fourier space, so these arrays should
  ! be complex arrays to match the velocity
  ! The size of the buffer array is 0:Nxm1,0:NZP-1
  real(rkind) ocpack(0:Nxm1, 0:Nzp - 1)
  real(rkind) icpack(0:Nxm1, 0:Nzp - 1)

  ! If we are using more than one processor, then we need to pass data

  if (NprocY > 1) then

    ! First, Pass data up the chain to higher ranked processes

    if (rankY == 0) then
      ! If we are the lowest ranked process, then we don't need to recieve
      ! data at the lower ghost cells. Instead, set NU_T=0 at the lower wall
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          nu_t(i, k, 1) = 0.d0
          nu_t(i, k, 2) = 0.d0
        end do
      end do

      ! Pass data up to the next process from GY(Nyp)
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ocpack(i, k) = nu_t(i, k, Nyp)
        end do
      end do
      ! Now, we have packed the data into a compact array, pass the data up
      call mpi_send(ocpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_precision &
                    , rankY + 1, 1, mpi_comm_y, ierror)

      ! End if RANK=0
    else if (rankY < NprocY - 1) then
      ! Here, we are one of the middle processes and we need to pass data
      ! up and recieve data from below
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ocpack(i, k) = nu_t(i, k, Nyp)
        end do
      end do
      ! Use MPI_SENDRECV since we need to recieve and send data
      call mpi_send(ocpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_precision &
                    , rankY + 1, 1, mpi_comm_y, ierror)

      call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_precision &
                    , rankY - 1, 1, mpi_comm_y, status, ierror)
      ! Now, unpack the data that we have recieved
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          nu_t(i, k, 1) = icpack(i, k)
        end do
      end do

    else
      ! Otherwise, we must be the uppermost process with RANK=Nprocs-1
      ! Here, we need to recieve data from below, but don't need to send data up
      call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_precision &
                    , rankY - 1, 1, mpi_comm_y, status, ierror)
      ! Unpack the data that we have recieved
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          nu_t(i, k, 1) = icpack(i, k)
        end do
      end do
    end if

    ! Now, we have hit the top process.  Set the BCs and pass data down

    if (rankY == NprocY - 1) then
      ! If we are the higest ranked process, then we don't need to recieve
      ! data at the upper ghost cells.
      ! Set NU_T=0 at the upper wall
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          nu_t(i, k, Nyp) = 0.d0
          nu_t(i, k, Nyp + 1) = 0.d0
        end do
      end do

      ! Now, send data down the chain
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ocpack(i, k) = nu_t(i, k, 2)
        end do
      end do
      ! Now, we have packed the data into a compact array, pass the data up
      call mpi_send(ocpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_precision &
                    , rankY - 1, 3, mpi_comm_y, ierror)
    else if (rankY > 0) then
      ! Here, we are one of the middle processes and we need to pass data
      ! down and recieve data from above us
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ocpack(i, k) = nu_t(i, k, 2)
        end do
      end do

      call mpi_send(ocpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_precision &
                    , rankY - 1, 3, mpi_comm_y, ierror)

      call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_precision &
                    , rankY + 1, 3, mpi_comm_y, status, ierror)
      ! Now, unpack the data that we have recieved
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          nu_t(i, k, Nyp + 1) = icpack(i, k)
        end do
      end do
    else
      ! Here, we must be the lowest process (RANK=0) and we need to recieve
      ! data from above
      call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_precision &
                    , rankY + 1, 3, mpi_comm_y, status, ierror)
      ! Unpack the data that we have recieved
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          nu_t(i, k, Nyp + 1) = icpack(i, k)
        end do
      end do
    end if

  else
    ! Here, NprocY=1, so we just need to set the boundary values
    ! Set NU_T=0 at the lower wall
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        nu_t(i, k, 1) = 0.d0
        nu_t(i, k, 2) = 0.d0
      end do
    end do
    ! Set NU_T=0 at the upper wall
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        nu_t(i, k, Nyp) = 0.d0
        nu_t(i, k, Nyp + 1) = 0.d0
      end do
    end do

  end if

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine ghost_les_th_mpi
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine is part of the MPI package for the LES subroutine
  ! Here, after calculating the SGS diffusivity, kappa_T on each core,
  ! We need to share the ghost cells between neighboring processors

  integer i, j, k, n

  ! Define the arrays that will be used for data packing.  This makes the
  ! communication between processes more efficient by only requiring one
  ! send and recieve.
  ! The communication will be done in Fourier space, so these arrays should
  ! be complex arrays to match the velocity
  ! The size of the buffer array is 0:Nxm1,0:NZP-1
  real(rkind) ocpack(0:Nxm1, 0:Nzp - 1)
  real(rkind) icpack(0:Nxm1, 0:Nzp - 1)

  ! If we are using more than one processor, then we need to pass data

  do n = 1, n_th
    ! First, Pass data up the chain to higher ranked processes

    if (NprocY > 1) then

      if (rankY == 0) then
        ! If we are the lowest ranked process, then we don't need to recieve
        ! data at the lower ghost cells. Instead, set NU_T=0 at the lower wall
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            kappa_t(i, k, 1, n) = 0.d0
            kappa_t(i, k, 2, n) = 0.d0
          end do
        end do

        ! Pass data up to the next process from GY(Nyp)
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            ocpack(i, k) = kappa_t(i, k, Nyp, n)
          end do
        end do
        ! Now, we have packed the data into a compact array, pass the data up
        call mpi_send(ocpack, (Nxm1 + 1) * (Nzp) &
                      , mpi_double_precision &
                      , rankY + 1, 1, mpi_comm_y, ierror)

        ! End if RANK=0
      else if (rankY < NprocY - 1) then
        ! Here, we are one of the middle processes and we need to pass data
        ! up and recieve data from below
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            ocpack(i, k) = kappa_t(i, k, Nyp, n)
          end do
        end do
        ! Use MPI_SENDRECV since we need to recieve and send data
        call mpi_send(ocpack, (Nxm1 + 1) * (Nzp) &
                      , mpi_double_precision &
                      , rankY + 1, 1, mpi_comm_y, ierror)

        call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                      , mpi_double_precision &
                      , rankY - 1, 1, mpi_comm_y, status, ierror)
        ! Now, unpack the data that we have recieved
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            kappa_t(i, k, 1, n) = icpack(i, k)
          end do
        end do

      else
        ! Otherwise, we must be the uppermost process with RANK=Nprocs-1
        ! Here, we need to recieve data from below, but don't need to send data up
        call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                      , mpi_double_precision &
                      , rankY - 1, 1, mpi_comm_y, status, ierror)
        ! Unpack the data that we have recieved
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            kappa_t(i, k, 1, n) = icpack(i, k)
          end do
        end do
      end if

      ! Now, we have hit the top process.  Set the BCs and pass data down

      if (rankY == NprocY - 1) then
        ! If we are the higest ranked process, then we don't need to recieve
        ! data at the upper ghost cells.
        ! Set NU_T=0 at the upper wall
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            kappa_t(i, k, Nyp, n) = 0.d0
            kappa_t(i, k, Nyp + 1, n) = 0.d0
          end do
        end do

        ! Now, send data down the chain
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            ocpack(i, k) = kappa_t(i, k, 2, n)
          end do
        end do
        ! Now, we have packed the data into a compact array, pass the data up
        call mpi_send(ocpack, (Nxm1 + 1) * (Nzp) &
                      , mpi_double_precision &
                      , rankY - 1, 3, mpi_comm_y, ierror)
      else if (rankY > 0) then
        ! Here, we are one of the middle processes and we need to pass data
        ! down and recieve data from above us
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            ocpack(i, k) = kappa_t(i, k, 2, n)
          end do
        end do

        call mpi_send(ocpack, (Nxm1 + 1) * (Nzp) &
                      , mpi_double_precision &
                      , rankY - 1, 3, mpi_comm_y, ierror)

        call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                      , mpi_double_precision &
                      , rankY + 1, 3, mpi_comm_y, status, ierror)
        ! Now, unpack the data that we have recieved
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            kappa_t(i, k, Nyp + 1, n) = icpack(i, k)
          end do
        end do
      else
        ! Here, we must be the lowest process (RANK=0) and we need to recieve
        ! data from above
        call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                      , mpi_double_precision &
                      , rankY + 1, 3, mpi_comm_y, status, ierror)
        ! Unpack the data that we have recieved
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            kappa_t(i, k, Nyp + 1, n) = icpack(i, k)
          end do
        end do
      end if

    else
      ! Here, NprocY=1, so we just need to set the boundary values
      ! Set NU_T=0 at the lower wall
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          kappa_t(i, k, 1, n) = 0.d0
          kappa_t(i, k, 2, n) = 0.d0
        end do
      end do
      ! Set NU_T=0 at the upper wall
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          kappa_t(i, k, Nyp, n) = 0.d0
          kappa_t(i, k, Nyp + 1, n) = 0.d0
        end do
      end do

    end if

  end do ! n = 1, N_th

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine ghost_CS2_mpi(CS2)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine is part of the MPI package for the LES subroutine
  ! to compute the strain
  ! Share down the CSij2 Nyp + 1 ghost cell between neighboring processors


  integer i, j, k, n
  complex(rkind), intent(inout) :: CS2(:,:,:) ! Can be CSij2 or CSIij2

  ! Define the arrays that will be used for data packing.  This makes the
  ! communication between processes more efficient by only requiring one
  ! send and recieve.
  ! The communication will be done in Fourier space, so these arrays should
  ! be complex arrays to match the velocity
  ! The size of the buffer array is 0:Nxm1,0:NZP-1
  complex(rkind) ocpack(0:Nxp - 1, 0:twoNkz)
  complex(rkind) icpack(0:Nxp - 1, 0:twoNkz)

  ! If we are using more than one processor, then we need to pass data

  if (NprocY > 1) then

    ! Set the BCs and pass data down

    if (rankY == NprocY - 1) then
      ! If we are the higest ranked process, then we don't need to recieve
      ! data at the upper ghost cells.

      ! Set CS2 (cu2(i, k, Nyp + 2) - cu2(i, k, Nyp + 1)) / dyf(Nyp + 1)
      ! We don't need CS2 at the boundaries (dictates CSij4' & 6' only)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          CS2(i, k, Nyp + 1) = -10000.d0
        end do
      end do

      ! Now, send data down the chain
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ocpack(i, k) = CS2(i, k, 2)
        end do
      end do
      ! Now, we have packed the data into a compact array, pass the data up
      call mpi_send(ocpack, Nxp * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY - 1, 3, mpi_comm_y, ierror)
    else if (rankY > 0) then
      ! Here, we are one of the middle processes and we need to pass data
      ! down and recieve data from above us
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          ocpack(i, k) = CS2(i, k, 2)
        end do
      end do

      call mpi_send(ocpack, Nxp * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY - 1, 3, mpi_comm_y, ierror)

      call mpi_recv(icpack, Nxp * (twoNkz + 1) &
                    , mpi_double_complex &
                    , rankY + 1, 3, mpi_comm_y, status, ierror)
      ! Now, unpack the data that we have recieved
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          CS2(i, k, Nyp + 1) = icpack(i, k)
        end do
      end do
    else
      ! Here, we must be the lowest process (RANK=0) and we need to receive
      ! data from above
      call mpi_recv(icpack, (Nxm1 + 1) * (Nzp) &
                    , mpi_double_complex &
                    , rankY + 1, 3, mpi_comm_y, status, ierror)
      ! Unpack the data that we have recieved
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          CS2(i, k, Nyp + 1) = icpack(i, k)
        end do
      end do
    end if

  else
    ! Here, NprocY=1, so we just need to set the boundary values
    ! Set CS2 (cu2(i, k, Nyp + 2) - cu2(i, k, Nyp + 1)) / dyf(Nyp + 1)
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        CS2(i, k, Nyp + 1) = -10000.d0
      end do
    end do

  end if

  return
end
