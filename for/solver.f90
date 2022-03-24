
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 subroutine user_rhs_chan_physical
   !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

   ! Here, you can add terms to the right hand side
   ! of the momentum and scalar equations.
   ! The right hand side forcing arrays, CF1, CF2, CF3, CFTH
   ! are in Fourier space.  The velocity and scalars are available
   ! in physical space.
   ! S1 is available as a working variable

    ! integer i, j, k, n
    !
    ! real(rkind) alpha
    !
    ! call fft_xz_to_physical(cf3, f3)
    !
    ! do j = 1, Nyp
    !   do k = 0, Nzp - 1
    !     do i = 0, Nxm1
    !       f3(i, k, j) = f3(i, k, j) + 0.001 * exp(-((gx(i) - Lx/2)**2 + (gyf(j) - Ly/2)**2)/(2.d0*0.01**2))*sin(0.7*time);
    !     end do
    !   end do
    ! end do
    !
    ! call fft_xz_to_fourier(f3, cf3)

    ! call slip_vel

    return
 end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 subroutine user_rhs_chan_fourier
   !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

   ! Here, you can add terms to the right hand side
   ! of the momentum and scalar equations.
   ! The right hand side forcing arrays, CF1, CF2, CF3, CFTH
   ! are in Fourier space.  The velocity and scalars are available
   ! in Fourier space.
   ! S1 is available as a working variable

   integer i, j, k, n

   ! Advection owing to thermal wind
   if ((flavor == 'Front') .and. (Ro_inv /= 0.d0)) then
     do n = 1, N_th
       ! Loop over all scalars

       ! Add thermal wind advection to the momentum equations
       do j = jstart, jend
         do k = 0, twoNkz
           do i = 0, Nxp - 1
             cf1(i, k, j) = cf1(i, k, j) &
                            - (dTHdX(n) * (gyf(j) - 0.5d0*Ly) * delta / Ro_inv) &
                            * cikz(k) * cu1(i, k, j) &
                            - (-1.d0 * dTHdZ(n) * (gyf(j) - 0.5d0*Ly) * delta / Ro_inv) &
                            * cikx(i) * cu1(i, k, j) &
                            - (-1.d0 * dTHdZ(n) * delta / Ro_inv) &
                            * 0.5d0 * (cu2(i, k, j) + cu2(i, k, j + 1))
             cf3(i, k, j) = cf3(i, k, j) &
                            - (dTHdX(n) * (gyf(j) - 0.5d0*Ly) * delta / Ro_inv) &
                            * cikz(k) * cu3(i, k, j) &
                            - (-1.d0 * dTHdZ(n) * (gyf(j) - 0.5d0*Ly) * delta / Ro_inv) &
                            * cikx(i) * cu3(i, k, j) &
                            - (dTHdX(n) * delta / Ro_inv) &
                            * 0.5d0 * (cu2(i, k, j) + cu2(i, k, j + 1))
           end do
         end do
       end do

       do j = 2, Nyp
         do k = 0, twoNkz
           do i = 0, Nxp - 1
             cf2(i, k, j) = cf2(i, k, j) &
                            - (dTHdX(n) * (gy(j) - 0.5d0*Ly) * delta / Ro_inv) &
                            * cikz(k) * cu2(i, k, j) &
                            - (-1.d0 * dTHdZ(n) * (gy(j) - 0.5d0*Ly) * delta / Ro_inv) &
                            * cikx(i) * cu2(i, k, j)
           end do
         end do
       end do

       ! Add advection by thermal wind to the scalar equations
       do j = jstart_th(n), jend_th(n)
         do k = 0, twoNkz
           do i = 0, Nxp - 1
             cfth(i, k, j, n) = cfth(i, k, j, n) &
                                - (delta / Ro_inv) * dTHdX(n) * (gyf(j) - 0.5d0*Ly) &
                                * cikz(k) * cth(i, k, j, n) &
                                - (delta / Ro_inv) * (-1.d0) * dTHdZ(n) * (gyf(j) - 0.5d0*Ly) &
                                * cikx(i) * cth(i, k, j, n)
           end do
         end do
       end do

       ! End do N_th
     end do

   end if

   ! Add sponge layer forcing
   ! do n = 1, N_th
   !   call sponge_th(n)
   ! end do
   ! call sponge_vel

   return
 end



 !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 subroutine sponge_th(n)
   !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
   ! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
   ! specified background state for the temperature field
   ! The intention is to allow an open boundary

   integer i, j, k, n
   real(rkind) L_sponge, L_bottom
   real(rkind) sponge_amp

   ! The following variables will store the background state
   real(rkind) th_0(-1:Nyp + 1)

   real(rkind) ri_b(0:Nyp + 1)

   ! This variable will hold the forcing rate
   real(rkind) sponge_sigma(0:Nyp + 1)

   ! Set the amplitude of the sponge
   sponge_amp = 0.005d0
   ! Set the top of the sponge layer in physical units
   L_sponge = -120.d0
   ! Set the bottom of the computational domain in physical units
   L_bottom = -140.d0
   do j = 0, Nyp + 1
     ! Quadratic damping at lower wall
     if (gyf(j) < L_sponge) then
       sponge_sigma(j) = sponge_amp * ((L_sponge - gyf(j)) &
                                       / (L_sponge - L_bottom))**2.d0
     else
       sponge_sigma(j) = 0.d0
     end if
   end do

   ! Set the profile for relaxing the mean TH
   do j = 0, Nyp + 1
     th_0(j) = th_BC_Ymin_c1(n) * gyf(j)
   end do

   ! For MLI latmix
   if (n == 1) then
     th_0(0) = 0.d0
     do j = 1, Nyp + 1
       ri_b(j) = 20.d0
       th_0(j) = th_0(j - 1) &
                 + dy(j) * ri_b(j) * (dTHdX(n))**2.d0 &
                 * (delta / Ro_inv)**2.d0
     end do
   else
     do j = 0, Nyp + 1
       th_0(j) = 0.d0
     end do
   end if

   ! Add damping to R-K terms
   ! Damp the perturbations towards 0
   do k = 0, twoNkz
     do i = 0, Nxp - 1
       if ((rankZ /= 0) .or. (i /= 0) .or. (k /= 0)) then
         do j = jstart_th(n), jend_th(n)
           cfth(i, k, j, n) = cfth(i, k, j, n) &
                              - sponge_sigma(j) * (cth(i, k, j, n) - 0.)
         end do
       end if
     end do
   end do
   ! Damp the mean gradient towards TH_0
   do j = jstart_th(n), jend_th(n)
     cfth(0, 0, j, n) = cfth(0, 0, j, n) - sponge_sigma(j) &
                        * (cth(0, 0, j, n) - th_0(j))
   end do

   return
 end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine slip_vel
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine adds advection by a slip velocity to some scalars

  integer i, j, k, n, j1_th(1:N_th), j2_th(1:N_th)
  real(rkind) w_s(0:Nyp + 1, 1:N_th)
  real(rkind), dimension(5) :: slip_def = (/ 0.d0, 0.d0, 0.00005d0, 0.0005d0, 0.005d0/)

  ! Set indices corresponding to start and end of GYF grid
  do n = 1, N_th
   if (rankY == NprocY - 1) then
     ! We are at the upper wall
     j1_th(n) = jstart_th(n)
     j2_th(n) = Nyp - 1
   else if (rankY == 0) then
     ! We are at the lower wall
     j1_th(n) = 2
     j2_th(n) = jend_th(n)
   else
     ! We are on a middle process
     j1_th(n) = jstart_th(n)
     j2_th(n) = jend_th(n)
   end if
  end do

  ! First, set the slip velocity
  do j = 0, Nyp + 1
   do n = 1, N_th
     w_s(j, n) = slip_def(n)
   end do
  end do

  if (rankY == NprocY - 1) then
   ! We are on a process at the top boundary
   ! Set the slip velocity to zero at GY(Nyp) (and ghost cells)
   do n = 1, N_th
     w_s(Nyp, n) = 0.d0
     w_s(Nyp + 1, n) = 0.d0
   end do
  else if (rankY == 0) then
   ! We are on a process at the bottom boundary
   ! Set the slip velocity to zero at GY(2) (and ghost cells)
   do n = 1, N_th
     w_s(0, n) = 0.d0
     w_s(1, n) = 0.d0
     w_s(2, n) = 0.d0
   end do
  end if

  do n = 1, N_th
   do j = j1_th(n), j2_th(n)
     do k = 0, Nzp - 1
       do i = 0, Nxm1
         ! Central differencing
         !              S1(I,K,J)=
         !     &     ((TH(I,K,J+1,N)*W_S(J+1,N) + TH(I,K,J,N)*W_S(J+1,N)
         !     &    -TH(I,K,J,N)*W_S(J,N)-TH(I,K,J-1,N)*W_S(J,n))/(2.d0*DYF(J)))
         ! Second order Upwinding
         !              S1(I,K,J)=(W_S(J+1,N)*TH(I,K,J,N)
         !     &               -W_S(J,N)*(TH(I,K,J,N)+TH(I,K,J-1,N))/2.d0)
         !     &                 /(GYF(j)-GY(j))
         ! First order upwinding
         s1(i, k, j) = (w_s(j + 1, n) * th(i, k, j, n) &
                        - w_s(j, n) * th(i, k, j - 1, n)) &
                       / (gyf(j) - gyf(j - 1))

         !              S1(I,K,J)=0.5d0*(W_S(J+1,N)+W_S(J,N))
         !     &              *(TH(I,K,J,N)-TH(I,K,J-1,N))/(GYF(J)-GYF(J-1))
       end do
     end do
   end do
   call fft_xz_to_fourier(s1, cs1)
   do j = j1_th(n), j2_th(n)
     do k = 0, twoNkz
       do i = 0, Nxp - 1
         cfth(i, k, j, n) = cfth(i, k, j, n) - cs1(i, k, j)
       end do
     end do
   end do
  end do

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine sponge_vel
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
  ! specified background state
  ! The intention is to allow an open boundary

  integer i, j, k

  real(rkind) L_sponge, L_bottom
  real(rkind) sponge_amp

  ! The following variables will store the background state
  real(rkind) u1_0(-1:Nyp + 1), u2_0(0:Nyp + 1), u3_0(-1:Nyp + 1)

  ! This variable will hold the forcing rate
  real(rkind) sponge_sigma(0:Nyp + 1)

  ! Set the amplitude of the sponge
  sponge_amp = 0.0001d0
  ! Set the top of the sponge layer in physical units
  L_sponge = -120.d0
  ! Set the bottom of the computational domain in physical units
  L_bottom = -140.d0
  do j = 0, Nyp + 1
   ! Quadratic damping at lower wall
   if (gyf(j) < L_sponge) then
     sponge_sigma(j) = sponge_amp * ((L_sponge - gyf(j)) &
                                     / (L_sponge - L_bottom))**2.d0
   else
     sponge_sigma(j) = 0.d0
   end if
  end do

  ! Set the background state
  ! Here, set the background to be geostrophic, with a linear temperature profile
  do j = 0, Nyp + 1
   u1_0(j) = 0.d0
   u3_0(j) = 0.d0
  end do
  do j = 0, Nyp + 1
   u2_0(j) = 0.d0
  end do

  ! Add damping function to explicit R-K
  do k = 0, twoNkz
   do i = 0, Nxp - 1 ! Nkx
     if ((i /= 0) .or. (k /= 0)) then
       do j = jstart, jend
         cf1(i, k, j) = cf1(i, k, j) - sponge_sigma(j) * (cu1(i, k, j) - 0.d0)
         cf3(i, k, j) = cf3(i, k, j) - sponge_sigma(j) * (cu3(i, k, j) - 0.d0)
       end do
       do j = 1, Nyp
         cf2(i, k, j) = cf2(i, k, j) - &
                        0.5 * (sponge_sigma(j) + sponge_sigma(j + 1)) * (cu2(i, k, j) - 0.d0)
       end do
     end if
   end do
  end do
  ! Damp mean flow
  do j = jstart, jend
   cf1(0, 0, j) = cf1(0, 0, j) - sponge_sigma(j) * (cu1(0, 0, j) - u1_0(j))
   cf3(0, 0, j) = cf3(0, 0, j) - sponge_sigma(j) * (cu3(0, 0, j) - u3_0(j))
  end do
  do j = 1, Nyp
   cf2(0, 0, j) = cf2(0, 0, j) - sponge_sigma(j) * (cu2(0, 0, j) - u2_0(j))
  end do

  return
end











!******************************************************************************|
! channel.f90, the channel-flow solvers for diablo.                  VERSION 1.0
! This solver was written by John R. Taylor.
!******************************************************************************|
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine rk_chan_1
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Main time-stepping algorithm for the channel-flow case.
  ! This algorithm uses Crank-Nicolson for viscous/diffusive terms and
  ! and 3rd order Runge-Kutta for the rest of the terms
  ! INPUTS  (in Fourier space):  CUi, CP, and (if k>1) CFi at (k-1)  (for i=1,2,3)
  ! OUTPUTS (in Fourier space):  CUi, CP, and (if k<3) CFi at (k)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|


  integer i, j, k, n, istart
  real(rkind) temp1, temp2, temp3, temp4, temp5, ubulk
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl,   matd,   matu, vec
  real(rkind), dimension(0:Nxp,0:Nyp+1) ::    matl_c, matd_c, matu_c
  complex(rkind), dimension(0:Nxp,0:Nyp+1) :: vec_c


  ! Define the constants that are used in the time-stepping
  ! For reference, see Numerical Renaissance
  temp1 = 0.5d0 * nu * h_bar(rk_step)
  temp2 = 0.5d0 * h_bar(rk_step)
  temp3 = zeta_bar(rk_step) * h_bar(rk_step)
  temp4 = h_bar(rk_step)
  temp5 = beta_bar(rk_step) * h_bar(rk_step)

  ! First, we will compute the explicit RHS terms and store in Ri
  ! Note, Momentum equation and hence the RHS is evaluated at the
  ! corresponding velocity points.

  ! Store the old velocity in the RHS vector
  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cr1(i, k, j) = cu1(i, k, j)
        cr3(i, k, j) = cu3(i, k, j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cr2(i, k, j) = cu2(i, k, j)
      end do
    end do
  end do

  ! Add the R-K term from the rk-1 step
  if (rk_step > 1) then
    do j = jstart, jend
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cr1(i, k, j) = cr1(i, k, j) + temp3 * cf1(i, k, j)
          cr3(i, k, j) = cr3(i, k, j) + temp3 * cf3(i, k, j)
        end do
      end do
    end do
    do j = 2, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cr2(i, k, j) = cr2(i, k, j) + temp3 * cf2(i, k, j)
        end do
      end do
    end do
  end if

  ! Take the y-derivative of the pressure at GY points in Fourier space
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = (cp(i, k, j) - cp(i, k, j - 1)) / dy(j)
      end do
    end do
  end do

  ! Add the pressure gradient to the RHS as explicit Euler
  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cr1(i, k, j) = cr1(i, k, j) - temp4 * (cikx(i) * cp(i, k, j))
        cr3(i, k, j) = cr3(i, k, j) - temp4 * (cikz(k) * cp(i, k, j))
      end do
    end do
  end do

  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cr2(i, k, j) = cr2(i, k, j) - temp4 * cs1(i, k, j)
      end do
    end do
  end do

  ! Here, add the constant, forcing pressure gradient
  ! There are two built-in ways of doing this
  ! f_type=1 -> Constant pressure gradient in the x-direction
  ! f_type=2 -> Oscillatory pressure gradient in the x-direction
  ! f_type=4 -> Oscillatory surface forcing on the top (to cu3)
  ! else -> No forcing added
  if (f_type == 1) then
    ! Add forcing for a constant pressure gradient
    do j = jstart, jend
      if (rankZ == 0) cr1(0, 0, j) = cr1(0, 0, j) - temp4 * px0
    end do
  else if (f_type == 2) then
    ! Oscillatory pressure gradient
    do j = jstart, jend
      if (rankZ == 0) cr1(0, 0, j) = cr1(0, 0, j) - &
                                       temp4 * (px0 + amp_omega0 * cos(omega0 * time))
    end do
    ! End if forcing type
  end if

  ! Now compute the term R-K term Ai
  ! Compile terms of Ai in CFi which will be saved for next time step
  ! First, store the horizontal viscous terms in CFi
  ! Note that Beta is the order of the Laplacian operator
  ! e.g. Beta=1 for second order viscosity/diffusivity
  ! or   Beta=2, etc. for hyper viscosity/diffusivity
  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf1(i, k, j) = -nu * (kx2(i) + kz2(k))**beta * cu1(i, k, j)
        cf3(i, k, j) = -nu * (kx2(i) + kz2(k))**beta * cu3(i, k, j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf2(i, k, j) = -nu * (kx2(i) + kz2(k))**beta * cu2(i, k, j)
      end do
    end do
  end do

  ! Add the terms owing to the system rotation (Coriolis terms)
  ! Assume that the flow is on an f-plane
  do k = 0, twoNkz
    do i = 0, Nxp - 1
      do j = jstart, jend
        cf1(i, k, j) = cf1(i, k, j) + (Ro_inv / delta) * cu3(i, k, j)
        cf3(i, k, j) = cf3(i, k, j) - (Ro_inv / delta) * cu1(i, k, j)
      end do
    end do
  end do

  ! Do for each scalar
  do n = 1, N_th

    ! If a scalar contributes to the denisty, RI is not equal to zero and
    ! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode in the
    ! y-direction, which corresponds to the plane-average.
    ! The plane averaged density balances the hydrostati! pressure
    do j = 2, Nyp
      do k = 1, twoNkz
        do i = 0, Nxp - 1
          ! Use second order interpolation
          cf2(i, k, j) = cf2(i, k, j) + grav_y * &
                         (cth(i, k, j, n) * dyf(j - 1) + &
                          cth(i, k, j - 1, n) * dyf(j)) / (2.d0 * dy(j))
        end do
      end do
      k = 0
      if (rankZ == 0) then
        istart = 1
      else
        istart = 0
      end if
      do i = istart, Nxp - 1
        cf2(i, k, j) = cf2(i, k, j) + grav_y * &
                       (cth(i, k, j, n) * dyf(j - 1) + &
                        cth(i, k, j - 1, n) * dyf(j)) / (2.d0 * dy(j))
      end do
    end do

    do j = jstart, jend
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cf1(i, k, j) = cf1(i, k, j) + grav_x * cth(i, k, j, n)
        end do
      end do
    end do

    do j = jstart, jend
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cf3(i, k, j) = cf3(i, k, j) + grav_z * cth(i, k, j, n)
        end do
      end do
    end do

    ! Now, compute the RHS vector for the scalar equations
    ! Since TH is defined at horizontal velocity points, the
    ! scalar update equation will be very similar to the horizontal
    ! velocity update.

    ! We will store the RHS scalar terms in CRTH, RTH
    ! The k-1 term for the R-K stepping is saved in FTH, CFTH

    ! First, build the RHS vector, use CRTH
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          crth(i, k, j, n) = cth(i, k, j, n)
        enddo
      end do
    end do
    ! Add term from k-2 step to free up CFTH variable
    if (rk_step > 1) then
      do j = jstart_th(n), jend_th(n)
        do k = 0, twoNkz
          do i = 0, Nxp - 1
            crth(i, k, j, n) = crth(i, k, j, n) + temp3 * cfth(i, k, j, n)
          end do
        end do
      end do
    end if

    ! Now compute the explicit R-K term Ai
    ! Compile terms of Ai in CFi which will be saved for next time step
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cfth(i, k, j, n) = -(nu / Pr(n)) * &
                              (kx2(i) + kz2(k))**beta * cth(i, k, j, n)
        end do
      end do
    end do

    ! Add advection acting on the background scalar gradient if present
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cfth(i, k, j, n) = cfth(i, k, j, n) - cu1(i, k, j) * dTHdX(n) &
                             - cu3(i, k, j) * dTHdZ(n)
        end do
      end do
    end do

    ! End do number of passive scalars (N_th)
  end do

  ! Optionally, add user forcing to the right hand side
  ! Here, we have U1, U2, U3, and TH in Fourier space
  call user_rhs_chan_fourier

  ! If we are considering an LES, then add the subgrid scale stress:
  ! Here, velocity and CFi should be in Fourier space
  ! The subgrid scale stress is added to CFi:   CFi=CFi - d/dx_i tau_ij

  if (use_LES .and. ((.not. create_new_flow) .or. (time_step > 100))) then
    ! If we have created new flow with random perturbations, wait for a
    ! spinup before applying the subgrid model for stability purposes
    ! In the process (les_chan), Ui is converted to physical space
    ! In the process (les_chan_th), th is converted to physical space

    call les_chan ! Calculates nu_t (and applies it explicitly in horizontal for RK)

    call les_chan_th ! Calculates kappa_t (and applies it explicitly in horizontal for RK)

  else

    if (use_LES .and. flag_save_LES) then
      ! Save out blanks for the LES Stats to keep up with flow stats
      flag_save_LES = .false.
      call save_stats_LES_OOL(.true.)
    end if

    ! If the subgrid model hasn't been called, then it is necessary to
    ! convert to physical space.
    call fft_xz_to_physical(cu1, u1)
    call fft_xz_to_physical(cu2, u2)
    call fft_xz_to_physical(cu3, u3)

    ! Transform THETA to physical space for computation of nonlinear terms
    ! Here pass the first location in memory of the array for scalar n
    do n = 1, N_th
      call fft_xz_to_physical(cth(:, :, :, n), th(:, :, :, n))
    end do

  end if



  ! Compute the nonlinear products in physical space, then transform
  ! back to Fourier space to compute the derivative.
  ! Here, we compute the horizontal derivatives of the nonlinear terms
  ! which will be treated with RKW3.
  ! Do terms one at a time to save on memory
  ! U1*U3
  do j = jstart, jend
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = u3(i, k, j) * u1(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_fourier(s1, cs1)

  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf1(i, k, j) = cf1(i, k, j) - cikz(k) * cs1(i, k, j)
        cf3(i, k, j) = cf3(i, k, j) - cikx(i) * cs1(i, k, j)
      end do
    end do
  end do

  ! U1*U1
  do j = jstart, jend
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = u1(i, k, j) * u1(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_fourier(s1, cs1)

  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf1(i, k, j) = cf1(i, k, j) - cikx(i) * cs1(i, k, j)
      end do
    end do
  end do

  ! U3*U3
  do j = jstart, jend
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = u3(i, k, j) * u3(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_fourier(s1, cs1)

  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf3(i, k, j) = cf3(i, k, j) - cikz(k) * cs1(i, k, j)
      end do
    end do
  end do

  ! U1*U2
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = ((dyf(j) * u1(i, k, j) &
                        + dyf(j - 1) * u1(i, k, j - 1)) / (2.*dy(j))) &
                      * u2(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_fourier(s1, cs1)

  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf2(i, k, j) = cf2(i, k, j) - cikx(i) * cs1(i, k, j)
      end do
    end do
  end do

  ! U3*U2
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = ((dyf(j) * u3(i, k, j) &
                        + dyf(j - 1) * u3(i, k, j - 1)) / (2.*dy(j))) &
                      * u2(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_fourier(s1, cs1)

  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf2(i, k, j) = cf2(i, k, j) - cikz(k) * cs1(i, k, j)
      end do
    end do
  end do

  ! Add the vertical derivative term
  do j = jstart, jend
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = &
          (u1(i, k, j + 1) * u2(i, k, j + 1) + u1(i, k, j) * u2(i, k, j + 1) &
           - u1(i, k, j) * u2(i, k, j) - u1(i, k, j - 1) * u2(i, k, j)) / (2.d0 * dyf(j))
      end do
    end do
  end do

  call fft_xz_to_fourier(s1, cs1)

  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf1(i, k, j) = cf1(i, k, j) - cs1(i, k, j)
      end do
    end do
  end do

  ! Add the vertical derivative term explicitly
  do j = jstart, jend
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = &
          (u3(i, k, j + 1) * u2(i, k, j + 1) + u3(i, k, j) * u2(i, k, j + 1) &
           - u3(i, k, j) * u2(i, k, j) - u3(i, k, j - 1) * u2(i, k, j)) / (2.d0 * dyf(j))
      end do
    end do
  end do

  call fft_xz_to_fourier(s1, cs1)

  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf3(i, k, j) = cf3(i, k, j) - cs1(i, k, j)
      end do
    end do
  end do

  ! Add the vertical derivative term explicitly
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = &
          (0.25d0 * (u2(i, k, j) + u2(i, k, j + 1))**2.d0 &
           - 0.25d0 * (u2(i, k, j) + u2(i, k, j - 1))**2.d0) / dy(j)
      end do
    end do
  end do

  call fft_xz_to_fourier(s1, cs1)

  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf2(i, k, j) = cf2(i, k, j) - cs1(i, k, j)
      end do
    end do
  end do

  ! -- At this point, we are done computing the nonlinear terms --

  ! Optionally, add user forcing to the right hand side
  ! Here, we have U1, U2, U3, and TH in physical space
  call user_rhs_chan_physical

  ! Finally, Add CFi to CRi
  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cr1(i, k, j) = cr1(i, k, j) + temp5 * cf1(i, k, j)
        cr3(i, k, j) = cr3(i, k, j) + temp5 * cf3(i, k, j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cr2(i, k, j) = cr2(i, k, j) + temp5 * cf2(i, k, j)
      end do
    end do
  end do

  ! Convert RHS terms to physical space
  call fft_xz_to_physical(cr1, r1)
  call fft_xz_to_physical(cr2, r2)
  call fft_xz_to_physical(cr3, r3)

  ! Compute the vertical viscous term in physical space and add to RHS
  ! This is the explicit part of the Crank-Nicolson term
  !   At this point, ri is the RK ~half-step for ui
  ! NU_V_SCALE is a coefficient defined as the ratio of the vertical
  ! to horizontal viscosity and can be used to add anisotropic viscosity
  do j = jstart, jend
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        r1(i, k, j) = r1(i, k, j) + temp1 * nu_v_scale * &
                      (((u1(i, k, j + 1) - u1(i, k, j)) / dy(j + 1) &
                        - (u1(i, k, j) - u1(i, k, j - 1)) / dy(j)) / dyf(j))
        r3(i, k, j) = r3(i, k, j) + temp1 * nu_v_scale * &
                      (((u3(i, k, j + 1) - u3(i, k, j)) / dy(j + 1) &
                        - (u3(i, k, j) - u3(i, k, j - 1)) / dy(j)) / dyf(j))
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        r2(i, k, j) = r2(i, k, j) + temp1 * nu_v_scale * &
                      (((u2(i, k, j + 1) - u2(i, k, j)) / dyf(j) &
                        - (u2(i, k, j) - u2(i, k, j - 1)) / dyf(j - 1)) / dy(j))
      end do
    end do
  end do

  ! If we are using a subgrid model, add the eddy viscosity term
  ! This is an added viscosity that will be treated just like the
  ! molecular viscosity with Crank-Nicolson for the vertical derivatives
  if (use_LES) then
    ! Note, nu_t is defined at GY points
    do j = jstart, jend
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          r1(i, k, j) = r1(i, k, j) + temp2 * &
                        ((nu_t(i, k, j + 1) * (u1(i, k, j + 1) - u1(i, k, j)) / dy(j + 1) &
                          - nu_t(i, k, j) * (u1(i, k, j) - u1(i, k, j - 1)) / dy(j)) &
                         / dyf(j))
          r3(i, k, j) = r3(i, k, j) + temp2 * &
                        ((nu_t(i, k, j + 1) * (u3(i, k, j + 1) - u3(i, k, j)) / dy(j + 1) &
                          - nu_t(i, k, j) * (u3(i, k, j) - u3(i, k, j - 1)) / dy(j)) &
                         / dyf(j))
        end do
      end do
    end do
    ! Here, interpolate nu_t to GYF points
    do j = 2, Nyp
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          r2(i, k, j) = r2(i, k, j) + temp2 * &
                        ((0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j + 1)) * (u2(i, k, j + 1) - u2(i, k, j)) &
                          / dyf(j) &
                          - 0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j - 1)) * (u2(i, k, j) - u2(i, k, j - 1)) &
                          / dyf(j - 1)) / dy(j))
        end do
      end do
    end do
  end if

  ! -- Here, we are done with computation of Velocity RHS, explicit terms --

  ! Now, build the explicit RHS terms for the passive scalar(s)

  do n = 1, N_th
    ! Do for each scalar:

    ! Compute the nonlinear terms that are present in the explicit term A
    ! U1*TH
    do j = jstart_th(n), jend_th(n)
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          s1(i, k, j) = th(i, k, j, n) * u1(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_fourier(s1, cs1)
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cfth(i, k, j, n) = cfth(i, k, j, n) - cikx(i) * cs1(i, k, j)
        end do
      end do
    end do
    ! U3*TH
    do j = jstart_th(n), jend_th(n)
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          s1(i, k, j) = th(i, k, j, n) * u3(i, k, j)
        end do
      end do
    end do
    call fft_xz_to_fourier(s1, cs1)
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cfth(i, k, j, n) = cfth(i, k, j, n) - cikz(k) * cs1(i, k, j)
        end do
      end do
    end do

    ! We are done with the horizontal derivatives of the nonlinear terms
    ! Add the vertical derivative term explicitly
    do j = jstart_th(n), jend_th(n)
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          s1(i, k, j) = &
            (th(i, k, j + 1, n) * u2(i, k, j + 1) + th(i, k, j, n) * u2(i, k, j + 1) &
             - th(i, k, j, n) * u2(i, k, j) - th(i, k, j - 1, n) * u2(i, k, j)) / (2.d0 * dyf(j))
        end do
      end do
    end do
    call fft_xz_to_fourier(s1, cs1)
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          cfth(i, k, j, n) = cfth(i, k, j, n) - cs1(i, k, j)
        end do
      end do
    end do

    ! Add CFTH to the RHS vector CRTH
    do j = jstart_th(n), jend_th(n)
      do k = 0, twoNkz
        do i = 0, Nxp - 1
          crth(i, k, j, n) = crth(i, k, j, n) + temp5 * cfth(i, k, j, n)
        end do
      end do
    end do
    ! Done with computation of RHS, explicit terms for the THETA equation
    ! Transform back to physical space

    call fft_xz_to_physical(crth(:, :, :, n), rth(:, :, :, n))

    ! Compute the Explicit part of the Crank-Nicolson terms for the TH equation
    ! First, the vertical derivative viscous term
    do j = jstart_th(n), jend_th(n)
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          rth(i, k, j, n) = rth(i, k, j, n) + (temp1 / Pr(n)) * nu_v_scale * ( &
                            ((th(i, k, j + 1, n) - th(i, k, j, n)) / dy(j + 1) &
                             - (th(i, k, j, n) - th(i, k, j - 1, n)) / dy(j)) / dyf(j))
        end do
      end do
    end do
    ! If we are using a subgrid model (LES) then add the eddy diffusivity here
    ! Note, KAPPA_T is defined at GY points
    if (use_LES) then
      do j = jstart_th(n), jend_th(n)
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            rth(i, k, j, n) = rth(i, k, j, n) + temp2 * ( &
                              (kappa_t(i, k, j + 1, n) * (th(i, k, j + 1, n) - th(i, k, j, n)) / dy(j + 1) &
                               - kappa_t(i, k, j, n) * (th(i, k, j, n) - th(i, k, j - 1, n)) / dy(j)) / dyf(j))
          end do
        end do
      end do
    end if

    ! -- Now, timestep the passive scalar equation --
    ! We solve the the passive scalar before the velocity so that
    ! it is advected with the velocity from the previous R-K step
    ! which we have already made divergence free

    ! Solve the implicit equation for THETA
    ! Note that the system size is Nyp+1, but only 1..Nyp are used

    ! Initialize the matrix used to store implicit coefficients
    matl = 0.
    matd = 1.
    matu = 0.
    vec = 0.

    ! Build implicit matrix
    ! Use quasi-second order interpolation for TH on GY points
    do k = 0, Nzp - 1
      do j = jstart_th(n), jend_th(n)
        do i = 0, Nxm1
          matl(i, j) = -(temp1 / Pr(n) * nu_v_scale) / (dy(j) * dyf(j))
          matd(i, j) = 1.+(temp1 / Pr(n) * nu_v_scale) / (dy(j + 1) * dyf(j)) &
                       + (temp1 / Pr(n) * nu_v_scale) / (dy(j) * dyf(j))
          matu(i, j) = -(temp1 / Pr(n) * nu_v_scale) / (dy(j + 1) * dyf(j))
          ! Define RHS vector
          vec(i, j) = rth(i, k, j, n)
        end do
      end do
      ! IF using a subgrid model (LES) then add the eddy diffusivity part implicitly
      if (use_LES) then
        do j = jstart_th(n), jend_th(n)
          do i = 0, Nxm1
            matl(i, j) = matl(i, j) - temp2 * kappa_t(i, k, j, n) &
                         / (dy(j) * dyf(j))
            matd(i, j) = matd(i, j) + temp2 * kappa_t(i, k, j + 1, n) &
                         / (dy(j + 1) * dyf(j)) &
                         + temp2 * kappa_t(i, k, j, n) &
                         / (dy(j) * dyf(j))
            matu(i, j) = matu(i, j) - temp2 * kappa_t(i, k, j + 1, n) &
                         / (dy(j + 1) * dyf(j))
          end do
        end do
      end if

      ! If we are using MPI, then solve the implicit system in separate forward
      ! and backward sweeps for efficiency
      call apply_BC_th_mpi(matl, matd, matu, vec, n)
      ! If we are using MPI, split the implicit solve into foward and
      ! backward sweeps for efficiency
      call thomas_forward_real_mpi(matl, matd, matu, vec, Nyp, Nx)
      call thomas_backward_real_mpi(matl, matd, matu, vec, Nyp, Nx)


      do j = jstart_th(n) - 1, jend_th(n) + 1
        do i = 0, Nxm1
          th(i, k, j, n) = vec(i, j)
        end do
      end do

      ! END do k
    end do

    ! End do number of passive scalars
  end do

  ! Initialize the matrix to zeros to be used for implicit solves
  ! Note that the system size is Nyp+1, but only 1..Nyp are used

  ! Initialize the matrix used to store implicit coefficients
  do j = 0, Nyp + 1
    do i = 0, Nxm1
      matl(i, j) = 0.
      matd(i, j) = 1.
      matu(i, j) = 0.
      vec(i, j) = 0.
    end do
  end do

  ! Build implicit matrix for U2
  do k = 0, Nzp - 1
    do j = 2, Nyp
      do i = 0, Nxm1
        matl(i, j) = -temp1 * nu_v_scale / (dyf(j - 1) * dy(j))
        matd(i, j) = 1.+temp1 * nu_v_scale / (dyf(j) * dy(j)) &
                     + temp1 * nu_v_scale / (dyf(j - 1) * dy(j))
        matu(i, j) = -temp1 * nu_v_scale / (dyf(j) * dy(j))
        vec(i, j) = r2(i, k, j)
      end do
    end do
    if (use_LES) then
      ! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
      do j = 2, Nyp
        do i = 0, Nxm1
          matl(i, j) = matl(i, j) &
                       - temp2 * 0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j - 1)) / (dyf(j - 1) * dy(j))
          matd(i, j) = matd(i, j) &
                       + temp2 * 0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j + 1)) / (dyf(j) * dy(j)) &
                       + temp2 * 0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j - 1)) / (dyf(j - 1) * dy(j))
          matu(i, j) = matu(i, j) &
                       - temp2 * 0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j + 1)) / (dyf(j) * dy(j))
        end do
      end do
    end if


    ! First, apply the boundary conditions
    call apply_BC_u2_mpi(matl, matd, matu, vec)
    ! If we are using MPI, split the implicit solve into forward and
    ! backward sweeps for efficiency
    call thomas_forward_real_mpi(matl, matd, matu, vec, Nyp, Nx)
    call thomas_backward_real_mpi(matl, matd, matu, vec, Nyp, Nx)

    do j = 1, Nyp + 1
      do i = 0, Nxm1
        u2(i, k, j) = vec(i, j)
      end do
    end do
    ! End do k
  end do

  ! Solve for U1
  ! Note, here the matrix will be indexed from 1...Nyp+1 corresponding to U1(0:Nyp)

  ! Initialize the matrix used to store implicit coefficients
  do j = 0, Nyp + 1
    do i = 0, Nxm1
      matl(i, j) = 0.
      matd(i, j) = 1.
      matu(i, j) = 0.
      vec(i, j) = 0.
    end do
  end do

  ! Build the implicit system of equations for U1
  do k = 0, Nzp - 1
    do j = jstart, jend
      do i = 0, Nxm1
        matl(i, j) = -temp1 * nu_v_scale / (dy(j) * dyf(j))
        matd(i, j) = 1.-temp1 * nu_v_scale * (-1./(dy(j + 1) * dyf(j)) &
                                              - 1./(dy(j) * dyf(j)))
        matu(i, j) = -temp1 * nu_v_scale / (dy(j + 1) * dyf(j))
        vec(i, j) = r1(i, k, j)
      end do
    end do
    ! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
    if (use_LES) then
      do j = jstart, jend
        do i = 0, Nxm1
          matl(i, j) = matl(i, j) - temp2 * nu_t(i, k, j) &
                       / (dy(j) * dyf(j))
          matd(i, j) = matd(i, j) + temp2 * nu_t(i, k, j + 1) &
                       / (dy(j + 1) * dyf(j)) &
                       + temp2 * nu_t(i, k, j) &
                       / (dy(j) * dyf(j))
          matu(i, j) = matu(i, j) - temp2 * nu_t(i, k, j + 1) &
                       / (dy(j + 1) * dyf(j))
        end do
      end do
    end if

    ! First, apply the boundary conditions
    call apply_BC_u1_mpi(matl, matd, matu, vec)
    ! If we are using MPI, split the implicit solve into forward and
    ! backward sweeps for efficiency
    call thomas_forward_real_mpi(matl, matd, matu, vec, Nyp, Nx)
    call thomas_backward_real_mpi(matl, matd, matu, vec, Nyp, Nx)


    do j = jstart - 1, jend + 1
      do i = 0, Nxm1
        u1(i, k, j) = vec(i, j)
      end do
    end do

    ! End do k
  end do

  ! Initialize the matrix used to store implicit coefficients
  do j = 0, Nyp + 1
    do i = 0, Nxm1
      matl(i, j) = 0.
      matd(i, j) = 1.
      matu(i, j) = 0.
      vec(i, j) = 0.
    end do
  end do

  ! Solve for U3
  ! Note, here the matrix will be indexed from 1...Nyp+1 corresponding to U1(0:Nyp)
  ! Build the implicit system of equations for U3
  do k = 0, Nzp - 1
    do j = jstart, jend
      do i = 0, Nxm1
        matl(i, j) = -temp1 * nu_v_scale / (dy(j) * dyf(j))
        matd(i, j) = 1.-temp1 * nu_v_scale * (-1./(dy(j + 1) * dyf(j)) &
                                              - 1./(dy(j) * dyf(j)))
        matu(i, j) = -temp1 * nu_v_scale / (dy(j + 1) * dyf(j))
        vec(i, j) = r3(i, k, j)
      end do
    end do
    ! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
    if (use_LES) then
      do j = jstart, jend
        do i = 0, Nxm1
          matl(i, j) = matl(i, j) - temp2 * nu_t(i, k, j) &
                       / (dy(j) * dyf(j))
          matd(i, j) = matd(i, j) + temp2 * nu_t(i, k, j + 1) &
                       / (dy(j + 1) * dyf(j)) &
                       + temp2 * nu_t(i, k, j) &
                       / (dy(j) * dyf(j))
          matu(i, j) = matu(i, j) - temp2 * nu_t(i, k, j + 1) &
                       / (dy(j + 1) * dyf(j))
        end do
      end do
    end if

    ! First, apply the boundary conditions
    call apply_BC_u3_mpi(matl, matd, matu, vec)
    ! If we are using MPI, split the implicit solve into forward and
    ! backward sweeps for efficiency
    call thomas_forward_real_mpi(matl, matd, matu, vec, Nyp, Nx)
    call thomas_backward_real_mpi(matl, matd, matu, vec, Nyp, Nx)


    do j = jstart - 1, jend + 1
      do i = 0, Nxm1
        u3(i, k, j) = vec(i, j)
      end do
    end do
    ! End do k
  end do

  ! -- Done getting U1hat, U2hat, U3hat at new RK Step --

  ! If we are on the final RK step, optionally update the timestep
  ! based on the CFL criteria. This won't affect the current timestep
  ! since the TEMP1, etc. variables have already been set using
  ! the current timestep
  if (variable_dt .and. (rk_step == 3) &
      .and. (mod(time_step, update_dt) == 0)) then
    call courant
  end if

  ! Transform TH and U to Fourier Space
  call fft_xz_to_fourier(u1, cu1)
  call fft_xz_to_fourier(u2, cu2)
  call fft_xz_to_fourier(u3, cu3)
  do n = 1, N_th
    call fft_xz_to_fourier(th(:, :, :, n), cth(:, :, :, n))
  end do

  ! Begin second step of the Fractional Step algorithm, making u divergence free
  ! The following subroutine projects Uhat onto divergence free space

  call rem_div_chan

  ! Now, phi is stored in CR1, use this to update the pressure field
  ! Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
  do j = jstart, jend
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cp(i, k, j) = cp(i, k, j) + cr1(i, k, j) / temp4
      end do
    end do
  end do

  ! Fix disparities at the boundary due to the thomas algorithm in parallel
  call ghost_chan_mpi

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine rk_chan_2
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Alternative time-stepping algorithm for the channel-flow case.
  ! This algorithm uses Crank-Nicolson for all viscous terms and
  ! third order Runge-Kutta for all nonlinear terms
  ! INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
  ! OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
  ! Each RK step, there are 11 FFT calls. 11 storage variables are used.
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

  integer i, j, k, n
  real(rkind) temp1, temp2, temp3, temp4, temp5, ubulk

  ! STOP -----------------------------------
  if (rank == 0) &
    write (*, '("RK_CHAN_2 not supported yet")')
  call mpi_finalize(ierror)
  stop
  ! ----------------------------------------

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine rem_div_chan
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Compute varphi, store in variable CR1.
  ! Solves for phi in computational space
  ! H_BAR has been absorbed into PHI, so we are solving for H_BAR*PHI

  integer i, j, k
  complex(rkind), dimension(0:Nxp,0:Nyp+1) :: vec_c
  real(rkind), dimension(0:Nxp,0:Nyp+1) ::    matl_c,matd_c,matu_c

  ! First, Initialize the matrix components
  do j = 0, Nyp + 1
    do i = 0, Nxp - 1
      matl_c(i, j) = 0.
      matd_c(i, j) = 1.
      matu_c(i, j) = 0.
      vec_c(i, j) = (0., 0.)
    end do
  end do

  ! The 2d FFT of Ui should have been taken and stored in CUi
  ! Solving for phi amounts to solving a tridiagonal system
  ! First, construct the system to be solved
  do k = 0, twoNkz
    do j = 1, Nyp
      do i = 0, Nxp - 1
        matl_c(i, j) = 1./(dy(j) * dyf(j))
        matd_c(i, j) = -kx2(i) - kz2(k) &
                       - 1./(dy(j + 1) * dyf(j)) - 1./(dy(j) * dyf(j))
        matu_c(i, j) = 1./(dy(j + 1) * dyf(j))
      end do
    end do

    ! Now, create the RHS vector
    do j = 1, Nyp
      do i = 0, Nxp - 1
        vec_c(i, j) = (cikx(i) * cu1(i, k, j) &
                       + (cu2(i, k, j + 1) - cu2(i, k, j)) / dyf(j) &
                       + cikz(k) * cu3(i, k, j))
      end do
    end do

    ! If we are using the MPI package...
    call apply_BC_rem_div_mpi(matl_c, matd_c, matu_c, vec_c, k)
    ! First, do all forward sweeps
    call thomas_forward_complex_mpi(matl_c, matd_c, matu_c, vec_c, Nyp, Nxp)
    ! Now, do the backward sweeps
    call thomas_backward_complex_mpi(matl_c, matd_c, matu_c, vec_c, Nyp, Nxp)


    do j = 1, Nyp
      do i = 0, Nxp - 1
        cr1(i, k, j) = vec_c(i, j)
      end do
    end do

  end do

  ! Now, Solve for CUi, the divergenceless velocity field
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu1(i, k, j) = cu1(i, k, j) - cikx(i) * cr1(i, k, j)
        cu3(i, k, j) = cu3(i, k, j) - cikz(k) * cr1(i, k, j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cu2(i, k, j) = cu2(i, k, j) - (cr1(i, k, j) &
                                       - cr1(i, k, j - 1)) / dy(j)
      end do
    end do
  end do

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine poisson_p_chan
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! We have CUi, need to compute CP.  Solve tridiagonal system exactly

  integer i, j, k, n
  complex(rkind), dimension(0:Nxp,0:Nyp+1) :: vec_c
  real(rkind), dimension(0:Nxp,0:Nyp+1) ::    matl_c,matd_c,matu_c

  if (flavor == 'Basic') then
    if (rank == 0) &
      write (*, '("Computing cp from cui divergence")')

  end if

  ! First, construct the RHS vector, (dui/dxj)(duj/dxi)
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cf1(i, k, j) = cikx(i) * cu1(i, k, j)
        cf2(i, k, j) = (cu2(i, k, j + 1) - cu2(i, k, j)) / dyf(j)
        cf3(i, k, j) = cikz(k) * cu3(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_physical(cf1, f1)
  call fft_xz_to_physical(cf2, f2)
  call fft_xz_to_physical(cf3, f3)

  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        f1(i, k, j) = f1(i, k, j)**2.
        f2(i, k, j) = f2(i, k, j)**2.
        f3(i, k, j) = f3(i, k, j)**2.
      end do
    end do
  end do

  call fft_xz_to_fourier(f1, cf1)
  call fft_xz_to_fourier(f2, cf2)
  call fft_xz_to_fourier(f3, cf3)

  ! Now we have the diagonal terms, add to the rhs term
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cs1(i, k, j) = cf1(i, k, j) + cf2(i, k, j) + cf3(i, k, j)
      end do
    end do
  end do

  ! Now get the first of the off-diagonal terms
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cf1(i, k, j) = (cu1(i, k, j + 1) - cu1(i, k, j - 1)) / (2.*dyf(j))
        cf2(i, k, j) = cikx(i) * 0.5 * (cu2(i, k, j) + cu2(i, k, j + 1))
      end do
    end do
  end do

  call fft_xz_to_physical(cf1, f1)
  call fft_xz_to_physical(cf2, f2)

  ! Compute product
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        f1(i, k, j) = 2.*f1(i, k, j) * f2(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_fourier(f1, cf1)

  ! Add to RHS term
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cs1(i, k, j) = cs1(i, k, j) + cf1(i, k, j)
      end do
    end do
  end do

  ! Now get the second of the off-diagonal terms
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cf1(i, k, j) = (cu3(i, k, j + 1) - cu3(i, k, j - 1)) / (2.*dyf(j))
        cf2(i, k, j) = cikz(k) * 0.5 * (cu2(i, k, j) + cu2(i, k, j + 1))
      end do
    end do
  end do

  ! Convert to Physical space
  call fft_xz_to_physical(cf1, f1)
  call fft_xz_to_physical(cf2, f2)

  ! Compute product
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        f1(i, k, j) = 2.*f1(i, k, j) * f2(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_fourier(f1, cf1)

  ! Add to RHS term
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cs1(i, k, j) = cs1(i, k, j) + cf1(i, k, j)
      end do
    end do
  end do

  ! Now get the third of the off-diagonal terms
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cf1(i, k, j) = cikz(k) * cu1(i, k, j)
        cf2(i, k, j) = cikx(i) * cu3(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_physical(cf1, f1)
  call fft_xz_to_physical(cf2, f2)

  ! Compute product
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        f1(i, k, j) = 2.*f1(i, k, j) * f2(i, k, j)
      end do
    end do
  end do

  call fft_xz_to_fourier(f1, cf1)

  ! Add to RHS term
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cs1(i, k, j) = cs1(i, k, j) + cf1(i, k, j)
      end do
    end do
  end do

  ! Finally, if the buoyancy force is active, then we need to add
  ! the contribution of the density to the pressure.  Note that the
  ! plane averaged density and the corresponding hydrostatic part of the
  ! pressure have been cancelled, so skip the 0,0 mode
  do n = 1, N_th
    do j = 2, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1 ! Nkx
          if ((rankZ /= 0) .or. (i /= 0) .or. (k /= 0)) then
            cs1(i, k, j) = cs1(i, k, j) + &
                           (cth(i, k, j + 1, n) - cth(i, k, j - 1, n)) / (gyf(j + 1) - gyf(j - 1))
          end if
        end do
      end do
    end do
  end do

  ! Now, the RHS term should be stored in cs1

  ! Construct the tridiagonal system in Fourier space to solve for CP
  ! First, zero the vectors
  do j = 0, Nyp + 1
    do i = 0, Nxp - 1 ! Nkx
      matl_c(i, j) = 0.d0
      matd_c(i, j) = 1.d0
      matu_c(i, j) = 0.d0
      vec_c(i, j) = (0., 0.)
    end do
  end do

  do k = 0, twoNkz
    do j = 2, Nyp
      do i = 0, Nxp - 1 ! Nkx
        matl_c(i, j) = 1./(dy(j) * dyf(j))
        matd_c(i, j) = -kx2(i) - kz2(k) - 1./(dy(j + 1) * dyf(j)) &
                       - 1./(dy(j) * dyf(j))
        matu_c(i, j) = 1./(dy(j + 1) * dyf(j))
        vec_c(i, j) = -1.*cs1(i, k, j)
      end do
    end do

    call apply_BC_poisson_mpi(matl_c, matd_c, matu_c, vec_c, k)
    ! First, do the forward sweeps
    call thomas_forward_complex_mpi(matl_c, matd_c, matu_c, vec_c, Nyp, Nxp)
    ! Now, do the backwared sweeps to put the solution in VEC_C
    call thomas_backward_complex_mpi(matl_c, matd_c, matu_c, vec_c, Nyp, Nxp)


    do j = 1, Nyp
      do i = 0, Nxp - 1 ! Nkx
        cp(i, k, j) = vec_c(i, j)
      end do
    end do
  end do

  return
end





!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine filter_chan(n)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies a filter to the highest wavenumbers
  ! It should be applied to the scalar in Fourier space
  ! The filter used is a sharpened raised cosine filter in the horizontal
  ! and a fourth order implicit compact filter in the vertical, with the
  ! parameter alpha determining the width of the vertical filtering window

  integer i, j, k, js, je, n

  ! Variables for horizontal filtering
  real(rkind) sigma(0:Nxp - 1, 0:twoNkz), sigma0

  ! Variables for vertical filtering
  real(rkind), parameter :: alpha = 0.d0

  ! Parameters for a larger stencil filter
  real(rkind) f_a, f_b, f_c

  ! Set the filtering constants for the horizontal direction
  do i = 0, Nxp - 1
    do k = 0, twoNkz
      sigma0 = 0.5d0 * (1.d0 + &
                        cos(sqrt((kx(i) * Lx * 1.d0 / float(Nx))**2.d0 &
                                 + (kz(k) * Lz * 1.d0 / float(Nz))**2.d0)))
      ! Apply a sharpened raised cosine filter
      sigma(i, k) = sigma0**4.d0 * (35.d0 - 84.d0 * sigma0 &
                                    + 70.d0 * sigma0**2.d0 - 20.d0 * sigma0**3.d0)
    end do
  end do

  ! Do the spectral filtering in the horizontal
  do k = 0, twoNkz
    do i = 0, Nxp - 1
      do j = jstart_th(n), jend_th(n)
        cth(i, k, j, n) = cth(i, k, j, n) * sigma(i, k)
      end do
    end do
  end do

  ! Filter the passive scalar, TH in the vertical direction
  ! Set the filtering constants
  !      f_a=(1.d0/8.d0)*(5.d0+6.d0*alpha)
  !      f_b=0.5d0*(1.d0+2.d0*alpha)
  !      f_c=(-1.d0/8.d0)*(1.d0-2.d0*alpha)
  ! First, zero the tridiagonal matrix components
  !      DO I=0,Nkx
  !        DO J=1,Nyp
  !          MATD_C(I,J)=1.d0
  !          MATL_C(I,J)=0.d0
  !          MATU_C(I,J)=0.d0
  !          VEC_C(I,J)=0.d0
  !        END DO
  !      END DO
  !      DO K=1,twoNkz
  !        DO I=1,Nkx
  ! Construct the centered difference terms
  !          DO J=2,Nyp-1
  !            MATL_C(I,J)=alpha
  !            MATD_C(I,J)=1.d0
  !            MATU_C(I,J)=alpha
  !            VEC_C(I,J)=f_a*CTH(I,K,J,n)
  !     &                +(f_b/2.d0)*(CTH(I,K,J+1,n)+CTH(I,K,J-1,n))
  !     &                +(f_c/2.d0)*(CTH(I,K,J+2,n)+CTH(I,K,J-2,n))
  !          END DO
  ! Now, construct the equations for the boundary nodes
  !          J=1
  !            MATL_C(I,J)=0.d0
  !            MATD_C(I,J)=1.d0
  !            MATU_C(I,J)=0.d0
  !            VEC_C(I,J)=CTH(I,K,J,n)
  !          J=Nyp
  !            MATL_C(I,J)=0.d0
  !            MATD_C(I,J)=1.d0
  !            MATU_C(I,J)=0.d0
  !            VEC_C(I,J)=CTH(I,K,J,n)
  !         END DO
  ! Now, solve the tridiagonal system
  !         CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,Nyp,Nkx)
  !         DO I=1,Nkx
  !           DO J=JSTART_TH(N),JEND_TH(N)
  !             CTH(I,K,J,n)=VEC_C(I,J)
  !           END DO
  !         END DO
  ! END DO K
  !       END DO

  return
end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine thomas_forward_real_mpi(a, b, c, g, iNyp, inx)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This version solves for one full row, then passes the data
  ! This subroutine performs the backward sweep of the Thomas algorithm
  ! Thomas algorithm solves Ax=b for tridiagonal A
  ! The RHS vector and solution are real
  ! Input lower, main, and upper diagonals, ld, md, ud, and rhs x
  ! Returns solution in x
  ! The indexing should be done by ROW, ie.
  ! [ b1  c1   0   0   0 ...
  ! [ a2  b2  c2   0   0 ...
  ! [  0  a3  b3   c3  0 ...

  integer i, j, inx, iNyp
  real(rkind), dimension(0:inx - 1, 0:iNyp + 1) :: a, b, c
  real(rkind), dimension(0:inx - 1, 0:iNyp + 1) :: g

  real(rkind) ocpack(0:inx - 1, 4), icpack(0:inx - 1, 4)

  if (rankY /= 0) then
    ! If we aren't the lowest process, then wait for data
    call mpi_recv(ocpack, 4 * inx, mpi_double_precision, rankY - 1, 12 &
                  , mpi_comm_y, status, ierror)
    ! Unpack the data
    do i = 0, inx - 1
      a(i, 1) = ocpack(i, 1)
      b(i, 1) = ocpack(i, 2)
      c(i, 1) = ocpack(i, 3)
      g(i, 1) = ocpack(i, 4)
    end do
    ! If we aren't the lowest process, start at J=2
    do j = 2, iNyp
      do i = 0, inx - 1
        a(i, j) = -a(i, j) / b(i, j - 1)
        b(i, j) = b(i, j) + a(i, j) * c(i, j - 1)
        g(i, j) = g(i, j) + a(i, j) * g(i, j - 1)
      end do
    end do
  else
    ! Here, we are the lowest process, start solving at J=1
    do j = 1, iNyp
      do i = 0, inx - 1
        a(i, j) = -a(i, j) / b(i, j - 1)
        b(i, j) = b(i, j) + a(i, j) * c(i, j - 1)
        g(i, j) = g(i, j) + a(i, j) * g(i, j - 1)
      end do
    end do
  end if

  if (rankY /= NprocY - 1) then
    do i = 0, inx - 1
      icpack(i, 1) = a(i, iNyp)
      icpack(i, 2) = b(i, iNyp)
      icpack(i, 3) = c(i, iNyp)
      icpack(i, 4) = g(i, iNyp)
    end do
    call mpi_send(icpack, 4 * inx, mpi_double_precision, rankY + 1, 12 &
                  , mpi_comm_y, ierror)
  else
    ! Here, we are at the upper process, so solve one more row containing
    ! the boundary conditions
    j = iNyp + 1
    do i = 0, inx - 1
      a(i, j) = -a(i, j) / b(i, j - 1)
      b(i, j) = b(i, j) + a(i, j) * c(i, j - 1)
      g(i, j) = g(i, j) + a(i, j) * g(i, j - 1)
    end do
  end if

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|------
subroutine thomas_forward_complex_mpi(a, b, c, g, iNyp, inx)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-----|
  ! This subroutine performs the backward sweep of the Thomas algorithm
  ! Thomas algorithm solves Ax=b for tridiagonal A
  ! The RHS vector and solution are real
  ! Input lower, main, and upper diagonals, ld, md, ud, and rhs x
  ! Returns solution in x
  ! The indexing should be done by ROW, ie.
  ! [ b1  c1   0   0   0 ...
  ! [ a2  b2  c2   0   0 ...
  ! [  0  a3  b3   c3  0 ...

  integer i, j, inx, iNyp
  real(rkind), dimension(0:inx, 0:iNyp + 1) :: a, b, c
  complex(rkind), dimension(0:inx, 0:iNyp + 1) :: g

  complex(rkind) ocpack(4), icpack(4)

  do i = 0, inx

    if (rankY /= 0) then
      ! If we aren't the lowest process, then wait for data
      call mpi_recv(ocpack, 4, mpi_double_complex, rankY - 1, 13 &
                    , mpi_comm_y, status, ierror)
      ! Unpack the data
      a(i, 1) = real(ocpack(1))
      b(i, 1) = real(ocpack(2))
      c(i, 1) = real(ocpack(3))
      g(i, 1) = ocpack(4)
    end if

    do j = 2, iNyp
      a(i, j) = -a(i, j) / b(i, j - 1)
      b(i, j) = b(i, j) + a(i, j) * c(i, j - 1)
      g(i, j) = g(i, j) + a(i, j) * g(i, j - 1)
    end do

    if (rankY /= NprocY - 1) then
      icpack(1) = a(i, iNyp)
      icpack(2) = b(i, iNyp)
      icpack(3) = c(i, iNyp)
      icpack(4) = g(i, iNyp)
      call mpi_send(icpack, 4, mpi_double_complex, rankY + 1, 13 &
                    , mpi_comm_y, ierror)
    end if

  end do

  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|------
subroutine thomas_backward_real_mpi(a, b, c, g, iNyp, inx)
  !----*|--.---------.---------.---------.---------.---------.---------.-|------
  ! This subroutine performs the backward sweep of the Thomas algorithm
  ! Thomas algorithm solves Ax=b for tridiagonal A
  ! The RHS vector and solution are real
  ! Input lower, main, and upper diagonals, ld, md, ud, and rhs x
  ! Returns solution in x
  ! The indexing should be done by ROW, ie.
  ! ￼￼￼C[b1c1 0 C[a2 b2 c2 C[0a3b3
  ! 0 0... 0 0... c30...

  integer i, j, inx, iNyp
  real(rkind), dimension(0:inx - 1, 0:iNyp + 1) :: a, b, c
  real(rkind), dimension(0:inx - 1, 0:iNyp + 1) :: g
  real(rkind) icpack(1), ocpack(1)
  do i = 0, inx - 1
    if (rankY /= NprocY - 1) then
      ! If we aren’t the highest process, then wait for data
      call mpi_recv(ocpack, 1, mpi_double_precision, rankY + 1, 10 &
                    , mpi_comm_y, status, ierror)
      g(i, iNyp + 1) = ocpack(1)
    else
      ! Else, if we are the highest process, compute the solution at j=INyp
      g(i, iNyp + 1) = g(i, iNyp + 1) / b(i, iNyp + 1)
    end if
    ! All processes solve from INyp..1
    do j = iNyp, 0, -1
      g(i, j) = (g(i, j) - c(i, j) * g(i, j + 1)) / b(i, j)
    end do
    if (rankY /= 0) then
      icpack(1) = g(i, 2)
      call mpi_send(icpack, 1, mpi_double_precision, rankY - 1, 10 &
                    , mpi_comm_y, ierror)
    end if
  end do
  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|------
subroutine thomas_backward_complex_mpi(a, b, c, g, iNyp, inx)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-----|
  ! This subroutine performs the backward sweep of the Thomas algorithm
  ! Thomas algorithm solves Ax=b for tridiagonal A
  ! The RHS vector and solution are real
  ! Input lower, main, and upper diagonals, ld, md, ud, and rhs x
  ! Returns solution in x
  ! The indexing should be done by ROW, ie.
  ! [ b1  c1   0   0   0 ...
  ! [ a2  b2  c2   0   0 ...
  ! [  0  a3  b3   c3  0 ...

  integer i, j, inx, iNyp
  real(rkind), dimension(0:inx, 0:iNyp + 1) :: a, b, c
  complex(rkind), dimension(0:inx, 0:iNyp + 1) :: g

  do i = 0, inx

    if (rankY /= NprocY - 1) then
      ! If we aren't the highest process, then wait for data
      call mpi_recv(g(i, iNyp + 1), 1, mpi_double_complex, rankY + 1, 11 &
                    , mpi_comm_y, status, ierror)
      j = iNyp
      g(i, j) = (g(i, j) - c(i, j) * g(i, j + 1)) / b(i, j)
    else
      ! Else, if we are the highest process, then compute the solution at j=INyp
      g(i, iNyp) = g(i, iNyp) / b(i, iNyp)
    end if

    do j = iNyp - 1, 1, -1
      g(i, j) = (g(i, j) - c(i, j) * g(i, j + 1)) / b(i, j)
    end do

    if (rankY /= 0) then
      call mpi_send(g(i, 2), 1, mpi_double_complex, rankY - 1, 11 &
                    , mpi_comm_y, ierror)
    end if

  end do

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|------
subroutine apply_BC_1_lower(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-----

  integer i
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::  matl, matd, matu, vec

  ! Bottom Wall:
  if (u_BC_Ymin == 0) then
   ! Dirichlet
   do i = 0, Nxm1
     matl(i, 0) = 0.
     matd(i, 0) = 1.
     matu(i, 0) = 0.
     vec(i, 0) = 0.

     matl(i, 1) = 0.
     matd(i, 1) = 1.
     matu(i, 1) = 0.
     vec(i, 1) = u_BC_Ymin_c1
   end do
  else
   ! Neumann
   do i = 0, Nxm1
     matl(i, 0) = 0.
     matd(i, 0) = 1.
     matu(i, 0) = 0.
     vec(i, 0) = 0.
   end do
   do i = 0, Nxm1
     matl(i, 1) = 0.
     matd(i, 1) = -1.
     matu(i, 1) = 1.
     vec(i, 1) = dy(2) * u_BC_Ymin_c1
   end do

  end if

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|----
subroutine apply_BC_1_upper(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|--

  integer i
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! Top wall
  if (u_BC_Ymax == 0) then
   ! Dirichlet
   do i = 0, Nxm1
     matl(i, Nyp + 1) = 0.
     matd(i, Nyp + 1) = 1.
     matu(i, Nyp + 1) = 0.
     vec(i, Nyp + 1) = 0.

     matl(i, Nyp) = 0.
     matd(i, Nyp) = 1.
     matu(i, Nyp) = 0.
     vec(i, Nyp) = u_BC_Ymax_c1
   end do
  else
   ! Neumann
   do i = 0, Nxm1
     matl(i, Nyp) = -1.
     matd(i, Nyp) = 1.
     matu(i, Nyp) = 0.
     vec(i, Nyp) = dy(Nyp) * u_BC_Ymax_c1
   end do
   do i = 0, Nxm1
     matl(i, Nyp + 1) = 0.
     matd(i, Nyp + 1) = 1.
     matu(i, Nyp + 1) = 0.
     vec(i, Nyp + 1) = 0.
   end do

  end if

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|---
subroutine apply_BC_2_lower(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|--

  integer i
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! Bottom Wall:
  if (v_BC_Ymin == 0) then
   ! Dirichlet
   do i = 0, Nxm1
     matl(i, 1) = 0.d0
     matd(i, 1) = 1.d0
     matu(i, 1) = 0.d0
     vec(i, 1) = v_BC_Ymin_c1

     matl(i, 2) = 0.d0
     matd(i, 2) = 1.d0
     matu(i, 2) = 0.d0
     vec(i, 2) = v_BC_Ymin_c1
   end do
  else if (v_BC_Ymin == 1) then
   ! Neumann
   do i = 0, Nxm1
     matd(i, 1) = -1.d0
     matu(i, 1) = 1.d0
     matl(i, 1) = 0.d0
     vec(i, 1) = dyf(1) * v_BC_Ymin_c1
   end do
  end if

  ! The following is only a placeholder, this row is used for U1 and U3
  do i = 0, Nxm1
    matl(i, 0) = 0.
    matd(i, 0) = 1.
    matu(i, 0) = 0.
    vec(i, 0) = 0.
  end do

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_2_upper(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|--

  integer i
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! Top wall
  if (v_BC_Ymax == 0) then
   ! Dirichlet
   do i = 0, Nxm1
     matl(i, Nyp + 1) = 0.
     matd(i, Nyp + 1) = 1.
     matu(i, Nyp + 1) = 0.
     vec(i, Nyp + 1) = v_BC_Ymax_c1

     matl(i, Nyp) = 0.
     matd(i, Nyp) = 1.
     matu(i, Nyp) = 0.
     vec(i, Nyp) = v_BC_Ymax_c1
   end do
  else if (v_BC_Ymax == 1) then
   ! Neumann
   do i = 0, Nxm1
     matl(i, Nyp + 1) = -1.
     matd(i, Nyp + 1) = 1.
     matu(i, Nyp + 1) = 0.
     vec(i, Nyp + 1) = dyf(Nyp) * v_BC_Ymax_c1
   end do
  end if
  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_3_lower(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|--

  integer i
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! Bottom Wall:
  if (w_BC_Ymin == 0) then
   ! Dirichlet
   do i = 0, Nxm1
     matl(i, 0) = 0.
     matd(i, 0) = 1.
     matu(i, 0) = 0.
     vec(i, 0) = 0.

     matl(i, 1) = 0.
     matd(i, 1) = 1.
     matu(i, 1) = 0.
     vec(i, 1) = w_BC_Ymin_c1
   end do
  else
   ! Neumann
   do i = 0, Nxm1
     matl(i, 0) = 0.
     matd(i, 0) = 1.
     matu(i, 0) = 0.
     vec(i, 0) = 0.
   end do
   do i = 0, Nxm1
     matl(i, 1) = 0.
     matd(i, 1) = -1.
     matu(i, 1) = 1.
     vec(i, 1) = dy(2) * w_BC_Ymin_c1
   end do

  end if

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_3_upper(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|--

  integer i
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! Top wall
  if (w_BC_Ymax == 0) then
   ! Dirichlet
   do i = 0, Nxm1
     matl(i, Nyp + 1) = 0.
     matd(i, Nyp + 1) = 1.
     matu(i, Nyp + 1) = 0.
     vec(i, Nyp + 1) = 0.

     matl(i, Nyp) = 0.
     matd(i, Nyp) = 1.
     matu(i, Nyp) = 0.
     vec(i, Nyp) = w_BC_Ymax_c1
   end do
  else
   ! Neumann
   if (f_type == 4)  w_BC_Ymax_c1_transient = w_BC_Ymax_c1 + amp_omega0 * sin(omega0 * (Ro_inv/delta * time - force_start) ) ! force_start is the phase to start _down-front_ forcing
   do i = 0, Nxm1
     matl(i, Nyp) = -1.
     matd(i, Nyp) = 1.
     matu(i, Nyp) = 0.
     vec(i, Nyp) = dy(Nyp) * w_BC_Ymax_c1_transient
   end do
   do i = 0, Nxm1
     matl(i, Nyp + 1) = 0.
     matd(i, Nyp + 1) = 1.
     matu(i, Nyp + 1) = 0.
     vec(i, Nyp + 1) = 0.
   end do

  end if

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_th_lower(matl, matd, matu, vec, n)
  !----*|--.---------.---------.---------.---------.---------.---------.-|--

  integer i, n
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! Bottom Wall:
  if (th_BC_Ymin(n) == 0) then
   ! Dirichlet
   do i = 0, Nxm1
     matl(i, 0) = 0.
     matd(i, 0) = 1.
     matu(i, 0) = 0.
     vec(i, 0) = 0.

     matl(i, 1) = 0.
     matd(i, 1) = 1.
     matu(i, 1) = 0.
     vec(i, 1) = th_BC_Ymin_c1(n)
   end do
  else
   ! Neumann
   ! NOTE: BC enforced at GY(2)
   do i = 0, Nxm1
     matl(i, 1) = 0.
     matd(i, 1) = -1.
     matu(i, 1) = 1.
     vec(i, 1) = dy(2) * th_BC_Ymin_c1(n)
   end do
   do i = 0, Nxm1
     matl(i, 0) = 0.
     matd(i, 0) = 1.
     matu(i, 0) = 0.
     vec(i, 0) = 0.
   end do
  end if
  return

end


!----*|--.---------.---------.---------.---------.---------.---------.-|--
subroutine apply_BC_th_upper(matl, matd, matu, vec, n)
  !----*|--.---------.---------.---------.---------.---------.---------.-|--

  integer i, n
  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! Top wall
  if (th_BC_Ymax(n) == 0) then
   ! Dirichlet
   do i = 0, Nxm1
     matl(i, Nyp + 1) = 0.
     matd(i, Nyp + 1) = 1.
     matu(i, Nyp + 1) = 0.
     vec(i, Nyp + 1) = 0.

     matl(i, Nyp) = 0.
     matd(i, Nyp) = 1.
     matu(i, Nyp) = 0.
     vec(i, Nyp) = th_BC_Ymax_c1(n)
   end do
  else
   ! Neumann
   ! NOTE: BC enforced at GY(Nyp)
   do i = 0, Nxm1
     matl(i, Nyp) = -1.
     matd(i, Nyp) = 1.
     matu(i, Nyp) = 0.
     vec(i, Nyp) = dy(Nyp) * th_BC_Ymax_c1(n)
   end do
   do i = 0, Nxm1
     matl(i, Nyp + 1) = 0.
     matd(i, Nyp + 1) = 1.
     matu(i, Nyp + 1) = 0.
     vec(i, Nyp + 1) = 0.
   end do
  end if
  return
end





!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_u1_mpi(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions to the
  ! velocity field prior to the implicit solve

  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec


  ! We first need to check to see which processor we are, if we are
  ! the upper or lowermost process, then apply boundary conditions
  if (rankY == 0) then
   ! If we have the lowest plane, apply the boundary conditions
   call apply_BC_1_lower(matl, matd, matu, vec)
  end if
  if (rankY == NprocY - 1) then
   ! If we have the highest plane, apply the boundary conditions
   call apply_BC_1_upper(matl, matd, matu, vec)
  end if
  return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_u2_mpi(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions to the
  ! velocity field prior to the implicit solve

  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! We first need to check to see which processor we are, if we are
  ! the upper or lowermost process, then apply boundary conditions
  if (rankY == 0) then
   ! If we have the lowest plane, apply the boundary conditions
   call apply_BC_2_lower(matl, matd, matu, vec)
  end if
  if (rankY == NprocY - 1) then
   ! If we have the highest plane, apply the boundary conditions
   call apply_BC_2_upper(matl, matd, matu, vec)
  end if
  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_u3_mpi(matl, matd, matu, vec)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions to the
  ! velocity field prior to the implicit solve

  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec

  ! We first need to check to see which processor we are, if we are
  ! the upper or lowermost process, then apply boundary conditions
  if (rankY == 0) then
   ! If we have the lowest plane, apply the boundary conditions
   call apply_BC_3_lower(matl, matd, matu, vec)
  end if
  if (rankY == NprocY - 1) then
   ! If we have the highest plane, apply the boundary conditions
   call apply_BC_3_upper(matl, matd, matu, vec)
  end if
  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_th_mpi(matl, matd, matu, vec, n)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions to the
  ! scalar fields prior to the implicit solve

  real(rkind), dimension(0:Nx-1,0:Nyp+1) ::   matl, matd, matu, vec
  integer n

  ! We first need to check to see which processor we are, if we are
  ! the upper or lowermost process, then apply boundary conditions
  if (rankY == 0) then
   ! If we have the lowest plane, apply the boundary conditions
   call apply_BC_th_lower(matl, matd, matu, vec, n)
  end if
  if (rankY == NprocY - 1) then
   ! If we have the upper plane, apply the boundary conditions
   call apply_BC_th_upper(matl, matd, matu, vec, n)
  end if
  return
end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_rem_div_mpi(matl_c, matd_c, matu_c, vec_c, k)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions for the Poisson Eq.

  complex(rkind), dimension(0:Nxp,0:Nyp+1) :: vec_c
  real(rkind), dimension(0:Nxp,0:Nyp+1) ::    matl_c,matd_c,matu_c
  integer i, k

  ! We first need to check to see which processor we are, if we are
  ! the upper or lowermost process, then apply boundary conditions
  ! Apply the boundary conditions
  if (rankY == 0) then
   do i = 0, Nxp - 1
     ! Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
     if ((k == 0) .and. (i == 0) .and. (rankZ == 0)) then
       ! Otherwise the matrix will be singular
       ! Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
       matl_c(i, 1) = 0.
       matd_c(i, 1) = 1.
       matu_c(i, 1) = 0.
       vec_c(i, 1) = (0., 0.)
     else
       ! Use Dirichlet boundary conditions, dp/dz=0 at walls
       matl_c(i, 1) = 0.
       matd_c(i, 1) = 1.
       matu_c(i, 1) = -1.
       vec_c(i, 1) = (0., 0.)
     end if
   end do
  end if
  ! Apply the boundary conditions
  if (rankY == NprocY - 1) then
   do i = 0, Nxp - 1
     matl_c(i, Nyp) = 1.
     matd_c(i, Nyp) = -1.
     matu_c(i, Nyp) = 0.
     vec_c(i, Nyp) = (0., 0.)
   end do
  end if

  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine apply_BC_poisson_mpi(matl_c, matd_c, matu_c, vec_c, k)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! This subroutine applies the boundary conditions for the Poisson Eq.

  complex(rkind), dimension(0:Nxp,0:Nyp+1) :: vec_c
  real(rkind), dimension(0:Nxp,0:Nyp+1) ::    matl_c,matd_c,matu_c
  integer i, k
  ! We first need to check to see which processor we are, if we are
  ! the upper or lowermost process, then apply boundary conditions
  ! Use dirichlet boundary condition at the lower wall to
  ! prevent the tridiagonal matrix from becomming singular for i,k=0
  if (rankY == 0) then
   do i = 0, Nxp - 1
     if ((i == 0) .and. (k == 0) .and. (rankZ == 0)) then
       matd_c(i, 1) = 1.
       matu_c(i, 1) = 0.
       vec_c(i, 1) = (0., 0.)
     else
       ! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
       matd_c(i, 1) = 1.
       matu_c(i, 1) = -1.
       vec_c(i, 1) = (0., 0.)
     end if
   end do
  end if
  ! Use dirichlet boundary condition at the lower wall to
  ! prevent the tridiagonal matrix from becomming singular for i,k=0
  if (rankY == NprocY - 1) then
   do i = 0, Nxp - 1
     matd_c(i, Nyp) = -1.
     matl_c(i, Nyp) = 1.
     vec_c(i, Nyp) = (0., 0.)
   end do
  end if

  return
end
