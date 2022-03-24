
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine save_stats_chan(movie,final)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Computes domain-integrated and horizontally-integrated (X-Z) stats

  character(len=35) fname
  character(len=20) gname
  logical movie,final
  integer i, j, k, n


  ! Scalar diagnostics
  real(rkind) thsum(0:Nyp + 1)
  ! Store/write 2D slices
  real(rkind) varxy(0:Nxm1, 1:Nyp), varzy(0:Nzp - 1, 1:Nyp), varxz(0:Nxm1, 0:Nzp - 1)

  ! HDF5 writing
  real(rkind) Diag(1:Nyp)
  real(rkind) DiagX(0:Nxp - 1)

  if (rank == 0) &
    write (*, '("Saving Flow Statistics for Time Step       " I10)')  time_step

  call mpi_barrier(mpi_comm_world, ierror)
  call apply_BC_vel_mpi_post ! Apply BCs FIRST (it screws up ghost cells...)
  call apply_BC_th_mpi_post
  call ghost_chan_mpi
  call ghost_chan_mpi_j0 ! Need the j = 0 boundary filled for les output


  if (rank == 0) write (*, '("Time    = " ES12.5 "       dt = " ES12.5)') time, dt ! Note: dt is the physical / CFL-constrained time-step

  ! Store FF CUi in cr1(), and keep PP Ui in u1()
  do j = 0, Nyp + 1
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cr1(i, k, j) = cu1(i, k, j)
        cr2(i, k, j) = cu2(i, k, j)
        cr3(i, k, j) = cu3(i, k, j)
        do n = 1, N_th
          crth(i, k, j, n) = cth(i, k, j, n)
        end do
      end do
    end do
  end do

  ! Convert to physical space
  call fft_xz_to_physical(cu1, u1)
  call fft_xz_to_physical(cu2, u2)
  call fft_xz_to_physical(cu3, u3)
  do n = 1, N_th
    call fft_xz_to_physical(cth(:, :, :, n), th(:, :, :, n))
  end do



  ! Compute ume(y), etc
  ! Also, computes Z-Averages ume_xy(x,y), etc, if not homogeneousX
  !   Otherwise, puts ume into ume_xy, etc
  call compute_averages(movie)


  !!! Dissipation Rate !!!
  call compute_TKE_diss(movie)
  call compute_MKE_diss(movie)
  if (use_LES) then
    call ghost_les_mpi ! Share nu_t
    call compute_TKE_diss_les
  end if

  !!! TKE / RMS Velocities !!!
  call compute_TKE(movie)

  !!! Production and Reynolds Stresses !!!
  call compute_TKE_Production(movie)


  !!! MKE (on GY Grid) !!!
  do j = 2, Nyp
    mke(j) = 0.d0
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        mke(j) = mke(j) &
                    +    (dyf(j - 1) * ume_xy(i, j)**2.d0 &
                         + dyf(j) * ume_xy(i, j - 1)**2.d0) / (2.d0 * dy(j)) &
                    +    (dyf(j - 1) * (wme_xy(i, j) + & ! Include the TWS!
                                    (1.d0 / (Ro_inv / delta)) * dTHdX(1) * (gyf(j) - 0.5d0*Ly))**2.d0 &
                         + dyf(j) * (wme_xy(i, j - 1) + &
                                         (1.d0 / (Ro_inv / delta)) * dTHdX(1) * (gyf(j - 1) - 0.5d0*Ly))**2.d0) / (2.d0 * dy(j)) &
                    +    vme_xy(i, j)**2.d0
      end do
    end do
  end do
  call mpi_allreduce(mpi_in_place, mke, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  mke = 0.5 * mke / float(Nx * Nz)


  !!! Gradient of _Mean_ Velocity !!!
  do j = 1, Nyp
    dudy(j) = (ume(j) - ume(j - 1)) / (gyf(j) - gyf(j - 1))
    dwdy(j) = (wme(j) - wme(j - 1)) / (gyf(j) - gyf(j - 1))
  end do


  !!! Mean Square Shear !!!
  do j = 1, Nyp
    shear(j) = 0.d0
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        f1(i, k, j) = ((u1(i, k, j + 1) - u1(i, k, j - 1)) / (2.d0 * dyf(j)))**2.d0 &
                    + ((u3(i, k, j + 1) - u3(i, k, j - 1)) / (2.d0 * dyf(j)))**2.d0
        shear(j) = shear(j) + f1(i, k, j)
      end do
    end do
  end do
  call mpi_allreduce(mpi_in_place, shear, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  shear = shear / float(Nx * Nz)

  gname = 'shear2_zstar'
  call Bin_Ystar_and_Write(gname, f1)



  call compute_Vorticity(movie)



  !!! Write Mean Stats f(y) !!!
  fname = 'mean.h5'
  gname = 'time'
  call WriteHDF5_real(fname, gname, time)

  if (rankZ == 0) then

    gname = 'gzf'
    Diag = gyf(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'ume'
    Diag = ume(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'wme'
    Diag = vme(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'vme'
    Diag = wme(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'mke'
    Diag = mke(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'urms'
    Diag = urms(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'wrms'
    Diag = vrms(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'vrms'
    Diag = wrms(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'uw'
    Diag = uv(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'uv'
    Diag = uw(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'wv'
    Diag = wv(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'uu_dudx'
    Diag = uu_dudx(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'ww_dwdz'
    Diag = vv_dvdy(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'vu_dvdx'
    Diag = wu_dwdx(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'wu_dwdx'
    Diag = vu_dvdx(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'uw_dudz'
    Diag = uv_dudy(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'vw_dvdz'
    Diag = wv_dwdy(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'dudz'
    Diag = dudy(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'dvdz'
    Diag = dwdy(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'cp'
    Diag = dble(cp(0, 0, 1:Nyp))
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'shear'
    Diag = shear(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'omega_x'
    Diag = omega_x(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'omega_z'
    Diag = omega_y(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

    gname = 'omega_y'
    Diag = omega_z(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)

  end if







  !!! Iterate through all TH Statistics !!!
  do n = 1, N_th
    ! Store FF CTH crth(), and keep PP in th() (Already done above)

    ! Compute the TH Gradients and store in CRi
    do j = 1, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1 ! Nkx
          ! Store gradients of TH(:,:,:,n) (if it is used) in CRi
          cr1(i, k, j) = cikx(i) * crth(i, k, j, n)
          cr2(i, k, j) = (crth(i, k, j + 1, 1) - crth(i, k, j - 1, 1)) / (gyf(j + 1) - gyf(j - 1))
          cr3(i, k, j) = cikz(k) * crth(i, k, j, n)
        end do
      end do
    end do
    ! Convert gradients to physical space
    call fft_xz_to_physical(cr1, r1)
    call fft_xz_to_physical(cr2, r2)
    call fft_xz_to_physical(cr3, r3)
    ! (Already have th in PP)


    !!! RMS TH !!!
    thvar_xy = 0.
    do j = 1, Nyp
      thsum(j) = 0.
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          thsum(j) = thsum(j) + (th(i, k, j, n) - thme_xy(i, j, n))**2.
          thvar_xy(i, j) = thvar_xy(i, j) + (th(i, k, j, n) - thme_xy(i, j, n))**2.
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, thsum, (Nyp + 2), &
                       mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

    thrms(:, n) = sqrt(thsum / float(Nx * Nz))
    thvar_xy = thvar_xy / float(Nz) ! Can't take sqrt, then sum next...

    if (n == 1 .and. movie .and. Nz > 1) then
      fname = 'mean_xz.h5'
      gname = 'thth_xz'
      call reduce_and_write_XYplane(fname, gname, thvar_xy, .false., movie)
    end if



    !!! TH Reynolds Stress, th*v -- i.e. Buoyancy Production !!!
    uvar_xy = 0.
    vvar_xy = 0.
    do j = 1, Nyp
      thsum(j) = 0.
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          thsum(j) = thsum(j) + (th(i, k, j, n) - thme_xy(i, j, n)) &
                     * (0.5 * (u2(i, k, j) + u2(i, k, j + 1)) &
                        - 0.5 * (vme_xy(i, j) + vme_xy(i, j + 1)))

          uvar_xy(i, j) = uvar_xy(i, j) + (th(i, k, j, 1) - thme_xy(i, j, n)) &
                                        * (u1(i, k, j) - ume_xy(i, j))

          f1(i, k, j) = (th(i, k, j, 1) - thme_xy(i, j, n)) &
                                    * 0.5d0 * ((u2(i, k, j) - vme_xy(i, j)) &
                                             + (u2(i, k, j + 1) - vme_xy(i, j + 1)))
          vvar_xy(i, j) = vvar_xy(i, j) + f1(i, k, j)
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, thsum, (Nyp + 2), &
                       mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

    thv(:, n) = thsum / float(Nx * Nz)
    uvar_xy = uvar_xy / float(Nz)
    vvar_xy = vvar_xy / float(Nz)


    if (n == 1 .and. movie .and. Nz > 1) then
      fname = 'mean_xz.h5'
      gname = 'thu_xz'
      call reduce_and_write_XYplane(fname, gname, uvar_xy, .false., movie)
      gname = 'thw_xz'
      call reduce_and_write_XYplane(fname, gname, vvar_xy, .false., movie)
    end if

    gname = 'thw_zstar'
    call Bin_Ystar_and_Write(gname, f1)





    !!! TH Mean Production, th_m*v_m (with _full_ buoyancy) !!!
    do j = 1, Nyp
      thsum(j) = 0.
        do i = 0, Nxm1
          thsum(j) = thsum(j) + (thme_xy(i, j, n) + dTHdX(n)*(gx(i) - 0.5*Lx)) &
                     * (0.5 * (vme_xy(i, j) + vme_xy(i, j + 1)))
        end do
    end do

    thv_m(:, n) = thsum / float(Nx)



    !!! Gradient of _Mean_ TH !!!
    do j = 1, Nyp
      dthdy(j, n) = (thme(j, n) - thme(j - 1, n)) / (gyf(j) - gyf(j - 1))
    end do


    !!! PE Dissipation (Chi!), grad(TH) \cdot grad(TH) !!!
    ! Store |grad b|^2 in r1
    vvar_xy = 0.
    do j = 1, Nyp
      thsum(j) = 0.d0
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          r1(i, k, j) =  (r1(i, k, j) + dTHdX(n))**2.d0 &
                       +  r2(i, k, j)**2.d0 &
                       + (r3(i, k, j) + dTHdZ(n))**2.d0
          vvar_xy(i, j) = vvar_xy(i, j) + r1(i, k, j)
          thsum(j)    = thsum(j) + r1(i, k, j)
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, thsum, (Nyp + 2), &
                       mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

    pe_diss(:, n) = thsum / float(Nx * Nz) ! NOT actually PE dissipation -- just (grad TH)^2
    vvar_xy = vvar_xy / float(Nz)

    if (n == 1 .and. movie .and. Nz > 1) then
      fname = 'mean_xz.h5'
      gname = 'chi_xz'
      call reduce_and_write_XYplane(fname, gname, vvar_xy, .false., movie)
    end if

    gname = 'chi_zstar'
    call Bin_Ystar_and_Write(gname, r1)



    if (n == 1) then
      call compute_BPE

      !!! Compute integrated y*u1 at the left boundary !!!
      do j = 1, Nyp
        u1y_left(j) = 0.d0
        do k = 0, Nzp - 1
          u1y_left(j) = u1y_left(j) + gyf(j) * u1(0, k, j)
        end do
      end do
      call mpi_allreduce(mpi_in_place, u1y_left, (Nyp + 2), &
                         mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

      u1y_left = u1y_left / float(Nz)
      if (rankZ == 0) then
        call integrate_y_var(u1y_left, u1y_left_b)
      end if

    end if



    !!! Write Movie TH Slices !!!
    if (movie) then
      if (rank == 0) &
        write (*, '("Saving Movie Slice Output")')

      fname = 'movie.h5'
      call mpi_barrier(mpi_comm_world, ierror)
      if (rankZ == rankzmovie) then
        do j = 1, Nyp
          do i = 0, Nxm1
            varxy(i, j) = th(i, NzMovie, j, n)
          end do
        end do
        write (gname,'("th", I0.1 "_xz")') n
        call WriteHDF5_XYplane(fname, gname, varxy)
      end if

      if (rankY == rankymovie) then
        do j = 0, Nzp - 1
          do i = 0, Nxm1
            varxz(i, j) = th(i, j, NyMovie, n)
          end do
        end do
        write (gname,'("th", I0.1 "_xy")') n
        call WriteHDF5_XZplane(fname, gname, varxz)
      end if

      do j = 1, Nyp
        do i = 0, Nzp - 1
          varzy(i, j) = th(NxMovie, i, j, n)
        end do
      end do
      write (gname,'("th", I0.1 "_yz")') n
      call WriteHDF5_ZYplane(fname, gname, varzy)
    end if

  end do ! Over passive scalars, n





  !!! Write Mean TH Stats f(y) !!!
  fname = 'mean.h5'
  if (rankZ == 0) then
    do n = 1, N_th
      Diag = thme(1:Nyp, n)
      write (gname,'("thme", I0.2)') n
      call WriteStatH5_Y(fname, gname, Diag)

      Diag = dthdy(1:Nyp, n)
      write (gname,'("dthdz", I0.2)') n
      call WriteStatH5_Y(fname, gname, Diag)

      Diag = thrms(1:Nyp, n)
      write (gname,'("thrms", I0.2)') n
      call WriteStatH5_Y(fname, gname, Diag)

      Diag = thv(1:Nyp, n) ! \bar{w'b'}
      write (gname,'("thw", I0.2)') n
      call WriteStatH5_Y(fname, gname, Diag)

      Diag = thv_m(1:Nyp, n) ! \bar{w}\bar{b} + BG b!
      write (gname,'("thw", I0.2, "_m")') n
      call WriteStatH5_Y(fname, gname, Diag)

      Diag = pe_diss(1:Nyp, n)
      write (gname,'("pe_diss", I0.2)') n
      call WriteStatH5_Y(fname, gname, Diag)
    end do
  end if
  gname = 'u1z_x0'
  call WriteHDF5_real(fname, gname, u1y_left_b)




  !!! Write All Movie Slices !!!
  if (movie) then
    fname = 'movie.h5'
    call mpi_barrier(mpi_comm_world, ierror)
    if (rankZ == rankzmovie) then
      do j = 1, Nyp
        do i = 0, Nxm1
          varxy(i, j) = u1(i, NzMovie, j)
        end do
      end do
      gname = 'u_xz'
      call WriteHDF5_XYplane(fname, gname, varxy)
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    if (rankZ == rankzmovie) then
      do j = 1, Nyp
        do i = 0, Nxm1
          ! Interpolate onto the GYF grid
          varxy(i, j) = 0.5 * (u2(i, NzMovie, j) + u2(i, NzMovie, j + 1))
        end do
      end do
      gname = 'w_xz'
      call WriteHDF5_XYplane(fname, gname, varxy)
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    if (rankZ == rankzmovie) then
      do j = 1, Nyp
        do i = 0, Nxm1
          varxy(i, j) = u3(i, NzMovie, j)
        end do
      end do
      gname = 'v_xz'
      call WriteHDF5_XYplane(fname, gname, varxy)
    end if

    if (use_LES) then
      call mpi_barrier(mpi_comm_world, ierror)
      if (rankZ == rankzmovie) then
        do j = 1, Nyp
          do i = 0, Nxm1
            varxy(i, j) = nu_t(i, NzMovie, j)
          end do
        end do
        gname = 'nu_t_xz'
        call WriteHDF5_XYplane(fname, gname, varxy)
      end if
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    if (rankY == rankymovie) then
      do j = 0, Nzp - 1
        do i = 0, Nxm1
          varxz(i, j) = u1(i, j, NyMovie)
        end do
      end do
      gname = 'u_xy'
      call WriteHDF5_XZplane(fname, gname, varxz)
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    if (rankY == rankymovie) then
      do j = 0, Nzp - 1
        do i = 0, Nxm1
          varxz(i, j) = 0.5 * (u2(i, j, NyMovie) + u2(i, j, NyMovie + 1))
        end do
      end do
      gname = 'w_xy'
      call WriteHDF5_XZplane(fname, gname, varxz)
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    if (rankY == rankymovie) then
      do j = 0, Nzp - 1
        do i = 0, Nxm1
          varxz(i, j) = u3(i, j, NyMovie)
        end do
      end do
      gname = 'v_xy'
      call WriteHDF5_XZplane(fname, gname, varxz)
    end if

    if (use_LES) then
      call mpi_barrier(mpi_comm_world, ierror)
      if (rankY == rankymovie) then
        do j = 0, Nzp - 1
          do i = 0, Nxm1
            varxz(i, j) = nu_t(i, j, NyMovie)
          end do
        end do
        gname = 'nu_t_xy'
        call WriteHDF5_XZplane(fname, gname, varxz)
      end if
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    do j = 1, Nyp
      do i = 0, Nzp - 1
        varzy(i, j) = u1(NxMovie, i, j)
      end do
    end do
    gname = 'u_yz'
    call WriteHDF5_ZYplane(fname, gname, varzy)

    call mpi_barrier(mpi_comm_world, ierror)
    do j = 1, Nyp
      do i = 0, Nzp - 1
        varzy(i, j) = 0.5 * (u2(NxMovie, i, j) + u2(NxMovie, i, j + 1))
      end do
    end do
    gname = 'w_yz'
    call WriteHDF5_ZYplane(fname, gname, varzy)

    call mpi_barrier(mpi_comm_world, ierror)
    do j = 1, Nyp
      do i = 0, Nzp - 1
        varzy(i, j) = u3(NxMovie, i, j)
      end do
    end do
    gname = 'v_yz'
    call WriteHDF5_ZYplane(fname, gname, varzy)

    call mpi_barrier(mpi_comm_world, ierror)
    if (use_LES) then
      do j = 1, Nyp
        do i = 0, Nzp - 1
          varzy(i, j) = nu_t(NxMovie, i, j)
        end do
      end do
      gname = 'nu_t_yz'
      call WriteHDF5_ZYplane(fname, gname, varzy)
    end if

    call mpi_barrier(mpi_comm_world, ierror)
  end if ! END IF MOVIE







  ! Convert velocity back to Fourier space
  call fft_xz_to_fourier(u1, cu1)
  call fft_xz_to_fourier(u2, cu2)
  call fft_xz_to_fourier(u3, cu3)
  do n = 1, N_th
    call fft_xz_to_fourier(th(:, :, :, n), cth(:, :, :, n))
  end do




  !!! Spectrum of FT(u)^2 (Y-Averaged) !!!

  do j = 0, Nyp + 1
    do i = 0, Nxp - 1 ! Nkx
      cuu1_yx(j, i) = cu1(i, 0, j) ! Include both horizontal mean & structure of meanY
    end do
  end do

  cuu1_yx = cuu1_yx*conjg(cuu1_yx)

  DiagX = 0.
  do i = 0, Nxp - 1
    do j = 2, Nyp
      DiagX(i) = DiagX(i) + 0.5 * (cuu1_yx(j, i) + cuu1_yx(j - 1, i)) * dy(j) ! Integrate cuu1 in y
    end do
  end do
  call mpi_allreduce(mpi_in_place, DiagX, Nxp, &
                     mpi_double_precision, mpi_sum, mpi_comm_y, ierror)
  DiagX = DiagX / Ly

  if (rankY == 0) then
    fname = 'mean.h5'
    gname = 'FTx_uu'
    call WriteStatH5_X(fname, gname, DiagX, Nxp) ! Resulting size is Nx/2
  end if



  call mpi_barrier(mpi_comm_world, ierror)
  return
end





!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine compute_averages(movie)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Calculates horizontal averages
  !   and also z-line averages (if homogeneousX == .false.)
  ! NOTE: Assumes Ui are PP

  character(len=35) fname
  character(len=20) gname
  logical movie
  integer i, j, k, n
  real(rkind) ubulk


  !!! Horizontally-Averaged Velocities ume(y) (and Broadcast), and Write Bulk U !!!
  if (rankZ == 0) then
    ume = dble(cr1(0, 0, :))
    vme = dble(cr2(0, 0, :)) ! Still on GY
    wme = dble(cr3(0, 0, :))
    do n = 1, N_th
      thme(:, n) = dble(crth(0, 0, :, n))
    end do
  end if

  call integrate_y_var(ume, ubulk)
  if (rank == 0) write (*,  '("U Bulk  = " ES26.18)') ubulk

  call mpi_bcast(ume, Nyp + 2, mpi_double_precision, 0, &
                 mpi_comm_z, ierror)
  call mpi_bcast(vme, Nyp + 2, mpi_double_precision, 0, &
                 mpi_comm_z, ierror)
  call mpi_bcast(wme, Nyp + 2, mpi_double_precision, 0, &
                 mpi_comm_z, ierror)
  call mpi_bcast(thme, (Nyp + 2) * N_th, &
                 mpi_double_precision, 0, mpi_comm_z, ierror)


  if (Nz > 1 .and. .not. homogeneousX) then ! If 3D and if x-direction is not homogeneous...

    !!! XY Mean Velocities and Buoyancy (and keep) !!!
    do j = 1, Nyp
      do i = 0, Nxm1
        ume_xy(i, j) = 0.d0
        var_xy(i, j) = 0.d0
        wme_xy(i, j) = 0.d0
        do k = 0, Nzp - 1
          ume_xy(i, j) = ume_xy(i, j) + u1(i, k, j)
          var_xy(i, j) = var_xy(i, j) + 0.5d0 * (u2(i, k, j) + u2(i, k, j + 1)) ! Mean on GYF (to write out)
          wme_xy(i, j) = wme_xy(i, j) + u3(i, k, j)
        end do
      end do
    end do

    do n = 1, N_th
      do j = 1, Nyp
        do i = 0, Nxm1
          thme_xy(i, j, n) = 0.d0
          do k = 0, Nzp - 1
            thme_xy(i, j, n) = thme_xy(i, j, n) + th(i, k, j, n)
          end do
        end do
      end do
    end do

    vme_xy = 0.d0
    do j = 1, Nyp + 1
      do i = 0, Nxm1
        do k = 0, Nzp - 1
          vme_xy(i, j) = vme_xy(i, j) + u2(i, k, j) ! Mean on GY (To use in code)
        end do
      end do
    end do

    ume_xy = ume_xy / float(Nz)
    vme_xy = vme_xy / float(Nz)
    var_xy = var_xy / float(Nz)
    wme_xy = wme_xy / float(Nz)
    thme_xy = thme_xy / float(Nz)

    ! Use allreduce so that each process in comm_z has the mean
    call mpi_allreduce(mpi_in_place, vme_xy, Nx * (Nyp + 1), &
                       mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

    fname = 'mean_xz.h5'

    gname = 'ume_xz'
    call reduce_and_write_XYplane(fname, gname, ume_xy, .true., movie)
    gname = 'wme_xz'
    call reduce_and_write_XYplane(fname, gname, var_xy, .false., movie)
    gname = 'vme_xz'
    call reduce_and_write_XYplane(fname, gname, wme_xy, .true., movie)
    gname = 'thme_xz'
    call reduce_and_write_XYplane(fname, gname, thme_xy(:,:,1), .true., movie)

    if (N_th > 1) then
      do n = 2, N_th
        call reduce_and_write_XYplane(fname, gname, thme_xy(:,:,n), .true., .false.)
      end do
    end if

  else
    ! Put ume into ume_xy, etc

    do j = 1, Nyp
      do i = 0, Nxm1
        ume_xy(i, j) = ume(j)
        wme_xy(i, j) = wme(j)
      end do
    end do
    do j = 1, Nyp + 1
      do i = 0, Nxm1
        vme_xy(i, j) = vme(j)
      end do
    end do
    do n = 1, N_th
      do j = 1, Nyp
        do i = 0, Nxm1
          thme_xy(i, j, n) = thme(j, n)
        end do
      end do
    end do

  end if

end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine compute_TKE_diss(movie)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Calculate the turbulent dissipation rate, epsilon
  ! Note that this is actually the pseudo-dissipation (see Pope, Turb. Flows)
  ! for an explanation
  ! Assumes that FF Ui are stored in cri(), and PP Ui are in ui()

  character(len=35) fname
  character(len=20) gname
  logical movie
  real(rkind) Diag(1:Nyp)
  real(rkind) varxy(0:Nxm1, 1:Nyp), varzy(0:Nzp - 1, 1:Nyp), varxz(0:Nxm1, 0:Nzp - 1)
  integer i, j, k
  integer k_start

  ! Store the 3D dissipation rate in F1
  f1 = 0.d0

  ! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
  ! epsilon will be calculated on the GY grid (2:Nyp)
  !  This is so that it remains conserved (as in the code)
  epsilon = 0.

  if (.not. homogeneousX) then
    k_start = 1 ! Skip the mean component!
  else
    k_start = 0
  end if


  ! Compute du/dx  NOTE: Remove mean if not homogeneousX (k_start)
  if (k_start == 1) cs1 = 0
  do j = 1, Nyp
    do k = k_start, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikx(i) * cr1(i, k, j)
      end do
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                   + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
        f1(i, k, j) = f1(i, k, j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                     + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
      end do
    end do
  end do


  ! Compute dv/dx  NOTE: Remove mean if not homogeneousX (k_start)
  if (k_start == 1) cs1 = 0
  do j = 2, Nyp
    do k = k_start, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikx(i) * cr2(i, k, j)
      end do
    end do
  end do
  ! Convert to physical space
  call fft_xz_to_physical(cs1, s1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do


  ! Compute du/dy at GY gridpoints  NOTE: Remove mean
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = ((u1(i, k, j) - ume_xy(i, j)) &
                       - (u1(i, k, j - 1) - ume_xy(i, j - 1))) &
                      / dy(j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do


  ! Compute dw/dx  NOTE: Remove mean if not homogeneousX (k_start)
  if (k_start == 1) cs1 = 0
  do j = 1, Nyp
    do k = k_start, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikx(i) * cr3(i, k, j)
      end do
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                   + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
        f1(i, k, j) = f1(i, k, j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                     + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
      end do
    end do
  end do


  ! Compute du/dz at GY gridpoints
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikz(k) * cr1(i, k, j)
      end do
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                   + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
        f1(i, k, j) = f1(i, k, j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                     + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
      end do
    end do
  end do


  ! Compute dv/dy at GY gridpoints  NOTE: Remove mean
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = ((u2(i, k, j + 1) - vme_xy(i, j + 1)) &
                       - (u2(i, k, j - 1) - vme_xy(i, j - 1))) &
                      / (gy(j + 1) - gy(j - 1))
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do


  ! Compute dw/dy at GY gridpoints  NOTE: Remove mean
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = ((u3(i, k, j) - wme_xy(i, j)) &
                       - (u3(i, k, j - 1) - wme_xy(i, j - 1))) &
                      / dy(j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do


  ! Compute dv/dz
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikz(k) * cr2(i, k, j)
      end do
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do


  ! Compute dw/dz
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikz(k) * cr3(i, k, j)
      end do
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                   + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
        f1(i, k, j) = f1(i, k, j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                     + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
      end do
    end do
  end do


  epsilon = nu * epsilon / float(Nx * Nz)
  f1 = nu * f1
  call mpi_allreduce(mpi_in_place, epsilon, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)




  ! Write mean / movie / mean_xy
  fname = 'mean.h5'
  if (rankZ == 0) then
    gname = 'epsilon'
    Diag = epsilon(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)
  end if

  gname = 'epsilon_zstar'
  call Bin_Ystar_and_Write(gname, f1)


  if (movie) then

    fname = 'movie.h5'
    call mpi_barrier(mpi_comm_world, ierror)
    if (rankZ == rankzmovie) then
      do i = 0, Nxm1
        do j = 1, Nyp
          varxy(i, j) = f1(i, NzMovie, j)
        end do
      end do
      gname = 'epsilon_xz'
      call WriteHDF5_XYplane(fname, gname, varxy)
    end if
    call mpi_barrier(mpi_comm_world, ierror)
    if (rankY == rankymovie) then
      do i = 0, Nxm1
        do j = 0, Nzp - 1
          varxz(i, j) = f1(i, j, NyMovie)
        end do
      end do
      gname = 'epsilon_xy'
      call WriteHDF5_XZplane(fname, gname, varxz)
    end if
    call mpi_barrier(mpi_comm_world, ierror)
    do i = 0, Nzp - 1
      do j = 1, Nyp
        varzy(i, j) = f1(NxMovie, i, j)
      end do
    end do
    gname = 'epsilon_yz'
    call WriteHDF5_ZYplane(fname, gname, varzy)


    ! Mean XY Slice
    do j = 1, Nyp
      do i = 0, Nxm1
        uvar_xy(i, j) = 0.d0
        do k = 0, Nzp - 1
          uvar_xy(i, j) = uvar_xy(i, j) + f1(i, k, j)
        end do
      end do
    end do
    uvar_xy = uvar_xy / float(Nz)

    fname = 'mean_xz.h5'
    gname = 'epsilon_xz'
    if (movie) &
      call reduce_and_write_XYplane(fname, gname, uvar_xy, .false., movie)

  end if

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine compute_MKE_diss(movie)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Calculate the mean dissipation rate, epsilon_m
  ! Note that this is actually the pseudo-dissipation (see Pope, Turb. Flows)
  ! Assumes that FF Ui are stored in cri(), and PP Ui are in ui()

  character(len=35) fname
  character(len=20) gname
  logical movie
  real(rkind) Diag(1:Nyp)
  real(rkind) varxy(0:Nxm1, 1:Nyp), tempxy(0:Nxm1, 1:Nyp)
  integer i, j

  ! Store the 2D dissipation rate in varxy
  varxy = 0.

  ! Compute the mean dissipation rate, epsilon_m=nu*<du_i/dx_j du_i/dx_j>
  ! epsilon will be calculated on the GY grid (2:Nyp)
  epsilon_m = 0.

  if (homogeneousX) then
    return
  end if


  ! Compute du/dx
  cs1 = 0
  do j = 1, Nyp
    do i = 0, Nxp - 1
      cs1(i, 0, j) = cikx(i) * cr1(i, 0, j)
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  do j = 1, Nyp
    do i = 0, Nxm1
      epsilon_m(j) = epsilon_m(j) + (dyf(j - 1) * s1(i, 0, j)**2.d0 &
                                 + dyf(j) * s1(i, 0, j - 1)**2.d0) / (2.d0 * dy(j))
      varxy(i, j) = varxy(i, j) + (dyf(j - 1) * s1(i, 0, j)**2.d0 &
                                   + dyf(j) * s1(i, 0, j - 1)**2.d0) / (2.d0 * dy(j))

      dudx_m(i, j) = s1(i, 0, j) ! Store for compute_TKE_Production
    end do
  end do


  ! Compute dv/dx
  cs1 = 0
  do j = 1, Nyp
    do i = 0, Nxp - 1
      cs1(i, 0, j) = cikx(i) * cr2(i, 0, j)
    end do
  end do
  ! Convert to physical space
  call fft_xz_to_physical(cs1, s1)
  do j = 1, Nyp
    do i = 0, Nxm1
      epsilon_m(j) = epsilon_m(j) + (s1(i, 0, j)**2.0)
      varxy(i, j) = varxy(i, j) + (s1(i, 0, j)**2.0)

      dvdx_m(i, j) = s1(i, 0, j) ! Store for compute_TKE_Production
    end do
  end do


  ! Compute du/dy at GY gridpoints
  do j = 2, Nyp
    do i = 0, Nxm1
      tempxy(i, j) = (ume_xy(i, j) - ume_xy(i, j - 1))  / dy(j)
    end do
  end do
  do j = 2, Nyp
    do i = 0, Nxm1
      epsilon_m(j) = epsilon_m(j) + (tempxy(i, j)**2.0)
      varxy(i, j) = varxy(i, j) + (tempxy(i, j)**2.0)
    end do
  end do


  ! Compute dw/dx
  cs1 = 0
  do j = 1, Nyp
    do i = 0, Nxp - 1
      cs1(i, 0, j) = cikx(i) * cr3(i, 0, j)
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  do j = 1, Nyp
    do i = 0, Nxm1
      epsilon_m(j) = epsilon_m(j) + (dyf(j - 1) * s1(i, 0, j)**2.d0 &
                                 + dyf(j) * s1(i, 0, j - 1)**2.d0) / (2.d0 * dy(j))
      varxy(i, j) = varxy(i, j) + (dyf(j - 1) * s1(i, 0, j)**2.d0 &
                                   + dyf(j) * s1(i, 0, j - 1)**2.d0) / (2.d0 * dy(j))

      dwdx_m(i, j) = s1(i, 0, j) ! Store for compute_TKE_Production
    end do
  end do


  ! Compute dv/dy at GY gridpoints
  do j = 2, Nyp
    do i = 0, Nxm1
      tempxy(i, j) = (vme_xy(i, j + 1) - vme_xy(i, j - 1)) &
                    / (gy(j + 1) - gy(j - 1))
    end do
  end do
  do j = 2, Nyp
    do i = 0, Nxm1
      epsilon_m(j) = epsilon_m(j) + (tempxy(i, j)**2.0)
      varxy(i, j) = varxy(i, j) + (tempxy(i, j)**2.0)
    end do
  end do


  ! Compute dw/dy at GY gridpoints
  do j = 2, Nyp
    do i = 0, Nxm1
      tempxy(i, j) = (wme_xy(i, j) - wme_xy(i, j - 1)) / dy(j) + &
                                (1.d0 / (Ro_inv / delta)) * dTHdX(1) ! Include TWS!
    end do
  end do
  do j = 2, Nyp
    do i = 0, Nxm1
      epsilon_m(j) = epsilon_m(j) + (tempxy(i, j)**2.0)
      varxy(i, j) = varxy(i, j) + (tempxy(i, j)**2.0)
    end do
  end do





  epsilon_m = nu * epsilon_m / float(Nx)
  varxy = nu * varxy


  ! Write mean / movie / mean_xy
  fname = 'mean.h5'
  if (rankZ == 0) then
    gname = 'epsilon_m'
    Diag = epsilon_m(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)
  end if

  if (movie) then
    fname = 'mean_xz.h5'
    gname = 'epsilon_m_xz'
    if (rankZ == 0) then
      call WriteHDF5_XYplane(fname, gname, varxy)
    end if
  end if

  return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine compute_TKE(movie)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Compute the RMS velocities
  !  Both the horizontal averages -- ume(y), etc
  !   and the z-averages -- ume_xz(x,y), etc
  ! Assumes PP Ui

  character(len=35) fname
  character(len=20) gname
  logical movie
  integer i, j, k



  !!! TKE f(y) & f(x,y) !!!
  ! urms and wrms are on the GYF grid
  ! vrms is defined on the GY grid
  uu_xy = 0.
  vv_xy = 0.
  ww_xy = 0.
  f1 = 0.

  do j = 1, Nyp
    urms(j) = 0.
    vrms(j) = 0.
    wrms(j) = 0.
    do k = 0, Nzp - 1
      do i = 0, Nxm1

        f1(i, k, j) = 0.5d0 * (u1(i, k, j) - ume_xy(i, j))**2. &
                    + 0.5d0 * (u2(i, k, j) - vme_xy(i, j))**2. &
                    + 0.5d0 * (u3(i, k, j) - wme_xy(i, j))**2.

        urms(j) = urms(j) + (u1(i, k, j) - ume_xy(i, j))**2.
        uu_xy(i, j) = uu_xy(i, j) + (u1(i, k, j) - ume_xy(i, j))**2.

        vrms(j) = vrms(j) + (u2(i, k, j) - vme_xy(i, j))**2.
        vv_xy(i, j) = vv_xy(i, j) + (u2(i, k, j) - vme_xy(i, j))**2.

        wrms(j) = wrms(j) + (u3(i, k, j) - wme_xy(i, j))**2.
        ww_xy(i, j) = ww_xy(i, j) + (u3(i, k, j) - wme_xy(i, j))**2.

      end do
    end do
  end do


  ! Communicate Horizontal Mean (Save them later)

  call mpi_allreduce(mpi_in_place, urms, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)
  call mpi_allreduce(mpi_in_place, vrms, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)
  call mpi_allreduce(mpi_in_place, wrms, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  urms = urms / float(Nx * Nz)
  vrms = vrms / float(Nx * Nz)
  wrms = wrms / float(Nx * Nz) ! Only take sqrt after summing across procs

  ! Get the bulk RMS value
  call integrate_y_var(urms, urms_b)
  call integrate_y_var(vrms, vrms_b)
  call integrate_y_var(wrms, wrms_b)
  ! Write out the bulk RMS Velocity
  if (rank == 0) then
    write (*,  '("<U_rms> = " ES26.18)') sqrt(urms_b)
    write (*,  '("<V_rms> = " ES26.18)') sqrt(wrms_b)
    write (*,  '("<W_rms> = " ES26.18)') sqrt(vrms_b)
  end if

  urms = sqrt(urms)
  vrms = sqrt(vrms)
  wrms = sqrt(wrms)



  ! Communicate Z-Mean

  uu_xy = uu_xy / float(Nz) ! Can't take sqrt, then sum next...
  vv_xy = vv_xy / float(Nz)
  ww_xy = ww_xy / float(Nz)

  if (Nz > 1) then

    ! Need these in rankZ = 0 to compute Production!
    fname = 'mean_xz.h5'
    gname = 'uu_xz'
    call reduce_and_write_XYplane(fname, gname, uu_xy, .false., movie)
    gname = 'ww_xz'
    call reduce_and_write_XYplane(fname, gname, vv_xy, .false., movie)
    gname = 'vv_xz'
    call reduce_and_write_XYplane(fname, gname, ww_xy, .false., movie)

  end if


  gname = 'TKE_zstar'
  call Bin_Ystar_and_Write(gname, f1)

  f1 = 0.5d0 * (u1**2. + u2**2. + u3**2.)
  gname = 'KE_zstar'
  call Bin_Ystar_and_Write(gname, f1)



end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine compute_TKE_Production(movie)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Compute the TKE Production terms, and also Reynolds Stresses
  !  Both the horizontal averages -- ume(y), etc
  !   and the z-averages -- ume_xz(x,y), etc
  ! Assumes PP Ui
  ! Requires dudx_m to be stored already from compute_MKE_diss

  character(len=35) fname
  character(len=20) gname
  logical movie
  integer i, j, k


  !!! Reynolds Stress !!!
  ! uv and wv are defined on the GY grid -- Then can compute that production conservatively!
  ! uw is defined on the GYF grid
  uvar_xy = 0. ! uv
  wvar_xy = 0. ! wv
  vvar_xy = 0. ! uw

  do j = 1, Nyp ! On GYF
    uw(j) = 0.
    do k = 0, Nzp - 1
      do i = 0, Nxm1

        uw(j) = uw(j) + (u1(i, k, j) - ume_xy(i, j)) &
                      * (u3(i, k, j) - wme_xy(i, j))
        vvar_xy(i, j) = vvar_xy(i, j) + (u1(i, k, j) - ume_xy(i, j)) &
                                      * (u3(i, k, j) - wme_xy(i, j))

      end do
    end do
  end do

  do j = 2, Nyp ! On GY
    uv(j) = 0.
    wv(j) = 0.
    do k = 0, Nzp - 1
      do i = 0, Nxm1

        uv(j) = uv(j) +   (dyf(j - 1) * (u1(i, k, j) - ume_xy(i, j)) + &
                            dyf(j) * (u1(i, k, j - 1) - ume_xy(i, j - 1))) &
                                    / (2.d0 * dy(j)) &
                              * (u2(i, k, j) - vme_xy(i, j))
        uvar_xy(i, j) = uvar_xy(i, j) +   (dyf(j - 1) * (u1(i, k, j) - ume_xy(i, j)) + &
                            dyf(j) * (u1(i, k, j - 1) - ume_xy(i, j - 1))) &
                                    / (2.d0 * dy(j)) &
                              * (u2(i, k, j) - vme_xy(i, j))

        wv(j) = wv(j) +   (dyf(j - 1) * (u3(i, k, j) - wme_xy(i, j)) + &
                            dyf(j) * (u3(i, k, j - 1) - wme_xy(i, j - 1))) &
                                    / (2.d0 * dy(j)) &
                              * (u2(i, k, j) - vme_xy(i, j))
        wvar_xy(i, j) = wvar_xy(i, j) +   (dyf(j - 1) * (u3(i, k, j) - wme_xy(i, j)) + &
                            dyf(j) * (u3(i, k, j - 1) - wme_xy(i, j - 1))) &
                                    / (2.d0 * dy(j)) &
                              * (u2(i, k, j) - vme_xy(i, j))

      end do
    end do
  end do



  ! Communicate Horizontal Mean (Save them later)

  call mpi_allreduce(mpi_in_place, uv, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)
  call mpi_allreduce(mpi_in_place, wv, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)
  call mpi_allreduce(mpi_in_place, uw, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)


  uv = uv / float(Nx * Nz)
  wv = wv / float(Nx * Nz)
  uw = uw / float(Nx * Nz)



  ! Communicate Z-Mean

  uvar_xy = uvar_xy / float(Nz)
  wvar_xy = wvar_xy / float(Nz)
  vvar_xy = vvar_xy / float(Nz)

  if (Nz > 1) then
    fname = 'mean_xz.h5'
    gname = 'uw_xz'
    call reduce_and_write_XYplane(fname, gname, uvar_xy, .false., movie)
    gname = 'wv_xz'
    call reduce_and_write_XYplane(fname, gname, wvar_xy, .false., movie)
    gname = 'uv_xz'
    call reduce_and_write_XYplane(fname, gname, vvar_xy, .false., movie)


    !!! Production Terms (Just on rankZ = 0) !!!
    ! Just multiply the Reynolds Stress with the mean_xy field...
    !   Only compute as f(y), since can get z-average by multiplication...
    ! Use the mean field gradients stored from compute_MKE_diss
    if (rankZ == 0) then

      do j = 1, Nyp
        uu_dudx(j) = 0.
        wu_dwdx(j) = 0.
        do i = 0, Nxm1

          ! On GYF
          uu_dudx(j) = uu_dudx(j) +  uu_xy(i, j) * dudx_m(i, j)
          wu_dwdx(j) = wu_dwdx(j) +  vvar_xy(i, j)     * dwdx_m(i, j)

        end do
      end do

      do j = 2, Nyp
        vu_dvdx(j) = 0.
        vv_dvdy(j) = 0.
        uv_dudy(j) = 0.
        wv_dwdy(j) = 0.
        do i = 0, Nxm1

          ! On GY
          vu_dvdx(j) = vu_dvdx(j) +  uvar_xy(i, j)     * dvdx_m(i, j)

          vv_dvdy(j) = vv_dvdy(j) +  (vme_xy(i, j + 1) - vme_xy(i, j - 1)) &
                                           / (gy(j + 1) - gy(j - 1)) &
                                    * vv_xy(i, j)
          uv_dudy(j) = uv_dudy(j) +  (ume_xy(i, j) - ume_xy(i, j - 1)) &
                                           / dy(j) &
                                    * uvar_xy(i, j)
          wv_dwdy(j) = wv_dwdy(j) +  ((wme_xy(i, j) - wme_xy(i, j - 1)) &
                                           / dy(j)  + &
                                              (1.d0 / (Ro_inv / delta)) * dTHdX(1) )  &  ! Include the TWS!
                                    * wvar_xy(i, j)

        end do
      end do

      uu_dudx = uu_dudx / float(Nx)
      wu_dwdx = wu_dwdx / float(Nx)
      vu_dvdx = vu_dvdx / float(Nx)
      vv_dvdy = vv_dvdy / float(Nx)
      uv_dudy = uv_dudy / float(Nx)
      wv_dwdy = wv_dwdy / float(Nx)

      ! Don't communicate, since they're already all on rankZ = 0!

    end if

  end if


end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine compute_Vorticity(movie)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Compute the RMS Vorticity and also Movie Slices of Vorticity
  ! Assumes PP Ui AND also FF CUi are stored in cr1(), etc

  character(len=35) fname
  character(len=20) gname
  logical movie
  integer i, j, k
  real(rkind) varxy(0:Nxm1, 1:Nyp), varzy(0:Nzp - 1, 1:Nyp), varxz(0:Nxm1, 0:Nzp - 1)



  !!! RMS Vorticity !!!
  ! X-component in FF space
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 !Nkx
        cs1(i, k, j) = cikz(k) * 0.5d0 * (cr2(i, k, j + 1) + cr2(i, k, j)) &
                       - (cr3(i, k, j + 1) - cr3(i, k, j - 1)) / (2.d0 * dyf(j))
      end do
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  ! RMS value
  do j = 1, Nyp
    omega_x(j) = 0.d0
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = s1(i, k, j) - (1.d0 / (Ro_inv / delta)) * dTHdX(1)  ! Include the TWS!
        omega_x(j) = omega_x(j) + s1(i, k, j)**2.d0
      end do
    end do
  end do
  call mpi_allreduce(mpi_in_place, omega_x, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  omega_x = sqrt(omega_x / float(Nx * Nz))



  ! Write Movie Slices for omega_x
  if (movie) then

    fname = 'movie.h5'
    call mpi_barrier(mpi_comm_world, ierror)
    if (rankZ == rankzmovie) then
      do j = 1, Nyp
        do i = 0, Nxm1
          varxy(i, j) = s1(i, NzMovie, j)
        end do
      end do
      write (gname,'("omegaX_xz")')
      call WriteHDF5_XYplane(fname, gname, varxy)
    end if

    if (rankY == rankymovie) then
      do j = 0, Nzp - 1
        do i = 0, Nxm1
          varxz(i, j) = s1(i, j, NyMovie)
        end do
      end do
      write (gname,'("omegaX_xy")')
      call WriteHDF5_XZplane(fname, gname, varxz)
    end if

    do j = 1, Nyp
      do i = 0, Nzp - 1
        varzy(i, j) = s1(NxMovie, i, j)
      end do
    end do
    write (gname,'("omegaX_yz")')
    call WriteHDF5_ZYplane(fname, gname, varzy)
  end if



  ! Y-component in FF space
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 !Nkx
        cs1(i, k, j) = cikx(i) * cr3(i, k, j) - cikz(k) * cr1(i, k, j)
      end do
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  ! RMS value
  do j = 1, Nyp
    omega_y(j) = 0.d0
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        omega_y(j) = omega_y(j) + s1(i, k, j)**2.d0
      end do
    end do
  end do
  call mpi_allreduce(mpi_in_place, omega_y, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  omega_y = sqrt(omega_y / float(Nx * Nz))



  ! Write Movie Slices for omega_y
  if (movie) then

    fname = 'movie.h5'
    call mpi_barrier(mpi_comm_world, ierror)
    if (rankZ == rankzmovie) then
      do j = 1, Nyp
        do i = 0, Nxm1
          varxy(i, j) = s1(i, NzMovie, j)
        end do
      end do
      write (gname,'("omegaZ_xz")')
      call WriteHDF5_XYplane(fname, gname, varxy)
    end if

    if (rankY == rankymovie) then
      do j = 0, Nzp - 1
        do i = 0, Nxm1
          varxz(i, j) = s1(i, j, NyMovie)
        end do
      end do
      write (gname,'("omegaZ_xy")')
      call WriteHDF5_XZplane(fname, gname, varxz)
    end if

    do j = 1, Nyp
      do i = 0, Nzp - 1
        varzy(i, j) = s1(NxMovie, i, j)
      end do
    end do
    write (gname,'("omegaZ_yz")')
    call WriteHDF5_ZYplane(fname, gname, varzy)
  end if




  ! Z-component in FF space
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        cs1(i, k, j) = (cr1(i, k, j + 1) - cr1(i, k, j - 1)) / (2.d0 * dyf(j)) &
                       - cikx(i) * 0.5d0 * (cr2(i, k, j + 1) + cr2(i, k, j))
      end do
    end do
  end do
  call fft_xz_to_physical(cs1, s1)
  ! RMS value
  do j = 1, Nyp
    omega_z(j) = 0.d0
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        omega_z(j) = omega_z(j) + s1(i, k, j)**2.d0
      end do
    end do
  end do
  call mpi_allreduce(mpi_in_place, omega_z, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  omega_z = sqrt(omega_z / float(Nx * Nz))



  ! Write Movie Slices for omega_z
  if (movie) then

    fname = 'movie.h5'
    call mpi_barrier(mpi_comm_world, ierror)
    if (rankZ == rankzmovie) then
      do j = 1, Nyp
        do i = 0, Nxm1
          varxy(i, j) = s1(i, NzMovie, j)
        end do
      end do
      write (gname,'("omegaY_xz")')
      call WriteHDF5_XYplane(fname, gname, varxy)
    end if

    if (rankY == rankymovie) then
      do j = 0, Nzp - 1
        do i = 0, Nxm1
          varxz(i, j) = s1(i, j, NyMovie)
        end do
      end do
      write (gname,'("omegaY_xy")')
      call WriteHDF5_XZplane(fname, gname, varxz)
    end if

    do j = 1, Nyp
      do i = 0, Nzp - 1
        varzy(i, j) = s1(NxMovie, i, j)
      end do
    end do
    write (gname,'("omegaY_yz")')
    call WriteHDF5_ZYplane(fname, gname, varzy)
  end if




end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine compute_BPE
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Compute the Background Potential Energy (BPE) -- Tseng & Ferziger 2001
  ! th(:,:,:,1) already in Physical Space
  ! Uses s1 as storage
  ! AFW 2020

  character(len=20) gname
  character(len=35) fname
  integer i, j, k, bin
  integer, parameter :: Nbin = Nx ! Useful for writing in parallel (vs 16*Ny)
  real(rkind) thmin, thmax, dTH, BPE
  real(rkind) PDF(0:Nbin - 1)
  real(rkind) Y_r(0:Nbin)
  real(rkind) th_bin(0:Nbin)
  real(rkind) DiagX(0:int(Nbin/NprocZ) - 1)

  ! Add background buoyancy gradient, store in s1
  !   Offset by Lx/2 to match the calculation of thv_m (Mean Buoyancy Production)
  s1 = 0. ! Need to clear the edges otherwise min/max value catches it!
  do j = 1, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = th(i, k, j, 1) + dTHdX(1) * (gx(i) - 0.5*Lx)
      end do
    end do
  end do


  if (homogeneousX) then ! Infinite Front, or the likes -- Use constant th bins

    ! Bounds of theta
    thmin = -dTHdX(1) * 0.5*Lx * 1.1
    thmax =  dTHdX(1) * 0.5*Lx * 1.1

    dTH = (thmax - thmin) / Nbin + 1.d-14

    do i = 0, Nbin
      th_bin(i) = thmin + i * dTH
    end do

    ! Compile the PDF in b
    PDF = 0.d0
    do j = jstart_th(1), jend_th(1) ! Avoid repeats at the vertical proc boundaries
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! Compute bindex
          bin = min(int((s1(i, k, j) - thmin) / dTH), Nbin - 1) ! Catch case when bin = Nbin...

          PDF(bin) = PDF(bin) + dyf(j)
        end do
      end do
    end do

  else ! Finite Front -- Use adaptive th bins with dTH ~ sech2(z/delta)

    ! Bounds of theta
    thmin = -Ro_inv
    thmax = +Ro_inv

    do i = 0, Nbin
      th_bin(i) = -(Ro_inv + 1.d-14) * cos(i * pi / Nbin)  ! (Ro_inv + 1.d-14) * tanh( (dble(i)/Nbin - 0.5d0) * Lx/(2.d0*delta) )
    end do

    ! Compile the PDF in b
    PDF = 0.d0
    do j = jstart_th(1), jend_th(1) ! Avoid repeats at the vertical proc boundaries
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! Do inverse transform to get bindex
          bin = min(int( acos( -min(max(s1(i, k, j),thmin),thmax) / (Ro_inv + 1.d-14) ) * Nbin / pi ), Nbin - 1) ! Catch case when bin = Nbin...

          PDF(bin) = PDF(bin) + dyf(j)
        end do
      end do
    end do


  endif



  call mpi_allreduce(mpi_in_place, PDF, Nbin, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, ierror)

  ! Enforce \int_B PDF dB = 1 exactly (small dyf/2 at BCs...)  vs. /(Ly * Nx * Nz)
  PDF = PDF / sum(PDF)

  ! Compute Y_r (corresponding to TH bin edges)
  Y_r(0) = 0.d0
  do i = 1, Nbin
    ! NOTE: PDF is the distribution function, NOT density (i.e. NOT divided by dTH, so don't need another dTH!)
    Y_r(i) = Y_r(i - 1) + PDF(i - 1) ! * (th_bin(i) - th_bin(i - 1))
  end do
  Y_r = Y_r * Ly

  ! Compute BPE
  BPE = 0.d0
  do i = 0, Nbin - 1  ! Integrate
    BPE = BPE - (0.5 * (th_bin(i + 1) + th_bin(i)) * 0.5 * (Y_r(i + 1) + Y_r(i))) * &
                    (Y_r(i + 1) - Y_r(i)) / Ly
  end do


  fname = 'mean.h5'
  gname = 'BPE'
  call WriteHDF5_real(fname, gname, BPE)

  ! Write out the entire PDF and extents to construct bins
  if (rankY == 0) then
    DiagX = PDF(rankZ * int(Nbin/NprocZ):(rankZ+1) * int(Nbin/NprocZ) - 1)
    gname = 'th1PDF'
    call WriteStatH5_X(fname, gname, DiagX, int(Nbin/NprocZ))

    DiagX = th_bin(rankZ * int(Nbin/NprocZ):(rankZ+1) * int(Nbin/NprocZ) - 1)
    gname = 'th1bin'
    call WriteStatH5_X(fname, gname, DiagX, int(Nbin/NprocZ))
  end if

end





!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine Bin_Ystar_and_Write(gname, field)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Bin field into z* Coordinates (Winters et al 1996)
  !   using Parallel PDF Sorting of Tseng & Ferziger 2001
  ! th(:,:,:,1) already in Physical Space
  ! field is to be binned
  ! Uses s1 as storage
  ! AFW 2021

  character(len=20) gname
  real(rkind), intent(in) :: field(:,:,:)

  character(len=35) fname
  integer i, j, k, bin
  integer, parameter :: Nbin = Nx ! Useful for writing in parallel (vs 16*Ny)
  real(rkind) thmin, thmax, dTH
  real(rkind) PDF(0:Nbin - 1)
  real(rkind) field_binned(0:Nbin - 1) ! For computing mean (field) on each z* surface
  real(rkind) DiagX(0:int(Nbin/NprocZ) - 1)

  ! Add background buoyancy gradient, store in s1
  !   Offset by Lx/2 to match the calculation of thv_m (Mean Buoyancy Production)
  s1 = 0. ! Need to clear the edges otherwise min/max value catches it!
  do j = 1, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = th(i, k, j, 1) + dTHdX(1) * (gx(i) - 0.5*Lx)
      end do
    end do
  end do


  if (homogeneousX) then ! Infinite Front, or the likes -- Use constant th bins

    ! Bounds of theta
    thmin = -dTHdX(1) * 0.5*Lx * 1.1
    thmax =  dTHdX(1) * 0.5*Lx * 1.1

    dTH = (thmax - thmin) / Nbin + 1.d-14

    ! Compile the PDF in b
    PDF = 0.d0
    field_binned = 0.d0
    do j = jstart_th(1), jend_th(1) ! Avoid repeats at the vertical proc boundaries
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! Compute bindex
          bin = min(int((s1(i, k, j) - thmin) / dTH), Nbin - 1) ! Catch case when bin = Nbin...

          PDF(bin) = PDF(bin) + dyf(j)
          field_binned(bin) = field_binned(bin) + field(i, k, j) * dyf(j)
        end do
      end do
    end do

  else ! Finite Front -- Use adaptive th bins with dTH ~ sech2(z/delta)

    ! Bounds of theta
    thmin = -Ro_inv
    thmax = +Ro_inv

    ! Compile the PDF in b
    PDF = 0.d0
    field_binned = 0.d0
    do j = jstart_th(1), jend_th(1) ! Avoid repeats at the vertical proc boundaries
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          ! Do inverse transform to get bindex
          bin = min(int( acos( -min(max(s1(i, k, j),thmin),thmax) / (Ro_inv + 1.d-14) ) * Nbin / pi ), Nbin - 1) ! Catch case when bin = Nbin...

          PDF(bin) = PDF(bin) + dyf(j)
          field_binned(bin) = field_binned(bin) + field(i, k, j) * dyf(j)
        end do
      end do
    end do

  endif


  call mpi_allreduce(mpi_in_place, PDF, Nbin, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, ierror)
  call mpi_allreduce(mpi_in_place, field_binned, Nbin, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, ierror)

  ! Divide by sum(dyf), i.e. PDF, to get a bin-average
  field_binned = field_binned / PDF

  fname = 'mean.h5'
  ! Write out the entire PDF and extents to construct bins
  if (rankY == 0) then
    DiagX = field_binned(rankZ * int(Nbin/NprocZ):(rankZ+1) * int(Nbin/NprocZ) - 1)
    call WriteStatH5_X(fname, gname, DiagX, int(Nbin/NprocZ))
  end if

end









!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine integrate_y_var(var, res)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Integrates a vector in Y

  integer j
  real(rkind) var(0:Nyp + 1), res

  res = 0.
  do j = 2, Nyp
    res = res + 0.5 * (var(j) + var(j - 1)) * dy(j)
  end do
  call mpi_allreduce(mpi_in_place, res, 1, &
                     mpi_double_precision, mpi_sum, mpi_comm_y, ierror)

  res = res / Ly ! To get average

end subroutine


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine integrate_z_var(var, res)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Integrates a full 3D cube across Z

  integer i, k, j
  real(rkind) var(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1), res(0:Nxm1, 1:Nyp)

  do i = 0, Nxm1
    do j = 1, Nyp
      res(i, j) = 0.
      do k = 0, Nzp
        res(i, j) = res(i, j) + var(i, k, j) * dz(k)
      end do
    end do
  end do
  call mpi_allreduce(mpi_in_place, res, Nx * Nyp, &
                     mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
  call mpi_barrier(mpi_comm_world, ierror)

  res = res / Lz ! Gives average

end subroutine



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine reduce_and_write_XYplane(fname, gname, res, allreduce, movie)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Communicates X-Y plane and sums across all Z processors.
  ! Sends the result to all processors if allreduce == .true.

  character(len=35) fname
  character(len=20) gname
  logical movie
  real(rkind) res(0:Nxm1, 1:Nyp)
  logical allreduce, writeHDF5

  if (allreduce) then
    call mpi_allreduce(mpi_in_place, res, Nx * Nyp, &
                       mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
  else
    ! if (.not. movie) return ! Why are we doing this otherwise?

    if (rankZ == 0) then
      call mpi_reduce(mpi_in_place, res, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    else ! Don't use mpi_in_place for other processes, except for allreduce...
      call mpi_reduce(res, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    end if
  end if

  call mpi_barrier(mpi_comm_world, ierror)

  if (movie .and. rankZ == 0) then
    call WriteHDF5_XYplane(fname, gname, res)
  end if

end subroutine











! ---- LES ----

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine compute_TKE_diss_les
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Calculate the componet of the SGS dissipation rate
  ! only includes the terms timestepped implicitly

  character(len=35) fname
  character(len=20) gname
  real(rkind) eps_sgs2(1:Nyp)
  real(rkind) Diag(Nyp)
  integer i, j, k

  ! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
  ! At *consistent* GY points!

  ! Compute the contribution at GYF first. Store in S1.
  do j = 1, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = u1(i, k, j) * &
                        ((nu_t(i, k, j + 1) * (u1(i, k, j + 1) - u1(i, k, j)) / dy(j + 1) &
                               - nu_t(i, k, j) * (u1(i, k, j) - u1(i, k, j - 1)) / dy(j)) &
                           / dyf(j)) &  ! At GYF
                     + u3(i, k, j) * &
                        ((nu_t(i, k, j + 1) * (u3(i, k, j + 1) - u3(i, k, j)) / dy(j + 1) &
                               - nu_t(i, k, j) * (u3(i, k, j) - u3(i, k, j - 1)) / dy(j)) &
                           / dyf(j))  ! At GYF
      end do
    end do
  end do

! Then, interpolate the u1 & u2 contribution onto GY
! so that it conserves the dissipation as in code
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        temp(i, k, j) = 0.5d0 * (s1(i, k, j) + s1(i, k, j - 1)) & ! u1 & u3 at GY
                      + u2(i, k, j) * &
                        ((0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j + 1)) * (u2(i, k, j + 1) - u2(i, k, j)) &
                            / dyf(j) &
                        - 0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j - 1)) * (u2(i, k, j) - u2(i, k, j - 1)) &
                            / dyf(j - 1)) / dy(j)) ! At GY
      end do
    end do
  end do

  ! Now calculate the horizontal average
  do j = 1, Nyp
    eps_sgs2(j) = 0.d0
    do i = 0, Nxm1
      do k = 0, Nzp - 1
        eps_sgs2(j) = eps_sgs2(j) + temp(i, k, j)
      end do
    end do
  end do

  call mpi_allreduce(mpi_in_place, eps_sgs2, Nyp &
                     , mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  fname = 'mean.h5'

  if (rankZ == 0) then
    gname = 'eps_sgs2'
    Diag = eps_sgs2(1:Nyp) / float(Nx * Nz)
    call WriteStatH5_Y(fname, gname, Diag)
  end if

end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine save_stats_LES_OOL(blank)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

  integer n
  character(len=35) fname
  character(len=20) gname
  logical blank
  real(rkind) :: Diag(1:Nyp)


  if (blank) then
    fname = 'mean.h5'

    if (rankZ == 0) then
      Diag = 0.d0
      gname = 'nu_sgs'
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'eps_sgs1'
      call WriteStatH5_Y(fname, gname, Diag)

      Diag = 0.d0
      gname = 'kappa_sgs'
      call WriteStatH5_Y(fname, gname, Diag)

    end if
  else
    ! Needed to write out LES Statistics without timestepping...
    ! DON'T run this except for when stopping the simulation!

    rk_step = 1
    flag_save_LES = .true.

    call les_chan
    call les_chan_th

    call fft_xz_to_fourier(u1, cu1)
    call fft_xz_to_fourier(u2, cu2)
    call fft_xz_to_fourier(u3, cu3)

    do n = 1, N_th
      call fft_xz_to_fourier(th(:, :, :, n), cth(:, :, :, n))
    end do

  end if

  return

end
