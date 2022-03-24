module fft
  !******************************************************************************|
  ! fft.f90, the FFT package for diablo.                               VERSION 3.0
  !  Updated for FFTW3, 4/10/21 AFW
  !
  ! This file isolates all calls to the FFTW package (available at: www.fftw.org)
  ! These wrapper routines were written by T. Bewley (spring 2001).
  !******************************************************************************|
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! The arrangement of the significant real numbers in the arrays (denoted by +)
  ! in physical space, in Fourier space, and in Fourier space after packing are
  ! shown below for the 2D (X-Z) plane.  The third direction (Y) is handled in
  ! an identical matter as the Z direction shown here.
  !
  !      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
  !      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
  ! Nz-1 ++++++++++++++++oo     -1  ++++++++++++oooooo         oooooooooooooooooo
  !      ++++++++++++++++oo     -2  ++++++++++++oooooo         oooooooooooooooooo
  !      ++++++++++++++++oo     -3  ++++++++++++oooooo         oooooooooooooooooo
  !      ++++++++++++++++oo         ++++++++++++oooooo         oooooooooooooooooo
  !      ++++++++++++++++oo    -Nkz ++++++++++++oooooo         oooooooooooooooooo
  !      ++++++++++++++++oo         oooooooooooooooooo     -1  ++++++++++++oooooo
  !      ++++++++++++++++oo         oooooooooooooooooo     -2  ++++++++++++oooooo
  !      ++++++++++++++++oo         oooooooooooooooooo     -3  ++++++++++++oooooo
  !      ++++++++++++++++oo         oooooooooooooooooo         ++++++++++++oooooo
  !      ++++++++++++++++oo         oooooooooooooooooo    -Nkz ++++++++++++oooooo
  !      ++++++++++++++++oo     Nkz ++++++++++++oooooo     Nkz ++++++++++++oooooo
  !      ++++++++++++++++oo         ++++++++++++oooooo         ++++++++++++oooooo
  !   3  ++++++++++++++++oo      3  ++++++++++++oooooo      3  ++++++++++++oooooo
  !   2  ++++++++++++++++oo      2  ++++++++++++oooooo      2  ++++++++++++oooooo
  !   1  ++++++++++++++++oo      1  ++++++++++++oooooo      1  ++++++++++++oooooo
  !   0  ++++++++++++++++oo      0  +o++++++++++oooooo      0  +o++++++++++oooooo
  !      ^^^^           ^           ^ ^ ^     ^                ^ ^ ^     ^
  !      0123           Nx-1        0 1 2     Nkx              0 1 2     Nkx
  !
  !       PHYSICAL SPACE              FOURIER SPACE         FOURIER SPACE (PACKED)
  !
  ! After the Real->Fourier transform, the significant coefficients are put next
  ! to each other in the array, so a loop such as
  !
  !        DO K=0,twoNkz           [where twoNkz = 2*Nkz = 2*(Nz/3) ]
  !          DO I=0,Nkx          [where  Nkx = Nx/3             ]
  !            CP(I,K,J)= ...
  !          END DO
  !        END DO
  !
  ! includes all the Fourier coefficients of interest.  The subsequent loops in
  ! Fourier space just work on these coefficients in the matrix.
  !
  ! Before a Fourier->Real transform, the significant coefficients are unpacked
  ! and the higher wavenumbers are SET TO ZERO before the inverse transform.
  ! This has the effect of doing the required dealiasing.
  !
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  use domain
  use iso_c_binding
  implicit none
  save
  include 'fftw3.f03'


  integer*8   fftw_x_to_p_plan, fftw_x_to_f_plan, &
              fftw_z_to_p_plan, fftw_z_to_f_plan

  complex(C_DOUBLE_COMPLEX), pointer, contiguous :: temp_fft(:,:,:)

contains


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine alloc_array3D(rdat,cdat)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|


    real(C_DOUBLE), pointer, contiguous, intent(out) :: rdat(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer, contiguous, intent(out) :: cdat(:,:,:)

    real(C_DOUBLE), pointer, contiguous :: r1_temp(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer, contiguous :: c1_temp(:,:,:)

    type(C_PTR) :: data

    integer num_el

    num_el = max(ceiling(0.5*(Nx+2)*(Nzp+2)*(Nyp+2)), &
                  (Nxp+1)*(Nz+2)*(Nyp+2))

    data = fftw_alloc_complex(int(num_el, C_SIZE_T))
    call c_f_pointer(data, r1_temp,  [(Nx+2),(Nzp+2),(Nyp+2)])
    call c_f_pointer(data, c1_temp, [(Nxp+1),(Nz+2),(Nyp+2)])


    rdat(0:Nx+1,0:Nzp+1,0:Nyp+1) => r1_temp
    cdat(0:Nxp,0:Nz+1,0:Nyp+1) => c1_temp

  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine alloc_array4D(rdat,cdat)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(C_DOUBLE), pointer, contiguous, intent(out) :: rdat(:,:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer, contiguous, intent(out) :: cdat(:,:,:,:)

    real(C_DOUBLE), pointer, contiguous :: r1_temp(:,:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer, contiguous :: c1_temp(:,:,:,:)

    type(C_PTR) :: data

    integer num_el

    num_el = N_th*max(ceiling(0.5*(Nx+2)*(Nzp+2)*(Nyp+2)), &
                  (Nxp+1)*(Nz+2)*(Nyp+2))

    data = fftw_alloc_complex(int(num_el, C_SIZE_T))
    call c_f_pointer(data, r1_temp,  [(Nx+2),(Nzp+2),(Nyp+2),N_th])
    call c_f_pointer(data, c1_temp, [(Nxp+1),(Nz+2),(Nyp+2),N_th])


    rdat(0:Nx+1,0:Nzp+1,0:Nyp+1,1:N_th) => r1_temp
    cdat(0:Nxp,0:Nz+1,0:Nyp+1,1:N_th) => c1_temp

  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine init_fft(u1,cu1)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    real(C_DOUBLE) :: u1(0:Nx+1,0:Nzp+1,0:Nyp+1)
    complex(C_DOUBLE_COMPLEX) :: cu1(0:Nxp,0:Nz+1,0:Nyp+1)
    integer i, j, k

    type(C_PTR) :: p_temp_fft
    complex(C_DOUBLE_COMPLEX), pointer, contiguous :: ttemp_fft(:,:,:)
    integer(C_SIZE_T) :: temp_fft_size
    integer(kind=mpi_address_kind) :: temp_fft_win_size
    integer :: disp_mem
    type(mpi_info) :: info


    integer howmany_n(2), howmany_stride(2), howmany_stride_c(2)

    if (verbosity > 2 .and. rank == 0) &
      write (*, '("Initializing FFTW3 package....")')

    pi = 4.*atan(1.0)
    ci = cmplx(0.0, 1.0)
    eps = 0.000000001



    ! Allocate the contiguous C2R FFT Temporary Working Array
    temp_fft_size = int((Nx/2+1) * (Nzp+2) * (Nyp+2), C_SIZE_T)
    p_temp_fft = fftw_alloc_complex(temp_fft_size)
    call c_f_pointer(p_temp_fft, ttemp_fft, [int(Nx/2+1), (Nzp+2), (Nyp+2)])
    temp_fft(0:Nx / 2, 0:Nzp + 1, 0:Nyp + 1) => ttemp_fft


#ifdef SHARED_MEMORY
    ! Make temp_fft a shared-memory access window across all of mpi_comm_z

    call mpi_info_create(info, ierror)
    call mpi_info_set(info, "same_size", "true", ierror)
    call mpi_info_set(info, "same_disp_unit", "true", ierror)
    temp_fft_win_size = int((Nx/2+1) * (Nzp+2) * (Nyp+2) * 16, mpi_address_kind)
    disp_mem = 16
    call mpi_win_create(temp_fft(0,0,0), temp_fft_win_size, &
                        disp_mem, info, mpi_comm_z, temp_fft_win, ierror)

#endif

    ! Plan the XY PP <--> FP
    howmany_n(1)      = Nzp
    howmany_stride(1) = Nx + 2
    howmany_n(2)      = Nyp + 2
    howmany_stride(2) = (Nzp + 2)*(Nx + 2)
    howmany_stride_c(1) = int(Nx/2 + 1)
    howmany_stride_c(2) = int(Nx/2 + 1)*(Nzp + 2)
    call dfftw_plan_guru_dft_r2c(fftw_x_to_f_plan, 1, Nx, 1, 1, &
                                  2, howmany_n, howmany_stride, howmany_stride_c, &
                                  u1(0,0,0), temp_fft(0,0,0), &
                                  fftw_patient + fftw_destroy_input)

    call dfftw_plan_guru_dft_c2r(fftw_x_to_p_plan, 1, Nx, 1, 1, &
                                  2, howmany_n, howmany_stride_c, howmany_stride, &
                                  temp_fft(0,0,0), u1(0,0,0), &
                                  fftw_patient + fftw_destroy_input)


    ! Plan the XY FP <--> FF
    howmany_n(1)      = Nxp
    howmany_stride(1) = 1
    howmany_n(2)      = Nyp + 2
    howmany_stride(2) = (Nxp + 1)*(Nz + 2)
    call dfftw_plan_guru_dft(fftw_z_to_f_plan, 1,  Nz, Nxp + 1, Nxp + 1,    &
                              2, howmany_n, howmany_stride, howmany_stride, &
                              cu1(0,0,0), cu1(0,0,0), &
                              fftw_forward, fftw_patient)

    call dfftw_plan_guru_dft(fftw_z_to_p_plan, 1,  Nz, Nxp + 1, Nxp + 1,    &
                              2, howmany_n, howmany_stride, howmany_stride, &
                              cu1(0,0,0), cu1(0,0,0), &
                              fftw_backward, fftw_patient)



    rNx = 1.0 * Nx
    do i = 0, Nxp - 1
      kx(i) = (i + Nxp * rankZ) * (2.*pi) / Lx
      kx2(i) = kx(i) * kx(i)
      cikx(i) = ci * kx(i)
    end do

    rNz = 1.0 * Nz
    do k = 0, Nkz
      kz(k) = k * (2.*pi) / Lz
    end do
    do k = 1, Nkz
      kz(twoNkz + 1 - k) = -k * (2.*pi) / Lz
    end do
    do k = 0, twoNkz
      kz2(k) = kz(k) * kz(k)
      cikz(k) = ci * kz(k)
    end do

    if (rank == 0) then
      write (*, '("Nkx   = " I10)') Nkx
      write (*, '("2Nkz  = " I10)') twoNkz
    end if

    if (verbosity > 2 .and. rank == 0) &
      write (*, '("FFTW3 package initialised.")')

    return
  end




  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine fft_xz_to_fourier(v, vv)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Optimised now always from 0: Nyp+1

    real(rkind) :: v(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    complex(rkind) :: vv(0:Nxp, 0:Nz + 1, 0:Nyp + 1)

    integer i, j, k, p

    ! FFT in X
    ! (Now computes all slices in one go, for efficiency)
    call dfftw_execute_dft_r2c(fftw_x_to_f_plan, v(0, 0, 0), temp_fft(0, 0, 0))
    do j = 0, Nyp + 1
      do k = 0, Nzp - 1
        do i = 0, Nkx
          temp_fft(i, k, j) = temp_fft(i, k, j) / Nx
        end do
        do i = Nkx + 1, Nx / 2
          temp_fft(i, k, j) = cmplx(0.d0, 0.d0)
        end do
      end do
    end do


#ifdef SHARED_MEMORY
    ! 1-point Communication within the node

    call mpi_barrier(mpi_comm_z, ierror)
    call mpi_win_lock_all(mpi_mode_nocheck, temp_fft_win, ierror)
    do p = 0, NprocZ - 1
      call mpi_get(vv(0, Nzp*p, 0), 1, type_cFF_full, p, &
                  int(rankZ * Nxp, mpi_address_kind), & ! Displacement in units of disp_mem! (NOT bytes!)
                  1, type_cFP_full, temp_fft_win, ierror)
    end do
    !call mpi_barrier(mpi_comm_z, ierror)
    call mpi_win_flush_local_all(temp_fft_win, ierror)
    call mpi_win_unlock_all(temp_fft_win, ierror)
    call mpi_barrier(mpi_comm_z, ierror)

#else

    call mpi_alltoall(temp_fft(0, 0, 0), 1, xy2zy_1, &
                       vv(0, 0, 0), 1, xy2zy_2, mpi_comm_z, ierror)

#endif

    ! FFT in Z
    ! (Now computes all slices in one go, for efficiency. #Guru.)
    call dfftw_execute_dft(fftw_z_to_f_plan, vv(0, 0, 0), vv(0, 0, 0))


    do j = 0, Nyp + 1
      do k = 0, Nkz
        do i = 0, Nxp - 1
          vv(i, k, j) = vv(i, k, j) / Nz
        end do
      end do
      ! PACK
      do k = 1, Nkz
        do i = 0, Nxp - 1
          vv(i, Nkz + k, j) = vv(i, Nz - 1 + k - Nkz, j) / Nz
        end do
      end do
    end do

  end subroutine


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine fft_xz_to_physical(vv, v)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! Optimised now always from 0: Nyp+1

    real(rkind) :: v(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    complex(rkind) :: vv(0:Nxp, 0:Nz + 1, 0:Nyp + 1)

    integer i, j, k, p

    ! FFT in Z
    ! UNPACK
    do j = 0, Nyp + 1
      do k = Nkz, 1, -1
        do i = 0, Nxp - 1
          vv(i, Nz - 1 + k - Nkz, j) = vv(i, Nkz + k, j)
        end do
      end do
      do k = Nkz + 1, Nz - Nkz - 1
        do i = 0, Nxp - 1
          vv(i, k, j) = cmplx(0.d0, 0.d0) ! De-alias in Z
        end do
      end do
    end do


    ! (Now computes all slices in one go, for efficiency. #Guru.)
    call dfftw_execute_dft(fftw_z_to_p_plan, vv(0, 0, 0), vv(0, 0, 0))



#ifdef SHARED_MEMORY
    ! 1-point Communication within the node

    !call mpi_barrier(mpi_comm_z, ierror)
    call mpi_win_lock_all(mpi_mode_nocheck, temp_fft_win, ierror)
    do p = 0, NprocZ - 1
      call mpi_put(vv(0, Nzp*p, 0), 1, type_cFF_full, p, &
                  int(rankZ * Nxp, mpi_address_kind), & ! Displacement in units of disp_mem! (NOT bytes!)
                  1, type_cFP_full, temp_fft_win, ierror)
    end do
    call mpi_barrier(mpi_comm_z, ierror)
    call mpi_win_flush_all(temp_fft_win, ierror)
    call mpi_win_unlock_all(temp_fft_win, ierror)
    call mpi_barrier(mpi_comm_z, ierror)

#else

    call mpi_alltoall(vv(0, 0, 0),  1, xy2zy_2, &
                      temp_fft(0, 0, 0), 1, xy2zy_1, mpi_comm_z, ierror)

#endif

    ! FFT in X
    do j = 0, Nyp + 1
      do k = 0, Nzp - 1
        do i = Nkx + 1, Nx / 2
          temp_fft(i, k, j) = cmplx(0.d0, 0.d0) ! De-alias in X
        end do
      end do
    end do
    ! (Now computes all slices in one go, for efficiency)
    call dfftw_execute_dft_c2r(fftw_x_to_p_plan, temp_fft(0, 0, 0), v(0, 0, 0))

  end subroutine


end module fft
