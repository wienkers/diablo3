module domain
  use hdf5
  use mpi_f08
  implicit none
  save

  integer  verbosity

  ! Specify data-types
  integer, parameter :: single_kind = kind(0.0)
  integer, parameter :: double_kind = kind(0.d0)
  integer, parameter :: rkind = double_kind

  ! Details of the Computational Domain
  ! (We hardwire these into the code so that the compiler may perform
  !  optimizations based on the grid size at compile time).
  integer num_per_dir
  integer :: Nx, Ny, Nz, N_th
  integer :: NprocY, NprocZ, Nprocs,  NprocShared
  include 'grid_def'
  include 'grid_mpi'
  integer, parameter :: Nyp = (Ny-1)/NprocY + 1
  integer, parameter :: Nxp = Nx/(2*NprocZ) ! Nkxp... Since only need half in Fourier Space for Reals
  integer, parameter :: Nzp = Nz/(  NprocZ)

  ! Derived Constants
  integer, parameter :: Nkx = Nx/3
  integer, parameter :: Nkz = Nz/3
  integer, parameter :: Nxm1 = Nx - 1
  integer, parameter :: Nypm1 = Nyp - 1
  integer, parameter :: Nzm1 = Nz - 1
  integer, parameter :: twoNkz = 2*Nkz ! Skips the de-aliased wavenumbers

  ! Grid
  real(rkind)   gx(0:Nx+1),   gy(0:Nyp+1),   gz(0:Nz+1), &
                dx(0:Nx+1),   dy(0:Nyp+1),   dz(0:Nz+1), &
                gxf(0:Nx+1),  gyf(0:Nyp+1),  gzf(0:Nz+1), &
                dxf(0:Nx+1),  dyf(0:Nyp+1),  dzf(0:Nz+1)
  integer       jstart, jend
  integer       jstart_th(1:N_th), jend_th(1:N_th)


  ! FFT Grid
  real(rkind)   kx(0:Nxp), kz(0:2*Nkz), &
                kx2(0:Nxp), kz2(0:2*Nkz), &
                pi, eps, rNx, rNyp, rNz
  complex(rkind)  cikx(0:Nxp), cikz(0:2*Nkz), ci



  ! MPI Details
  logical :: use_mpi
  integer :: rank, rankY, rankZ, rankShared
  type(mpi_comm) :: mpi_comm_y, mpi_comm_z, mpi_comm_shared
  type(mpi_status) :: status
  integer :: ierror
  integer :: joff
  type(mpi_datatype) :: type1, xy2zy_1, xy2zy_2
  type(mpi_datatype) :: type_cFP_full, type_cFF_full
  type(mpi_win) :: temp_fft_win

  ! Movie
  real(rkind)   XcMovie, YcMovie, ZcMovie
  integer       NxMovie, NyMovie, NzMovie, rankyMovie, rankzMovie


  ! Physical Domain
  real(rkind)         Lx, Ly, Lz

  ! BCs & Values
  integer   u_BC_Xmin, v_BC_Xmin, w_BC_Xmin, th_BC_Xmin(1:N_th)
  integer   u_BC_Xmax, v_BC_Xmax, w_BC_Xmax, th_BC_Xmax(1:N_th)
  integer   u_BC_Ymin, v_BC_Ymin, w_BC_Ymin, th_BC_Ymin(1:N_th)
  integer   u_BC_Ymax, v_BC_Ymax, w_BC_Ymax, th_BC_Ymax(1:N_th)
  integer   u_BC_Zmin, v_BC_Zmin, w_BC_Zmin, th_BC_Zmin(1:N_th)
  integer   u_BC_Zmax, v_BC_Zmax, w_BC_Zmax, th_BC_Zmax(1:N_th)

  real(rkind)   u_BC_Xmin_c1, v_BC_Xmin_c1, w_BC_Xmin_c1
  real(rkind)   u_BC_Ymin_c1, v_BC_Ymin_c1, w_BC_Ymin_c1
  real(rkind)   u_BC_Zmin_c1, v_BC_Zmin_c1, w_BC_Zmin_c1
  real(rkind)   th_BC_Xmin_c1(1:N_th), th_BC_Ymin_c1(1:N_th), th_BC_Zmin_c1(1:N_th)
  real(rkind)   u_BC_Xmax_c1, v_BC_Xmax_c1, w_BC_Xmax_c1
  real(rkind)   u_BC_Ymax_c1, v_BC_Ymax_c1, w_BC_Ymax_c1
  real(rkind)   u_BC_Zmax_c1, v_BC_Zmax_c1, w_BC_Zmax_c1
  real(rkind)   th_BC_Xmax_c1(1:N_th), th_BC_Ymax_c1(1:N_th), th_BC_Zmax_c1(1:N_th)



contains

  !----*|--.---------.---------.---------.---------.---------.---------.-|------
  subroutine init_mpi
    !----*|--.---------.---------.---------.---------.---------.---------.-|-----|
    ! This subroutine initializes all MPI variables

    integer bl(2)
    integer(kind=mpi_address_kind) :: disp(2)
    type(mpi_datatype) :: types(2)

    integer iprocs

    type(mpi_comm) :: comm_cart
    integer dims(2)
    logical mflag(2), perdim(2)

    integer, dimension(3) :: sizes, subsizes, starts


    integer i, j, k

    character(len=55) fname


    call mpi_init(ierror)
    call mpi_comm_size(mpi_comm_world, iprocs, ierror)
    call mpi_comm_rank(mpi_comm_world, rank, ierror)

    if (Nprocs /= iprocs) then
      if (rank == 0) write (*, '("Error: Compiled with ", I10 " cores, running with ", I10, " cores.")') Nprocs, iprocs
      call mpi_finalize(ierror)
      stop
    end if

    if (mod(Nprocs, NprocY) /= 0) then
      if (rank == 0) write (*, '("Error: Nprocs is Not a Multiple of NprocZ!")')
      call mpi_finalize(ierror)
      stop
    end if

#ifdef SHARED_MEMORY
    ! Optimise communication by forcing all of mpi_comm_z on a single node...

    call mpi_comm_split_type(mpi_comm_world, mpi_comm_type_shared, 0, mpi_info_null, mpi_comm_shared, ierror)
    call mpi_comm_rank(mpi_comm_shared, rankShared, ierror)
    call mpi_comm_size(mpi_comm_shared, NprocShared, ierror)
    if (NprocShared < NprocZ) then
      if (rank == 0) write (*, '("Error: NprocShared must be larger than NprocY!")')
      call mpi_finalize(ierror)
      stop
    end if
    if (mod(NprocShared, NprocZ) /= 0) then
      if (rank == 0) write (*, '("Error: NprocShared is Not a Multiple of NprocY!")')
      call mpi_finalize(ierror)
      stop
    end if

    ! Only need NprocZ, so adjust the rankShared => rankZ
    rankZ = mod(rankShared, NprocZ)

    ! Split by rankZ, so that we can get the Nnodes & colour for sorting
    call mpi_comm_split(mpi_comm_world, rankZ, 0, mpi_comm_y)
    call mpi_comm_rank(mpi_comm_y, rankY, ierror)
    call mpi_comm_split(mpi_comm_world, rankY, 0, mpi_comm_z)

#else
    ! All Distributed Communications (diablo2)

    dims(2) = NprocY
    dims(1) = NprocZ
    mflag(:) = .false.

    call mpi_cart_create(mpi_comm_world, 2, dims, mflag, .false., &
                         comm_cart, ierror)
    perdim = (/ .false. , .true.  /)
    call mpi_cart_sub(comm_cart, perdim, mpi_comm_y, ierror)
    perdim = (/ .true.  , .false. /)
    call mpi_cart_sub(comm_cart, perdim, mpi_comm_z, ierror)

    call mpi_comm_rank(mpi_comm_y, rankY, ierror)
    call mpi_comm_rank(mpi_comm_z, rankZ, ierror)

#endif



    if (rank == 0) &
      write (*, '("MPI Initialised with ", I10, " processors")') Nprocs
    call mpi_barrier(mpi_comm_world, ierror)
    if (verbosity > 2 .and. rank == 0) &
      write (*, '("Rank, rankZ, rankY: " 3I10)') rank, rankY, rankZ


    ! ------------------------------
    !  Define datatypes
    ! ------------------------------

    ! The full local 3D sub-array of FP to be sent to a particular rank:  temp_fft(,,j)
    sizes = (/ Nx/2+1, Nzp+2, Nyp+2 /)
    subsizes = (/ Nxp, Nzp,   Nyp+2 /)
    starts = (/ 0, 0, 0 /)
    call mpi_type_create_subarray(3, sizes, subsizes, starts, mpi_order_fortran, &
                                  mpi_double_complex, type_cFP_full, ierror)
    call mpi_type_commit(type_cFP_full, ierror)


    ! The full local 3D sub-array of FF to be sent to a particular rank:  e.g. cu1(,,j)
    sizes = (/ Nxp+1, Nz+2, Nyp+2 /)
    subsizes = (/ Nxp, Nzp, Nyp+2 /)
    starts = (/ 0, 0, 0 /)
    call mpi_type_create_subarray(3, sizes, subsizes, starts, mpi_order_fortran, &
                                  mpi_double_complex, type_cFF_full, ierror)
    call mpi_type_commit(type_cFF_full, ierror)


#ifdef SHARED_MEMORY

    if (rank == 0) &
      write (*, '("Optimising diablo for Shared Memory X-Y Communication...")')

#else

    if (rank == 0) &
      write (*, '("Using Only Distributed Memory MPI")')

    ! Need to translate each one with mpi upper bound for mpi_alltoall

    bl(1:2) = (/1, 1/)
    disp(1:2) = (/0, Nxp * 16/) ! Spacing between consecutive block starts in temp_fft
    types = (/type_cFP_full, mpi_ub/)

    call mpi_type_create_struct(2, bl, disp, types, xy2zy_1, ierror)
    call mpi_type_commit(xy2zy_1, ierror)
    call mpi_type_free(type_cFP_full, ierror)


    bl(1:2) = (/1, 1/)
    disp(1:2) = (/0, Nzp * (Nxp + 1) * 16/) ! Spacing between consecutive block starts in e.g. cu1
    types = (/type_cFF_full, mpi_ub/)

    call mpi_type_create_struct(2, bl, disp, types, xy2zy_2, ierror)
    call mpi_type_commit(xy2zy_2, ierror)
    call mpi_type_free(type_cFF_full, ierror)


#endif



    return
  end




  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine create_grid_chan
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    character(len=55) fname
    integer i, j, k

    if (rank == 0) &
      write (*, '("Fourier in X")')
    do i = 0, Nx
      gx(i) = (i * Lx) / Nx
      dx(i) = Lx / Nx
      if (verbosity > 3 .and. rank == 0) &
        write (*, *) 'GX(', i, ') = ', gx(i)
    end do
    if (rank == 0) &
      write (*, '("Fourier in Y")')
    do k = 0, Nz
      gz(k) = (k * Lz) / Nz
      dz(k) = Lz / Nz
      if (rank == 0 .and. verbosity > 3) &
        write (*, *) 'GY(', k, ') = ', gz(k)
    end do
    if (rank == 0) &
      write (*, '("Finite-difference in Z")')

    fname = 'grid.h5'
    call mpi_barrier(mpi_comm_world, ierror)
    call ReadGridHDF5(fname, 2)


    ! Define grid spacing
    do j = 1, Nyp + 1
      dy(j) = (gyf(j) - gyf(j - 1))
    end do
    do j = 1, Nyp
      dyf(j) = (gy(j + 1) - gy(j))
    end do
    ! dyf(Nyp + 1) = dyf(Nyp)

    ! Communicate  dyf(0)  and  dyf(Nyp + 1)
    if (rankY == 0) then
      dyf(0) = dyf(1)
    else
      call mpi_send(dyf(2), 1, mpi_double_precision, rankY - 1, &
                    100 + rankY, mpi_comm_y, ierror)
      call mpi_recv(dyf(0), 1, mpi_double_precision, rankY - 1, &
                    110 + rankY - 1, mpi_comm_y, status, ierror)
    end if
    if (rankY == NprocY - 1) then
      dyf(Nyp + 1) = dyf(Nyp)
    else
      call mpi_send(dyf(Nyp - 1), 1, mpi_double_precision, rankY + 1, &
                    110 + rankY, mpi_comm_y, ierror)
      call mpi_recv(dyf(Nyp + 1), 1, mpi_double_precision, rankY + 1, &
                    100 + rankY + 1, mpi_comm_y, status, ierror)
    end if

    return
  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|------
  subroutine init_chan_mpi
    !----*|--.---------.---------.---------.---------.---------.---------.-|-----|
    ! Initialize any constants here

    integer j, n

    if (rank == 0) then
      write (*, *)
      write (*, *) '             ****** In Init_Chan_MPI ******'
      write (*, *)
    end if


    ! Set the upper and lower bounds for timestepping
    if (rankY == 0) then
      jend = Nyp
      if (u_BC_Ymin == 0) then
        jstart = 2
      else if (u_BC_Ymin == 1) then
        jstart = 1
      else
        jstart = 2
      end if
      ! Now, set the indexing for the scalar equations
      do n = 1, N_th
        jend_th(n) = Nyp
        if (th_BC_Ymin(n) == 0) then
          jstart_th(n) = 2
        else if (th_BC_Ymin(n) == 1) then
          jstart_th(n) = 1
        else
          jstart_th(n) = 2
        end if
      end do
    else if (rankY == NprocY - 1) then
      jstart = 2
      if (u_BC_Ymax == 0) then
        jend = Nyp - 1
      else if (u_BC_Ymax == 1) then
        jend = Nyp
      else
        jend = Nyp - 1
      end if

      ! Set the upper and lower limits of timestepping of the scalar equations
      do n = 1, N_th
        jstart_th(n) = 2
        if (th_BC_Ymax(n) == 0) then
          jend_th(n) = Nyp - 1
        else if (th_BC_Ymax(n) == 1) then
          jend_th(n) = Nyp
        else
          jend_th(n) = Nyp - 1
        end if
      end do

    else
      ! Here, we are on a middle process
      jstart = 2
      jend = Nyp
      do n = 1, N_th
        jstart_th(n) = 2
        jend_th(n) = Nyp
      end do
    end if

    return
  end




  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine init_chan_movie
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    integer i, j, k

    ! Get the indices
    Nxmovie = int(XcMovie * Nx / Lx)
    NyMovie = nint(YcMovie * (((Nyp - 1) * Nprocs) / Ly))
    NzMovie = int(ZcMovie * Nz / Lz)

    rankzMovie = int(NzMovie / Nzp)
    NzMovie = NzMovie - rankzMovie * Nzp

    rankyMovie = -1
    if (gyf(jstart) <= YcMovie .and. gyf(jend + 1) > YcMovie) then
      rankyMovie = rankY
      i = 1
      do while (.not. &
                (gyf(i) <= YcMovie .and. gyf(i + 1) > YcMovie))
        i = i + 1
      end do
      NyMovie = i;
    end if

    if (rankY == rankyMovie .and. rankZ == rankzMovie) then
      write (*, '("Movie Parameters")')
      write (*, '("Rank = " I10)') rank
      write (*, '("  Xc = " ES26.18)') gx(Nxmovie)
      write (*, '("  Yc = " ES26.18)') gz(rankzMovie * Nzp + NzMovie)
      write (*, '("  Zc = " ES26.18)') gyf(NyMovie)
    end if

    return
  end






  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine ReadGridHDF5(fname, coord)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    use hdf5

    character(len=55) fname

    ! Identifiers
    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: filspace_id, memspace_id
    ! Identifiers
    integer(hid_t) :: gid, selspace_id
    integer(hid_t) :: plist_id_w, plist_id_d

    ! Dataset names
    integer           :: coord
    character(len=10) :: cname
    ! Dimensions in the memory and in the file
    integer(hsize_t), dimension(1)  :: dimsm, dimsf
    integer(hsize_t), dimension(1)  :: idims, imaxd

    integer(hsize_t), dimension(1)  :: count, offset
    integer(hsize_t), dimension(1)  :: stride, block, offset_m

    integer           :: Error, rHDF5, ith

    ! Initialize interface
    call h5open_f(Error)

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(h5p_file_access_f, plist_id_d, Error)
    call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world%mpi_val, &
                            mpi_info_null%mpi_val, Error)

    ! Create the file collectively
    call h5fopen_f(trim(fname), h5f_acc_rdonly_f, &
                   file_id, Error, access_prp=plist_id_d)
    call h5pclose_f(plist_id_d, Error)


    select case (coord)
    case (1)
      write (*, '("Error 235454. Not implemented yet!")')
      stop
    case (2)
      cname = 'y'

      dimsm = Nyp + 2
      dimsf = (Nyp - 1) * NprocY + 1

      ! Stride and count for number of rows and columns in each dimension
      stride = 1
      count = 1

      ! Offset determined by the rank of a processor
      block = Nyp + 1
      offset = rankY * (Nyp - 1)
    case (3)
      write (*, '("Error 235455. Not implemented yet!")')
      stop
    end select

    call h5gopen_f(file_id, "/grids", gid, Error)

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id_w, Error)
    call h5pset_dxpl_mpio_f(plist_id_w, h5fd_mpio_collective_f, &
                            Error)

    ! Dataspace in memory
    rHDF5 = 1
    call h5screate_simple_f(rHDF5, dimsm, memspace_id, &
                            Error)

    call h5dopen_f(gid, trim(cname), dset_id, Error)
    call h5dget_space_f(dset_id, filspace_id, Error)

    ! Check for the dimensions
    call h5sget_simple_extent_dims_f(filspace_id, idims, imaxd, Error)
    if (idims(1) - 1 /= dimsf(1)) then
      if (rank == 0) then
        write (*, '("Grid file and program do not match.")')
        write (*, *) '   gridfile (', trim(cname), '): ', idims - 1
        write (*, *) '   program     : ', dimsf
      end if
      call mpi_finalize(ierror)
      stop
    end if

    ! Select hyperslab in the file.
    call h5sselect_hyperslab_f(filspace_id, h5s_select_set_f, &
                               offset, count, Error, stride, block)

    offset_m = 1
    call h5sselect_hyperslab_f(memspace_id, h5s_select_set_f, &
                               offset_m, count, Error, stride, block)

    ! Write the dataset collectively
    call h5dread_f(dset_id, h5t_native_double, &
                   gy, &
                   dimsm, Error, file_space_id=filspace_id, &
                   mem_space_id=memspace_id) !, xfer_prp = plist_id_w)

    ! Close dateset
    call h5sclose_f(filspace_id, Error)
    call h5dclose_f(dset_id, Error)

    ! Close the dataspace for the memory
    call h5sclose_f(memspace_id, Error)

    ! Close the properties for the reading
    call h5pclose_f(plist_id_w, Error)

    ! Close groups
    call h5gclose_f(gid, Error)
    call h5fclose_f(file_id, Error)
    call h5close_f(Error)

    ! Calculate the GYF in the interior
    do ith = 1, Nyp
      gyf(ith) = 0.5 * (gy(ith) + gy(ith + 1))
    end do

    ! ###############################
    !    Get the outer ghost cells
    ! ###############################

    ! in the lower part of the domain ...
    if (rankY == 0) then
      gyf(0) = 2.0 * gyf(1) - gyf(2)
    else
      call mpi_send(gyf(2), 1, mpi_double_precision, rankY - 1, &
                    100 + rankY, mpi_comm_y, ierror)
      call mpi_recv(gyf(0), 1, mpi_double_precision, rankY - 1, &
                    110 + rankY - 1, mpi_comm_y, status, ierror)
    end if

    ! in the lower part of the domain ...
    if (rankY == NprocY - 1) then
      gyf(Nyp + 1) = 2.d0 * gyf(Nyp) - gyf(Nyp - 1)
    else
      call mpi_send(gyf(Nyp - 1), 1, mpi_double_precision, rankY + 1, &
                    110 + rankY, mpi_comm_y, ierror)
      call mpi_recv(gyf(Nyp + 1), 1, mpi_double_precision, rankY + 1, &
                    100 + rankY + 1, mpi_comm_y, status, ierror)
    end if


  end subroutine




end module domain
