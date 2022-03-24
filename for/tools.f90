module tools
  use domain
  use parameters
  implicit none
  save

contains

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine get_minimum_mpi(val)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) val, vmin

    call mpi_allreduce(val, vmin, 1, mpi_double_precision, &
                       mpi_min, mpi_comm_world, ierror)

    val = vmin

    return

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine get_maximum_mpi(val)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) val, vmax

    call mpi_allreduce(val, vmax, 1, mpi_double_precision, &
                       mpi_max, mpi_comm_world, ierror)

    val = vmax

    return

  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine end_run(flag)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    logical flag, file_exists

    flag = .false.
    ! Check for the time
    call wall_time(end_wall_time)
    if (end_wall_time - start_wall_time > wall_time_limit) then
      if (rank == 0) &
        write (*, '("STOP because of wall-time hit!")')
      flag = .true.
    end if

    if (time >= time_limit) then
      if (rank == 0) &
        write (*, '("STOP because of simulation end-time hit!")')
      flag = .true.
    end if

    inquire (file="stop.now", exist=file_exists)
    if (file_exists) then
      if (rank == 0) &
        write (*, '("STOP because of stop.now file!")')
      flag = .true.
    end if

    return
  end



  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine end_run_mpi(flag)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    logical flag

    if (rank == 0) then
      call end_run(flag)
    end if
    call mpi_bcast(flag, 1, mpi_logical, 0, mpi_comm_world, ierror)

  end





  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine wall_time(wt)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    !
    ! Return wall-clock time as seconds after Jan. 1, 2016.
    ! Support for leap year is not included anymore.
    !
    ! By using a 'save' statement, the wall-time after the first
    ! call to the subroutine could be computed, but that is not
    ! intended with the present subroutine (e.g. the history file)
    !
    implicit none

    real(kind(0.d0)) wt
    integer val(8), i, shift, day

    integer mon(12, 2)
    data mon/ &
      31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
      31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
    !
    ! Get current date and time
    ! val(1) : year
    ! val(2) : month
    ! val(3) : day
    ! val(4) : difference to GMT
    ! val(5) : hour
    ! val(6) : minute
    ! val(7) : second
    ! val(8) : 1/1000 second
    !
    call date_and_time(values=val)
    !
    ! Determine leap year
    !
    if (mod(val(1), 4) == 0) then
      if (mod(val(1), 100) == 0) then
        if (mod(val(1), 400) == 0) then
          shift = 2
        else
          shift = 1
        end if
      else
        shift = 2
      end if
    else
      shift = 1
    end if
    !
    ! Construct day of the year
    !
    day = val(3) - 1
    do i = 1, val(2) - 1
      day = day + mon(i, shift)
    end do
    !
    ! And compute wall-clock time
    !
    wt = (val(1) - 2016) * 365 * 86400 + &
         day * 86400 + val(5) * 3600 + val(6) * 60 + val(7) + dble(val(8) / 1000.d0)

  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine time_string(cdt)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    !
    ! Construct string in the format '19-DEC-2005 22:47:06'
    !
    implicit none

    integer i

    integer val(8)
    character(len=20) cdt
    character(len=3) monc

    call date_and_time(values=val)

    if (val(2) == 1) then
      monc = 'JAN'
    else if (val(2) == 2) then
      monc = 'FEB'
    else if (val(2) == 3) then
      monc = 'MAR'
    else if (val(2) == 4) then
      monc = 'APR'
    else if (val(2) == 5) then
      monc = 'MAY'
    else if (val(2) == 6) then
      monc = 'JUN'
    else if (val(2) == 7) then
      monc = 'JUL'
    else if (val(2) == 8) then
      monc = 'AUG'
    else if (val(2) == 9) then
      monc = 'SEP'
    else if (val(2) == 10) then
      monc = 'OCT'
    else if (val(2) == 11) then
      monc = 'NOV'
    else if (val(2) == 12) then
      monc = 'DEC'
    else
      monc = 'XXX'
    end if

    write (cdt, '(i2,a1,a3,a1,i4,a1,i2,a1,i2,a1,i2)') &
      val(3), '-', monc, '-', val(1), ' ', val(5), ':', val(6), ':', val(7)
    do i = 1, 2
      if (cdt(i:i) == ' ') then
        cdt(i:i) = '0'
      end if
    end do
    do i = 13, 20
      if (cdt(i:i) == ' ') then
        cdt(i:i) = '0'
      end if
    end do

  end subroutine time_string



  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine swapzy(in, out)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) in(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    real(rkind) out(1:Nx, 1:Nyp, 1:Nzp)
    integer x, z, y

    out = 0.d0
    do x = 0, Nx - 1
      do y = 1, Nyp
        do z = 0, Nzp - 1
          out(x + 1, y, z + 1) = in(x, z, y)
        end do
      end do
    end do

  end subroutine

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine swapyz(in, out)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) out(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    real(rkind) in(1:Nx, 1:Nyp, 1:Nzp)
    integer x, z, y

    out = 0.d0
    do x = 0, Nx - 1
      do y = 1, Nyp
        do z = 0, Nzp - 1
          out(x, z, y) = in(x + 1, y, z + 1)
        end do
      end do
    end do

  end subroutine


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine gf2g(var)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) var(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    integer x, z, y
    type(mpi_datatype) :: xzbox

    ! Define new data type
    call mpi_type_vector(Nzp, Nx, Nx + 2, mpi_double_precision, &
                         xzbox, ierror)
    call mpi_type_commit(xzbox, ierror)

    if (rankY /= NprocY - 1) then
      call mpi_recv(var(0, 0, Nyp + 1), 1, xzbox, rankY + 1, 101 + rankY, &
                    mpi_comm_y, status, ierror)
    else
      do x = 0, Nx - 1
        do z = 0, Nzp - 1
          var(x, z, Nyp + 1) = var(x, z, Nyp)
        end do
      end do
    end if

    do x = 0, Nx - 1
      do z = 0, Nzp - 1
        do y = Nyp, 1, -1
          var(x, z, y) = 2 * var(x, z, y) - var(x, z, y + 1)
        end do
      end do
    end do

    if (rankY /= 0) call mpi_send(var(0, 0, 2), 1, &
                                    xzbox, rankY - 1, 100 + rankY, &
                                    mpi_comm_y, ierror)

    ! Impose the values at the boundary as prescribed in the
    ! code in order to have zero mass flux
    if (rankY == 0) var(:, :, 1) = -var(:, :, 2)
    if (rankY == NprocY - 1) var(:, :, Nyp + 1) = -var(:, :, Nyp)

  end subroutine


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine g2gf(var)
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    real(rkind) var(0:Nx + 1, 0:Nzp + 1, 0:Nyp + 1)
    integer x, z, y

    do x = 0, Nx - 1
      do z = 0, Nzp - 1
        if (rankY == 0) var(x, z, 1) = var(x, z, 2)
        if (rankY == NprocY - 1) var(x, z, Nyp + 1) = var(x, z, Nyp)
        do y = 1, Nyp
          var(x, z, y) = 0.5 * (var(x, z, y) + var(x, z, y + 1))
        end do
      end do
    end do

  end subroutine






end module tools
