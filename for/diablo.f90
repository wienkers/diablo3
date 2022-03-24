!******************************************************************************|
! diablo.f90                                                        VERSION 3.0
! v2.0 JRT, 2012
! v3.0 AFW, 2021
!
! This Fortran code computes incompressible flow in a channel.
!
! Primative variables (u,v,w,p) are used, and continuity is enforced with a
! fractional step algorithm.
!
! SPATIAL DERIVATIVES:
!   The 1 & 3 (X/Z) directions are periodic and handled spectrally
!   The 2 (Y) direction is taken to be bounded by walls and handled with
!   momentum- and energy-conserving second-order central finite differences.
!
! TIME ADVANCEMENT
!   Two main approaches are implemented:
!     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
!     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
!
! A few simple high-performance programming constructs are used:
!   -> The inner 2 loops are broken out in such a way as to enable out-of-order
!      execution of these loops as much as possible, thereby leveraging
!      vector and superscalar CPU architectures.
!   -> The outer loops are fairly long (including as many operations as
!      possible inside on a single J plane of data) in order to make effective
!      use of cache.
!
!******************************************************************************|
!
! This code is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation; either version 2 of the License, or (at your
! option) any later version. This code is distributed in the hope that it
! will be useful, but WITHOUT any WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details. You should have received a
! copy of the GNU General Public License along with this code; if not,
! write to the Free Software Foundation, Inc., 59 Temple Place - Suite
! 330, Boston, MA 02111-1307, USA.
!
!******************************************************************************|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
program diablo
  use domain
  use parameters
  use flow
  use tools

  integer n
  logical flag

  call init_parameters
  call init_flow

  ! Initialize start_wall_time for run timing
  call wall_time(start_wall_time)

  first_time = .true.

  do
    time_step = time_step + 1

    if (verbosity > 2 .and. rank == 0) &
      write (*, '("Now beginning time step ", I10)') time_step

    do rk_step = 1, 3
      if (time_ad_meth == 1) call rk_chan_1
      if (time_ad_meth == 2) call rk_chan_2
    end do
    time = time + delta_t
    first_time = .false.

    ! Optionally apply a filter to the scalar field
    do n = 1, N_th
      if (filter_th(n) &
          .and. (mod(time_step, filter_int(n)) == 0)) then
        call filter_chan(n)
      end if
    end do

    call end_run_mpi(flag)


    ! Save statistics to an output file
    if (time >= save_stats_time) then

      flag_save_LES = .true.
      if (time >= save_movie_time) then
        call save_stats(.true.,.false.)
        save_movie_time = save_stats_time + save_movie_dt - save_stats_dt*1.d-5 ! Increment from stats_time
      else
        call save_stats(.false.,.false.)
      end if
      if (flag .and. use_LES) then ! Won't get into the next RK to save LES
        call save_stats_LES_OOL(.false.)
      end if
      save_stats_time = save_stats_time + save_stats_dt

      ! Timing Diagnostics
      call wall_time(end_wall_time)
      if (rank == 0) then
        write (*,'("Wall Seconds per Stats Output: ", ES15.3)') &
                  (end_wall_time - previous_wall_time)
        write (*,'("Wall Seconds per Iteration: ", ES18.3)') &
                  (end_wall_time - previous_wall_time) / float(time_step - previous_time_step)
        write (*,'("Wall Seconds per Simulation Time: ", ES12.3)') &
                  (end_wall_time - previous_wall_time) / save_stats_dt

        write (*, *)

      end if
      call wall_time(previous_wall_time)
      previous_time_step = time_step

    end if

    ! Save entire flow to a file (Only really a restart file if end.h5)
    if (time >= save_flow_time) then
      save_flow_time = save_flow_time + save_flow_dt
      call save_flow(.false.)
    end if

    if (flag) then ! We're done running
      exit
    end if

  end do

  ! Calculate and display the runtime for the simulation
  call wall_time(end_wall_time)
  if (rank == 0) then
    write (*, '("Elapsed Wall Time (sec): ", ES11.3)') end_wall_time - start_wall_time
    write (*, '("Wall Seconds per Iteration: ", ES12.5)') (end_wall_time - start_wall_time) / time_step
  end if

  call save_flow(.true.)
  !call deallocate_all
  call mpi_finalize(ierror)

end
