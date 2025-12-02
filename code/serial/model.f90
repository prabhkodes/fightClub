program atmosphere_model
  use mpi
  use calculation_types, only : wp
  use module_physics, only : dt, oldstat, newstat, flux, tend, ref
  use module_physics, only : init, finalize
  use module_physics, only : rungekutta, total_mass_energy
  use module_output, only : create_output, write_record, close_output
  use dimensions , only : sim_time, output_freq, nx 
  use iodir, only : stdout
  use parallel_timer
  use parallel_parameters                            
  implicit none

  real(wp) :: etime
  real(wp) :: ptime
  real(wp) :: output_counter
  real(wp) :: pctime
  real(wp) :: mass0, te0
  real(wp) :: mass1, te1
  integer(8) :: t1, t2, rate

  integer :: dt_values(8) 
  real(8) :: final_duration
  integer(8) :: pp


  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, ierr) 

  call parallel_setup(nx)

  ! --- DEBUG: PRINT DOMAIN COVERAGE ---
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  do pp = 0, n_procs-1
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     if (my_rank == pp) then
        write(stdout, '(A,I4,A,I4,A,I5,A,I5,A,I5)') &
           "Rank ", my_rank, ": NX_Local=", nx_local, &
           " Starts=", i_beg_local, " Ends=", i_end_local
     end if
  end do
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if (my_rank == 0) then
      write(stdout, *) "Expected Global NX:", nx
  end if
  ! ------------------------------------

  if (my_rank == 0) then
      call date_and_time(values=dt_values)
      write(stdout, *) 'SIMPLE ATMOSPHERIC MODEL STARTING.'
      write(stdout, '(a,I2.2,a,I2.2,a,I2.2,a,I3.3)') &
        "Wall Clock Start: ", dt_values(5), ":", dt_values(6), ":", dt_values(7), ".", dt_values(8)
  end if

  call init(etime,output_counter,dt)
  call total_mass_energy(mass0,te0)

  ! --- DEBUG PRINT 1: INITIAL STATE ---
  if (my_rank == 0) then
     write(stdout, '(a, E25.16)') "DEBUG: Initial Total Mass: ", mass0
  end if
  ! ------------------------------------
  call create_output( )
  call write_record(oldstat,ref,etime)

  call system_clock(t1, rate) 

  ptime = int(sim_time/10.0)
  do while (etime < sim_time)

    if (etime + dt > sim_time) dt = sim_time - etime

    call rungekutta(oldstat,newstat,flux,tend,dt)

    if ( mod(etime,ptime) < dt ) then
      pctime = (etime/sim_time)*100.0_wp
      if (my_rank == 0) then
         write(stdout,'(1x,a,i2,a)') 'TIME PERCENT : ', int(pctime), '%'
      end if
    end if

    etime = etime + dt
    output_counter = output_counter + dt

    if (output_counter >= output_freq) then
      output_counter = output_counter - output_freq
      call write_record(oldstat,ref,etime)
    end if

  end do

  call total_mass_energy(mass1,te1)
  ! --- DEBUG PRINT 1: INITIAL STATE ---
  if (my_rank == 0) then
     write(stdout, '(a, E25.16)') "DEBUG: Final Total Mass: ", mass1
  end if
  call close_output( )

  if (my_rank == 0) then
      write(stdout,*) "----------------- Atmosphere check ----------------"
      write(stdout,*) "Fractional Delta Mass  : ", (mass1-mass0)/mass0
      write(stdout,*) "Fractional Delta Energy: ", (te1-te0)/te0
      write(stdout,*) "---------------------------------------------------"
  end if

  call finalize()
  call print_timing_results()
  
  call system_clock(t2)

  if (my_rank == 0) then
      call date_and_time(values=dt_values)
      write(stdout,*) "SIMPLE ATMOSPHERIC MODEL RUN COMPLETED."
      
      write(stdout, '(a,I2.2,a,I2.2,a,I2.2,a,I3.3)') &
        "Wall Clock End:   ", dt_values(5), ":", dt_values(6), ":", dt_values(7), ".", dt_values(8)
      
      final_duration = dble(t2-t1)/dble(rate)
      write(stdout, '(a,F18.6,a)') "USED CPU TIME:    ", final_duration, " seconds"
  end if

  call MPI_FINALIZE(ierr)

end program atmosphere_model