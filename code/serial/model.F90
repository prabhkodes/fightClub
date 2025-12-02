program atmosphere_model
  use mpi
  use calculation_types, only : wp
  use module_physics, only : dt, oldstat, newstat, flux, tend, ref
  use module_physics, only : init, finalize
  use module_physics, only : rungekutta, total_mass_energy
  use module_output, only : create_output, write_record, close_output
  use dimensions , only : sim_time, output_freq, set_dimensions
  use dimensions , only : default_nx, default_sim_time, default_output_freq
  use iodir, only : stdout
  use parallel_timer
#ifdef USE_OPENACC
  use openacc
#endif
#ifdef USE_OPENMP
  use omp_lib
#endif
  implicit none

  real(wp) :: etime
  real(wp) :: ptime
  real(wp) :: output_counter
  real(wp) :: pctime
  real(wp) :: mass0, te0
  real(wp) :: mass1, te1
  integer(8) :: t1, t2, rate

  integer :: ierr
  integer :: my_rank
  integer :: dt_values(8)
  integer :: arg_count
  integer :: arg_status
  integer :: nx_cli
  real(wp) :: sim_time_cli
  real(wp) :: output_freq_cli
  character(len=64) :: arg
  real(8) :: final_duration
  integer :: num_procs


  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

#ifdef USE_OPENACC
  ! Assign each MPI rank to a different GPU (round-robin)
  block
    integer :: num_gpus, my_gpu
    num_gpus = acc_get_num_devices(acc_device_nvidia)
    if (num_gpus > 0) then
      my_gpu = mod(my_rank, num_gpus)
      call acc_set_device_num(my_gpu, acc_device_nvidia)
    end if
  end block
#endif

  if (my_rank == 0) then
    write(stdout,*) "================= Execution Info =================="
    write(stdout,'(a,i4)') "  Number of MPI tasks: ", num_procs
#ifdef USE_OPENMP
    write(stdout,'(a,i4)') "  Number of OpenMP threads: ", omp_get_max_threads()
    write(stdout,*) "  OpenMP: ENABLED"
#else
    write(stdout,*) "  OpenMP: DISABLED"
#endif
#ifdef USE_OPENACC
    write(stdout,'(a,i4)') "  Number of GPUs available: ", acc_get_num_devices(acc_device_nvidia)
    write(stdout,*) "  OpenACC: ENABLED"
#else
    write(stdout,*) "  OpenACC: DISABLED"
#endif
    write(stdout,*) "==================================================="
  end if

#ifdef USE_OPENACC
  ! Print GPU assignment for each MPI rank
  write(stdout,'(a,i4,a,i2)') "  MPI Rank ", my_rank, " -> GPU ", acc_get_device_num(acc_device_nvidia)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

  nx_cli = default_nx
  sim_time_cli = default_sim_time
  output_freq_cli = default_output_freq

  arg_count = command_argument_count()

  if (arg_count >= 1) then
    call get_command_argument(1, arg)
    read(arg, *, iostat=arg_status) nx_cli
    if (arg_status /= 0) nx_cli = default_nx
  end if

  if (arg_count >= 2) then
    call get_command_argument(2, arg)
    read(arg, *, iostat=arg_status) sim_time_cli
    if (arg_status /= 0) sim_time_cli = default_sim_time
  end if

  if (arg_count >= 3) then
    call get_command_argument(3, arg)
    read(arg, *, iostat=arg_status) output_freq_cli
    if (arg_status /= 0) output_freq_cli = default_output_freq
  end if

  call set_dimensions(nx_cli, sim_time_cli, output_freq_cli)

  ! --- PRINT START TIME (Rank 0 Only) ---
  if (my_rank == 0) then
      call date_and_time(values=dt_values)
      write(stdout, *) 'SIMPLE ATMOSPHERIC MODEL STARTING.'
      ! Format: HH:MM:SS.Milliseconds
      write(stdout, '(a,I2.2,a,I2.2,a,I2.2,a,I3.3)') &
        "Wall Clock Start: ", dt_values(5), ":", dt_values(6), ":", dt_values(7), ".", dt_values(8)
  end if

  call init(etime,output_counter,dt)
  call total_mass_energy(mass0,te0)
  call create_output( )
  call write_record(oldstat,ref,etime)

#ifdef USE_OPENACC
  !$acc enter data copyin(oldstat%mem, newstat%mem, flux%mem, tend%mem) &
  !$acc              copyin(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
#endif

  ! Get initial tick count and the clock rate (ticks per second)
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
#ifdef USE_OPENACC
      !$acc update self(oldstat%mem)
#endif
      call write_record(oldstat,ref,etime)
    end if

  end do

  call total_mass_energy(mass1,te1)
  call close_output( )

#ifdef USE_OPENACC
  !$acc exit data delete(oldstat%mem, newstat%mem, flux%mem, tend%mem) &
  !$acc             delete(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
#endif

  if (my_rank == 0) then
      write(stdout,*) "----------------- Atmosphere check ----------------"
      write(stdout,*) "Fractional Delta Mass  : ", (mass1-mass0)/mass0
      write(stdout,*) "Fractional Delta Energy: ", (te1-te0)/te0
      write(stdout,*) "---------------------------------------------------"
  end if

  call finalize()
  call print_timing_results()
  
  ! Get final tick count
  call system_clock(t2)

  ! --- PRINT END TIME (Rank 0 Only) ---
  if (my_rank == 0) then
      call date_and_time(values=dt_values)
      write(stdout,*) "SIMPLE ATMOSPHERIC MODEL RUN COMPLETED."
      
      ! 1. Print Wall Clock with Milliseconds
      write(stdout, '(a,I2.2,a,I2.2,a,I2.2,a,I3.3)') &
        "Wall Clock End:   ", dt_values(5), ":", dt_values(6), ":", dt_values(7), ".", dt_values(8)
      
      ! 2. Print Duration with Microsecond Precision
      final_duration = dble(t2-t1)/dble(rate)
      write(stdout, '(a,F18.6,a)') "USED CPU TIME:    ", final_duration, " seconds"
  end if

  call MPI_FINALIZE(ierr)

end program atmosphere_model
