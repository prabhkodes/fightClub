module parallel_timer
    use mpi
    implicit none
    
    ! --- Internal Storage ---
    integer, parameter :: MAX_TIMERS = 100
    integer, parameter :: LABEL_LEN = 64
    
    type :: timer_store_t
        character(len=LABEL_LEN) :: label = ""
        real(8) :: total_time = 0.0d0
        integer :: calls = 0
        logical :: active = .false.
    end type timer_store_t
    
    ! Global storage for accumulated times
    type(timer_store_t), save :: database(MAX_TIMERS)
    integer, save :: num_timers = 0

    ! --- The Scoped Timer Object ---
    type :: CTimer
        integer :: db_index = 0
        real(8) :: start_time = 0.0d0
    contains
        final :: stop_timer ! This acts as the Destructor (~CTimer)
    end type CTimer

    ! Constructor interface to look like C++
    interface CTimer
        module procedure new_timer
    end interface

contains

    ! --- Constructor: CTimer t("Label") ---
    function new_timer(label) result(t)
        character(len=*), intent(in) :: label
        type(CTimer) :: t
        integer :: i
        
        t%start_time = MPI_Wtime()
        
        ! Find or create entry in database
        t%db_index = -1
        do i = 1, num_timers
            if (trim(database(i)%label) == trim(label)) then
                t%db_index = i
                exit
            end if
        end do
        
        ! If new, add it
        if (t%db_index == -1) then
            if (num_timers < MAX_TIMERS) then
                num_timers = num_timers + 1
                t%db_index = num_timers
                database(num_timers)%label = label
                database(num_timers)%active = .true.
            else
                print *, "Timer Error: Max timers exceeded"
            end if
        end if
    end function new_timer

    ! --- Destructor: Automatically called at end of scope ---
    impure subroutine stop_timer(this)
        type(CTimer), intent(inout) :: this
        real(8) :: end_time, duration
        
        ! Check if index is valid (prevents issues with temporary copies)
        if (this%db_index > 0) then
            end_time = MPI_Wtime()
            duration = (end_time - this%start_time) * 1.0d6 ! Convert to microseconds
            
            database(this%db_index)%total_time = database(this%db_index)%total_time + duration
            database(this%db_index)%calls = database(this%db_index)%calls + 1
            
            ! Reset index so it doesn't trigger again on deallocation
            this%db_index = 0 
        end if
    end subroutine stop_timer

    ! --- Reporting ---
    subroutine print_timing_results()
        integer :: i, rank, ierr, world_size
        real(8) :: local_t, max_t, min_t, sum_t
        integer :: calls_sum
        
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)

        if (rank == 0) then
            print *, repeat("-", 100)
            print *, "TIMING STATISTICS (Microseconds)"
            print "(A30, A15, A15, A15, A10)", "Function", "Max Time", "Min Time", "Avg Time", "Calls"
            print *, repeat("-", 100)
        end if

        do i = 1, num_timers
            local_t = database(i)%total_time
            
            ! Reduce data to get Min/Max/Avg across all ranks
            call MPI_Reduce(local_t, max_t, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(local_t, min_t, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(local_t, sum_t, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            
            ! Just taking calls from rank 0 for simplicity, or sum them up
            call MPI_Reduce(database(i)%calls, calls_sum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            if (rank == 0) then
                print "(A30, F15.2, F15.2, F15.2, I10)", &
                      trim(database(i)%label), max_t, min_t, sum_t/world_size, calls_sum
            end if
        end do
        
        if (rank == 0) print *, repeat("-", 100)
    end subroutine print_timing_results

end module parallel_timer