module calculation_types
  use netcdf, only : nf90_real, nf90_double
  use iso_fortran_env
  implicit none
  public
  integer, parameter :: wp = real64
  integer, parameter :: iowp = nf90_double
  integer, parameter :: ip = int32
end module calculation_types

module physical_constants
  use calculation_types, only : wp
  implicit none
  public
  real(wp), parameter :: pi = 3.14159265358979323846264338327_wp
  real(wp), parameter :: hpi = 0.5_wp * pi
  real(wp), parameter :: grav = 9.80665_wp
  real(wp), parameter :: boltzk = 1.3806490e-23_wp
  real(wp), parameter :: navgdr = 6.02214076e23_wp
  real(wp), parameter :: rgasmol = navgdr*boltzk
  real(wp), parameter :: amd = 28.96454_wp
  real(wp), parameter :: amw = 18.01528_wp
  real(wp), parameter :: rgas = (rgasmol/amd)*1000.0_wp
  real(wp), parameter :: cp = 3.5_wp*rgas
  real(wp), parameter :: cv = 2.5_wp*rgas
  real(wp), parameter :: rd = rgas
  real(wp), parameter :: p0 = 1.e5_wp
  real(wp), parameter :: t0 = 298.0_wp
  real(wp), parameter :: cdocv = 3.5_wp/2.5_wp
  real(wp), parameter :: cvocd = 2.5_wp/3.5_wp
  real(wp), parameter :: rdocp = rd/cp
  real(wp), parameter :: rdocv = rd/cv
  real(wp), parameter :: c0 = (rd**cdocv)*p0**(-rdocv)
  real(wp), parameter :: theta0 = 300.0_wp
  real(wp), parameter :: exner0 = 1.0_wp
end module physical_constants

module physical_parameters
  use calculation_types, only : wp
  implicit none
  public
  real(wp), parameter :: xlen = 20000.0_wp
  real(wp), parameter :: zlen = 10000.0_wp
  real(wp), parameter :: hxlen = 0.5_wp*xlen
  real(wp), parameter :: p1 = xlen/10.0_wp
  real(wp), parameter :: p2 = 4.0_wp*p1
  real(wp), parameter :: p3 = 2.0_wp*p1
  real(wp), parameter :: mnt_width = xlen/8.0_wp
  real(wp), parameter :: hv_beta = 0.05_wp
  real(wp), parameter :: cfl = 1.5_wp
  real(wp), parameter :: max_speed = 450.0_wp
  real(wp), parameter :: x0 = mnt_width
  real(wp), parameter :: z0 = 1000.0_wp
  real(wp), parameter :: xrad = 500.0_wp
  real(wp), parameter :: zrad = 500.0_wp
  real(wp), parameter :: amp = 0.01_wp
end module physical_parameters

module indexing
  implicit none
  public
  integer, parameter :: NVARS = 4
  integer, parameter :: STEN_SIZE = 4
  integer, parameter :: I_DENS  = 1
  integer, parameter :: I_UMOM  = 2
  integer, parameter :: I_WMOM  = 3
  integer, parameter :: I_RHOT  = 4
  integer, parameter :: DIR_X = 1
  integer, parameter :: DIR_Z = 2
end module indexing

module iodir
  use iso_fortran_env
  public
  integer, public, parameter :: stdin = input_unit
  integer, public, parameter :: stdout = output_unit
  integer, public, parameter :: stderr = error_unit
end module iodir

module parallel_parameters
  implicit none
  public
  integer, parameter :: i_beg = 1
  integer, parameter :: k_beg = 1
  integer, parameter :: hs = 2
end module parallel_parameters

module legendre_quadrature
  use calculation_types, only : wp
  implicit none
  public
  integer , parameter :: nqpoints = 3
  real(wp), parameter :: qpoints(nqpoints) =    &
      [ 0.112701665379258311482073460022E0_wp , &
        0.500000000000000000000000000000E0_wp , &
        0.887298334620741688517926539980E0_wp ]
  real(wp), parameter :: qweights(nqpoints) =   &
      [ 0.277777777777777777777777777779E0_wp , &
        0.444444444444444444444444444444E0_wp , &
        0.277777777777777777777777777779E0_wp ]
end module legendre_quadrature

module dimensions
  use calculation_types, only : wp
  use physical_parameters, only : zlen, xlen
  use indexing
  implicit none
  public
  integer, parameter :: default_nx = 100
  real(wp), parameter :: default_sim_time = 1000.0_wp
  real(wp), parameter :: default_output_freq = 10.0_wp

  integer :: nx = default_nx
  integer :: nz = int(default_nx * zlen/xlen)
  real(wp) :: sim_time = default_sim_time
  real(wp) :: output_freq = default_output_freq

  ! MPI domain decomposition variables
  integer :: nprocs = 1           ! Total number of MPI processes
  integer :: myrank = 0           ! Rank of this process
  integer :: nx_local             ! Local nx for this process
  integer :: nx_global            ! Global nx (original nx)
  integer :: i_start_global       ! Global starting index for this process
  integer :: i_end_global         ! Global ending index for this process
  integer :: left_rank            ! Left neighbor MPI rank
  integer :: right_rank           ! Right neighbor MPI rank

contains

  subroutine set_dimensions(nx_in, sim_time_in, output_freq_in)
    integer, intent(in) :: nx_in
    real(wp), intent(in) :: sim_time_in
    real(wp), intent(in) :: output_freq_in

    nx = nx_in
    nx_global = nx_in
    sim_time = sim_time_in
    output_freq = output_freq_in
    nz = int(nx * zlen/xlen)
  end subroutine set_dimensions

  subroutine setup_domain_decomposition(rank, num_procs)
    implicit none
    integer, intent(in) :: rank, num_procs
    integer :: base_nx, remainder

    myrank = rank
    nprocs = num_procs
    nx_global = nx

    ! Divide nx among processes
    base_nx = nx_global / nprocs
    remainder = mod(nx_global, nprocs)

    ! Distribute remainder among first processes
    if (myrank < remainder) then
      nx_local = base_nx + 1
      i_start_global = myrank * (base_nx + 1) + 1
    else
      nx_local = base_nx
      i_start_global = remainder * (base_nx + 1) + (myrank - remainder) * base_nx + 1
    end if
    i_end_global = i_start_global + nx_local - 1

    ! Update nx to local value for this process
    nx = nx_local

    ! Neighbors (periodic in x)
    left_rank = mod(myrank - 1 + nprocs, nprocs)
    right_rank = mod(myrank + 1, nprocs)

  end subroutine setup_domain_decomposition

end module dimensions
