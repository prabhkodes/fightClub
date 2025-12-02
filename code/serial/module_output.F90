module module_output
  use mpi
  use calculation_types, only : wp, iowp
  use parallel_parameters, only : i_beg, k_beg
  ! We need i_start_global and i_end_global to know where to paste data
  use dimensions, only : nx, nz, nx_global, nprocs, myrank, i_start_global, i_end_global
  use module_types, only : atmospheric_state, reference_state
  use iodir, only : stderr
  use netcdf

  implicit none

  private

  public :: create_output
  public :: write_record
  public :: close_output

  ! Local arrays for calculation
  real(wp), allocatable, dimension(:,:) :: dens
  real(wp), allocatable, dimension(:,:) :: uwnd
  real(wp), allocatable, dimension(:,:) :: wwnd
  real(wp), allocatable, dimension(:,:) :: theta

  ! Global arrays (Allocated ONLY on Rank 0)
  real(wp), allocatable, dimension(:,:) :: global_dens
  real(wp), allocatable, dimension(:,:) :: global_uwnd
  real(wp), allocatable, dimension(:,:) :: global_wwnd
  real(wp), allocatable, dimension(:,:) :: global_theta

  ! NetCDF IDs
  integer :: ncid
  integer :: dens_varid, uwnd_varid, wwnd_varid, theta_varid, t_varid
  integer :: rec_out

  contains

  subroutine create_output
    implicit none
    integer :: t_dimid, x_dimid, z_dimid
    
    ! 1. Allocate LOCAL arrays
    allocate(dens(nx,nz))
    allocate(uwnd(nx,nz))
    allocate(wwnd(nx,nz))
    allocate(theta(nx,nz))

    ! 2. Setup Global Arrays and File (RANK 0 ONLY)
    if (myrank == 0) then
      allocate(global_dens(nx_global, nz))
      allocate(global_uwnd(nx_global, nz))
      allocate(global_wwnd(nx_global, nz))
      allocate(global_theta(nx_global, nz))

      ! Create ONE single file
      call ncwrap(nf90_create('output.nc', nf90_clobber, ncid), __LINE__)

      ! Define dimensions using GLOBAL size
      call ncwrap(nf90_def_dim(ncid, 'time', nf90_unlimited, t_dimid), __LINE__)
      call ncwrap(nf90_def_dim(ncid, 'x', nx_global, x_dimid), __LINE__) 
      call ncwrap(nf90_def_dim(ncid, 'z', nz, z_dimid), __LINE__)

      call ncwrap(nf90_def_var(ncid,'time',iowp,[t_dimid],t_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'rho',iowp, [x_dimid,z_dimid,t_dimid],dens_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'u',iowp,   [x_dimid,z_dimid,t_dimid],uwnd_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'w',iowp,   [x_dimid,z_dimid,t_dimid],wwnd_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'theta',iowp,[x_dimid,z_dimid,t_dimid],theta_varid), __LINE__)
      
      call ncwrap(nf90_enddef(ncid), __LINE__)
    end if

    rec_out = 1
  end subroutine create_output

  subroutine write_record(atmostat,ref,etime)
    implicit none
    type(atmospheric_state), intent(in) :: atmostat
    type(reference_state), intent(in) :: ref
    real(wp), intent(in) :: etime
    integer :: i, k, src, ierr
    integer :: r_nx, r_istart
    integer, dimension(1) :: st1, ct1
    integer, dimension(3) :: st3, ct3
    real(wp), dimension(1) :: etimearr
    real(wp), allocatable, dimension(:,:) :: temp_buf
    integer :: status(MPI_STATUS_SIZE)

    ! --- 1. Calculate Local Fields ---
    do k = 1, nz
      do i = 1, nx
        dens(i,k) = atmostat%dens(i,k)
        uwnd(i,k) = atmostat%umom(i,k)/(ref%density(k)+dens(i,k))
        wwnd(i,k) = atmostat%wmom(i,k)/(ref%density(k)+dens(i,k))
        theta(i,k) = (atmostat%rhot(i,k) + ref%denstheta(k)) / &
            (ref%density(k) + dens(i,k)) - ref%denstheta(k)/ref%density(k)
      end do
    end do

    ! --- 2. Smart Gather to Rank 0 ---
    if (myrank == 0) then
       ! A. Copy MY OWN data to the global array first
       ! We use our known global indices to place it correctly
       global_dens(i_start_global:i_end_global, :) = dens(:,:)
       global_uwnd(i_start_global:i_end_global, :) = uwnd(:,:)
       global_wwnd(i_start_global:i_end_global, :) = wwnd(:,:)
       global_theta(i_start_global:i_end_global, :) = theta(:,:)

       ! B. Receive data from all other processors
       do src = 1, nprocs - 1
          ! Receive the remote nx and start index so we know where to paste
          call MPI_Recv(r_nx, 1, MPI_INTEGER, src, 100, MPI_COMM_WORLD, status, ierr)
          call MPI_Recv(r_istart, 1, MPI_INTEGER, src, 101, MPI_COMM_WORLD, status, ierr)
          
          ! Allocate a temp buffer of the correct size
          allocate(temp_buf(r_nx, nz))

          ! Recv Density
          call MPI_Recv(temp_buf, r_nx*nz, MPI_DOUBLE_PRECISION, src, 200, MPI_COMM_WORLD, status, ierr)
          global_dens(r_istart : r_istart + r_nx - 1, :) = temp_buf(:,:)
          
          ! Recv U-Wind
          call MPI_Recv(temp_buf, r_nx*nz, MPI_DOUBLE_PRECISION, src, 201, MPI_COMM_WORLD, status, ierr)
          global_uwnd(r_istart : r_istart + r_nx - 1, :) = temp_buf(:,:)

          ! Recv W-Wind
          call MPI_Recv(temp_buf, r_nx*nz, MPI_DOUBLE_PRECISION, src, 202, MPI_COMM_WORLD, status, ierr)
          global_wwnd(r_istart : r_istart + r_nx - 1, :) = temp_buf(:,:)

          ! Recv Theta
          call MPI_Recv(temp_buf, r_nx*nz, MPI_DOUBLE_PRECISION, src, 203, MPI_COMM_WORLD, status, ierr)
          global_theta(r_istart : r_istart + r_nx - 1, :) = temp_buf(:,:)

          deallocate(temp_buf)
       end do

       ! --- 3. Write to File (Only Rank 0) ---
       st1 = [ rec_out ]
       ct1 = [ 1 ]
       etimearr(1) = etime
       call ncwrap(nf90_put_var(ncid, t_varid, etimearr, st1, ct1), __LINE__)

       st3 = [ 1, 1, rec_out ]
       ct3 = [ nx_global, nz, 1 ]
       call ncwrap(nf90_put_var(ncid, dens_varid, global_dens, st3, ct3), __LINE__)
       call ncwrap(nf90_put_var(ncid, uwnd_varid, global_uwnd, st3, ct3), __LINE__)
       call ncwrap(nf90_put_var(ncid, wwnd_varid, global_wwnd, st3, ct3), __LINE__)
       call ncwrap(nf90_put_var(ncid, theta_varid, global_theta, st3, ct3), __LINE__)

    else
       ! --- Workers Send Data to Rank 0 ---
       ! Send my size and start position
       call MPI_Send(nx, 1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, ierr)
       call MPI_Send(i_start_global, 1, MPI_INTEGER, 0, 101, MPI_COMM_WORLD, ierr)
       
       ! Send the actual data arrays
       call MPI_Send(dens, nx*nz, MPI_DOUBLE_PRECISION, 0, 200, MPI_COMM_WORLD, ierr)
       call MPI_Send(uwnd, nx*nz, MPI_DOUBLE_PRECISION, 0, 201, MPI_COMM_WORLD, ierr)
       call MPI_Send(wwnd, nx*nz, MPI_DOUBLE_PRECISION, 0, 202, MPI_COMM_WORLD, ierr)
       call MPI_Send(theta, nx*nz, MPI_DOUBLE_PRECISION, 0, 203, MPI_COMM_WORLD, ierr)
    end if

    rec_out = rec_out + 1
  end subroutine write_record

  subroutine close_output
    implicit none
    ! Deallocate Local
    if ( allocated(dens) ) deallocate(dens, uwnd, wwnd, theta)
    
    if (myrank == 0) then
      ! Deallocate Global
      if ( allocated(global_dens) ) deallocate(global_dens, global_uwnd, global_wwnd, global_theta)
      call ncwrap(nf90_close(ncid), __LINE__)
    end if
  end subroutine close_output

  subroutine ncwrap(ierr,line)
    implicit none
    integer, intent(in) :: ierr
    integer, intent(in) :: line
    if (ierr /= nf90_noerr) then
      write(stderr,*) 'NetCDF Error (Rank ', myrank, ') at line: ', line
      write(stderr,*) nf90_strerror(ierr)
      stop
    end if
  end subroutine ncwrap

end module module_output