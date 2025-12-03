module module_output
  use mpi
  use calculation_types, only : wp, iowp
  use parallel_parameters, only : i_beg, k_beg
  use dimensions, only : nx, nz, nx_global, nprocs, myrank, i_start_global
  use module_types, only : atmospheric_state, reference_state
  use iodir, only : stderr
  use netcdf

  implicit none

  private

  public :: create_output
  public :: write_record
  public :: close_output

  real(wp), allocatable, dimension(:,:) :: dens
  real(wp), allocatable, dimension(:,:) :: uwnd
  real(wp), allocatable, dimension(:,:) :: wwnd
  real(wp), allocatable, dimension(:,:) :: theta

  ! Global arrays for gathering data (only on rank 0)
  real(wp), allocatable, dimension(:,:) :: global_dens
  real(wp), allocatable, dimension(:,:) :: global_uwnd
  real(wp), allocatable, dimension(:,:) :: global_wwnd
  real(wp), allocatable, dimension(:,:) :: global_theta

  ! Arrays for MPI_Gatherv
  integer, allocatable, dimension(:) :: recvcounts
  integer, allocatable, dimension(:) :: displs

  integer :: ncid
  integer :: dens_varid, uwnd_varid, wwnd_varid, theta_varid, t_varid
  integer :: rec_out

  contains

  subroutine create_output
    implicit none
    integer :: t_dimid, x_dimid, z_dimid
    integer :: i, ierr, base_nx, remainder

    ! Allocate local arrays
    allocate(dens(nx,nz))
    allocate(uwnd(nx,nz))
    allocate(wwnd(nx,nz))
    allocate(theta(nx,nz))

    ! Setup for MPI_Gatherv
    allocate(recvcounts(nprocs))
    allocate(displs(nprocs))

    ! Calculate receive counts and displacements for each process
    base_nx = nx_global / nprocs
    remainder = mod(nx_global, nprocs)
    
    displs(1) = 0
    do i = 1, nprocs
      if (i-1 < remainder) then
        recvcounts(i) = (base_nx + 1) * nz
      else
        recvcounts(i) = base_nx * nz
      end if
      if (i > 1) then
        displs(i) = displs(i-1) + recvcounts(i-1)
      end if
    end do

    ! Only rank 0 creates file and allocates global arrays
    if (myrank == 0) then
      allocate(global_dens(nx_global,nz))
      allocate(global_uwnd(nx_global,nz))
      allocate(global_wwnd(nx_global,nz))
      allocate(global_theta(nx_global,nz))

      call ncwrap(nf90_create('output.nc',nf90_clobber,ncid), __LINE__)
      call ncwrap(nf90_def_dim(ncid,'time',nf90_unlimited,t_dimid), __LINE__)
      call ncwrap(nf90_def_dim(ncid,'x',nx_global,x_dimid), __LINE__)
      call ncwrap(nf90_def_dim(ncid,'z',nz,z_dimid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'time',iowp,[t_dimid],t_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'rho',iowp, &
          [x_dimid,z_dimid,t_dimid],dens_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'u',iowp, &
          [x_dimid,z_dimid,t_dimid],uwnd_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'w',iowp, &
          [x_dimid,z_dimid,t_dimid],wwnd_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'theta',iowp, &
          [x_dimid,z_dimid,t_dimid],theta_varid), __LINE__)
      call ncwrap(nf90_enddef(ncid), __LINE__)
    end if

    rec_out = 1
  end subroutine create_output

  subroutine write_record(atmostat,ref,etime)
    implicit none
    type(atmospheric_state), intent(in) :: atmostat
    type(reference_state), intent(in) :: ref
    real(wp), intent(in) :: etime
    integer :: i, k, ierr
    integer, dimension(1) :: st1, ct1
    integer, dimension(3) :: st3, ct3
    real(wp), dimension(1) :: etimearr

    ! Calculate local fields
    do k = 1, nz
      do i = 1, nx
        dens(i,k) = atmostat%dens(i,k)
        uwnd(i,k) = atmostat%umom(i,k)/(ref%density(k)+dens(i,k))
        wwnd(i,k) = atmostat%wmom(i,k)/(ref%density(k)+dens(i,k))
        theta(i,k) = (atmostat%rhot(i,k) + ref%denstheta(k)) / &
            (ref%density(k) + dens(i,k)) - ref%denstheta(k)/ref%density(k)
      end do
    end do

    ! Gather all data to rank 0
    call MPI_Gatherv(dens, nx*nz, MPI_DOUBLE_PRECISION, &
                     global_dens, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                     0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(uwnd, nx*nz, MPI_DOUBLE_PRECISION, &
                     global_uwnd, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                     0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(wwnd, nx*nz, MPI_DOUBLE_PRECISION, &
                     global_wwnd, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                     0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(theta, nx*nz, MPI_DOUBLE_PRECISION, &
                     global_theta, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                     0, MPI_COMM_WORLD, ierr)

    ! Only rank 0 writes to file
    if (myrank == 0) then
      st3 = [ 1, 1, rec_out ]
      ct3 = [ nx_global, nz, 1 ]
      call ncwrap(nf90_put_var(ncid,dens_varid,global_dens,st3,ct3), __LINE__)
      call ncwrap(nf90_put_var(ncid,uwnd_varid,global_uwnd,st3,ct3), __LINE__)
      call ncwrap(nf90_put_var(ncid,wwnd_varid,global_wwnd,st3,ct3), __LINE__)
      call ncwrap(nf90_put_var(ncid,theta_varid,global_theta,st3,ct3), __LINE__)

      st1 = [ rec_out ]
      ct1 = [ 1 ]
      etimearr(1) = etime
      call ncwrap(nf90_put_var(ncid,t_varid,etimearr,st1,ct1), __LINE__)
    end if

    rec_out = rec_out + 1
  end subroutine write_record

  subroutine close_output
    implicit none
    if ( allocated(dens) ) then
      deallocate(dens)
      deallocate(uwnd)
      deallocate(wwnd)
      deallocate(theta)
    end if
    if ( allocated(recvcounts) ) then
      deallocate(recvcounts)
      deallocate(displs)
    end if
    if (myrank == 0) then
      if ( allocated(global_dens) ) then
        deallocate(global_dens)
        deallocate(global_uwnd)
        deallocate(global_wwnd)
        deallocate(global_theta)
      end if
      call ncwrap(nf90_close(ncid), __LINE__)
    end if
  end subroutine close_output

  subroutine ncwrap(ierr,line)
    implicit none
    integer, intent(in) :: ierr
    integer, intent(in) :: line
    if (ierr /= nf90_noerr) then
      write(stderr,*) 'NetCDF Error at line: ', line
      write(stderr,*) nf90_strerror(ierr)
      stop
    end if
  end subroutine ncwrap

end module module_output
