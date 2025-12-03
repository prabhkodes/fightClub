# P1.8 - Best Practices in Scientific Software Development
## 2-D Euler Equations
First rule of fight club is we don't talk about fight club.

## Leonardo
### Project Path
```
/leonardo/pub/userexternal/jgordill/fightClub
```
### Modules
#### Open MP
```
module purge

# Compiler
module load gcc/12.2.0

# MPI
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2

# NetCDF Fortran compile with MPI + GCC
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22
```

### Open ACC
```
module purge

module load nvhpc/24.5
module load hpcx-mpi/2.19

module load netcdf-fortran/4.6.1--hpcx-mpi--2.19--nvhpc--24.5

module load binutils/2.42
```
