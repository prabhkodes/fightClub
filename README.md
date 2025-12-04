# P1.8 - Best Practices in Scientific Software Development
## 2-D Euler Equations
First rule of fight club is we don't talk about fight club.

## Directory structure
- `code` contains the source code, scripts to run scaling tests in Leonardo
  - `serial` contains the main source code
    - `build` contains the generated files and executable when building the CMake project
      - `doc/html` contains the generated documentation
  - `results` contains results and plots
- `perf_results` contains performance counter stats for running the program
- `.github/workflows/ci.yaml` contains GitHub actions for automated integration testing
- `requirements.txt` contains Python packages to run the output file comparison tests
- `Dockerfile` is used to build the image and run the container to build and run the program
- `runenv.sh` contains the commands to initialize the containerized environment

## Source code structure
- `model.F90` contains the main routine
- `module_output.F90` contains the routines relevant for output file generation in parallel
- `module_parameters.F90` contains parameters for domain decomposition and solvers, as well as physical constants
- `module_physics.F90` takes care of initial and boundary conditions, numerical solution and calculation of mass and energy budgets
- `module_types.F90` calculates the atmospheric state type, flux and tendency, as well as handle halo exchange
- `parallel_timer.F90` is used to time the execution of several routines among all ranks and generates a summary on total time, max time, average time and number of calls.

## How to build and run
Make sure your system has the following packages installed: `cmake`, `graphviz`, `doxygen`, `gfortran`, `libmpich-dev`, `libnetcdff-dev`, `python3`, `python3-pip`, `build-essential`, CUDA device drivers (optional, to offload computation to GPU) and should be able to be found by CMake.

1. Change working directory
2. Build the CMake project
3. Test the program if it runs properly
4. Generate the documentation (optional)

```bash
cd fightClub/serial
mkdir build
cd build
cmake .. -DUSE_OPENACC=ON # if compiling on a system with MPI+GPU. otherwise don't include this flag and it will only compile with MPI+OpenMP flags
make -j 4
make test
make doc
```

If you want to manually run the executable file, you can do the following:

```bash
cd serial/build
mpirun -n 4 ./model 100 1000 10 # No of grid points in x-axis, number of timesteps, output frequency
```

## Running on Leonardo
You can see the slurm scripts to see how the program was built and run on the cluster:
- `cpu_model.sh`, `gpu_model.sh`
- `serial/batch.sh`, `serial/gpu_batch.sh`, `serial/gpu.sh`

### Project Path
```
/leonardo/pub/userexternal/jgordill/fightClub
```
### Modules
#### MPI+OpenMP
```
module purge

# Compiler
module load gcc/12.2.0

# MPI
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2

# NetCDF Fortran
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22
```

### MPI+OpenACC
```
module purge

# Compiler
module load nvhpc/24.5

# MPI
module load hpcx-mpi/2.19

# NetCDF Fortran
module load netcdf-fortran/4.6.1--hpcx-mpi--2.19--nvhpc--24.5
module load binutils/2.42
```
