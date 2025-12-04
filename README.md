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

## Running the program
### Option A
Make sure your system has the following packages installed: `cmake`, `graphviz`, `doxygen`, `gfortran`, `libmpich-dev`, `libnetcdff-dev`, `python3`, `python3-pip`, `build-essential` and should be able to be found by CMake. These are the package names in Ubuntu 22.02. You may install CUDA device drivers so you offload computation to GPU, if available. If you do this, skip Option B and continue below.

### Option B
If you want to quickly install all of packages in an isolated environment, you can simply run `runenv.sh` and continue below. This assumes you have Docker engine installed on your machine

### Build and run
When you have all the packages installed:
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

### Sample output
Parallel, with compiler optimization flags on Ubuntu 22.02 with GNU compiler
```bash
root@91d7785e2597:/app/build# OMP_NUM_THREADS=2 mpirun -n 4 ./model 100 1000 0
 ================= Execution Info ==================
  Number of MPI tasks:    4
  Number of OpenMP threads:    2
   OpenMP: ENABLED
   OpenACC: DISABLED
 ===================================================
 --------------- Domain Decomposition --------------
  Global nx:    100
  Local nx per process: ~    25
 ---------------------------------------------------
 SIMPLE ATMOSPHERIC MODEL STARTING.
Wall Clock Start: 12:06:53.808
 INITIALIZING MODEL STATUS.
 nx_global  :          100
 nx_local   :           25
 nz         :           50
 dx         :    200.00000000000000     
 dz         :    200.00000000000000     
 dt         :   0.66666666666666663     
 final time :    1000.0000000000000     
 MODEL STATUS INITIALIZED.
 TIME PERCENT :  0%
 TIME PERCENT : 10%
 TIME PERCENT : 20%
 TIME PERCENT : 30%
 TIME PERCENT : 40%
 TIME PERCENT : 50%
 TIME PERCENT : 60%
 TIME PERCENT : 70%
 TIME PERCENT : 80%
 TIME PERCENT : 90%
 ----------------- Atmosphere check ----------------
 Fractional Delta Mass  :   -6.0578657789430107E-015
 Fractional Delta Energy:    1.0051785904252513E-004
 ---------------------------------------------------
 ----------------------------------------------------------------------------------------------------
 PARALLEL TIMING STATISTICS (Microseconds)
                      Function      Max Total       Max Excl      Avg Total     Calls
 ----------------------------------------------------------------------------------------------------
                          INIT        3778.42         685.63        3758.13         4
          Computation: thermal        3118.97        1623.71        3089.61     57228
Computation: hydrostatic_const        1587.39        1587.39        1555.97     57228
Computation: total_mass_energy         281.50         104.25         217.03         8
            MPI: Communication       73227.01       73227.01       72265.45     18020
       Computation: rungekutta      415414.62        1440.36      415298.07      6004
             Computation: step      414268.06      343668.06      414068.18     36024
 ----------------------------------------------------------------------------------------------------
 * Max Total/Excl: The slowest rank for that function.
 * Avg Total: Average wall time across all ranks.
 SIMPLE ATMOSPHERIC MODEL RUN COMPLETED.
Wall Clock End:   12:06:54.228
USED CPU TIME:              0.416140 seconds
```

Parallel, with compiler optimization flags and GPU offloading on Leonardo
```bash
[jrayo000@lrdn0322 build]$ OMP_NUM_THREADS=1 mpirun -n 4 ./model 100 1000 10
 ================= Execution Info ==================
  Number of MPI tasks:    4
   OpenMP: DISABLED
  Number of GPUs available:    4
   OpenACC: ENABLED
 ===================================================
  MPI Rank    0 -> GPU  0
  MPI Rank    1 -> GPU  1
  MPI Rank    3 -> GPU  3
  MPI Rank    2 -> GPU  2
 --------------- Domain Decomposition --------------
  Global nx:    100
  Local nx per process: ~    25
 ---------------------------------------------------
 SIMPLE ATMOSPHERIC MODEL STARTING.
Wall Clock Start: 13:16:30.984
 INITIALIZING MODEL STATUS.
 nx_global  :           100
 nx_local   :            25
 nz         :            50
 dx         :     200.0000000000000     
 dz         :     200.0000000000000     
 dt         :    0.6666666666666666     
 final time :     1000.000000000000     
 MODEL STATUS INITIALIZED.
 TIME PERCENT :  0%
 TIME PERCENT : 10%
 TIME PERCENT : 20%
 TIME PERCENT : 30%
 TIME PERCENT : 40%
 TIME PERCENT : 50%
 TIME PERCENT : 60%
 TIME PERCENT : 70%
 TIME PERCENT : 80%
 TIME PERCENT : 90%
 ----------------- Atmosphere check ----------------
 Fractional Delta Mass  :   -1.9541502512719523E-016
 Fractional Delta Energy:    1.0051785903963127E-004
 ---------------------------------------------------
 ----------------------------------------------------------------------------------------------------
 PARALLEL TIMING STATISTICS (Microseconds)
                      Function      Max Total       Max Excl      Avg Total     Calls
 ----------------------------------------------------------------------------------------------------
                          INIT        7003.29         709.00        6988.19         4
          Computation: thermal        6304.66        2306.16        6291.93     57228
Computation: hydrostatic_const        4019.37        4019.37        4004.93     57228
Computation: total_mass_energy         764.49         150.51         752.23         8
            MPI: Communication       23811.84       23811.84       22877.72     18020
       Computation: rungekutta      764785.78        1006.35      764769.74      6004
             Computation: step      763792.72      742497.44      763774.35     36024
 ----------------------------------------------------------------------------------------------------
 * Max Total/Excl: The slowest rank for that function.
 * Avg Total: Average wall time across all ranks.
 SIMPLE ATMOSPHERIC MODEL RUN COMPLETED.
Wall Clock End:   13:16:31.759
USED CPU TIME:              0.766376 seconds
```

Serial, with compiler optimization flags on Ubuntu 22.02 with GNU compiler
```bash
root@3f0e2852e6ec:/app# OMP_NUM_THREADS=1 ./model 100 1000 10
 SIMPLE ATMOSPHERIC MODEL STARTING.
 INITIALIZING MODEL STATUS.
 nx         :          100
 nz         :           50
 dx         :    200.00000000000000     
 dz         :    200.00000000000000     
 dt         :   0.66666666666666663     
 final time :    1000.0000000000000     
 MODEL STATUS INITIALIZED.
 TIME PERCENT :  0%
 TIME PERCENT : 10%
 TIME PERCENT : 20%
 TIME PERCENT : 30%
 TIME PERCENT : 40%
 TIME PERCENT : 50%
 TIME PERCENT : 60%
 TIME PERCENT : 70%
 TIME PERCENT : 80%
 TIME PERCENT : 90%
 ----------------- Atmosphere check ----------------
 Fractional Delta Mass  :    2.7748933568062396E-014
 Fractional Delta Energy:    1.0051785903549762E-004
 ---------------------------------------------------
 SIMPLE ATMOSPHERIC MODEL RUN COMPLETED.
 USED CPU TIME:   0.91237400000000002
```

Serial, with compiler optimization flags on Leonardo
```bash
[jrayo000@lrdn2570 serial]$ OMP_NUM_THREADS=1 ./model 500 1000 10
 SIMPLE ATMOSPHERIC MODEL STARTING.
 INITIALIZING MODEL STATUS.
 nx         :          500
 nz         :          250
 dx         :    40.000000000000000     
 dz         :    40.000000000000000     
 dt         :   0.13333333333333333     
 final time :    1000.0000000000000     
 MODEL STATUS INITIALIZED.
 TIME PERCENT :  0%
 TIME PERCENT : 10%
 TIME PERCENT : 20%
 TIME PERCENT : 30%
 TIME PERCENT : 40%
 TIME PERCENT : 50%
 TIME PERCENT : 60%
 TIME PERCENT : 70%
 TIME PERCENT : 80%
 TIME PERCENT : 90%
 ----------------- Atmosphere check ----------------
 Fractional Delta Mass  :   -1.9541502512719538E-016
 Fractional Delta Energy:    1.2056948351433392E-004
 ---------------------------------------------------
 SIMPLE ATMOSPHERIC MODEL RUN COMPLETED.
 USED CPU TIME:    343.14903285899999
```

## Running on Leonardo
You can see the slurm scripts to see how the program was built and run on the cluster:
- `cpu_model.sh`, `gpu_model.sh`
- `serial/batch.sh`, `serial/gpu_batch.sh`, `serial/gpu.sh`

Sample command to upload your files to the cluster. This assumes that you have SSH key acquired from `step`:
```bash
rsync -arvzP code leo:/leonardo_scratch/large/userexternal/jrayo000
```

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
