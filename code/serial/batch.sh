#!/bin/bash
#SBATCH -A ICT25_MHPC
#SBATCH -p dcgp_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2             
#SBATCH --cpus-per-task=56             
#SBATCH --time=00:05:00
#SBATCH --job-name=openmp_only_node
#SBATCH --output=logs/%x_%j.out

# 1. Load Modules
module purge
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22 

# 2. Setup Environment Variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

NODES=${SLURM_NNODES:-1}

echo "--- Building Model ---"
make clean || exit 1
make -j 4 USE_OPENACC=0 || exit 1


echo "--- Running Simulation ---"
echo "Nodes = ${NODES}"
echo "Tasks / Node = 1 (Pure OpenMP)"
echo "CPUs / Task = ${SLURM_CPUS_PER_TASK} (OMP Threads)"

# 3. Execution 

mpirun -np 2 ./model

#4. END ba ba ba

echo "--- Ending Simulation ---"
