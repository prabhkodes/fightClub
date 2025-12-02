#!/bin/bash
#SBATCH -A ICT25_MHPC
#SBATCH -p dcgp_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:05:00
#SBATCH --job-name=openmp_only_model
#SBATCH --output=logs/%x_%j.out


OMP_NUM_THREADS=112 


module purge
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22 


export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

NODES=${SLURM_NNODES:-1}



echo "--- Building Model ---"
make clean || exit 1
make -j 4 || exit 1
echo "--- Running Pure OpenMP Simulation ---"
echo "Threads = , Nodes = "

./model