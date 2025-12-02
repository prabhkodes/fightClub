#!/bin/bash
#SBATCH --job-name=2dEulerEquations
#SBATCH --account=ICT25_MHPC_0
#SBATCH --partition=boost_usr_prod    # Booster (GPU)
#SBATCH --qos=boost_qos_dbg           # Cola de debug (máx ~30 min)
#SBATCH --time=00:10:00               # 10 minutos para probar
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4           # 4 MPI ranks
#SBATCH --cpus-per-task=8             # 8 hilos CPU por rank
#SBATCH --gres=gpu:1                  # 1 GPU
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module purge
module load nvhpc/24.5
module load hpcx-mpi/2.19
module load netcdf-fortran/4.6.1--hpcx-mpi--2.19--nvhpc--24.5
module load binutils/2.42

echo "--- Building Model ---"
make clean || exit 1
make -j 4 || exit 1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "--- Running Simulation ---"
echo "Nodes = ${NODES}"
echo "Tasks / Node = 1 (Pure OpenMP)"
echo "CPUs / Task = ${SLURM_CPUS_PER_TASK} (OMP Threads)"

./model

echo "--- Ending Simulation ---"
