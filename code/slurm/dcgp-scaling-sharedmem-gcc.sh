#!/bin/bash

#SBATCH --job-name=2dEulerEquations    # Job name
#SBATCH --exclusive                 # grab the whole node
#SBATCH --hint=multithread
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1          # max 112
#SBATCH --mem=0                     # all memory on the node
#SBATCH --mem-bind=local
#SBATCH --distribution=block:block
#SBATCH --time=00:30:00             # Time limit hrs:min:sec
#SBATCH --account=ICT25_MHPC
#SBATCH --partition=dcgp_usr_prod
##SBATCH --qos=dcgp_qos_dbg
#SBATCH --error=slurm-%x-%j.err
#SBATCH --output=slurm-%x-%j.out

# prepare
PROJ_DIR=$SCRATCH/code/serial
OUTPUT_DIR=$SCRATCH/code/output
# BUILD_DIR=$SCRATCH/build
# SRC_DIR=$PROJ_DIR/src
# INCLUDE_DIR=$PROJ_DIR/include
TARGET=model

# load modules
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22

# compile
# mpiicpx -std=c++17 -O3 -march=icelake-server -mtune=icelake-server \
#     -DNDEBUG -funroll-loops -ffast-math -fp-model=fast -vec-threshold \
#     -fno-alias -qopt-mem-layout-trans=3 -qopt-zmm-usage=high -ipo -qopenmp \
#     $SRC_DIR/$TARGET.cpp -o $BUILD_DIR/$TARGET -I$INCLUDE_DIR \
#     -lmkl_rt -liomp5

# ---- PINNING & THREAD ENV ----------------------------------------------------
# Pin more aggressively:
# We would therefore like consecutive threads to be bound close together, as is done with KMP_AFFINITY=compact,
# so that communication overhead, cache line invalidation overhead, and page thrashing are minimized.
export OMP_PLACES=cores           # no need for this if u specify KMP_AFFINITY
export OMP_PROC_BIND=close        # no need for this if u specify KMP_AFFINITY
export OMP_DYNAMIC=FALSE
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} # 14 threads per core/rank

OUTPUT=time_${SLURM_NNODES}_${SLURM_JOB_ID}.txt

echo "--- Building Model ---"
cd $PROJ_DIR
make clean
make -j 4

echo "--- Running Simulation ---"
echo "Nodes = ${SLURM_NNODES}"
echo "Tasks / Node = ${SLURM_NTASKS}"
echo "CPUs / Task = ${SLURM_CPUS_PER_TASK} (OMP Threads)"

srun --cpu-bind=cores $PROJ_DIR/$TARGET 896 1000 0

echo "--- Ending Simulation ---"