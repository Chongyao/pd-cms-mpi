#!/bin/bash

# --- Slurm Resource Configuration ---
# SBATCH --job-name=mp_cms 
# SBATCH --partition=amd_512
# SBATCH --nodes=9               # Number of nodes requested
# SBATCH --ntasks-per-node=1       # x MPI ranks per node,  <!!!notice: ntasks-per-node*nodes = number of ranks!!!> , If memory is sufficient, the more individual node processes there are, the better.
# SBATCH --cpus-per-task=128        # x OpenMP threads per MPI rank (ntasks-per-node*cpus-per-task = 64 cores utilized per node)
# SBATCH --time=4:00:00           # Maximum running time (12 hours)
# SBATCH --output=log_%j.out  
# SBATCH --error=log_%j.err   

# --- Environment Setup ---
module purge
source /public3/home/m6s001980/env.sh

# --- Data Dir ---
DATA_DIR="/public3/home/m6s001980/Project/CMS/Data/PCB_Coarse"


# --- Parallelism Configuration ---
# Thread counts must match the value set by --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK  
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# AMD Specific Optimization: Bind threads to cores to minimize overhead (NUMA optimization)
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Print run information for debugging
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks per Node: $SLURM_NTASKS_PER_NODE"
echo "Threads per Task: $OMP_NUM_THREADS"
echo "Total MPI Ranks: $SLURM_NTASKS"
echo "Start Time: $(date)"
echo "=========================================="


# No need to hardcode -np, mpirun/srun automatically detects the total ranks ($SLURM_NTASKS=ntasks-per-node*nodes)


mpirun -n 9 ./script/profile.sh ./build/bin/step3_mpi_loader ../zcy-cms-mpi-working-dir/PCB_Coarse 300 0 1e-8

echo "Job end at $(date)"
