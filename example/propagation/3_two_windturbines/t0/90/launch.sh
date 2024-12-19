#!/bin/bash
#SBATCH --job-name=c0t0h90
#SBATCH --output=out.out # output messages go here
#SBATCH --error=err.err    # error messages go here
#SBATCH --mail-user=jules.colas@ecl17.ec-lyon.fr
#SBATCH --mail-type=ALL
#SBATCH --partition=haswell # partition name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
module purge
module load  HDF5/1.10.1-intel-2018a
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_AFFINITY=granularity=fine,compact,1,0
export OMP_STACKSIZE=1g
ulimit -s unlimited
time ./PE_2D
