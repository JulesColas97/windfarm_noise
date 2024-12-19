#!/bin/bash
#SBATCH --job-name=PE0
#SBATCH --output=out.out # output messages go here
#SBATCH --error=err.err    # error messages go here
#SBATCH --mail-user=jules.colas@ecl17.ec-lyon.fr
#SBATCH --mail-type=ALL
#SBATCH --partition=test # partition name
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH --time=01:00:00
module purge
module load  HDF5/1.10.1-intel-2018a
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_AFFINITY=granularity=fine,compact,1,0
export OMP_STACKSIZE=1g
ulimit -s unlimited
cd t0/50/ 
time ./PE_2D >out.out & 
cd ../../ 
cd t0/90/ 
time ./PE_2D >out.out & 
cd ../../ 
cd t0/130/ 
time ./PE_2D >out.out & 
cd ../../ 
cd t1/50/ 
time ./PE_2D >out.out & 
cd ../../ 
cd t1/90/ 
time ./PE_2D >out.out & 
cd ../../ 
cd t1/130/ 
time ./PE_2D >out.out & 
cd ../../ 
wait
