#!/bin/bash


#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time=12:00:00


module purge
module load gcc
module load mvapich2


for i in 500 1000 2000 4000 8000 16000 32000 128000
do
srun ../NBodyNew $i 1 0.01 25 3 1000
done

