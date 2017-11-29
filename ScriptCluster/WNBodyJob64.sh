#!/bin/bash


#SBATCH --ntasks 64
#SBATCH --cpus-per-task 1
#SBATCH --time=4:00:00


module purge
module load gcc
module load mvapich2



for i in 6400 12800
do
srun ../NBodyNew $i 1 0.01 25 3 1000
done

