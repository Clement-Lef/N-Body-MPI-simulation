#!/bin/bash


#SBATCH --ntasks 32
#SBATCH --cpus-per-task 1
#SBATCH --time=6:00:00


module purge
module load gcc
module load mvapich2



for i in 3200 6400
do
srun ../NBodyNew $i 1 0.01 25 3 1000
done

