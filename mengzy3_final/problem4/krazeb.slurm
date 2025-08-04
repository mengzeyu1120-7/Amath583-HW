#!/bin/bash
#SBATCH --job-name=krazeb_test
#SBATCH --account=amath
#SBATCH --partition=gpu-rtx6k
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:45:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

module load gcc
module load intel/oneAPI/2023.2.1

source cmpl.krazeb

mpirun -np 1  ./xkrazeb 1 1 48 48 4096 >> krazeb_results.txt 2>&1
mpirun -np 2  ./xkrazeb 1 2 48 48 4096 >> krazeb_results.txt 2>&1
mpirun -np 4  ./xkrazeb 2 2 48 48 4096 >> krazeb_results.txt 2>&1
mpirun -np 8  ./xkrazeb 2 4 48 48 4096 >> krazeb_results.txt 2>&1
mpirun -np 16 ./xkrazeb 4 4 48 48 4096 >> krazeb_results.txt 2>&1
