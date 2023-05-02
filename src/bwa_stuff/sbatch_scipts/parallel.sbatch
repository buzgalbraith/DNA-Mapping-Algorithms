#!/bin/bash
#SBATCH --job-name=paralell
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=2GB
#SBATCH --time=01:00:00
#SBATCH --output=parallel.out

module purge
module load python/intel/3.8.6
module load openmpi/intel/4.0.5


cd ..
cd ./optimized_parallel
python ./bwa_optimization_call.py 5 2 3 5 10 20 30 40 50 100 150 200 250