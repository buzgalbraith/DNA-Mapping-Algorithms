#!/bin/bash 
#SBATCH --nodes=3                        # requests 5 compute servers
#SBATCH --ntasks-per-node=2              # runs 2 tasks on each server
#SBATCH --cpus-per-task=1                # uses 1 compute core per task
#SBATCH --time=1:00:00 
#SBATCH --mem=2GB
#SBATCH --job-name=naive
#SBATCH --output=naive.out

module purge
module load python/intel/3.8.6
cd ..
python ./naive_bwa.py 2 3 1 5 
