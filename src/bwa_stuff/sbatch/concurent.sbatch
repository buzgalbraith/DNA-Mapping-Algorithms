#!/bin/bash 
#SBATCH --nodes=3                        # requests 5 compute servers
#SBATCH --ntasks-per-node=2              # runs 2 tasks on each server
#SBATCH --cpus-per-task=1                # uses 1 compute core per task
#SBATCH --time=1:00:00 
#SBATCH --mem=2GB
#SBATCH --job-name=concurent_optimze
#SBATCH --output=concurent_optimze.out

module purge
module load python/intel/3.8.6
cd ..
python ./optmized_concurent_bwa.py 2 3 5 10 20 30 40 50 100 150 200 250
