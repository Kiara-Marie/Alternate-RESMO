#!/bin/bash
#SBATCH --time=06:00:00 
#SBATCH --job-name="resmoparallel3"
#SBATCH --mem=1000 
#SBATCH --account=def-edgrant
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks=32
module load matlab/2018a
matlab -nodisplay -singleCompThread -r "start_sim_plasma"
