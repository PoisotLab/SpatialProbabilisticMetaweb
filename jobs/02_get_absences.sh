#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --array=1-3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=248G
#SBATCH --time=04:00:00
#SBATCH --job-name=02_get_absences
#SBATCH --output=jobs/job_02_get_absences-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/scratch/2022-SpatialProbabilisticMetaweb
julia --project --threads=63 -e 'JOBARRAY = true; CAN = true; include("02_get_absences.jl")'
