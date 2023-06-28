#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --array=1-3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --job-name=02_get_absences
#SBATCH --output=jobs/job_02_get_absences-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
julia --project --threads=63 -e 'JOBARRAY = true; CAN = true; WITHIN_RADIUS = true; quiet = true; include("02_get_absences.jl")'
