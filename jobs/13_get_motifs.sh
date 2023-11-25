#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --array=1-500
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=248G
#SBATCH --time=02:00:00
#SBATCH --job-name=13_get_S4
#SBATCH --output=jobs/out/%x-%J.out

module load StdEnv/2020
module load julia/1.9.1

cd $HOME/projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
julia --project --threads=63 -e 'MOTIF = :S4; TOTAL_JOBS = 500; JOBARRAY = true; CAN = true; quiet = true; include("13_get_motifs.jl")'
