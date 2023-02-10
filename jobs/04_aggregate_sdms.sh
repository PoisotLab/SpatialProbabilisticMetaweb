#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=48G
#SBATCH --time=00:17:00
#SBATCH --job-name=04_aggregate_sdms
#SBATCH --output=jobs/job_04_aggregate_sdms-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/scratch/2022-SpatialProbabilisticMetaweb
julia --project --threads=63 -e 'CAN = true; include("04_aggregate_sdms.jl")'
