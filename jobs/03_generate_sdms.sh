#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=00:30:00
#SBATCH --job-name=03_generate_sdms
#SBATCH --output=jobs/job_03_generate_sdms-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/scratch/2022-SpatialProbabilisticMetaweb
julia --project --threads=31 03_generate_sdms.jl
