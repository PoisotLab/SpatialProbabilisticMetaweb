#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=01:30:00
#SBATCH --job-name=03_generate_sdms
#SBATCH --output=jobs/job_03_generate_sdms-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/scratch/2022-SpatialProbabilisticMetaweb
julia --project --threads=63 -e 'CAN = true; quiet = true; include("03_generate_sdms.jl")'
