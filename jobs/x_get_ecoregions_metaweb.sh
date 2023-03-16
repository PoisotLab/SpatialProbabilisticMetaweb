#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=450G
#SBATCH --time=04:30:00
#SBATCH --job-name=x_get_ecoregions_metaweb.jl
#SBATCH --output=jobs/job_x_get_ecoregions_metaweb-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/scratch/2022-SpatialProbabilisticMetaweb
julia --project --threads=63 -e 'CAN = true; quiet = true; include("x_get_ecoregions_measures.jl")'