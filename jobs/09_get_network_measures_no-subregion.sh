#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --time=04:00:00
#SBATCH --job-name=09_get_network_measures_no-subregion.jl
#SBATCH --output=jobs/job_09_get_network_measures_no-subregion-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/scratch/2022-SpatialProbabilisticMetaweb
julia --project --threads=63 -e 'CAN = true; include("09_get_network_measures_no-subregion.jl")'
