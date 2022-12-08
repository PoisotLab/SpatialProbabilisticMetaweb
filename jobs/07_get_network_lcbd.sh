#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=02:30:00
#SBATCH --job-name=07_get_network_lcbd.jl
#SBATCH --output=jobs/job_07_get_network_lcbd-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/scratch/2022-SpatialProbabilisticMetaweb
julia --project --threads=63 07_get_network_lcbd.jl
