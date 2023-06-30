#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=140G
#SBATCH --time=02:00:00
#SBATCH --job-name=05_assemble_networks.sh
#SBATCH --output=jobs/job_05_assemble_networks-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
julia --project --threads=63 -e 'CAN = true; quiet = true; include("05_assemble_networks.jl")'
