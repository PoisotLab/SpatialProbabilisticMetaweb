#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=248G
#SBATCH --time=08:00:00
#SBATCH --job-name=13_get_motifs
#SBATCH --output=jobs/out/%x-%J.out

module load StdEnv/2020
module load julia/1.9.1

cd $HOME/projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
julia --project --threads=63 -e 'quiet = true; include("09_get_network_measures.jl")'