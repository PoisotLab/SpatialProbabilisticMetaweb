#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=00:60:00
#SBATCH --job-name=test1
#SBATCH --output=test1-%J.out
module load StdEnv/2020 julia/1.8.1
julia projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/03_generate_sdms1.jl