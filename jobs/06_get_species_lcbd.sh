#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=160G
#SBATCH --time=01:00:00
#SBATCH --job-name=06_get_species_lcbd
#SBATCH --output=jobs/out/job_06_get_species_lcbd-%J.out

module load StdEnv/2020
module load julia/1.9.1

cd $HOME/projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
julia --project --threads=63 -e 'CAN = true; quiet = true; include("06_get_species_lcbd.jl")'
