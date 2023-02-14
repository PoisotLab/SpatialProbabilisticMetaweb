#!/bin/bash
#SBATCH --account=ctb-tpoisot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=130G
#SBATCH --time=00:25:00
#SBATCH --job-name=06_get_species_lcbd
#SBATCH --output=jobs/job_06_get_species_lcbd-%J.out

module load StdEnv/2020
module load julia/1.8.1

cd $HOME/scratch/2022-SpatialProbabilisticMetaweb
julia --project --threads=63 -e 'CAN = true; include("06_get_species_lcbd.jl")'
