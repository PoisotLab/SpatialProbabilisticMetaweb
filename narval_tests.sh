salloc --ntasks=1 --mem=4000M --time=00:05:00 --account=ctb-tpoisot

scp ./03_generate_sdms1.jl narval:projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
scp ./narval_test1.sh narval:projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/

sbatch projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/narval_test1.sh # 10807379

sq -u $USER