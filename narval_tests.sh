# Interactive session
salloc --ntasks=1 --mem=4000M --time=00:05:00 --account=ctb-tpoisot

# Copy the scripts to the cluster
scp ./03_generate_sdms.jl narval:projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
scp ./03_generate_sdms1.jl narval:projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
scp ./03_generate_sdms3.jl narval:projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/
scp narval_test1.sh narval:
scp narval_test2.sh narval:
scp narval_test3.sh narval:

# Run the job
sbatch narval_test1.sh
sbatch narval_test2.sh
sbatch narval_test3.sh
sbatch narval_test4.sh

# Check status
sq -u $USER

# Check job efficiency
seff 10807379 # test1, 4GB, 1 CPU, 30 min => OOM kill
seff 10817420 # test1, 8GB, 1 CPU, 30 min => 5GB used, timeout after 80 species
seff 10818154 # test1, 8GB, 1 CPU, 1h => completed, 58 min, 5 GB used
seff 10818453 # test2, 4GB, 16 CPU, 1h => timeout, 4 nodes, 4 nodes/core, 6.09% CPU (so not threaded?), 8.4% MEM (5GB/62GB)
seff 10843641
seff 10818463 # test3, 4GB, 16 CPU, 1h, threads => OOM kill (??), 3 nodes, 4 nodes/core,  6.23% of RAM (3.9 GB/62 GB)
seff 10843604 # test3, 8GB, 32 CPU, 15 min, threads ==> OOM kill, 8 nodes, 4 cores/node, 10 min CPU (3%), 7.79 GB of RAM (3.12% of 250.0)
seff 10843722 # test4, 1 node, 32 CPU, 8GB, 15 min
seff 10843978
seff 10844590 # test4, took 16 min

# Check the output files
cat test1-10807379.out
cat test1-10817420.out
cat test1-10818154.out
cat test1-10818453.out
cat test1-10818463.out
cat test4-10844590.out

# Check the results
ls projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/data/sdms/
rm projects/def-tpoisot/2022-SpatialProbabilisticMetaweb/data/sdms/*.tif