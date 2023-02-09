# Command to time & evaluate memory used
/usr/bin/time -v ls |& grep -E "Maximum resident|wall clock"

# Test on the package loading script first
/usr/bin/time -v julia --project A0_required.jl |& grep -E "Maximum resident|wall clock"

# Now the script with 9 regions
/usr/bin/time -v julia --project ./test_09_get_network_measures_9-subregions.jl |& grep -E "Maximum resident|wall clock"
# Elapsed (wall clock) time (h:mm:ss or m:ss): 4:11.83
# Maximum resident set size (kbytes): 13107472 => 13,107472 GB

# With 1 region
/usr/bin/time -v julia --project ./test_09_get_network_measures_1-subregions.jl |& grep -E "Maximum resident|wall clock"
# Elapsed (wall clock) time (h:mm:ss or m:ss): 4:07.84
# Maximum resident set size (kbytes): 12170496 => 12,170496 GB

# Without the subregion division
/usr/bin/time -v julia --project ./test_09_get_network_measures_0-subregions.jl |& grep -E "Maximum resident|wall clock"
# Elapsed (wall clock) time (h:mm:ss or m:ss): 4:03.68
# Maximum resident set size (kbytes): 13255796 => 13,255796 GB