# Packages required for ./03_generate_sdms.jl
using SimpleSDMLayers
using EvoTrees
using StatsBase
using StatsPlots
using ProgressMeter
using DataFrames
import CSV

# Install the same version as in the Manifest for these packages only
# To minimize the number of packages installed on the clusters for now
using Pkg
Pkg.add(name="SimpleSDMLayers", version="0.8.3")
Pkg.add(name="EvoTrees", version="0.10.0")
Pkg.add(name="StatsBase", version="0.33.21")
# Pkg.add(name="StatsPlots", version="0.14.33") # maybe not for now
Pkg.add(name="ProgressMeter", version="1.7.2")
Pkg.add(name="DataFrames", version="1.3.4")
Pkg.add(name="CSV", version="0.10.4")