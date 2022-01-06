import Pkg; Pkg.activate(".")

import CSV
using Chain
using DataFrames
using Distributions
using EcologicalNetworks
using ProgressMeter
using SimpleSDMLayers
using Random
using Statistics
using StatsPlots

default(; dpi=200)

include("A1_LCBD.jl")
include("A2_colors.jl")