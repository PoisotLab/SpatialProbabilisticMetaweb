import Pkg; Pkg.activate(".")

import CSV
using Chain
using DataFrames
using DataFramesMeta
using Distributions
using EcologicalNetworks
using ProgressMeter
using SimpleSDMLayers
using Random
using Statistics
using StatsPlots
using Base.Threads: @threads

default(; dpi=200)

include("A1_LCBD.jl")
include("A2_colors.jl")