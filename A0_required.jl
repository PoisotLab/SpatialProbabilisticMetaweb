import Pkg; Pkg.activate(".")

import CSV
using Chain
using DataFrames
using DataFramesMeta
using Distributions
using EcologicalNetworks
using EvoTrees
using ProgressMeter
using Random
using Shapefile
using SimpleSDMLayers
using SparseArrays
using Statistics
using StatsBase
using StatsPlots
using Base.Threads: @threads
using JLD2, CodecZlib
import ColorBlendModes
using Plots.PlotMeasures

default(; dpi=200)

include("A1_LCBD.jl")
include("A2_colors.jl")
include("A3_bivariate_values.jl")
include("A4_shapefile.jl")