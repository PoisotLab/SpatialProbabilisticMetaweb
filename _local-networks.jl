include("A0_required.jl")

# Run 05_assemble_networks.jl to define networks_thr and reference_layer

# Copy the previous network
networks = copy(networks_thr)

# Extract a single site
local_site = networks[1, :, :, :]

# Extract local interaction matrix
local_mat = dropdims(mean(local_site; dims=ndims(local_site)), dims=ndims(local_site))

# Transform into network
local_network = UnipartiteProbabilisticNetwork(local_mat, species(P))

# Transform into layer
local_layer = SimpleSDMResponse([fill(local_network, 2, 2); nothing nothing])
local_layer.grid

# Network operations
connectance(local_network)
# connectance.(local_layer) nope
connectance.(collect(local_layer)) # but not a layer

# Test layer operations
reference_layer.grid
(reference_layer * 10).grid
mean(reference_layer)
abs(reference_layer).grid
broadcast(abs, reference_layer).grid
broadcast(connectance, local_layer).grid # here we go

## Create a zero type

# What we did for Normal distribution
Base.zero(::Type{Normal{T}}) where T = Normal(zero(T), zero(T))
Normal{Float64} # works by default
zero(Normal{Float64}) # works because defined previously
Normal(0.0, 0.0) # this is what is defined by the zero

# What we need for UnipartiteProbabilisticNetwork
UnipartiteProbabilisticNetwork{Float64} # works by default
# zero(UnipartiteProbabilisticNetwork{Float64}) # need to be defined
UnipartiteProbabilisticNetwork(local_mat)
UnipartiteProbabilisticNetwork(zeros(eltype(local_mat), size(local_mat)))

# Attempt to create zero type
Base.zero(::Type{UnipartiteProbabilisticNetwork{T}}) where T = UnipartiteProbabilisticNetwork(zeros(T, (2,2)))
zero(UnipartiteProbabilisticNetwork{Float64})
Base.zero(::Type{UnipartiteProbabilisticNetwork{T, String}}) where T = zero(UnipartiteProbabilisticNetwork{T})
zero(UnipartiteProbabilisticNetwork{Float64, String})

# Why do I need the zero type??
# Do we need it in similar?
similar(local_layer)
similar(local_layer).grid
similar(local_layer).grid[1] # yeah ok

# What happens inside similar?
layer = local_layer
TC = eltype(layer)
emptygrid = convert(Matrix{Union{Nothing,TC}}, zeros(TC, size(layer))) # here's where we need the zero
emptygrid[findall(isnothing, layer.grid)] .= nothing
SimpleSDMResponse(emptygrid, layer.left, layer.right, layer.bottom, layer.top)

# Do we need elsewhere?
# convert?
# convert(UnipartiteProbabilisticNetwork{Float32, String}, local_layer)
# convert(SimpleSDMResponse{UnipartiteProbabilisticNetwork{Float32, String}}, local_layer)
# Doesn't work and we won't need it anyways

# similar (when changing the type)
similar(local_layer, Float32) # ok here it's needed

# broadcastable?
# connectance.(local_layer)
# abs.(reference_layer)
# Doesn't work, is this the half-baked broadcast solution?

# multiple layer operations (min, max, +, etc.)

# mask (requires similar)
# also sliding, radius, etc.

# Delete the zeros for now
Base.delete_method(@which zero(UnipartiteProbabilisticNetwork{Float64, String}))
Base.delete_method(@which zero(UnipartiteProbabilisticNetwork{Float64}))