## Attempt to compute omnivory

# Test on single network
N = layer[1]

# Omnivory attempt
# try omnivory(N) catch err; println(err) end
try omnivory(N) catch err; showerror(stdout, err) end

# Let's look inside omnivory
# methods(omnivory)
# function omnivory(N::T) where {T <: UnipartiteNetwork}
function omnivory(N::T) where {T <: Union{UnipartiteNetwork, UnipartiteProbabilisticNetwork}}
    OI = Dict([s => 0.0 for s in species(N)])

    TL = fractional_trophic_level(N) # failing here
    k = degree(N; dims=1)

    for sp_i in species(N)

        # Species with no interaction have an omnivory index of 0
        k[sp_i] > 0 || continue

        # For every species, we set its initial omnivory to 0
        oi = 0.0
        for (j, sp_j) in enumerate(species(N))
            # Then for every species it consumes, we ha
            tl_diff = (TL[sp_j] - (TL[sp_i]-1.0)).^2.0
            corr = N[sp_i,sp_j]/k[sp_i]
            oi += tl_diff * corr
        end
        OI[sp_i] = oi
    end

    OI
    return OI
end

# How about fractional trophic level
try fractional_trophic_level(N) catch err; showerror(stdout, err) end
# Looking inside
# function fractional_trophic_level(N::T) where {T<:UnipartiteNetwork}
function fractional_trophic_level(N::T) where {T<:Union{UnipartiteNetwork, UnipartiteProbabilisticNetwork}}
    Y = nodiagonal(N)
    producers = keys(filter(spedeg -> spedeg.second == 0, degree(Y; dims=1)))
    sp = shortest_path(Y) # not working
    prod_id = findall(isequal(0), vec(sum(sp; dims=2)))
    return Dict(zip(species(Y; dims=1), maximum(sp[:,prod_id]; dims=2).+1))
end

# And how about shortest_path
try shortest_path(N) catch err; showerror(stdout, err) end
# function shortest_path(N::UnipartiteNetwork; nmax::Int64=50)
function shortest_path(N::T; nmax::Int64=50) where {T <: Union{UnipartiteNetwork, UnipartiteProbabilisticNetwork}}
    # We will have a matrix of the same size at the adjacency matrix
    D = EcologicalNetworks.spzeros(Int64, size(N.edges)...)
    D[findall(!iszero, N.edges)] .= 1
    # for i in 2:nmax
    for i in 2:50
      P = number_of_paths(N, n=i)
      D[(P .> 0).*(D .== 0)] .= i
    end
    return D
end

tmp_res = omnivory(N)
tmp_df = DataFrame(
    species = collect(keys(tmp_res)),
    omnivory = collect(values(tmp_res)),
    degree = collect(values(degree(N, dims=1)))
)
println(sort(tmp_df, [:omnivory, :degree], rev=true))

# Let's try with motifs instead
unipartitemotifs() # all motifs
unipartitemotifs().S2 # omnivory according to Stouffer et al. (2007)
find_motif(unipartitemotifs().S2, unipartitemotifs().S2)
find_motif(unipartitemotifs().S2, unipartitemotifs().S2) |> length # motif is there once
@time N_motifs = find_motif(N, unipartitemotifs().S2)
N_motifs
length(N_motifs)
expected_motif_count(N_motifs)
@time expected_motif_count(find_motif(N, unipartitemotifs().S2))[1]
@time [expected_motif_count(find_motif(N, unipartitemotifs().S2))[1] for N in layer.grid[1:4]]

# Motifs in a loop
test = zeros(4)
@time for (i, N) in enumerate(layer.grid[1:4])
    test[i] = expected_motif_count(find_motif(N, unipartitemotifs().S2))[1]
end
test

# Motifs in parallel
@time @threads for N in layer.grid[1:4]
    expected_motif_count(find_motif(N, unipartitemotifs().S2))[1]
end

# Motifs in parallel on layer
test = similar(layer, Float64)
@time @threads for k in keys(test)[1:4]
    test[k] = expected_motif_count(find_motif(layer[k], unipartitemotifs().S2))[1]
end

# Faster with function?
_omnivory(N) = expected_motif_count(find_motif(N, unipartitemotifs().S2))[1]
@time _omnivory(N)

# Motifs in parallel on layer for 1% of values
test = similar(layer, Float64)
@time @threads for k in keys(test)[1:Int64(ceil(length(keys(test))/1000))]
    test[k] = _omnivory(layer[k])
end

# Motifs in parallel on layer for 1% of values
test = similar(layer, Float64)
@time @threads for k in keys(test)[1:Int64(ceil(length(keys(test))/1000))]
    test[k] = _omnivory(layer[k])
end

## Example MWE

# Check the omnivory motif
unipartitemotifs().S2 |> adjacency

# Create networks
_mat = [0 0.7 0.3; 0 0 0.5; 0 0 0]
_mat = [0 0.5 0.5; 0 0 0.5; 0 0 0]
Nprob = UnipartiteProbabilisticNetwork(_mat)
Nbin = UnipartiteNetwork(_mat .> 0.0)

# Test omnivory function
omnivory(Nbin)
omnivory(Nprob)

# Test motifs
expected_motif_count(find_motif(Nprob, unipartitemotifs().S2))[1]
find_motif(Nbin, unipartitemotifs().S2)

# Compare omnivory vs motifs
function omnivory_single(Nbin)
    _res = omnivory(Nbin)
    _res_exp = sum(collect(values(_res)))^(length(_res))
    return _res_exp
end
_res1 = omnivory_single(Nbin)
_res2 = expected_motif_count(find_motif(Nprob, unipartitemotifs().S2))[1]

# Test on larger network
_res1 = omnivory_single(N)
_res2 = expected_motif_count(find_motif(N, unipartitemotifs().S2))[1]

# Final MWE
_mat = [0 0.7 0.3; 0 0 0.5; 0 0 0];
Nprob = UnipartiteProbabilisticNetwork(_mat)
expected_motif_count(find_motif(Nprob, unipartitemotifs().S2))[1]
