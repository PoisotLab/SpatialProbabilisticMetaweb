# Checking network entropy
N = layer[1]
EcologicalNetworks.entropy(N)
H = broadcast(EcologicalNetworks.entropy, layer)
plot(H; c=:cividis, title="Entropy")

# Network conversion
make_joint_distribution(N) |> x -> convert(Array{Float64}, x)
# Indices
EcologicalNetworks.entropy(N)
EcologicalNetworks.entropy(N; dims=1)
EcologicalNetworks.entropy(N; dims=2)
conditional_entropy(N, 1)
conditional_entropy(N, 2)
conditional_entropy(N, 3) # ???
mutual_information(N)
variation_information(N)
potential_information(N)
diff_entropy_uniform(N)
# Decomposition
information_decomposition(N)
# Effective interactions
# convert2effective

sum(N)
rand(N)
degree(N)
specificity(N)
species(N)
interactions(N)
[interaction for interaction in N]

all_hp_data = filter(x -> occursin("Hadfield", x.Reference), web_of_life())
ids = getfield.(all_hp_data, :ID)
networks = convert.(BinaryNetwork, web_of_life.(ids))
N = networks[1]
P1 = null2(N)
R1 = rand(P1, 9)

βos(layer[1], layer[2])