## Mean function

# Create test matrices
test = rand(Bool, (3,2,2,2))
minitest = test[1, :, :, :]

# Custom function to calculate sum on 4th dimension
function handmade_sum(test::Array{Bool, 4})
    ts = size(test)[1:3]
    by_site2 = zeros(Int64, (ts..., 1))
    for i in 1:ts[1]
        for j in 1:ts[2]
            for k in 1:ts[3]
                by_site2[i, j, k, 1] = sum(test[i, j, k, :])
            end
        end
    end
    return by_site2
end

# Is the custom sum function equivalent
@time test_sum_hand = handmade_sum(test)
@time test_sum = sum(test; dims=4)
test_sum_hand == test_sum # yep
# and it's pretty fast for such a small task

# What about for the whole network?
@time sum_networks_hand = handmade_sum(networks_thr); # 105 sec
@time sum_networks = sum(networks_thr; dims=4); # 4.5 sec
sum_networks_hand == sum_networks_hand # it's the same
# But the speed difference is huge!!

# Is mapslices comparable?
@time sum(minitest, dims=3)
@time mapslices(sum, minitest, dims=3)
# For a small one yes

# Is mapslices comparable for something larger?
@time sum(test; dims=4) # instant, 8 alloc
@time mapslices(sum, test; dims=4) # instant, but 175 alloc
# @time mapslices(sum, networks_thr; dims=4) # so long I stopped it along the way

# Does threading help in any way?
tmp_nothread = zeros(Int64, (2, 2, 1))
tmp_thread = zeros(Int64, (2, 2, 1))
@time for i in 1:size(minitest)[1]
    for j in 1:size(minitest)[2]
        tmp_nothread[i, j, 1] = sum(minitest[i, j, :])
    end
end
@time Threads.@threads for i in 1:size(minitest)[1]
    for j in 1:size(minitest)[2]
        tmp_thread[i, j, 1] = sum(minitest[i, j, :])
    end
end
tmp_nothread
tmp_thread
tmp_nothread == tmp_thread
# Slower for mini matrix of course

# How about a bigger one?
tmp_nothread = zeros(Int64, size(test)[1:3]..., 1)
tmp_thread = zeros(Int64, size(test)[1:3]..., 1)
@time for i in 1:size(test)[1]
    for j in 1:size(test)[2]
        for k in 1:size(test)[3]
            tmp_nothread[i, j, k, 1] = sum(test[i, j, k, :])
        end
    end
end
@time Threads.@threads for i in 1:size(test)[1]
    for j in 1:size(test)[2]
        for k in 1:size(test)[3]
            tmp_nothread[i, j, k, 1] = sum(test[i, j, k, :])
        end
    end
end
tmp_nothread
tmp_thread
tmp_nothread == tmp_thread
# Still much slower

## Attempt at sqrt function

# Generate values
n = 10
values = rand(0:1, n, n)
values_sum =  sum(values; dims=2)
values_means = mean(values, dims=2)

# Custom std function
function handmade_std(values)
    n = size(values, 1)
    values_sum =  sum(values; dims=2)
    values_means = mean(values, dims=2)
    sqrt.((values_sum .* (1 .- values_means).^2 .+ (n .- values_sum) .* (0 .- values_means).^2)./(n-1))
end

# Test things
@time sq1 = sqrt.((values_sum .* (1 .- values_means).^2 .+ (n .- values_sum) .* (0 .- values_means).^2)./(n-1))
@time sq2 = map(std, eachrow(values))
@time sq3 = handmade_std(values)
sq1 == sq2
isapprox(sq1, sq2)
isapprox(sq3, sq2)