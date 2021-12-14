
test = rand(Bool, (3,2,2,2))
minitest = test[1, :, :, :]

function handmade_sum(test::Array{Bool, 4})
    ts = size(test)[1:3]
    by_site2 = zeros(Int64, (ts..., 1))
    for i in 1:ts[1]
        for j in 1:ts[2]
            for k in 1:[3]
                by_site2[i, j, k, 1] = sum(test[i, j, k, :])
            end
        end
    end
    return by_site2
end

@time test_sum_hand = handmade_sum(test)
@time test_sum = sum(test; dims=4)

test_sum_hand == test_sum

@time sum_networks_hand = handmade_sum(networks) # 105 sec
@time sum_networks = sum(networks; dims=4)

sum_networks_hand == sum_networks_hand


@time sum(minitest, dims=3)
@time mapslices(sum, minitest, dims=3)

@time mapslices(sum, test; dims=4)
@time mapslices(sum, networks; dims=4)

tmp = zeros(Int64, (2, 2, 1))
@time Threads.@threads for i in size(minitest)[1]
    for j in size(minitest)[2]
        tmp[i, j, 1] = sum(minitest[i, j, :])
    end
end

test

## Attempt at sqrt function
n = 10
values = rand(0:1, n, n)
values_sum =  sum(values; dims=2)
values_means = mean(values, dims=2)

function handmade_std(values)
    n = size(values, 1)
    values_sum =  sum(values; dims=2)
    values_means = mean(values, dims=2)
    sqrt.((values_sum .* (1 .- values_means).^2 .+ (n .- values_sum) .* (0 .- values_means).^2)./(n-1))
end

sq1 = sqrt.((values_sum .* (1 .- values_means).^2 .+ (n .- values_sum) .* (0 .- values_means).^2)./(n-1))
sq2 = map(std, eachrow(values))
sq3 = handmade_std(values)
sq1 == sq2
isapprox(sq1, sq2)
isapprox(sq3, sq2)