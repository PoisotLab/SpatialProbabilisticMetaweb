function LCBD(Y)
    S = (Y .- mean(Y; dims=1)).^2.0
    SStotal = sum(S)
    BDtotal = SStotal / (size(Y,1)-1)
    SSj = sum(S; dims=1)
    SCBDj = SSj ./ SStotal
    SSi = sum(S; dims=2)
    LCBDi = SSi ./ SStotal
    return LCBDi, SCBDj, BDtotal
end

function hellinger(Y::AbstractMatrix{T}) where {T <: Number}
    yi = sum(Y; dims=2)
    return sqrt.(Y ./ yi)
end