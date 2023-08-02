import StatsBase.sample # to avoid conflict with fit and predict
using Random
using RCall
using Printf
using StatsBase

"""
    div_into_folds(N::Int; K = 10, seed = 1234)

Equally divide `1:N` into `K` folds with random seed `seed`. Specially,

- If `seed` is negative, it is a non-random division, where the `i`-th fold would be the `i`-th equidistant range.
- If `seed = 0`, it is a non-random division, where each fold consists of equidistant indexes.
"""
@inline function div_into_folds(N::Int; K::Int = 10, seed = 1234)
    if seed > 0
        idxes = sample(MersenneTwister(seed), 1:N, N, replace = false)
    elseif seed == 0
        @assert N % K == 0
        n = Int(N / K)
        return [collect(i:K:N) for i in 1:K]
    else
        idxes = collect(1:N)
    end
    # maximum quota per fold
    n = Int(ceil(N/K))
    # number folds for the maximum quota
    k = N - (n-1)*K
    # number fols for n-1 quota: K-k
    folds = Array{Array{Int, 1}, 1}(undef, K)
    for i = 1:k
        folds[i] = idxes[collect(n*(i-1)+1:n*i)]
    end
    for i = 1:K-k
        folds[k+i] = idxes[collect((n-1)*(i-1)+1:(n-1)*i) .+ n*k]
    end
    return folds
end

"""
    coverage_prob(CIs::AbstractMatrix, y0::AbstractVector)

Calculate coverage probability given `n x 2` CI matrix `CIs` and true vector `y0` of size `n`.
"""
function coverage_prob(CIs::AbstractMatrix, y0::AbstractVector)
    @assert size(CIs, 1) == length(y0)
    @assert size(CIs, 2) == 2
    return sum((CIs[:, 1] .<= y0) .* (CIs[:, 2] .>= y0)) / length(y0)
end

"""
    conf_band_width(CIs::AbstractMatrix)

Calculate width of confidence bands.
"""
function conf_band_width(CIs::AbstractMatrix)
    len = CIs[:, 2] - CIs[:, 1]
    return mean(len)
end
