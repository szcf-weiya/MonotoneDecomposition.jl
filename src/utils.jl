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

"""
    cv_one_se_rule(μs::AbstractVector{T}, σs::AbstractVector{T}; small_is_simple = true)
    cv_one_se_rule(μs::AbstractMatrix{T}, σs::AbstractMatrix{T}; small_is_simple = [true, true])
    cv_one_se_rule2(μs::AbstractMatrix{T}, σs::AbstractMatrix{T}; small_is_simple = [true, true])

Return the index of parameter(s) (1dim or 2-dim) that minimize the CV error with one standard error rule.

For 2-dim parameters, `cv_one_se_rule2` adopts a grid search for `μ+σ` while `cv_one_se_rule` searchs after fixing one optimal parameter. 
The potential drawback of `cv_one_se_rule2` is that we might fail to determine the simplest model when both parameters are away from the optimal parameters.
So we recommend `cv_one_se_rule`.
"""
function cv_one_se_rule(μs::AbstractVector{T}, σs::AbstractVector{T}; small_is_simple = true) where T <: AbstractFloat
    iopt = argmin(μs)
    err_plus_se = μs[iopt] + σs[iopt]
    if small_is_simple
        ii = iopt-1:-1:1
    else
        ii = iopt+1:length(μs)
    end
    for i = ii
        if μs[i] < err_plus_se
            iopt = i
        # else
        #     break
        end
    end
    return iopt
end

function cv_one_se_rule(μs::AbstractMatrix{T}, σs::AbstractMatrix{T}; small_is_simple = [true, true]) where T <: AbstractFloat
    ind = argmin(μs)
    iopt, jopt = ind[1], ind[2]
    jopt1 = cv_one_se_rule(μs[iopt,:], σs[iopt,:], small_is_simple = small_is_simple[2])

    iopt1 = cv_one_se_rule(μs[:,jopt], σs[:,jopt], small_is_simple = small_is_simple[1])

    @debug "without 1se rule: iopt = $iopt, jopt = $jopt"
    # determine the minimum one: compare (iopt1, jopt) vs (iopt, jopt1), take the larger one since both are in 1se
    if μs[iopt1, jopt] < μs[iopt, jopt1]
        return iopt, jopt1
    else
        return iopt1, jopt
    end
end

function cv_one_se_rule2(μs::AbstractMatrix{T}, σs::AbstractMatrix{T}; small_is_simple = [true, true]) where T <: AbstractFloat
    ind = argmin(μs)
    iopt, jopt = ind[1], ind[2]
    err_plus_se = μs[iopt, jopt] + σs[iopt, jopt]
    ii = ifelse(small_is_simple[1], iopt:-1:1, iopt:size(μs, 1))
    jj = ifelse(small_is_simple[2], jopt:-1:1, jopt:size(μs, 2))
    @debug ii
    @debug jj
    best_err = 0.0
    for i = ii
        for j = jj
            if μs[i, j] < err_plus_se
                if μs[i, j] > best_err # pick the one with largest error (different from univariate case, WARN: no a single direction to a simpler model)
                    best_err = μs[i, j]
                    iopt = i
                    jopt = j
                end
            end
        end
    end
    return iopt, jopt
end
