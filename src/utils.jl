import StatsBase.sample # to avoid conflict with fit and predict
using Random
using RCall
using StatsBase

"""
    div_into_folds(N::Int; K = 10, seed = 1234)

Equally divide `1:N` into `K` folds with random seed `seed`. If `seed` is negative, it is a non-random division, where the `i`-th fold would be the `i`-th equidistant range.
"""
@inline function div_into_folds(N::Int; K::Int = 10, seed = 1234)
    if seed >= 0
        idxes = sample(MersenneTwister(seed), 1:N, N, replace = false)
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
    stratify_div_into_folds(x; K, percent)

Perform a stratified fold division. The top `percent` and the remaining perform `div_into_folds(x; K)` respectively. 
"""
function stratify_div_into_folds(x::Vector{Float64}; K::Int = 10, percent = 0.2)
    # from small to large
    idx = sortperm(x)
    n = size(x, 1)
    n1 = ceil(Int, percent * n)
    n2 = n - n1
    idx1 = idx[end-n1+1:end]
    folds1 = div_into_folds(n1; K = K)
    folds2 = div_into_folds(n2; K = K)
    for i = 1:K
        folds2[i] .= idx[1:end-n1][folds2[i]]
        append!(folds2[i], idx[end-n1+1:end][folds1[i]])
    end
    return folds2
end

"""
    pick_knots(x::AbstractVector{T})

Partial code for picking knots in R's `smooth.spline`.

The source code of `smooth.spline` can be directly accessed via typing `smooth.spline` is an R session. Note that there might be different in different R versions. The code is adapted based on R 3.6.3.
"""
# refer to `smooth.spline`
function pick_knots(x::AbstractVector{T}; tol = 1e-6 * iqr(x), all_knots = false, scaled = false) where T <: AbstractFloat
    xx = round.(Int, (x .- mean(x)) / tol )
    # https://stackoverflow.com/questions/50899973/indices-of-unique-elements-of-vector-in-julia
    # unique index (Noted in techNotes)
    ud = unique(i -> xx[i], 1:length(xx))
    ux = sort(x[ud])
    idx0 = sortperm(x[ud])
    nx = length(ux)
    if all_knots
        nknots = nx
    else
        nknots = Int(rcopy(R".nknots.smspl($nx)"))
    end
    idx = round.(Int, range(1, nx, length = nknots))
    if scaled
        rx = (ux[end] - ux[1])
        mx = ux[1]
        ux = (ux .- mx) ./ rx
        return ux[idx], mx, rx, (1:length(x))[ud][idx0][idx], (1:length(x))[ud][idx0]
    else
        # to keep the same output
        return ux[idx], nothing, nothing, nothing, nothing 
    end
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