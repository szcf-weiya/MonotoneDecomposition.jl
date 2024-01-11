using GaussianProcesses
using Statistics
using LinearAlgebra
using Distributions
using LaTeXStrings
sigmoid(x::Float64; a = 5.0) = 1 / (1 + exp(-a*x))

"""
    gen_data(n, σ, f::Union{Function, String}; xmin = -1, xmax = 1, k = 10)

Generate `n` data points `(xi, yi)` from curve `f` with noise level `σ`, i.e., `yi = f(xi) + N(0, σ^2)`.

For `f`,

- if `f` is a `Function`, just take `y = f(x)`
- if `f = "MLP"`, it will be a simple neural network with one layer.
- otherwise, it accepts the string with format `KernelName_Para[_OtherPara]` representing some Gaussian Processes, including
    - `SE`, `Mat12`, `Mat32`, `Mat52`, `Para`: the length scale parameter `ℓ`
    - `Poly`: `Para` is the degree parameter `p`
    - `RQ`: `Para` is `ℓ` and `OtherPara` is `α`

It returns four vectors, `x, y, x0, y0`, where

- `x, y`: pair points of length `n`.
- `x0, y0`: true curve without noise, represented by `k*n` points.
"""
function gen_data(n::Int, σ::Real, f::Union{Function, String}; k=10, xmin = -1, xmax = 1)
    x0 = sort(rand(k*n-(k-1)) * (xmax - xmin) .+ xmin)
    if isa(f, Function)
        y0 = f.(x0)
    elseif occursin("MLP", f)
        M = parse(Int, split(f, '_')[2])
        w = rand(M)
        w = w / sum(w)
        α = randn(M)
        y0 = [sum(w .* sigmoid.(α .+ 5z)) for z in x0]
    else
        y0 = gp(x0, kernel = f)
    end
    x = x0[1:k:end]
    y = y0[1:k:end] + randn(n) * σ
    return x, y, x0, y0
    #return x, y, x0[2:2:end], y0[2:2:end] # to reduce evaluation cost
end

# mutable struct MLP1{T <: AbstractFloat}
#     w::AbstractVector{T}
#     α::AbstractVector{T}
#     β::AbstractVector{T}
# end

# function mlp1(m::MLP1, x::Real)
#     return sum(m.w * sigmoid.(m.α + m.β * x))
# end

"""
    gp(x; K)

Generate a random Gaussian vector with mean zero and covariance matrix `Σij = K(xi, xj)`.

The candidates of kernel `K` include SE, Mat12, Mat32, Mat52.

See also: <https://stats.hohoweiya.xyz/2021/12/13/GP/>
"""
function gp(x::AbstractVector{T}; kernel = "SE_0.1") where T <: AbstractFloat
    strs = split(kernel, "_")
    str_k = strs[1]
    str_ℓ = strs[2]
    # https://stackoverflow.com/questions/34016768/julia-invoke-a-function-by-a-given-string
    f = getfield(GaussianProcesses, ifelse(str_k in ["Poly", "Periodic"], Symbol(str_k), Symbol("$(str_k)Iso")))
    ℓ = parse(Float64, str_ℓ)
    # if occursin("SE", kernel)
    #     Σ = SEIso(log(ℓ), 0)
    # elseif occursin("Mat12", kernel)
    #     # Mat12_0.1
    #     Σ = Mat12Iso(log(ℓ), 0)
    # elseif occursin("Mat32", kernel)
    #     Σ = Mat32Iso(log(ℓ), 0)
    # elseif occursin("Mat52", kernel)
    # end
    # Cannot use ifelse, since it would also evaluate othe conditions
    # Σ = ifelse(str_k == "Poly", f(0.0, 0.0, ℓ),  # ℓ is the degree in that case
    #                             ifelse(str_k == "RQ", f(log(ℓ), 0.0, 1/2), 
    #                                                   f(log(ℓ), 0.0)
    #                                   )
    #           )
    Σ = SEIso(0, 0)
    if str_k == "Poly"
        # ℓ is the degree in that case
        Σ = f(0.0, 0.0, Int.(ℓ))
    elseif str_k in ["RQ", "Periodic"]
        if length(strs) == 2
            Σ = f(log(ℓ), 0.0, log(1/2))
        else
            α = parse(Float64, strs[3])
            Σ = f(log(ℓ), 0.0, log(α))
        end
    else
        Σ = f(log(ℓ), 0.0)
    end
    K = cov(Σ, reshape(x, 1, :))
    ϵ = sqrt(eps())
    return rand(MvNormal(K + ϵ * 1.0I))
end

function gen_kern(x::AbstractVector{T}, ℓ::Real = 1) where T <: AbstractFloat
    return gen_kern(x, x, ℓ)
end

function gen_kern(x::AbstractVector, y::AbstractVector, ℓ::Real = 1)
    n1 = length(x)
    n2 = length(y)
    K = zeros(n1, n2)
    for i = 1:n1
        for j = 1:n2
            K[i, j] = exp( -(x[i] - y[j])^2/(2ℓ^2) )
            # K[i, j] = exp( -abs(y[j] - x[i])/2 )
        end
    end
    return K
end
