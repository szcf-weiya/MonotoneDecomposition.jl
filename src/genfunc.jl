using GaussianProcesses
using Statistics
using LinearAlgebra
using Distributions
using LaTeXStrings
sigmoid(x::Float64) = 1 / (1 + exp(-5x))

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

function gen_data_GP(n::Int, m::Int; ℓ = 1, n0 = 5, prediction = false, σ = 0.0)
    # x = rand(n) * 2 .- 1
    x0 = rand(n0) * 2 .- 1
    K0 = gen_kern(x0)
    f = rand(MvNormal(K0)) 
    x = range(-1, 1, length = n)
    K1 = gen_kern(x)
    K01 = gen_kern(x0, x)
    if !prediction
        # if without prediction, return here
        # return x0, f
        return x, rand(MvNormal(K1 + sqrt(eps()) * 1.0I ), m)
    else
        # by prediction
        iK0 = inv(K0 + σ^2 * 1.0I)
        μ = K01' * iK0 * f
        Σ = K1 - K01' * iK0 * K01
        println(eigen(Σ).values)
        # dist = MvNormal(K)
        dist = MvNormal(μ, Symmetric(Σ))
        return x0, f, x, rand(dist, m)
    end
end

# Square Exponential
function sq_exp(x::AbstractVector{T}; ℓ = 1) where T <: AbstractFloat
    K = gen_kern(x, ℓ)
    ϵ = sqrt(eps())
    return rand(MvNormal(K + ϵ * 1.0I))
end

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

function plot_functions()
    n = 100
    Random.seed!(1)
    x = sort(rand(n)) * 2 .- 1
    figsize = (400, 400)
    fig1 = plot(x, x .^ 3, label = L"x^3", ls = :solid, size = figsize, legend = :top)
    plot!(fig1, x, x .^ 2, label = L"x^2", ls = :dash)
    plot!(fig1, x, exp.(x) .- 1, label = L"\exp(x)-1", ls = :dot)
    plot!(fig1, x, sigmoid.(x), label = L"1/(1+\exp(-5x))", ls = :dashdot)
    plot!(fig1, x, sin.(2π * x), label = L"\sin(2\pi x)", ls = :dashdotdot)
    # plot(x, gp(x, kernel = "SE_0.1"))
    fig2 = plot(x, gp(x, kernel = "SE_1"), label = "SE_1", ls = :solid, size = figsize, legend = :topright)
    plot!(fig2, x, gp(x, kernel = "SE_0.1"), label = "SE_0.1", ls = :dash)
    plot!(fig2, x, gp(x, kernel = "SE_0.5"), label = "SE_0.5", ls = :dot)
    fig3 = plot(x, gp(x, kernel = "Mat12_1"), label = "Mat12", ls = :solid, size = figsize, legend = :topleft)
    plot!(fig3, x, gp(x, kernel = "Mat32_1"), label = "Mat32", ls = :dash)
    plot!(fig3, x, gp(x, kernel = "Mat52_1"), label = "Mat52", ls = :dot)
    plot!(fig3, x, gp(x, kernel = "RQ_1_0.5"), label = "RQ", ls = :dashdot)
    plot!(fig3, x, gp(x, kernel = "Periodic_1_1"), ls = :dashdotdot, label = "Periodic")
    plot(fig1, fig2, fig3, layout = (1, 3), size = (3*figsize[1], figsize[2]) )
    savefig("~/PGitHub/overleaf/MonoDecomp/figs/funcs.pdf")
end