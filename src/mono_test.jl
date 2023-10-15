using Plots
using LaTeXTables
using QuadGK
using ProgressMeter
using LaTeXStrings

using RCall


function gen_data_bowman(;n = 50, a = 0, σ = 0.1, plt = false, kw...)
    xs = range(0, 1, length = n)
    f(x, a) = 1 + x - a * exp(-1/2 * (x-0.5)^2 / 0.1^2)
    ys0 = [f(x, a) for x in xs] 
    ys = ys0 + randn(n) * σ
    if plt
        plot(xs, ys0; legend = false, title = "n = $n, a = $a, σ = $σ", kw...)
        scatter!(xs, ys)
    else
        return xs, ys
    end
end

function demo_data(; figfolder = "/tmp" # "../notes"
                   )
    figs = Plots.Plot[]
    figsize = (400, 400)
    for a in [0.0, 0.15, 0.25, 0.45]
        push!(figs, gen_data_bowman(a = a, plt = true, size = figsize))
    end
    plot(figs..., layout = (2, 2), size = (2 * figsize[1], 2 * figsize[2]))
    savefig(joinpath(figfolder, "ex-bowman.pdf"))
    figs2 = gen_data_ghosal(plt = true)
    plot(figs2..., layout = (2, 2), size = (2 * figsize[1], 2 * figsize[2]))
    savefig(joinpath(figfolder, "ex-ghosal.pdf"))
end

function gen_mono_data(; n = 100, σ = 0.1, equidistant = false)
    if equidistant
        x = range(0, 1, length = n)
    else
        x = rand(n)
    end
    m1 = x + randn(n) * σ
    m2 = x .^3 + randn(n) * σ
    m3 = cbrt.(x) + randn(n) * σ
    m4 = exp.(x) .- 1 + randn(n) * σ
    m5 = 1 ./ (1 .+ exp.(-x)) + randn(n) * σ
    return x, m1, m2, m3, m4, m5
end

function gen_desc_data(; n = 100, σ = 0.1)
    x, m1, m2, m3, m4, m5 = gen_mono_data(; n = n, σ = σ)
    return x, -m1, -m2, -m3, -m4, -m5
end

function gen_data_ghosal(; n = 100, σ = 0.1, plt = false)
    xs = rand(n)
    m10 = zeros(n)
    m20 = xs .* (1 .- xs)
    m30 = xs .+ 0.415 * exp.(-50 * xs.^2)
    m40 = (10 * (xs .- 0.5).^3 - exp.(-100 * (xs .- 0.25).^2) ) .* (xs .< 0.5) + (0.1 * (xs .- 0.5) - exp.(-100 * (xs .- 0.25).^2)) .* (xs .>= 0.5)
    m1 = m10 + randn(n) * σ
    m2 = m20 + randn(n) * σ
    m3 = m30 + randn(n) * σ
    m4 = m40 + randn(n) * σ
    if plt
        ind = sortperm(xs)
        p1 = plot(xs[ind], m10[ind], legend = false)
        scatter!(p1, xs, m1)
        p2 = plot(xs[ind], m20[ind], legend = false)
        scatter!(p2, xs, m2)
        p3 = plot(xs[ind], m30[ind], legend = false)
        scatter!(p3, xs, m3)
        p4 = plot(xs[ind], m40[ind], legend = false)
        scatter!(p4, xs, m4)
        return [p1, p2, p3, p4]
    else
        return xs, m1, m2, m3, m4
    end
end

function ghosal_λ()
    # -1 < x < 1
    k(x) = 0.75 * (1 - x^2)
    dk(x) = -1.5x
    d2k(x) = -1.5
    # ∫z<=x k(x) dz
    K(x) = 0.75(x+1) - 0.25(x^3 + 1)
    denom = quadgk(x-> (2K(x) - 1)^2 * k(x) ^ 2, -1, 1 )[1]
    num1 = 6 * quadgk(x-> (2K(x) - 1) * k(x)^2 * dk(x), -1, 1)[1]
    num2 = quadgk(x-> (2K(x)-1)^2 * k(x) * d2k(x), -1, 1)[1]
    return -(num1 + num2) / denom
end


# const LAMBDA = ghosal_λ() # 9.974576271186443
const LAMBDA = 9.974576271186443

function ghosal_critical_value(a = 0.05, b = 0.95, α = 0.05)
    cs = Dict()
    for n in [50, 100, 200, 500]
        hn = 0.5 * n^(-1/5)
        c1 = sqrt(2log((b-a)/hn))
        c2 = log(-log(1-α)) - log(sqrt(LAMBDA) / 2π)
        c = c1 - c2 / c1
        cs[n] = c
        # println("n = $n, bw = $hn, c = $c")
    end
    return cs
end

const C_GHOSAL = ghosal_critical_value()

function ghosal_Un(t, X, Y, kn)
    n = length(X)
    u = 0
    for i = 1:n
        for j = i+1:n
            u += sign(Y[j] - Y[i]) * sign(X[i] - X[j]) * kn(X[i] - t) * kn(X[j] - t)
        end
    end
    return u * 2 / n / (n - 1)
end

function ghosal_σn(t, X, kn)
    n = length(X)
    σ2 = 0
    for i = 1:n
        for j = 1:n
            for k = 1:n
                if i == j == k
                    continue
                end
                σ2 += sign(X[i] - X[j]) * sign(X[i] - X[k]) * kn(X[j] - t) * kn(X[k] - t) * kn(X[i] - t)^2
            end
        end
    end
    return sqrt(σ2 * 4 / 3n / (n-1) / (n-2))
end

function ghosal_S1n(X, Y, c; a = 0.05, b = 0.95, α = 0.05)
    k = x -> 0.75(1-x^2) * (-1 < x < 1)
    n = length(X)
    hn = 0.5 * n^(-1/5)
    nt = round(Int, 1 / hn * 2)
    kn(x) = k(x / hn) / hn
    S(t) = sqrt(n) * ghosal_Un(t, X, Y, kn) / ghosal_σn(t, X, kn)
    # Question: how to calculate the supremum
    ts = range(a, b, length = nt)
    S1n = maximum(S.(ts))
    return S1n > c
end

function ghosal_S1n(X, Y; a = 0.05, b = 0.95, α = 0.05)
    k = x -> 0.75(1-x^2) * (-1 < x < 1)
    n = length(X)
    hn = 0.5 * n^(-1/5)
    nt = round(Int, 1 / hn * 2)
    kn(x) = k(x / hn) / hn
    S(t) = sqrt(n) * ghosal_Un(t, X, Y, kn) / ghosal_σn(t, X, kn)
    # Question: how to calculate the supremum
    ts = range(a, b, length = nt)
    S1n = maximum(S.(ts))
    an = sqrt(2log((b-a)/hn))
    bn = an + log(sqrt(LAMBDA) / 2π) / an
    return 1-exp(-exp(-an * (S1n - bn)))
end


ghosal(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real = ghosal_S1n(x, y)

function meyer(x::AbstractVector{T}, y::AbstractVector{T}, nsim = 100, k = 6) where T <: Real
    meyer_rfile = joinpath(@__DIR__, "testmonotonicity.R")
    R"source($meyer_rfile)" # not found if put outside the function when wrapped into a package
    n = length(x)
    w = ones(n)
    return rcopy(R"montest($x, $y, 0, $w, $k, $nsim, 1)$pval")
end

function bowman(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real
    pval = 1
    try
        pval = rcopy(R"sm::sm.monotonicity($x, $y, display='none')$p")
        if isnothing(pval)
            pval = 1
        end
    catch
        pval = 1
    end
    return pval
end

# run single experiment to make comparisons on bowan's curves
function single_test_compare_bowman(;
        as = [0, 0.15, 0.25, 0.45],
        ns = [50, 100, 200],
        σs = [0.001, 0.01, 0.1],
        nrep = 500, kw...
    )
    # proposed, ghosal, meyer
    pvals = zeros(length(as), length(ns), length(σs), 5)
    for (j, a) in enumerate(as)
        for (k, n) in enumerate(ns)
            for (l, σ) in enumerate(σs)
                x, y = gen_data_bowman(a = a, n = n, σ = σ)
                pvals[j, k, l, 1] = meyer(x, y)
                pvals[j, k, l, 2] = ghosal(x, y)
                pvals[j, k, l, 3] = bowman(x, y)
                pvals[j, k, l, 4] = mono_test_bootstrap_cs(x, y; nrep = nrep, kw...)[1]
                pvals[j, k, l, 5] = mono_test_bootstrap_ss(x, y; nrep = nrep, kw...)[1]
            end
        end
    end
    return pvals
end

# run single experiment to make comparisons on ghosal's curves
function single_test_compare_ghosal(;
        ns = [50, 100, 200],
        σs = [0.001, 0.01, 0.1],
        nrep = 500, kw...
    )
    # proposed, ghosal, meyer, sm
    pvals = zeros(length(ns), length(σs), 4, 5)
    for (i, n) in enumerate(ns)
        for (k, σ) in enumerate(σs)
            x, m1, m2, m3, m4 = gen_data_ghosal(n = n, σ = σ)
            for (j, y) in enumerate([m1, m2, m3, m4])
                pvals[i, k, j, 1] = meyer(x, y)
                pvals[i, k, j, 2] = ghosal(x, y)
                pvals[i, k, j, 3] = bowman(x, y)
                pvals[i, k, j, 4] = mono_test_bootstrap_cs(x, y; nrep = nrep, kw...)[1]
                pvals[i, k, j, 5] = mono_test_bootstrap_ss(x, y; nrep = nrep, kw...)[1]
            end
        end
    end
    return pvals
end

function single_test_compare_mono(;
        ns = [50, 100, 200],
        σs = [0.001, 0.01, 0.1],
        equidistant = false,
        nrep = 500, kw...
    )
    # proposed, ghosal, meyer, sm
    pvals = zeros(length(ns), length(σs), 5, 5)
    for (i, n) in enumerate(ns)
        for (k, σ) in enumerate(σs)
            x, m1, m2, m3, m4, m5 = gen_mono_data(n = n, σ = σ, equidistant = equidistant)
            for (j, y) in enumerate([m1, m2, m3, m4, m5])
                pvals[i, k, j, 1] = meyer(x, y)
                pvals[i, k, j, 2] = ghosal(x, y)
                pvals[i, k, j, 3] = bowman(x, y)
                pvals[i, k, j, 4] = mono_test_bootstrap_cs(x, y; nrep = nrep, kw...)[1]
                pvals[i, k, j, 5] = mono_test_bootstrap_ss(x, y; nrep = nrep, kw...)[1]
            end
        end
    end
    return pvals
end

function mono_test_bootstrap_cs(x::AbstractVector{T}, y::AbstractVector{T}; nrep = 100, 
                                            μs = 10.0 .^ (-6:0.1:2), Js = 4:50, fixJ = true,
                                            nblock = -1, # wild bootstrap
                                            one_se_rule = false, 
                                            one_se_rule_pre = false, 
                                            figname = nothing,
                                            nfold = 10, nfold_pre = 10,
                                            use_GI = true,
                                            h0_mono = false, # increasing or decreasing
                                            kw...)::Tuple{T, MonoDecomp{T}} where T <: Real
    D, μ = cv_mono_decomp_cs(x, y, ss = μs, one_se_rule = one_se_rule, fixJ = fixJ, Js = Js, one_se_rule_pre = one_se_rule_pre, figname = figname, nfold = nfold, nfold_pre = nfold_pre, use_GI = use_GI)
    @debug D.γdown
    J = D.workspace.J
    # scatter(x, y)
    # plot!(x, res.yhat)
    c = mean(D.yhat) / 2
    error = y .- D.yhat
    @debug maximum(max.(error))
    # σ = std(err)
    tobs = var(D.γdown)
    flag_incr = true
    if h0_mono
        tobs1 = var(D.γup)
        if tobs1 < tobs
            flag_incr = false
            tobs = tobs1
        end
    end
    ts = zeros(nrep)
    # println("address of D.w: ", pointer_from_objref(D.workspace))
    yhat = D.workspace.B * D.γup .+ c
    if !flag_incr
        yhat = D.workspace.B * D.γdown .+ c
    end
    for i = 1:nrep
        yi = construct_bootstrap_y(y, error, yhat, nblock = nblock)
        try
            Di = mono_decomp_cs(x, yi, s = μ, s_is_μ = true, J = J, workspace = D.workspace, use_GI = use_GI)
            # println("address of Di.w: ", pointer_from_objref(Di.workspace))
            if h0_mono
                ts[i] = min(var(Di.γdown), var(Di.γup))
            else
                ts[i] = var(Di.γdown)
            end
        catch e
            @warn "due to error $e in optimization, assign test statistic as Inf"
            ts[i] = Inf
        end
    end
    @debug ts
    @debug tobs
    pval = sum(ts .> tobs) / nrep
    return pval, D
end

function mono_test_bootstrap_sup(x::AbstractVector{T}, y::AbstractVector{T}; 
                                                            nrep = 100, nμ = 10, 
                                                            nfold = 2, fixJ = true, 
                                                            nblock = 10,
                                                            kw...) where T <: AbstractFloat
    n = length(y)
    μs0 = 10.0 .^ (-6:0.5:6)
    D1, μ0 = cv_mono_decomp_cs(x, y, ss = μs0, one_se_rule = true, fixJ = fixJ, nfold = nfold)
    μ1 = D1.μ
    J = D1.workspace.J
    # μ0 < μ1
    if μ1 == μ0
        μ1 = maximum(μs0)
        μ0 = minimum(μs0)
    end
    # μ0 is without 1se rule
    μl = μ0 - 0.75 * (μ1 - μ0)
    if μl < 0
        μl = μ0 * 0.25
    end
    μr = μ1 #μ1 + 0.75 * (μ1 - μ0)
    μs = vcat(range(μl, μr, length = nμ), μ1, μ0)
    # μs = range(μ0, μ1, length = nμ)
    pval = Float64[]
    for (k, μ) in enumerate(μs)
        D = mono_decomp_cs(x, y, s = μ, s_is_μ = true, J = J)
        error = y - D.yhat
        tobs = var(D.γdown) #/ var(y - D.yhat)
        ts = zeros(nrep)
        c = mean(D.yhat) / 2
        for i = 1:nrep
            yi = construct_bootstrap_y(y, error, D.workspace.B, D.γup, c, nblock = nblock)
            Di = mono_decomp_cs(x, yi, s = μ, s_is_μ = true, J = J)
            # ts[i] = var(Di.γdown) / var(y - Di.yhat)
            ts[i] = var(Di.γdown) 
        end
        # pval[k] = sum(ts .> tobs) / nrep
        append!(pval, sum(ts .> tobs) / nrep)
    end
    # return pval
    # must include one
    # if length(pval) == 0
    #     @warn "$ρ is too small, and no mono decomp satisfies the condition"
    # end
    @debug "pval = $pval"
    return [maximum(pval) < 0.05, pval[end] < 0.05, pval[end-1] < 0.05]
end

maxgap(x::AbstractVector{T}) where T <: Real = maximum(x) - minimum(x)

## aim for hete error, but if we can change different md decomposition method such that the md is homo, then no need (paper#12). 
function block_bootstrap_idx(n::Int; nblock = 10)
    # suppose x is sorted
    blocks = div_into_folds(n, K = nblock, seed = -1)
    idx = vcat([sample(z, length(z)) for z in blocks]...)
    ## x might not be sorted, idx is for sorted x
    ## sort(x)[idx] = x[?] => x = sort(x)[idx][invperm(?)]
    ## idx for sorted x, what is the index for x?
    ## sort(x)[idx]
    ## note that sort(x) = x[sortperm(x)]
    return idx
end

function block_bootstrap_idx(e::AbstractVector; nblock = 10)
    n = length(e)
    blocks = div_into_folds(n, K = nblock, seed = -1)
    idxs = [sample(b, length(b)) for b in blocks]
    μb = vcat([mean(e[b]) * ones(length(b)) for b in blocks]...)
    μi = vcat([mean(e[b]) * ones(length(b)) for b in idxs]...)
    idx = vcat(idxs...)
    return idx, μb, μi
end

## yhat = B*γ .+ c
function construct_bootstrap_y(y::AbstractVector{T}, e::AbstractVector{T}, yhat::AbstractVector{T}; nblock = -1, σe = std(e), debias_mean_yi = true) where T <: AbstractFloat
    n = length(y)
    if nblock > 0
        idx, μb, μi = block_bootstrap_idx(e; nblock = nblock)
        ei = e[idx] - μi + μb
    elseif nblock == 0
        ei = randn(n) * σe
        ei = ei .- mean(ei) .+ mean(e)
    elseif nblock == -1 # wild bootstrap
        ei = randn(n) .* e
    elseif nblock == -2
        ei = ifelse.(randn(n) .> 0, 1, -1) .* e
    else
        idx = sample(1:n, n)
        ei = e[idx] .- mean(e[idx]) .+ mean(e)
    end
    yi = yhat .+ ei
    if debias_mean_yi
        yi .= yi .- mean(yi) .+ mean(y) 
    end
    return yi
end

function construct_bootstrap_y(y::AbstractVector{T}, e::AbstractVector{T}, B::AbstractMatrix{T}, γ::AbstractVector{T}, c::T; nblock = -1, σe = std(e), debias_mean_yi = true) where T <: AbstractFloat
    return construct_bootstrap_y(y, e, B * γ .+ c, nblock = nblock, σe = σe, debias_mean_yi = debias_mean_yi)
end

function mono_test_bootstrap_supss(x::AbstractVector{T}, y::AbstractVector{T}; 
                                    nrep = 100, nμ = 10, nfold = 2, seed = rand(UInt64),
                                    opstat::Union{String, Function} = var,
                                    md_method = "single_lambda",
                                    tol = 1e-4,
                                    nblock = -1,
                                    use_σ_from_ss = false,
                                    debias_mean_yi = true,
                                    kw...
                                    ) where T <: Real
    # for block index
    idx = sortperm(x)
    x = x[idx]
    y = y[idx]
    n = length(y)
    res, μ0, μs0, errs, σerrs, yhat, yhatnew = cv_mono_decomp_ss(x, y; one_se_rule = true, nfold = nfold, seed = seed, method = md_method, tol = tol, kw...)
    σe0 = std(y - yhat)
    μ1 = res.μ
    # μ0 < μ1
    # μ0 is not with 1se rule
    if μ1 == μ0
        μ1 = maximum(μs0)
        μ0 = minimum(μs0)
    end
    μl = μ0 - 0.75 * (μ1 - μ0)
    if μl < 0
        # μl = μ0 * 0.25
        μl = min(cbrt(eps()), 0.25μ0)
    end
    μr = μ1#μ1 + 0.75 * (μ1 - μ0)
    μs = vcat(range(μl, μr, length = nμ), μ1, μ0)
    @debug "perform monotone test: μ_min = $μ0, μ_1se = $μ1"
    @debug "construct composite hypothesis in the range [$μl, $μr]"
    # μs = range(μ0, μ1, length = nμ)
    pval = Float64[]
    for (k, μ) in enumerate(μs)
        D = mono_decomp_ss(res.workspace, x, y, res.λ, μ)
        error = y - D.yhat
        σe = std(error)
        ts = zeros(nrep)
        c = mean(D.yhat) / 2
        if isa(opstat, Function)
            tobs = opstat(D.γdown)
        else
            tobs = sum((D.γdown .- c).^2)
        end
        @debug "σe = $σe"
        for i = 1:nrep
            yi = construct_bootstrap_y(y, error, D.workspace.B, D.γup, c, nblock = nblock, 
                                        σe = ifelse(use_σ_from_ss, σe0, σe),
                                        debias_mean_yi = debias_mean_yi)
            try
                Di = mono_decomp_ss(res.workspace, x, yi, res.λ, μ, strict = true)                
                if isa(opstat, Function)
                    ts[i] = opstat(Di.γdown)
                else
                    ts[i] = sum((Di.γdown .- c).^2)
                end
            catch
                @warn "due to error in optimization, assign test statistic as Inf"
                # TODO: is Inf reasonable?
                # alternatively, exclude Inf when calculting p-value
                ts[i] = Inf
            end
            # ts[i] = var(Di.γdown) / var(y - Di.yhat)
            # savefig(scatter(x, yi), "/tmp/toy-$i.png")
            # ts[i] = sum((Di.γdown .- c).^2)
        end
        # if k == 2
        #     fig = histogram(ts[.!isinf.(ts)])
        #     vline!(fig, [tobs])
        #     savefig(fig, "/tmp/mu.png")
        # end
        # println(ts)
        # println(tobs)
        # pval[k] = sum(ts .> tobs) / nrep
        append!(pval, sum(ts .> tobs) / nrep)
    end
    # return pval
    # must include one
    # if length(pval) == 0
    #     @warn "$ρ is too small, and no mono decomp satisfies the condition"
    # end
    @debug "pval = $pval"
    return [maximum(pval) < 0.05, pval[end] < 0.05, pval[end-1] < 0.05]
end

"""
    mono_test_bootstrap_ss(x, y)

Perform monotonicity test after monotone decomposition with smoothing splines.
"""
function mono_test_bootstrap_ss(x::AbstractVector{T}, y::AbstractVector{T}; nrep = 100, 
                                                                one_se_rule = false,
                                                                one_se_rule_pre = false,
                                                                nfold = 5,
                                                                seed = rand(UInt64),
                                                                md_method = "double_grid",
                                                                rλs=10.0 .^ (0:0),
                                                                nblock = -1, # wild bootstrap
                                                                kw...)::Tuple{T, MonoDecomp{T}} where T <: Real
    D, _ = cv_mono_decomp_ss(x, y; one_se_rule = one_se_rule, 
            one_se_rule_pre = one_se_rule_pre,
            rλs = rλs,
            nfold = nfold, seed = seed, method = md_method, kw...)
    error = y - D.yhat
    c = mean(D.yhat) / 2
    error = y - D.yhat
    tobs = var(D.γdown)
    ts = zeros(nrep)
    yhat = D.workspace.B * D.γup .+ c
    for i = 1:nrep
        yi = construct_bootstrap_y(y, error, yhat, nblock = nblock)
        try
            Di = mono_decomp_ss(D.workspace, x, yi, D.λ, D.μ, strict = true)
            ts[i] = var(Di.γdown)
        catch e
            @warn "due to error $e in optimization, assign test statistic as Inf"
            ts[i] = Inf
        end
    end
    pval = sum(ts .> tobs) / nrep
    return pval, D
end

function mono_test(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real
    # x, y = gen_data(a = 0.15, σ = 0.025)
    res, _ = cv_mono_decomp_cs(x, y, Js = 4:20, ss = 10.0 .^ (-6:0.5:-1), fixJ = false)
    # scatter(x, y)
    # plot!(x, res.yhat)
    c = mean(res.yhat) / 2
    n = length(y)
    σhat2 = var(y - res.yhat) / 4n
    tobs = sum((res.γdown .- c).^2 / σhat2)
    pval = 1 - cdf(Chisq(res.workspace.J), tobs)
    # println("pval = $pval, tobs = $tobs, critical = $(quantile(Chisq(J), 0.95))")
    return pval < 0.05
end

# summary experiment results
function summary_mono_test(resfile::String, task = "typeI_error")
    texname = resfile[1:end-4] * ".tex"
    res0 = deserialize(resfile)
    res = [x .< 0.05 for x in res0]
    methods = ["Meyer", "Ghosal", "Bowman",
            #    "MD (CS) sup", "MD (CS) min", "MD (CS) 1se", 
                "MDCS",
                "MDSS"
            #    "MD (SS) sup", "MD (SS) min", "MD (SS) 1se"
               ]
    nmethod = length(methods)
    μ = mean(res)
    @assert size(μ, 4) == nmethod
    A = Array{Matrix, 1}(undef, nmethod)
    name_σs = ["\$\\sigma = 0.001\$", "\$\\sigma = 0.01\$", "\$\\sigma = 0.1\$"]
    if task == "typeI_error"
        name_curves = ["\$x\$", "\$x^3\$", "\$x^{1/3}\$", "\$e^x\$", "\$1/(1+e^{-x})\$"]
    elseif task == "bowman"
        name_curves = "a = " .* string.([0,0.15,0.25,0.45])
    elseif task == "ghosal"
        name_curves = "m" .* string.(1:4)
    end
    for i in 1:nmethod
        if task == "bowman"
            A[i] = hcat([μ[:, :, j, i] for j = 1:3]...)
        else
            A[i] = hcat([μ[:, j, :, i]' for j = 1:3]...)
        end
    end
    print2tex(A, methods, name_σs, 
                    subcolnames=["n = 50", "100", "200"], 
                    subrownames = name_curves, 
                    colnames_of_rownames = ["Methods", "Curves"], 
                    file = texname, format="raw")
    if task == "typeI_error"
        mkpath(resfile[1:end-4])
        for i = 1:3
            for j = 1:3
                for k = 1:5
                    plot([histogram([x[i, j, k, ℓ] for x in res0], xlims = (0, 1.1), bins = 0:0.05:1.1, normalize = :probability, legend = false, title = methods[ℓ], margin=5Plots.mm) for ℓ = 1:5]..., layout = (1, 5), size = (2000, 400))
                    savefig(joinpath(resfile[1:end-4], "pval_$i$j$k.pdf"))
                end
            end
        end
    end
end
