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

function ghosal_S1n(X, Y, c, a = 0.05, b = 0.95, α = 0.05)
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

ghosal(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real = ghosal_S1n(x, y, C_GHOSAL[length(x)])

function meyer(x::AbstractVector{T}, y::AbstractVector{T}, nsim = 100, k = 6) where T <: Real
    meyer_rfile = joinpath(@__DIR__, "testmonotonicity.R")
    R"source($meyer_rfile)" # not found if put outside the function when wrapped into a package
    n = length(x)
    w = ones(n)
    return rcopy(R"montest($x, $y, 0, $w, $k, $nsim, 1)$pval") < 0.05
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
    return pval < 0.05
end

# run single experiment to make comparisons on bowan's curves
function single_test_compare_bowman(;
        as = [0, 0.15, 0.25, 0.45],
        ns = [50, 100, 200],
        # ns = [200]
        # ns = [100]
        # σs = [0.001, 0.01, 0.025, 0.05, 0.1] 
        σs = [0.001, 0.01, 0.1],
        only_proposed = false,
        nrep = 500, kw...
        # σs = [0.001, 0.1]     
    )
    # proposed, ghosal, meyer
    props = zeros(length(as), length(ns), length(σs), 9)
    for (j, a) in enumerate(as)
        for (k, n) in enumerate(ns)
            for (l, σ) in enumerate(σs)
                x, y = gen_data_bowman(a = a, n = n, σ = σ)
                # props[j, k, l, 1] = mono_test(x, y)
                # props[j, k, l, 1:5] .= mono_test_bootstrap(x, y, data_obey_H0 = false, one_se_rule=true)
                # props[j, k, l, 6:10] .= mono_test_bootstrap_ss(x, y, data_obey_H0 = false, one_se_rule=true)
                # props[j, k, l, 11:15] .= mono_test_bootstrap(x, y, data_obey_H0 = true, one_se_rule=true)
                # props[j, k, l, 16:20] .= mono_test_bootstrap_ss(x, y, data_obey_H0 = true, one_se_rule=true)
                # # props[j, k, l, 1] = mono_test_bootstrap(x, y)
                # # props[j, k, l, 2] = mono_test_bootstrap_ss(x, y) # might be too slow
                # props[j, k, l, 21] = meyer(x, y)
                # props[j, k, l, 22] = ghosal_S1n(x, y, C_GHOSAL[n])
                # props[j, k, l, 23] = bowman(x, y)
                if !only_proposed
                    props[j, k, l, 1] = meyer(x, y)
                    props[j, k, l, 2] = ghosal(x, y)
                    props[j, k, l, 3] = bowman(x, y)
                end
                props[j, k, l, 4:6] .= mono_test_bootstrap_sup(x, y; nrep = nrep, nμ = 5, kw...)
                props[j, k, l, 7:9] .= mono_test_bootstrap_supss(x, y; nrep = nrep, nμ = 5, kw...)
            end
        end
    end
    return props
end

# run single experiment to make comparisons on ghosal's curves
function single_test_compare_ghosal(;
        ns = [50, 100, 200],
        # ns = [200]
        σs = [0.001, 0.01, 0.1],
        # σs = [0.001]
        only_proposed = false,
        nrep = 500, kw...
    )
    # proposed, ghosal, meyer, sm
    props = zeros(length(ns), length(σs), 4, 9)
    for (i, n) in enumerate(ns)
        for (k, σ) in enumerate(σs)
            x, m1, m2, m3, m4 = gen_data_ghosal(n = n, σ = σ)
            for (j, y) in enumerate([m1, m2, m3, m4])
            # for (j, y) in enumerate([m1, m2])
                # props[i, j, 1] = mono_test(x, y)
                # props[i, k, j, 1:5] .= mono_test_bootstrap(x, y, data_obey_H0 = false, one_se_rule=true)
                # props[i, k, j, 6:10] .= mono_test_bootstrap_ss(x, y, data_obey_H0 = false, one_se_rule=true)
                # props[i, k, j, 11:15] .= mono_test_bootstrap(x, y, data_obey_H0 = true, one_se_rule=true)
                # props[i, k, j, 16:20] .= mono_test_bootstrap_ss(x, y, data_obey_H0 = true, one_se_rule=true)
                # pvals[j, k, l, 1] = mono_test_bootstrap(x, y) # might be too slow
                # props[i, k, j, 21] = meyer(x, y)
                # props[i, k, j, 22] = ghosal_S1n(x, y, C_GHOSAL[n])
                # props[i, k, j, 23] = bowman(x, y)
                # props[i, k, j, 24] = mono_test_bootstrap_sup(x, y)
                # props[i, k, j, 25] = mono_test_bootstrap_supss(x, y)
                if !only_proposed
                    props[i, k, j, 1] = meyer(x, y)
                    props[i, k, j, 2] = ghosal(x, y)
                    props[i, k, j, 3] = bowman(x, y)
                end
                props[i, k, j, 4:6] .= mono_test_bootstrap_sup(x, y; nrep = nrep, nμ = 5, kw...)
                props[i, k, j, 7:9] .= mono_test_bootstrap_supss(x, y; nrep = nrep, nμ = 5, kw...)
            end
        end
    end
    return props
end

function single_test_compare_desc(;
        # ns = [50, 100, 200]
        ns = [200],
        # σs = [0.025, 0.05, 0.1]
        # σs = 10 .^ (-5:0.5:-1)
        # σs = 10 .^ (-3:0.5:-1)
        σs = [0.001, 0.01, 0.1],
        nrep = 100
    )
    # proposed, ghosal, meyer, sm
    props = zeros(length(ns), length(σs), 5, 15)
    for (i, n) in enumerate(ns)
        for (k, σ) in enumerate(σs)
            x, m1, m2, m3, m4, m5 = gen_desc_data(n = n, σ = σ)
            for (j, y) in enumerate([m1, m2, m3, m4, m5])
                props[i, k, j, 1:3] = mono_test_bootstrap(x, y, data_obey_H0 = false, nrep = nrep)[1:3]
                props[i, k, j, 4:6] = mono_test_bootstrap_ss(x, y, data_obey_H0 = false, nrep = nrep)[1:3]
                props[i, k, j, 7:9] = mono_test_bootstrap(x, y, data_obey_H0 = true, nrep = nrep)[1:3]
                props[i, k, j, 10:12] = mono_test_bootstrap_ss(x, y, data_obey_H0 = true, nrep = nrep)[1:3]
                props[i, k, j, 13] = meyer(x, y)
                props[i, k, j, 14] = ghosal(x, y)
                props[i, k, j, 15] = bowman(x, y)
            end
        end
    end
    return props
end

function single_test_compare_mono(;
        ns = [50, 100, 200],
        # ns = [200]
        # σs = [0.025, 0.05, 0.1]
        # σs = 10 .^ (-5:0.5:-1)
        # σs = 10 .^ (-3:0.5:-1)
        σs = [0.001, 0.01, 0.1],
        only_proposed = false,
        equidistant = false,
        nrep = 500, kw...
    )
    # proposed, ghosal, meyer, sm
    props = zeros(length(ns), length(σs), 5, 9)
    for (i, n) in enumerate(ns)
        for (k, σ) in enumerate(σs)
            x, m1, m2, m3, m4, m5 = gen_mono_data(n = n, σ = σ, equidistant = equidistant)
            for (j, y) in enumerate([m1, m2, m3, m4, m5])
                # props[i, k, j, 1:5] = mono_test_bootstrap(x, y, data_obey_H0 = false, one_se_rule=true)
                # props[i, k, j, 6:10] = mono_test_bootstrap_ss(x, y, data_obey_H0 = false, one_se_rule=true)
                # props[i, k, j, 11:15] = mono_test_bootstrap(x, y, data_obey_H0 = true, one_se_rule=true)
                # props[i, k, j, 16:20] = mono_test_bootstrap_ss(x, y, data_obey_H0 = true, one_se_rule=true)
                # props[i, k, j, 21] = meyer(x, y)
                # props[i, k, j, 22] = ghosal_S1n(x, y, C_GHOSAL[n])
                # props[i, k, j, 23] = bowman(x, y)
                if !only_proposed
                    props[i, k, j, 1] = meyer(x, y)
                    props[i, k, j, 2] = ghosal(x, y)
                    props[i, k, j, 3] = bowman(x, y)
                end
                props[i, k, j, 4:6] .= mono_test_bootstrap_sup(x, y; nrep = nrep, nμ = 5, kw...)
                props[i, k, j, 7:9] .= mono_test_bootstrap_supss(x, y; nrep = nrep, nμ = 5, kw...)
            end
        end
    end
    return props
end

function mono_test_bootstrap(x::AbstractVector{T}, y::AbstractVector{T}; nrep = 100, debug = false, data_obey_H0 = false, one_se_rule = false) where T <: Real
    # x, y = gen_data(a = 0.15, σ = 0.025)
    # J, μ, res = cvMBspl(x, y, x, Js = 4:30, ss = 10.0 .^ (-6:0.5:-1))
    n = length(y)
    # J, μ, res = cvMBspl(x, y, x, Js = 4:30, ss = 10.0 .^ (-6:0.5:0), one_se_rule = true)
    res, μ, μs = cv_mono_decomp_cs(x, y, ss = 10.0 .^ (-6:0.5:6), one_se_rule = one_se_rule, fixJ = true)
    J = res.workspace.J
    # scatter(x, y)
    # plot!(x, res.yhat)
    c = mean(res.yhat) / 2
    error = y - res.yhat
    # σ = std(err)
    tobs = sum((res.γdown .- c).^2)
    # tobs = var(res.γdown)
    ts = zeros(nrep)
    ts2 = zeros(nrep)
    ts3 = zeros(nrep)
    ts4 = zeros(nrep)
    ts5 = zeros(nrep)
    tobs3 = tobs
    # tobs4 = maximum(res.γdown) - minimum(res.γdown)
    tobs4 = var(res.γdown)
    tobs5 = var(res.γdown) / var(y - res.yhat)
    # tobs_up = var(res.γup)
    tobs_up = sum((res.γup .- c).^2)
    ts_up = zeros(nrep)
    if debug
        γsdown = zeros(J, nrep+1)
        γsdown[:, 1] = res.γdown
    end
    for i = 1:nrep
        idx = sample(1:n, n)
        if data_obey_H0
            yi = res.workspace.B * res.γup .+ c + error[idx]
        else
            yi = res.yhat + error[idx]
        end
        # if tobs < tobs_up
        #     yi = res.B * res.γup .+ c + error[idx]
        # else
        #     yi = res.B * res.γdown .+ c + error[idx]
        # end
        yi = yi .- mean(yi) .+ mean(y) # Note that mean(yi) = mean(res.yhat) / 2 + c + mean(error[idx])
        # so it is equivalent to -mean(err[idx]) + mean(err[idx])
        # yi = res.B * res.γup .+ c + randn(n) * σ
        Di = mono_decomp_cs(x, yi, s = μ, s_is_μ = true, J = J)
        # _, _, Di = cvMBspl(x, y, x, Js = J, ss = range(0.1μ, 2μ, length=10))
        # ts[i] = sum((Di.γdown .- c).^2) / var(y - Di.yhat)
        # if data_obey_H0
        #     ts[i] = sum((Di.γdown - res.γdown ).^2)
        # else
        #     ts[i] = sum((Di.γdown .- mean(Di.yhat) / 2 - (res.γdown .- c) ).^2)
        # end
        ts[i] = sum((Di.γdown .- c ).^2)
        ts2[i] = sum((Di.γdown .- mean(Di.yhat) / 2 - (res.γdown .- c) ).^2)
        # ts3[i] = abs( sum((Di.γdown .- mean(Di.yhat) / 2).^2) - tobs3 )
        # ts4[i] = abs(maximum(Di.γdown) - minimum(Di.γdown) - tobs4)
        # ts5[i] = abs(var(Di.γdown) - tobs5)
        # ts3[i] = sum((Di.γdown .- mean(Di.yhat) / 2).^2) - tobs3
        # ts4[i] = maximum(Di.γdown) - minimum(Di.γdown) - tobs4
        # ts5[i] = var(Di.γdown) - tobs5
        ts3[i] = sum((Di.γdown .- mean(Di.yhat) / 2).^2)
        # ts4[i] = maximum(Di.γdown) - minimum(Di.γdown)
        ts4[i] = var(Di.γdown)
        ts5[i] = var(Di.γdown) / var(y - Di.yhat)
        # ts3[i] = sum((Di.γup + Di.γdown - res.γup .- mean(Di.yhat) / 2).^2)
        # ts4[i] = maximum(Di.γup + Di.γdown - res.γup) - minimum(Di.γup + Di.γdown - res.γup)
        # ts5[i] = var(Di.γup + Di.γdown - res.γup)
        # ts[i] = var(Di.γdown)

        # yi_up = error[idx]
        # _, _, Di = monoBspl(J)(x, yi_up, x, s = μ, s_is_μ = true)
        # ts_up[i] = var(Di.γup)
        # ts_up[i] = sum((Di.γup .- c).^2)
        # ti_up = sum((Di.γup .- mean(Di.yhat)/2 ).^2)
        # ti_down = sum((Di.γdown .- mean(Di.yhat)/2 ).^2)
        # if ti_down < ti_up
        #     if tobs < tobs_up
        #         ts2[i] = sum((Di.γdown .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        #     else
        #         ts2[i] = sum((Di.γdown .- mean(Di.yhat) / 2 - (res.γup .- c)).^2)
        #     end
        # else
        #     if tobs < tobs_up
        #         ts2[i] = sum((Di.γup .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        #     else
        #         ts2[i] = sum((Di.γup .- mean(Di.yhat) / 2 - (res.γup .- c)).^2)
        #     end
        # end
        if debug
            γsdown[:, i+1] = Di.γdown
        end
    end
    if debug
        return ts, ts_up, tobs, tobs_up, γsdown
    end
    # pval = sum( (ts .> tobs) .| (ts_up .< tobs_up) ) / nrep
    # pval = 1 - sum( (ts .< tobs) .* (ts_up .> tobs) ) / nrep
    # pval = 1 - sum( (ts .< tobs) .* (ts_up .> ts) ) / nrep
    # println("tobs = $tobs, tobs_up = $tobs_up")
    # if tobs < tobs_up
    #     pval = 1 - sum( (ts .< tobs) ) / nrep
    # else
    #     pval = 1 - sum( (ts_up .< tobs_up) ) / nrep
    # end
    # pval = 1 - sum( (min.(ts, ts_up) .< min(tobs, tobs_up) ) ) / nrep
    # pval = sum(ts2 .> min(tobs, tobs_up)) / nrep
    pval = sum(ts .> tobs) / nrep
    pval2 = sum(ts2 .> tobs) / nrep
    pval3 = sum(ts3 .> tobs3) / nrep
    pval4 = sum(ts4 .> tobs4) / nrep
    pval5 = sum(ts5 .> tobs5) / nrep
    return [pval < 0.05, pval2 < 0.05, pval3 < 0.05, pval4 < 0.05, pval5 < 0.05]
end

function mono_test_bootstrap_sup(x::AbstractVector{T}, y::AbstractVector{T}; nrep = 100, nμ = 10, nfold = 2, fixJ = true, kw...) where T <: AbstractFloat
    n = length(y)
    D1, μ0, μs0 = cv_mono_decomp_cs(x, y, ss = 10.0 .^ (-6:0.5:6), one_se_rule = true, fixJ = fixJ, nfold = nfold)
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
        err = norm(error)
        tobs = var(D.γdown) #/ var(y - D.yhat)
        ts = zeros(nrep)
        c = mean(D.yhat) / 2
        for i = 1:nrep
            idx = sample(1:n, n)
            yi = D.workspace.B * D.γup .+ c + error[idx]
            yi = yi .- mean(yi) .+ mean(y)
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

function mono_test_bootstrap_supss(x::AbstractVector{T}, y::AbstractVector{T}; 
                                    nrep = 100, nμ = 10, nfold = 2, seed = rand(UInt64),
                                    opstat::Function = var,
                                    md_method = "single_lambda",
                                    tol = 1e-4,
                                    nblock = 10,
                                    kw...
                                    ) where T <: Real
    # for block index
    idx = sortperm(x)
    x = x[idx]
    y = y[idx]
    n = length(y)
    res, μ0, μs0 = cv_mono_decomp_ss(x, y; one_se_rule = true, nfold = nfold, seed = seed, method = md_method, tol = tol, kw...)
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
        ts = zeros(nrep)
        c = mean(D.yhat) / 2
        tobs = opstat(D.γdown)
        # tobs = sum((D.γdown .- c).^2)
        for i = 1:nrep
            # idx = sample(1:n, n)
            idx = block_bootstrap_idx(n; nblock = nblock)
            yi = res.workspace.B * D.γup .+ c + error[idx]
            yi = yi .- mean(yi) .+ mean(y)
            Di = mono_decomp_ss(res.workspace, x, yi, res.λ, μ)
            # ts[i] = var(Di.γdown) / var(y - Di.yhat)
            ts[i] = opstat(Di.γdown)
            # ts[i] = sum((Di.γdown .- c).^2)
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

function mono_test_bootstrap_ss(x::AbstractVector{T}, y::AbstractVector{T}; nrep = 100, debug = false, data_obey_H0 = false, one_se_rule = false) where T <: Real
    n = length(y)
    res, _ = cv_mono_decomp_ss(x, y, one_se_rule = one_se_rule)
    # res, workspace = cvMSS(x, y)
    # J, μ, res = cvMBspl_fixJ(x, y, ss = 10.0 .^ (-6:0.2:0), one_se_rule = true, nfold = n)
    # scatter(x, y)
    # plot!(x, res.yhat)
    c = mean(res.yhat) / 2
    error = y - res.yhat
    σ = std(error)
    # println("c = $c, μup = ", mean(workspace.B * res.γup), "μdown = ", mean(workspace.B * res.γdown))
    # ϵ = mean(workspace.B * res.γup - workspace.B * res.γdown) / 2
    # tobs = sum((res.γdown .+ ϵ  .- c).^2)
    # tobs = sum((workspace.B * res.γdown .+ ϵ  .- c).^2)
    tobs = sum((res.γdown .- c).^2)
    # tobs = var(res.γdown)
    ts = zeros(nrep)
    ts2 = zeros(nrep)
    ts3 = zeros(nrep)
    ts4 = zeros(nrep)
    ts5 = zeros(nrep)
    tobs3 = tobs
    tobs4 = maximum(res.γdown) - minimum(res.γdown)
    tobs5 = var(res.γdown)
    # tobs5 = var(res.γdown) / var(y - res.yhat)
    # tobs_up = sum((res.γup .- ϵ .- c).^2)
    # tobs_up = sum((workspace.B * res.γup .- ϵ .- c).^2)
    tobs_up = sum((res.γup .- c).^2)
    println("λ = ", res.λ, " μ = ", res.μ, " t1 = ", tobs)
    # tobs_up = var(res.γup)
    ts_up = zeros(nrep)
    if debug
        γsdown = zeros(length(res.γdown), nrep+1)
        γsdown[:, 1] = res.γdown
    end
    for i = 1:nrep
        idx = sample(1:n, n)
        if data_obey_H0
            yi = res.workspace.B * res.γup .+ c + error[idx]
        else
            yi = res.yhat + error[idx]
        end
        # # if tobs < tobs_up
        # yi = workspace.B * res.γup .+ c + error[idx]
        # if var(res.γdown) < var(res.γup)
        # # if rand() < 0.5
        #     yi = workspace.B * res.γup .+ c + error[idx]
        # else
        #     yi = workspace.B * res.γdown .+ c + error[idx]
        # end
        yi = yi .- mean(yi) .+ mean(y) # Note that mean(yi) = mean(res.yhat) / 2 + c + mean(error[idx])
        # so it is equivalent to -mean(err[idx]) + mean(err[idx])
        # yi = res.B * res.γup .+ c + randn(n) * σ
        Di = mono_decomp_ss(res.workspace, x, yi, res.λ, res.μ)
        # ϵi = mean(workspace.B * Di.γup - workspace.B * Di.γdown) / 2
        # _, _, Di = cvMBspl(x, y, x, Js = J, ss = range(0.1μ, 2μ, length=10))
        # ts[i] = sum((Di.γdown .- c).^2) / var(y - Di.yhat)
        # ts[i] = sum((Di.γdown .- c).^2)
        # ts[i] = sum((Di.γdown .+ ϵi .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        ts[i] = sum((Di.γdown .- c).^2)
        ts2[i] = sum((Di.γdown .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        # ts3[i] = sum((Di.γdown - res.γdown ).^2)
        # ts3[i] = abs( sum((Di.γdown .- mean(Di.yhat) / 2).^2) - tobs3 )
        # ts4[i] = abs(maximum(Di.γdown) - minimum(Di.γdown) - tobs4)
        # ts5[i] = abs(var(Di.γdown) - tobs5)
        ts3[i] = sum((Di.γdown .- mean(Di.yhat) / 2).^2)
        ts4[i] = maximum(Di.γdown) - minimum(Di.γdown)
        # ts4[i] = var(Di.γdown) 
        ts5[i] = var(Di.γdown) #/ var(y - Di.yhat)
        # ts3[i] = sum((Di.γhat - res.γup .- mean(Di.yhat) / 2).^2)
        # ts4[i] = maximum(Di.γhat - res.γup) - minimum(Di.γhat - res.γup)
        # ts5[i] = var(Di.γhat - res.γup)
        # ts3[i] = sum((Di.γdown .- mean(Di.yhat) / 2).^2) - tobs3
        # ts4[i] = maximum(Di.γdown) - minimum(Di.γdown) - tobs4
        # ts5[i] = abs(var(Di.γdown) - tobs5)
        # if data_obey_H0
        #     ts[i] = sum((Di.γdown .- c).^2)
        # else
        #     ts[i] = sum((Di.γdown .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        # end
        # ts[i] = sum((Di.γdown .- mean(Di.γdown) - (res.γdown .- mean(res.γdown))).^2)
        # ts[i] = var(Di.γdown)

        # yi_up = error[idx]
        # _, _, Di = monoBspl(J)(x, yi_up, x, s = μ, s_is_μ = true)
        # ts_up[i] = sum((Di.γup .- c).^2)
        # ts_up[i] = sum((Di.γup .- ϵi .- mean(Di.yhat) / 2 - (res.γup .- c)).^2)
        # ts_up[i] = sum((Di.γup .- mean(Di.yhat) / 2 - (res.γup .- c)).^2)
        # ts_up[i] = sum((Di.γup .- mean(Di.γup) - (res.γup .- mean(res.γup))).^2)
        # ts_up[i] = var(Di.γdown)
        # if var(Di.γdown) < var(Di.γup) # !!!NB: the variance is not exactly the same as ‖ γ - c ‖ since c is the mean of yhat instead of the mean of γ
        ## TODO: discuss the difference between var and c
        # ti_up = sum((Di.γup .- mean(Di.yhat)/2 ).^2)
        # ti_down = sum((Di.γdown .- mean(Di.yhat)/2 ).^2)
        # ti_up = sum((workspace.B * Di.γup .- ϵi .- mean(Di.yhat)/2 ).^2)
        # ti_down = sum((workspace.B * Di.γdown .+ ϵi .- mean(Di.yhat)/2 ).^2)
        # ti_up = var(Di.γup)
        # ti_down = var(Di.γdown)
        # ts2[i] = min(var(Di.γup), var(Di.γdown))
        # if ti_down < ti_up
        # if var(Di.γdown) < var(Di.γup)
        #     # if tobs < tobs_up
        #     if var(res.γdown) < var(res.γup)
        #         # ts2[i] = sum((Di.γdown .+ ϵi .- mean(Di.yhat) / 2 - (res.γdown .+ ϵ .- c)).^2)
        #         ts2[i] = sum((Di.γdown .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        #         # ts2[i] = sum((Di.γhat - res.γup .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        #         # ts2[i] = sum((workspace.B * Di.γdown .+ ϵi .- mean(Di.yhat) / 2 - (workspace.B * res.γdown .+ ϵ .- c)).^2)
        #         # ts2[i] = sum((Di.γdown - res.γdown).^2)
        #         # ts2[i] = sum((Di.γdown .- mean(Di.γdown) - (res.γdown .- mean(res.γdown))).^2)
        #     else
        #         # ts2[i] = sum((Di.γdown .+ ϵi .- mean(Di.yhat) / 2 - (res.γup .- ϵ .- c)).^2)
        #         ts2[i] = sum((Di.γdown .- mean(Di.yhat) / 2 - (res.γup .- c)).^2)
        #         # ts2[i] = sum((Di.γhat - res.γup .- mean(Di.yhat) / 2 - (res.γup .- c)).^2)
        #         # ts2[i] = sum((workspace.B * Di.γdown .+ ϵi .- mean(Di.yhat) / 2 - (workspace.B * res.γup .- ϵ .- c)).^2)
        #         # ts2[i] = sum((Di.γdown - res.γup).^2)
        #         # ts2[i] = sum((Di.γdown .- mean(Di.γdown) - (res.γup .- mean(res.γup))).^2)
        #     end
        # else
        #     if tobs < tobs_up
        #         # ts2[i] = sum((Di.γup .- ϵi .- mean(Di.yhat) / 2 - (res.γdown .+ ϵ .- c)).^2)
        #         ts2[i] = sum((Di.γup .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        #         # ts2[i] = sum((Di.γhat - res.γdown .- mean(Di.yhat) / 2 - (res.γdown .- c)).^2)
        #         # ts2[i] = sum((workspace.B * Di.γup .- ϵi .- mean(Di.yhat) / 2 - (workspace.B * res.γdown .+ ϵ .- c)).^2)
        #         # ts2[i] = sum((Di.γup - res.γdown).^2)
        #         # ts2[i] = sum((Di.γup .- mean(Di.γup) - (res.γdown .- mean(res.γdown))).^2)
        #     else
        #         # ts2[i] = sum((Di.γup .- ϵi .- mean(Di.yhat) / 2 - (res.γup .- ϵ .- c)).^2)
        #         ts2[i] = sum((Di.γup .- mean(Di.yhat) / 2 - (res.γup .- c)).^2)
        #         # ts2[i] = sum((Di.γhat - res.γdown .- mean(Di.yhat) / 2 - (res.γup .- c)).^2)
        #         # ts2[i] = sum((workspace.B * Di.γup .- ϵi .- mean(Di.yhat) / 2 - (workspace.B * res.γup .- ϵ .- c)).^2)
        #         # ts2[i] = sum((Di.γup - res.γup).^2)
        #         # ts2[i] = sum((Di.γup .- mean(Di.yup) - (res.γup .- mean(res.γup))).^2)
        #     end
        # end
        if debug
            γsdown[:, i+1] = Di.γdown
        end
    end
    if debug
        return γsdown
    end
    # pval = sum( (ts .> tobs) .| (ts_up .< tobs_up) ) / nrep
    # pval = 1 - sum( (ts .< tobs) .* (ts_up .> tobs) ) / nrep
    # pval = 1 - sum( (ts .< tobs) .* (ts_up .> ts) ) / nrep
    # if tobs < tobs_up
    #     pval = 1 - sum( (ts .< tobs) ) / nrep
    # else
    #     pval = 1 - sum( (ts_up .< tobs_up) ) / nrep
    # end
    # return ts2, tobs, tobs_up
    # pval = 1 - sum( (min.(ts, ts_up) .< min(tobs, tobs_up) ) ) / nrep
    # pval = sum(ts2 .> min(tobs, tobs_up)) / nrep
    # pval1 = sum(ts_up .> tobs_up) / nrep
    # return ts, ts2, ts3, ts4, ts5, tobs3, tobs4, tobs5
    pval = sum(ts .> tobs) / nrep
    pval2 = sum(ts2 .> tobs) / nrep
    pval3 = sum(ts3 .> tobs3) / nrep
    pval4 = sum(ts4 .> tobs4) / nrep
    pval5 = sum(ts5 .> tobs5) / nrep
    # return pval, pval1, pval2
    return [pval < 0.05, pval2 < 0.05, pval3 < 0.05, pval4 < 0.05, pval5 < 0.05]
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
    res = deserialize(resfile)
    methods = ["Meyer", "Ghosal", "Bowman",
               "MD (CS) sup", "MD (CS) min", "MD (CS) 1se", 
               "MD (SS) sup", "MD (SS) min", "MD (SS) 1se"]
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
end
