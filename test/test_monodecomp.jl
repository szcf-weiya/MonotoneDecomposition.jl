module TestMonoDecomp 

import Test: Test, finish
# Do not omit the following two lines, otherwise, it throws no such function.
using Test: DefaultTestSet, Broken
using Test: parse_testset_args
using Test
using Plots
using Random
using RCall
using LinearAlgebra
using StatsBase
using ProgressMeter
using MonotoneDecomposition

ecos()

"""
Skip a testset: (https://discourse.julialang.org/t/skipping-a-whole-testset/65006/4)

Use `@testset_skip` to replace `@testset` for some tests which should be skipped.

Usage
-----
Replace `@testset` with `@testset "reason"` where `"reason"` is a string saying why the
test should be skipped (which should come before the description string, if that is
present).
"""
macro testset_skip(args...)
    isempty(args) && error("No arguments to @testset_skip")
    length(args) < 2 && error("First argument to @testset_skip giving reason for "
                              * "skipping is required")

    skip_reason = args[1]

    desc, testsettype, options = parse_testset_args(args[2:end-1])

    ex = quote
        # record the reason for the skip in the description, and mark the tests as
        # broken, but don't run tests
        local ts = DefaultTestSet(string($desc, " - ", $skip_reason))
        push!(ts.results, Broken(:skipped, "skipped tests"))
        local ret = finish(ts)
        ret
    end

    return ex
end


@testset "mean(yup) == mean(ydown)" begin
    for σ in [0.001, 0.01, 0.1]
        x, y, x0, y0 = gen_data(100, σ, "SE_0.1")
        res, _ = cv_mono_decomp_ss(x, y)
        yup, ydown = predict(res.workspace, x0, res.γup, res.γdown)
        d = abs(mean(yup) - mean(ydown))
        println("σ = $σ, MDSS: | mean(yup) - mean(ydown) | = ", d)
        # @test d < cbrt(eps())

        res, _ = cv_mono_decomp_cs(x, y)
        yup, ydown = predict(res.workspace, x0, res.γup, res.γdown)
        d = abs(mean(yup) - mean(ydown))
        println("σ = $σ, MDCS: | mean(yup) - mean(ydown) | = ", d)
        # @test d < cbrt(eps())
    end
end

@testset "one se rule" begin
    x, y, x0, y0 = gen_data(100, 0.001, "SE_0.1")
    res = cv_mono_decomp_ss(x, y, one_se_rule = true, figname = "/tmp/p.png")
end

@testset "gcv in ss" begin
    f = x -> x^3
    n = 50
    σ = 0.5
    x, y, xnew = gen_data(n, σ, f)
    # ss
    yhat, yhatnew, Ω, λopt, spl = smooth_spline(x, y, x)

    λs = range(λopt/10, 2λopt, length = 10)
    gcvs = Float64[]
    rgcvs = Float64[]
    for λ in λs
        spl = R"smooth.spline($x, $y, lambda = $λ)"
        yhat = rcopy(R"predict($spl, $x)$y")
        knots = rcopy(R"$spl$fit$knot")[4:end-3]
        bbasis = R"fda::create.bspline.basis(breaks = $knots, norder = 4)"
        B = rcopy(R"predict($bbasis, $knots)")
        Ω = rcopy(R"fda::eval.penalty($bbasis, 2)")
        gcv = norm(yhat - y)^2 / (1 - tr(B * inv(B' * B + λ * Ω) * B')/n )^2
        rgcv = rcopy(R"$spl$cv.crit")
        append!(gcvs, gcv)
        append!(rgcvs, rgcv)
    end
    fig1 = plot(λs, gcvs)
    vline!(fig1, [λopt])
    fig2 = plot(λs, rgcvs)
    vline!(fig2, [λopt])
    plot(fig1, fig2, layout = (2, 1))
end

@testset "shrinkage on ss" begin
    f = x -> x^3
    # f = exp
    # f = x -> 1/(1+exp(-x))
    # f = x -> x^2
    # f = x -> sin(2π*x)
    n = 100
    ks = 0.5:0.01:1.5
    nrep = 10
    σs = 0.1:0.1:2.0
    nσ = length(σs)
    res = zeros(length(ks), nrep, nσ)
    @showprogress for j = 1:nrep
        for l = 1:nσ
            σ = σs[l]
            x, y, xnew = gen_data(n, σ, f)
            yhat, yhatnew, Ω, λopt = smooth_spline(x, y, xnew)
            for (i, k) in enumerate(ks)
                res[i, j, l] = norm(yhatnew * k - f.(xnew))
            end
        end
    end
    μres = mean(res, dims = 2)[:, 1, :]
    figs = Plots.Plot[]
    for i = 1:nσ
        σ = σs[i]
        plt = plot(ks, μres[:, i], title = "shrinkage (σ = $σ)")
        vline!(plt, [1.0])
        push!(figs, plt)
    end
    # save_grid_plots(figs)
end

@testset "cross-validation for monotone decomposition with smoothing splines" begin
    f = x -> x^3
    n = 100
    σ = 0.1
    x, y, x0, y0 = gen_data(n, σ, f)
    yhat, yhatnew, Ω, λopt = smooth_spline(x, y, x0)
    @testset "fix λ as λopt in ss" begin
        D, _ = cvfit(x, y, 0.05:0.05:0.99, λopt, figname = nothing)
        [norm(predict(D, x0) - y0), norm(yhatnew - y0)]
        # plot([x, y], [x0, y0], workspace, D, yhatnew, digits = 5)
        all(yhat[2:end] - yhat[1:end-1] .>= 0)
    end

    @testset "set μmax" begin
        D, _ = cvfit(x, y, 1.0, range(λopt/2, λopt*2, length = 3), nμ = 10, figname = nothing)
        @test all(D.γup[2:end] - D.γup[1:end-1] .>= 0)
    end 

    @testset "divide search region" begin
        D, _ = cvfit(x, y, 10.0, range(λopt/2, λopt*2, length = 3), nμ = 10, figname = nothing, ρ = 0.11)
        @test all(D.γup[2:end] - D.γup[1:end-1] .>= 0)        
    end

    @testset "extend search region" begin
        D, _ = cvfit(x, y, 0.01, range(λopt/2, λopt*2, length = 3), nμ = 10, figname = nothing, ρ = 0.11)
        @test all(D.γup[2:end] - D.γup[1:end-1] .>= 0)
    end

    @testset "fix μ, vary λ" begin
        D, _ = cvfit(x, y, 0.1, λopt, nλ = 10, figname = nothing)
        @test all(D.γup[2:end] - D.γup[1:end-1] .>= 0)
    end

    @testset "extend search region" begin
        D, _ = cvfit(x, y, 0.0, λopt*0.9, nλ = 10, figname = nothing, rλ = 0.1, ρ = 0.05)
        @test all(D.γup[2:end] - D.γup[1:end-1] .>= 0)
    end
end

@testset "fix lambda in ss" begin
    @testset "y = x^2" begin
        f = x -> x^2
        n = 100
        σ = 0.5
        x = sort(rand(n) * 2 .- 1)
        y = f.(x) + σ * randn(n)
            
        yhat, yhatnew, Ω, λ = smooth_spline(x, y, x, keep_stuff = false)
        err_ss = norm(yhat - f.(x))
        err_md = Float64[]
        # !! Ω + I can avoid positive definite issue, and interesting, it can let err_md achieves err_ss
        # ss = 2.0 .^ (-4:10)
        s0 = norm(yhat)
        ss = range(s0/4, s0, length = 10)
        for s in ss
            res = mono_decomp_ss(x, y, λ = λ, s = s, s_is_μ = false)
            append!(err_md, norm(res.yhat - f.(x)))
        end
        @test (minimum(err_md) - err_ss) / err_ss < 0.15 # negative is even better
        # plot(log2.(ss), err_md, xlab = "2^s")
        # hline!(log2.(ss), [err_ss])
        # Omega is for original ss, if in cross validation, the dimension does not match.
        # fig, obj_val, yup, ydown, γup, γdown, df, J = cvfit(x, y, ss = 0.5:0.5:5, λ = λ, Ω = Ω + I, nfold = 10, one_se_rule = false)
        # fig, obj_val, yup, ydown, γup, γdown, df, J = cvfit(x, y, ss = 0.5:0.5:5, λ = λ, nfold = 10, one_se_rule = false, figname = "/tmp/optim_cv.png")
    end

    @testset "y = x^3" begin
        f = x -> x^3
        n = 50
        σ = 0.5
        x, y, xnew = gen_data(n, σ, f)
        # ss
        yhat, yhatnew, Ω, λ = smooth_spline(x, y, x)
        figs = Plots.Plot[]
        errs = Float64[]
        crit = Float64[]
        crit1 = Float64[]
        crit2 = Float64[]
        # for s = 1000:1000:10000
        s0 = norm(yhat)
        tol = sqrt(eps())
        ss = range(s0/4, 2s0, length = 10)
        for s in ss
            res = mono_decomp_ss(x, y, λ = λ, s = s, s_is_μ = false)
            # push!(figs, fig)
            append!(errs, norm(res.yhat - f.(x)))
            rss = norm(res.yhat - y)^2
            trS = tr(res.workspace.B * inv(res.workspace.B' * res.workspace.B + λ * res.workspace.L * res.workspace.L') * res.workspace.B')
            dfup = res.workspace.J - sum(abs.(diff(res.γup) .< tol))
            dfdown = res.workspace.J - sum(abs.(diff(res.γdown) .< tol))
            gcv = log(rss) - 2log(1-trS/n)
            df = (dfup + dfdown) / 2
            append!(crit1, gcv)
            append!(crit2, 2df / n)
            append!(crit, gcv + 2df / n)
        end
        # save_grid_plots(figs, 1, 10)
        err = norm(f.(x) - yhat)
        findfirst(errs .< err)
        @test (minimum(errs) - err) / err < 0.1
    end
end

@testset "verify proposition (y = x^3)" begin
    f = x -> x^3
    n = 100
    σ = 0.1
    J = 10
    x, y, xnew = gen_data(n, σ, f)
    @testset "sum MSE(x_i) = Jσ^2" begin
        mses = zeros(100)
        for i = 1:100
            # x, y, xnew = gen_data(n, σ, f)
            y = f.(x) + randn(n) * σ # fix x
            yhat, _, spl = cubic_spline(J)(x, y, x)
            mses[i] = norm(yhat - x.^3)^2
        end
        # mean(mses) ≈ σ^2 * J
        @test abs(mean(mses) - σ^2 * J) / (J * σ^2) < 0.15
    end
    @testset "check solution" begin
        x, y, xnew = gen_data(n, σ, x -> x + 3)
        # if J = 4, is always D.γup[2] == D.γup[3]??
        μ = 0.01
        D = mono_decomp_cs(x, y, s = μ, J = 5, s_is_μ = true)
        γlm = inv(D.workspace.B' * D.workspace.B) * D.workspace.B' * y
        γlm / (μ + 1) + (μ - 1) / (μ + 1) * D.γdown
        @testset "γ1 < γ2 < ... < γk = ... = γℓ = ... = γJ" begin
            x, y, xnew = gen_data(100, 0.01, x -> x^3)
            μ = 0.01
            D = mono_decomp_cs(x, y, s = μ, J = 5, s_is_μ = true)
            G = zeros(4, 5) # γ1 < γ2 = γ3 < γ4 < γ5
            G[1, 1] = 1
            G[2, 2] = G[2, 3] = 1
            G[3, 4] = G[4, 5] = 1
            inv(G * D.workspace.B' * D.workspace.B * G') * G * D.workspace.B' * y / (μ + 1) .+ (μ - 1)/(μ+1) * D.γdown[1]
        end

        @testset "multiple equal" begin
            # γ1 < ... < γk1 = ... = γk2 < ... < γk3 = ... = γk4 < ... < γJ
            # x, y, xnew = gen_data(100, 0.01, x -> sin(2π * x))
            g(x, a = 3) = ((x + 1/a) ^3 - 1/a^3) * (x < 0) + ((x - 1/a) ^3 + 1/a^3) * (x > 0)
            x, y, xnew = gen_data(100, 0.001, x -> g(x, 2))
            μ = 0.01
            D = mono_decomp_cs(x, y, s = μ, J = 20, s_is_μ = true)
        end
        # k1 = 5, k2 = 8, k3 = 13, k4 = 16
        G = zeros(20 - 3 - 3, 20)
        for i = 1:4
            G[i, i] = 1
        end
        G[5, 5:8] .= 1#1/4
        for i = 9:12
            G[i-3, i] = 1
        end
        G[10, 13:16] .= 1#1/4
        for i = 17:20
            G[i-6, i] = 1
        end
        γup1 = G' * inv(G * D.workspace.B' * D.workspace.B * G') * G * D.workspace.B' * y / (μ + 1) .+ (μ - 1)/(μ+1) * D.γdown[1]
    end
    @testset "check shrinkage" begin
        for k = 0.1:0.1:0.9
            i = 0
            while true
                i += 1
                # println("i = $i")
                x, y, xnew = gen_data(n, σ, f)
                yhat, _, spl = cubic_spline(J)(x, y, x)
                B = rcopy(spl.H)
                γ = spl.β
                if all(γ[2:end] - γ[1:end-1] .> 0)
                    # since the condition has been normalized by norm(B)
                    res = mono_decomp_cs(x, y, J = J, s = k * norm(B * γ), s_is_μ = false)
                    ratio = res.γup ./ spl.β
                    println("k = $k, ratio = $ratio")
                    # @test abs.(mean(ratio) - k) / k < 0.1
                    @test abs.(ratio[1] - k) / k < 0.3
                    @test abs.(ratio[end] - k) / k < 0.3
                    break
                else
                    continue
                end
            end
        end
    end
    @testset "MSE(x_0) individually" begin
        yhat, _, spl = cubic_spline(J)(x, y, x)
        B = rcopy(spl.H)
        γ = spl.β
    
        # all(γ[2:end] - γ[1:end-1] .> 0)
        invBB = inv(B' * B)
        sopts = zeros(n)
        lows = zeros(n)
        for i = 1:n
            a = (γ' * B[i, :])^2 
            b = σ^2 * B[i, :]' * invBB * B[i, :]
            sopts[i] = a / (a + b)
            lows[i] = (a - b) / (a + b)
        end
        sopts *= norm(B * γ)
        lows *= norm(B * γ)            
    end

    @testset "MSE(̃f) < MSE(̂f)" begin
        err_md = zeros(100)
        err_ss = zeros(100)
        a = sum(f.(x) .^ 2)
        b = J * σ^2
        kstar = a / (a + b)
        @testset "check optimal MSE" begin
            for i = 1:100
                y = f.(x) + randn(n) * σ # fix x
                yhat, _, spl = cubic_spline(J)(x, y, x)
                B = rcopy(spl.H)
                γ = spl.β
                sstar = kstar * norm(B * γ)    
                res = mono_decomp_cs(x, y, J = J, s = sstar, s_is_μ = false)
                err_md[i] = norm(res.yhat - f.(x))
                err_ss[i] = norm(yhat - f.(x))
            end
            scatter(err_ss, err_md, xlab = "Bspl", ylab = "monodecomp")
            Plots.abline!(1, 0)
            @test mean(err_md) < mean(err_ss)                    
        end
        kstar1 = (a - b) / (a + b)
        @testset "check different ks" begin
            ks = range(kstar1 - 10(1-kstar1), 1, length = 10)
            flags = zeros(Bool, 10)
            for (j, k) in enumerate(ks)
                println("k = $k")
                for i = 1:100
                    y = f.(x) + randn(n) * σ # fix x
                    yhat, _, spl = cubic_spline(J)(x, y, x)
                    B = rcopy(spl.H)
                    γ = spl.β
                    res = mono_decomp_cs(x, y, J = J, s = k * norm(B * γ) / norm(B), s_is_μ = false)
                    err_md[i] = norm(res.yhat - f.(x))
                    err_ss[i] = norm(yhat - f.(x))
                end
                flags[j] = mean(err_md) < mean(err_ss)
            end
            # it turns out when k smaller than kstar1, it can also have err_md < err_ss
        end
    end
    @testset "fix lambda" begin
        @testset "y = x^3" begin
            f = x -> x^3
            n = 100
            σ = 0.5
            x, y, xnew = gen_data(n, σ, f)
            yhat, yhatnew, Ω, λ = smooth_spline(x, y, x)
            # ks = 0.5:0.05:1
            ks = 0.9:0.01:1
            # ks = 0.99:0.001:1
            nk = length(ks)
            err_ss = norm(yhat - f.(x))
            err_md = zeros(nk, 2)
            for (j, k) in enumerate(ks)
                s = norm(yhat) * k
                res1 = mono_decomp_ss(x, y, λ = λ / k, s = s, s_is_μ = false)
                res2 = mono_decomp_ss(x, y, λ = λ, s = s, s_is_μ = false)
                err_md[j, 1] = norm(res1.yhat - f.(x))
                err_md[j, 2] = norm(res2.yhat - f.(x))
            end
        end
    end
end
function sol_c(y, yhat, k; nc = 100)
    err = zeros(nc)
    cs = range(minimum(yhat), maximum(yhat), length = nc)
    for (i, c) in enumerate(cs)
        err[i] = abs.(mean(y - k * (max.(yhat .- c, 0) .+ c)) - (1-k)c)
    end
    ind = argmin(err)
    return cs[ind], err
end

@testset "null space of B" begin
    f = x -> x^3
    n = 100
    σ = 0.1
    J = 20#10
    x, y, xnew = gen_data(n, σ, f)
    yhat, _, spl = cubic_spline(J)(x, y, x)
    B = rcopy(spl.H)

    B1 = zeros(J-1, n)
    B1[1:J-2, :] .= B[:, 1:J-2]'
    B1[J-1, :] .= B[:, J-1] + B[:, J]
    B1 * B
    γhat = inv(B' * B) * B' * y
    l = nullspace(B1 * B)[:]
    λ = (γhat[end-1] - γhat[end]) / (l[end-1] - l[end])
    γd = (γhat - λ * l)/2

    c = sum([B[:, i]' * B * γhat for i = 1:J]) / 2n
end

@testset "GI solver for MDCS" begin
    tol = 1e-4
    x, y, _ = gen_data(100, 0.1, x->x^3)
    D1 = mono_decomp_cs(x, y, s = 1.0, J = 10, use_GI = true)
    D2 = mono_decomp_cs(x, y, s = 1.0, J = 10, use_GI = false)
    @test all(abs.(D1.γdown - D2.γdown) .< tol)
    @test all(abs.(D1.γup - D2.γup) .< tol)    

    cvD1 = cv_mono_decomp_cs(x, y, ss = [0.1, 1.0], Js = 10:10, use_GI = true, nfold = 100) # LOOCV avoid cv randomness
    cvD2 = cv_mono_decomp_cs(x, y, ss = [0.1, 1.0], Js = 10:10, use_GI = false, nfold = 100)
    @test all(abs.(cvD1[3] - cvD2[3]) .< tol)
    @test all(abs.(cvD1[4] - cvD2[4]) .< tol)
    @test all(abs.(cvD1[1].γdown -cvD2[1].γdown) .< tol)
    @test all(abs.(cvD1[1].γup -cvD2[1].γup) .< tol)
end

@testset "verify proposition (y = x^2)" begin
    f = x -> x^2
    n = 100
    σ = 0.1
    J = 20#10
    @testset "check shrinkage" begin
        for k = 0.5:0.05:0.9
            x, y, xnew = gen_data(n, σ, f)
            yhat, _, spl = cubic_spline(J)(x, y, x)
            B = rcopy(spl.H)
            γ = spl.β
            γ1, γ2 = mono_decomp(γ)
            c, _ = sol_c(y, yhat, k)
            # s = k * norm(yhat .- c)
            s = k * norm(max.(yhat .- c, 0))
            res = mono_decomp_cs(x, y, J = J, s = s, s_is_μ = false)
            # theoretical
            γ_theo = k * max.(γ .- c, 0) .+ c/2
            γ_pred = max.(res.γup, res.γdown)
            err = norm(γ_theo - γ_pred) / norm(γ_theo)
            println("k = $k, err = $err")
        end
    end
    @testset "verify approximation" begin
        k = 0.9
        nrep = 100
        γs = zeros(J, nrep)
        μ_theo = zeros(J, nrep)
        for i = 1:nrep
            x, y, xnew = gen_data(n, σ, f)
            yhat, _, spl = cubic_spline(J)(x, y, x)
            B = rcopy(spl.H)
            γ = spl.β
            c, _ = sol_c(y, yhat, k)
            # s = k * norm(yhat .- c)
            s = k * norm(max.(yhat .- c, 0))
            res = mono_decomp_cs(x, y, J = J, s = s, s_is_μ = false)
            γs[:, i] = max.(res.γup, res.γdown)
            Ψ = 1 ./ γs[:, i]
            Ψ[Ψ .== Inf] .= 0
            W = 1 ./ (1 .+ Ψ * c)
            μ_theo[:, i] = k * W .* γ .+ c
            # σ_theo[i] = σ^2 * diagm(W) * inv(B' * B) * diagm(W)
        end
    end
end
@testset "explore" begin
    @testset "y = x^2 check the coef" begin
        f = x -> x^2
        n = 100
        σ = 0.1
        x = sort(rand(n) * 2 .- 1)
        y = f.(x) + σ * randn(n)
        J = 10
        yhat, _, spl = cubic_spline(J)(x, y, x)
        B = rcopy(spl.H)
        γ = spl.β
        plot(γ) # check if unique peak
        # the optimal s at norm(f(x)) is also useful
        k = 0.7
        γ1, γ2 = mono_decomp(γ)
        # correct γ1 and γ2 (specifically)
        # c = (γ1[1] + γ2[end])
        # c = 0.13754534234793778
        # c = 0.2 
        # treat it as a tuning parameter
        # c = mean(y - k * yhat) / (1-k)
        c, _ = sol_c(y, yhat, k)
        # ~~s = k * norm(yhat .- c)~~
        # ~~seems that B(B'B)^{-1}B1 = 1, then it equals to k*norm(yhat-c)~~
        # ~~s = norm(k * yhat + (1-k) * B * inv(B' * B) * B' * ones(n) * c .- c) ~~
        # it seems that B * ones(J) = ones(J), then it agains equal to k*norm(yhat.-c)
        s = k * norm(max.(yhat .- c, 0))
        res = mono_decomp_cs(x, y, J = J, s = s, s_is_μ = false)
        γ1 .-= (γ1[1] - c/2)
        γ2 .-= (γ2[end] - c/2)
        # ~~offset = (1-k) / 2 * inv(B' * B) * B' * ones(n) * c~~
        offset = (1-k)/2 * c
        γup2 = k * γ1 .+ offset
        γdown2 = k * γ2 .+ offset
        γ_theo = k * max.(γ .- c, 0) .+ c/2
        γ_pred = max.(res.γup, res.γdown)
        err = norm(γ_theo - γ_pred) / norm(γ_theo)
        plot(γ, label = "γhat", legend = :top, lw = 1.5, title = "err = $err")
        # plot!(γ1, label = "γ1")
        # plot!(γ2, label = "γ2")
        plot!(res.γup, label = "γup (k=$k)", ls = :dash, lw = 2)
        plot!(res.γdown, label = "γdown (k=$k)", ls = :dash, lw = 2)
        # plot!(γup2, label = "γup2 (k=$k)", lw = 2)
        # plot!(γdown2, label = "γdown2 (k=$k)", lw = 2)
        # plot!(max.(γup, γdown), label = "max(γup, γdown)", lw = 1)
        # plot!(k * γ .+ (1/2 -k)c, label = "γu or γd", lw = 2)
        plot!(k * max.(γ .- c, 0) .+ c/2, label = "γu or γd", lw = 2)
    end
end    

@testset "benchmarking experiments" begin
    benchmarking(jplot = false, nrep=1, σs = [0.1], competitor = "ss_single_lambda", nfold = 2, one_se_rule=false, resfolder = "/tmp")
    MonotoneDecomposition.summary(σs=[0.1], curves=["x^3"], nrep=1, resfolder = "/tmp")
    @test isfile("/tmp/tmp.tex")
    @testset "cubic splines" begin
        res = benchmarking_cs(fixJ = false, figname_cv="/tmp/cv.png", figname_fit="/tmp/fit.png", dataseed = 1234)
        @test isfile("/tmp/cv_Jmu.sil")
        savefig(cvplot("/tmp/cv_Jmu.sil", "bspl2"), "/tmp/cv2.png")
        @test isfile("/tmp/cv2.png")
        
        @test res[2] < res[4]
        
        res = benchmarking_cs(fixJ = true, dataseed = 1234)
        @test res[2] < res[4]
    end
    @testset "smoothing splines" begin
        for one_se_rule in [false, true]
            for method in ["single_lambda", "fix_ratio", "grid_search"]#, "iter_search"]
                res = benchmarking_ss(100, nothing, method = method, one_se_rule = one_se_rule, 
                                        one_se_rule_pre = one_se_rule,
                                        nfold = 5, nfold_pre = 5,
                                        dataseed = 1,
                                        snr = 0.1,
                                        nλ = 2, # only for grid_search or iter_search
                                        nk = 50 # only for fix_ratio
                                    )
                if one_se_rule
                    @test res[2] > 0
                else
                    @test res[2] < res[4]
                end
            end
        end
    end
end

end