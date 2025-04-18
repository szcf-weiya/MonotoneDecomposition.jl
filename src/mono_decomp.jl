using StatsBase
import StatsBase.predict
using Statistics
using LinearAlgebra
using RCall
using Plots
import Plots.plot
using JuMP
using ECOS
using LaTeXTables
using Serialization
using ProgressMeter
using Gurobi
using GoldfarbIdnaniSolver
using MonotoneSplines
import MonotoneSplines.build_model
import MonotoneSplines.pick_knots
# using Distributed
# using SharedArrays

OPTIMIZER = ECOS.Optimizer
OK_STATUS = [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

# https://www.gurobi.com/documentation/current/refman/threads.html#parameter:Threads
function gurobi(nthread::Int = 0)
    GRB_ENV = Gurobi.Env()
    # disable print of `set parameter ..`
    # find in https://github.com/jump-dev/Gurobi.jl/blob/master/src/gen91/libgrb_api.jl
    # and inspired by https://support.gurobi.com/hc/en-us/community/posts/4412624836753-Do-not-print-Set-parameter-Username-in-console
    GRBsetintparam(GRB_ENV, GRB_INT_PAR_OUTPUTFLAG, 0)
    GRBsetintparam(GRB_ENV, GRB_INT_PAR_THREADS, nthread)
    # GRBsetdblparam(GRB_ENV, GRB_DBL_PAR_OPTIMALITYTOL, 1e-9)
    # GRBsetdblparam(GRB_ENV, GRB_DBL_PAR_FEASIBILITYTOL, 1e-9)
    # GRBsetdblparam(GRB_ENV, GRB_DBL_PAR_BARCONVTOL, 1e-9)
    # GRBsetdblparam(GRB_ENV, GRB_DBL_PAR_BARQCPCONVTOL, 1e-9)
    global OPTIMIZER = () -> Gurobi.Optimizer(GRB_ENV)
end

function ecos()
    global OPTIMIZER = ECOS.Optimizer
end

abstract type WorkSpace end

# difference compared to MonoDecomp
# workspace can be shared with different tuning parameters
# MonoDecomp records the specific results with particular parameters
mutable struct WorkSpaceSS <: WorkSpace
    evaluated::Bool
    J::Int
    B::Matrix{Float64}
    L::Matrix{Float64}
    H::Matrix{Int}
    mx::Real
    rx::Real
    idx::Vector{Int}
    idx0::Vector{Int}
    Bend::Matrix{Float64}
    Bendd::Matrix{Float64}
    knots::Vector{Float64}
    # Bnew::AbstractMatrix{Float64} # not for cross-validation
    WorkSpaceSS() = new()
    # WorkSpaceSS(e, j, b, l, h, m, r, i, i0, be, bd, k) = new(e, j, b, l, h, m, r, i, i0, be, bd, k)
end

mutable struct WorkSpaceCS <: WorkSpace
    evaluated::Bool
    J::Int
    B::Matrix{Float64}
    rB::RObject
    H::Matrix{Int}
    W::Matrix{Float64} # [B'; -B'] * [B -B]
    V::Matrix{Float64} # [B'; B'] * [B B]
    WorkSpaceCS() = new()
end

struct MonoDecomp{T <: AbstractFloat}
    γup::Vector{T}
    γdown::Vector{T}
    γhat::Vector{T}
    yhat::Vector{T}
    λ::T
    μ::T
    workspace::WorkSpace
end

"""
    mono_decomp(y::AbstractVector)

Perform monotone decomposition on vector `y`, and return `yup`, `ydown`.
"""
function mono_decomp(y::AbstractVector{T}; tol = 1e-10) where T <: AbstractFloat
    n = length(y)
    yup = zeros(n)
    ydown = zeros(n)
    for i = 1:n
        if i == 1
            yup[i] = y[i]
        else
            if y[i] > y[i-1]
                yup[i] = y[i] - y[i-1] + yup[i-1]
                ydown[i] = ydown[i-1]
            else
                yup[i] = yup[i-1]
                ydown[i] = y[i] - y[i-1] + ydown[i-1]
            end
        end
    end
    μup = mean(yup)
    μdown = mean(ydown)
    μ = (μup + μdown) / 2
    yup = yup .- μup .+ μ
    ydown = ydown .- μdown .+ μ
    err = norm(y - yup - ydown)
    if err > tol
        @warn "err $err larger than tol $tol"
    end
    return yup, ydown
end

"""
    recover(Σ)

Recover matrix from the vector-stored `Σ`.
"""
function recover(Σ::AbstractVector{T}) where T <: AbstractFloat
    n = length(Σ)
    p = round(Int, n/4)
    res = zeros(p, p)
    for i = 1:p
        res[i, i] = Σ[i]
    end
    for i = 1:p-1
        res[i, i+1] = Σ[p + i]
        res[i+1, i] = Σ[p + i]
    end
    for i = 1:p-2
        res[i, i+2] = Σ[2p + i]
        res[i+2, i] = Σ[2p + i]
    end
    for i = 1:p-3
        res[i, i+3] = Σ[3p + i]
        res[i+3, i] = Σ[3p + i]
    end
    return res
end

"""
    smooth_spline(x::AbstractVector, y::AbstractVector, xnew::AbstractVector)

Perform smoothing spline on `(x, y)`, and make predictions on `xnew`.

Returns: `yhat`, `ynewhat`,....
"""
function smooth_spline(x::AbstractVector{T}, y::AbstractVector{T}, xnew::AbstractVector{T}; keep_stuff = false, design_matrix = false, LOOCV = false) where T <: AbstractFloat
    spl = R"smooth.spline($x, $y, keep.stuff = $keep_stuff, cv = $LOOCV)"
    Σ = nothing
    if keep_stuff
        Σ = recover(rcopy(R"$spl$auxM$Sigma"))
    end
    B = nothing
    if design_matrix
        knots = rcopy(R"$spl$fit$knot")[4:end-3]
        bbasis = R"fda::create.bspline.basis(breaks = $knots, norder = 4)"
        B = rcopy(R"predict($bbasis, $knots)")
    end
    λ = rcopy(R"$spl$lambda")
    coef = rcopy(R"$spl$fit$coef")
    return rcopy(R"predict($spl, $x)$y"), rcopy(R"predict($spl, $xnew)$y"), Σ, λ, spl, B, coef
end

function cv_smooth_spline(x::AbstractVector{T}, y::AbstractVector{T}, λs::AbstractVector{T}; nfold = 5, seed = rand(UInt64), one_se_rule = false, figname  = nothing) where T <: Real
    n = length(x)
    folds = div_into_folds(n, K = nfold, seed = seed)
    nλ = length(λs)
    errs = zeros(nfold, nλ)
    for k in 1:nfold
        test_idx = folds[k]
        train_idx = setdiff(1:n, test_idx)
        for (i, lam) in enumerate(λs)
            spl = R"smooth.spline($(x[train_idx]), $(y[train_idx]), lambda = $lam)"
            yhat = rcopy(R"predict($spl, $(x[test_idx]))$y")
            # for LOOCV, yhat is a scalar, while y[test_idx] is a vector, use `.-` can help
            errs[k, i] = norm(yhat .- y[test_idx])^2 / length(test_idx) 
        end
    end
    μerr = dropdims(mean(errs, dims = 1), dims = 1)
    σerr = dropdims(std(errs, dims = 1), dims = 1) / sqrt(nfold)
    if one_se_rule
        ind = cv_one_se_rule(μerr, σerr, small_is_simple = false)
    else
        ind = argmin(μerr)
    end
    if ind == nλ
        @warn "the optimal is on the right boundary of λs"
    elseif ind == 1
        @warn "the optimal is on the left boundary of λs"
    end
    optlam = λs[ind]
    if !isnothing(figname)
        # save to sil
        serialize(figname[1:end-4] * ".sil", [μerr, σerr, λs, nfold, ind])
        savefig(cvplot(μerr, σerr, λs, nfold = nfold, ind0 = ind, lbl = "\\lambda"), figname)
    end
    spl = R"smooth.spline($x, $y, lambda = $optlam, keep.stuff = TRUE)"
    Σ = recover(rcopy(R"$spl$auxM$Sigma"))
    knots = rcopy(R"$spl$fit$knot")[4:end-3]
    bbasis = R"fda::create.bspline.basis(breaks = $knots, norder = 4)"
    B = rcopy(R"predict($bbasis, $knots)")
    λ = rcopy(R"$spl$lambda")
    coef = rcopy(R"$spl$fit$coef")
    return rcopy(R"predict($spl, $x)$y"), Σ, λ, spl, B, coef, μerr
end

"""
    mono_decomp_cs(x::AbstractVector, y::AbstractVector)

Monotone Decomposition with Cubic B-splines by solving an optimization problem.
"""
function mono_decomp_cs(x::AbstractVector{T}, y::AbstractVector{T}; 
                    s::T = 1.0, s_is_μ = true,
                    J = 4,
                    workspace = nothing,
                    use_GI = false,
                    )::MonoDecomp{T} where T <: AbstractFloat
    if isnothing(workspace) || !workspace.evaluated
        workspace = WorkSpaceCS()
        workspace.J = J
        build_model!(workspace, x, use_GI = use_GI)
    end
    if use_GI
        _γhat = try
            _optim(y, workspace, s)
        catch e
            @warn "GI solver failed due to $e; use default solver"
            _optim(y, workspace, [s])[:]
        end
    else
        if s_is_μ
            _γhat = _optim(y, workspace, [s])[:]
        else
            _γhat = _optim(y, workspace.J, workspace.B, s, workspace.H)
        end    
    end
    γup = _γhat[1:J]
    γdown = _γhat[J+1:2J]
    γhat = γup .+ γdown
    yhat = workspace.B * γhat
    return MonoDecomp(γup, γdown, γhat, yhat, 0.0, s, workspace)
end

MDCS = mono_decomp_cs

function cubic_spline(x::AbstractVector{T}, y::AbstractVector{T}, xnew::AbstractVector{T}; J = 4) where T <: AbstractFloat
    spl = MonotoneSplines.bspline(x, y, J)
    return MonotoneSplines.predict(spl, x), MonotoneSplines.predict(spl, xnew), spl
end

mono_decomp_cs(J) = (x, y; kw...) -> mono_decomp_cs(x, y; J = J, kw...)
cubic_spline(J) = (x, y, xnew) -> cubic_spline(x, y, xnew; J = J)

"""
    cv_cubic_spline(x::AbstractVector, y::AbstractVector, xnew::AbstractVector)

B-spline fitting with `nfold` CV-tuned `J` from `Js`.
"""
function cv_cubic_spline(x::AbstractVector{T}, y::AbstractVector{T}, xnew::AbstractVector{T}; 
                    nfold = 10, Js = 4:50, one_se_rule = false, figname = nothing) where T <: AbstractFloat
    n = length(x)
    folds = div_into_folds(n, K = nfold, seed = rand(UInt64))
    errs = zeros(nfold, length(Js))
    for k = 1:nfold
        test_idx = folds[k]
        train_idx = setdiff(1:n, test_idx)
        for (j, J) in enumerate(Js)
            _, yhat = cubic_spline(J)(x[train_idx], y[train_idx], x[test_idx])
            errs[k, j] = norm(yhat .- y[test_idx])^2 / length(test_idx)
        end
    end
    μerr = dropdims(mean(errs, dims = 1), dims = 1)
    σerr = dropdims(std(errs, dims = 1), dims = 1) / sqrt(nfold)
    if one_se_rule
        ind = cv_one_se_rule(μerr, σerr, small_is_simple = true)
    else
        ind = argmin(μerr)
    end
    Jopt = Js[ind]
    @debug "one_se_rule = $one_se_rule: Jopt = $Jopt"
    if !isnothing(figname)
        silname = figname[1:end-4] * "_J.sil"
        serialize(silname, [μerr, σerr, Js, nfold, ind])
        savefig(cvplot(μerr, σerr, 1.0 .* Js, nfold = nfold, ind0 = ind, lbl = "J"), figname)
    end
    yhat, yhatnew = cubic_spline(Jopt)(x, y, xnew)
    return Jopt, yhat, yhatnew
end

"""
    cv_mono_decomp_cs(x::AbstractVector, y::AbstractVector, xnew::AbstractVector; )
    cv_mono_decomp_cs(x::AbstractVector, y::AbstractVector; fixJ = true)

Cross-validation for Monotone Decomposition with Cubic B-splines. Parameters `J` and `s` (`μ` if `s_is_μ`) are tuned by cross-validation.

- if `fixJ == true`, then `J` is CV-tuned by the corresponding cubic B-spline fitting
- if `fixJ == false`, then both `J` and `s` would be tuned by cross-validation.

# Arguments

- `figname`: if not `nothing`, then the CV erro figure will be saved to the given name (can include the path)
"""
function cv_mono_decomp_cs(x::AbstractVector{T}, y::AbstractVector{T}, xnew::AbstractVector{T}; 
                                nfold = 10, Js = 4:50, 
                                ss = 10.0 .^ (-6:0.1:-1), 
                                s_is_μ = true, figname = nothing, 
                                one_se_rule = false, use_GI = false) where T <: AbstractFloat
    n = length(x)
    folds = div_into_folds(n, K = nfold, seed = rand(UInt64))
    errs = zeros(nfold, length(Js), length(ss))
    n2 = length(folds[1])
    n1 = n - n2
    # @assert n % nfold == 0
    for k = 1:nfold
        test_idx = folds[k]
        train_idx = setdiff(1:n, test_idx)
        for (j, J) in enumerate(Js)
            workspace = WorkSpaceCS()
            workspace.J = J
            build_model!(workspace, x[train_idx], use_GI = use_GI)
            # workspace = nothing
            @assert s_is_μ
            if use_GI
                γhats = hcat([_optim(y[train_idx], workspace, μ) for μ in ss]...)
            else
                γhats = _optim(y[train_idx], workspace, ss)
            end
            ynews = predict(workspace, x[test_idx], γhats[1:J, :] .+ γhats[J+1:2J, :])
            for (i, s) in enumerate(ss)
                errs[k, j, i] = norm(ynews[:, i] .- y[test_idx])^2 / length(test_idx)
            end
        end
    end
    μerr = dropdims(mean(errs, dims = 1), dims = 1)
    σerr = dropdims(std(errs, dims = 1), dims=1) / sqrt(nfold)
    if one_se_rule
        ind = cv_one_se_rule(μerr, σerr, small_is_simple = [true, false])
    else
        ind = argmin(μerr)
    end
    Jopt = Js[ind[1]]
    sopt = ss[ind[2]]
    @debug "Jopt = $Jopt, μopt = $sopt"
    if !isnothing(figname)
        silname = figname[1:end-4] * "_Jmu.sil"
        serialize(silname, [μerr, σerr, Jopt, sopt, Js, ss, nfold])
        savefig(cvplot(μerr, σerr, 1.0 .* Js, ss, lbl = ["J", "μ"], nfold = nfold, ind0 = ind), figname)
    end
    D = mono_decomp_cs(Jopt)(x, y; s = sopt, use_GI = use_GI)
    return D, ss[argmin(μerr)[2]], μerr, σerr  # use for recording without 1se
end

"""
    cv_mono_decomp_cs(x::AbstractVector, y::AbstractVector)

Cross-validation for monotone decomposition with cubic B-splines when the fixed `J` is CV-tuned by the corresponding cubic B-spline fitting method.
"""
function cv_mono_decomp_cs(x::AbstractVector{T}, y::AbstractVector{T}; ss = 10.0 .^ (-6:0.1:-1), 
                                                            figname = nothing, 
                                                            nfold = 10, 
                                                            nfold_pre = 10,
                                                            fixJ = true,
                                                            x0 = x,
                                                            Js = 4:50,
                                                            one_se_rule = false,
                                                            use_GI = false, # currently only for single μ
                                                            one_se_rule_pre = false)::Tuple{MonoDecomp{T}, T, Array{T}, Array{T}} where T <: AbstractFloat
    if length(Js) == 1
        fixJ = true # only one J can be selected
    end
    if fixJ
        if length(Js) == 1
            J = Js[1]
            if length(ss) == 1
                D = mono_decomp_cs(J)(x, y; s = ss[1], use_GI = use_GI)
                return D, ss[1], [0.0], [0.0]
            end    
        else
            J, _ = cv_cubic_spline(x, y, x0, one_se_rule = one_se_rule_pre, nfold = nfold_pre, Js = Js,
                                                figname = isnothing(figname) ? figname : figname[1:end-4] * "_bspl.png")
        end
        return cv_mono_decomp_cs(x, y, x0, Js = J:J, ss = ss, figname = figname, nfold = nfold, one_se_rule = one_se_rule, use_GI = use_GI)
    else
        return cv_mono_decomp_cs(x, y, x0, ss = ss, Js = Js, figname = figname, nfold = nfold, one_se_rule = one_se_rule, use_GI = use_GI)
    end
end

"""
    cv_mono_decomp_ss(x::AbstractVector, y::AbstractVector)

Cross Validation for Monotone Decomposition with Smoothing Splines. With `λ` tuned by smoothing spline, and then perform golden search for `μ`.

# Returns

- `D`: a `MonoDecomp` object.
- `workspace`: workspace contained some intermediate results
- `μmin`: the parameter `μ` that achieve the smallest CV error
- `μs`: the investigated parameter `μ`

# Example

```julia
x, y, x0, y0 = gen_data(100, 0.001, "SE_0.1")
res, workspace = cv_mono_decomp_ss(x, y, one_se_rule = true, figname = "/tmp/p.png", tol=1e-3)
yup = workspace.B * res.γup
ydown = workspace.B * res.γdown
scatter(x, y)
scatter!(x, yup)
scatter!(x, ydown)
```
"""
function cv_mono_decomp_ss(x::AbstractVector{T}, y::AbstractVector{T}; figname = nothing, 
                                                            verbose = true,
                                                            nfold = 10, nfold_pre = 10,
                                                            tol = 1e-5,
                                                            tol_boundary = 1e-1,
                                                            x0 = x, # test point of x
                                                            method = "single_lambda", # fix_ratio, iter_search, grid_search
                                                            kmin = 0.05, kmax = 1 - 1e-8, nk = 100, #used in fix_ratio
                                                            multi_fix_ratio = false,
                                                            nλ = 10, rλ = 0.5, # used in grid_search
                                                            μrange = nothing,
                                                            ngrid_μ = 20, # used in double grid
                                                            rλs = nothing,
                                                            rel_tol = 1e-1, maxiter = 10, # iter_search
                                                            k_magnitude = 1,
                                                            seed = rand(UInt64),
                                                            prop_nknots = 1.0,
                                                            minmaxμ = 1e-1,
                                                            include_boundary = false,
                                                            same_J_after_CV = false,
                                                            λs_in_ss = 10.0 .^ (-16:0.25:2),
                                                            search_μ_in_log_scale = false,
                                                            rerun_check = false,
                                                            one_se_rule = false, 
                                                            one_se_rule_pre = false, 
                                                            use_r_ss = false,
                                                            kw...) where T <: AbstractFloat
    if use_r_ss
        yhat, yhatnew, Ω, λ, spl, B, γss = smooth_spline(x, y, x0, design_matrix = true, keep_stuff = true, LOOCV = true)
    else
        yhat, Ω, λ, spl, B, γss = cv_smooth_spline(x, y, λs_in_ss, nfold = nfold_pre, seed = seed, 
                                                one_se_rule = one_se_rule_pre, 
                                                figname = isnothing(figname) ? figname : figname[1:end-4] * "_ss.png")
    end
    @debug "λopt in ss = $λ"
    yhatnew = rcopy(R"predict($spl, $x0)$y")
    γup, γdown = mono_decomp(γss)
    s0 = norm(B * (γup .- γdown))
    s_residual = norm(y - yhat)^2
    s_smoothness = λ * (γup .+ γdown)' * Ω * (γup .+ γdown)
    # s_discrepancy = [max(eps(), min(s_residual, s_smoothness)) / 10^k_magnitude, 
    #                             max(s_residual, s_smoothness) * 10^k_magnitude]
    small_tol = sqrt(eps())
    s_discrepancy = [small_tol, max(s_residual, s_smoothness, small_tol) * 10^k_magnitude]
    if isnothing(μrange)
        if s0 == 0
            μrange = [small_tol, small_tol * 10^k_magnitude]
        else
            μrange = s_discrepancy / s0^2
        end
    end
    verbose && @info "μrange for compatible terms: $μrange"
    verbose && @info "s0 = $s0, s_residual = $s_residual, s_smoothness = $s_smoothness, s_discrepancy = $s_discrepancy"
    if μrange[2] < minmaxμ # cbrt of eps()
        μrange[2] = minmaxμ
    end
    verbose && @info "μrange: $μrange"
    if method == "single_lambda"
        verbose && @info "Smoothing Splines with fixed λ"
        D, μs, errs, σerrs = cvfit_gss(x, y, μrange, λ, figname = figname, nfold = nfold, tol = tol, seed = seed, prop_nknots = prop_nknots, include_boundary = include_boundary)
    elseif method == "fix_ratio"
        verbose && @info "Smoothing Splines with fixratio strategy"
        #ks = range(kmin, kmax, length = nk)
        #ks = exp.(range(log.(kmin), log.(kmax), length = nk))
        ks = 1 ./ (exp.(range(log.(1 ./ kmax .- 1), log.(1 ./ kmin .- 1), length = nk)) .+ 1)
        D, μs, errs, σerrs = cvfit(x, y, ks, λ, figname = figname, nfold = nfold, seed = seed, prop_nknots = prop_nknots, multi_fix_ratio = multi_fix_ratio)
    elseif method == "iter_search" #88
        verbose && @info "Smoothing Splines with iter-search: λ -> μ -> λ -> ... -> μ"
        iter = 0
        λ0 = λ # make a backup
        λrange = [min(0.1*λ0, eps()), max(10λ0, 1e-6)]
        while true
            iter += 1
            ## tune mu given lambda
            @debug "tune mu given lambda = $λ"
            # D1, workspace1 = cvfit(x, y, μmax * r, [λ], nfold = nfold, figname = figname, nμ = nμ, ρ = ρ)
            ## if needed, perform one se rule on the last iteration given lambda
            D1, μs, errs, σerrs = cvfit_gss(x, y, μrange, λ, nfold = nfold, 
                                            figname = isnothing(figname) ? figname : figname[1:end-4] * "$iter-mu.png", 
                                            tol = tol, seed = seed, prop_nknots = prop_nknots, same_J_after_CV = same_J_after_CV)
            μrange = [min(μrange[1], 0.5 * D1.μ), max(μrange[2], 2 * D1.μ)]
            # since D is not defined in the first iteration, so use `if..else`, and hence cannot use `ifelse`
            if iter == 1
                err_μ = 1.0
            else
                err_μ = abs(D1.μ - D.μ) / D1.μ
            end
            ## re-tune lambda given mu
            @debug "tune lambda given mu = $(D1.μ)"
            # D, workspace = cvfit(x, y, D1.μ, λ, nfold = nfold, figname = figname, nλ = nλ, ρ = ρ)
            D, _ = cvfit_gss(x, y, λrange, D1.μ, nfold = nfold, 
                            figname = isnothing(figname) ? figname : figname[1:end-4] * "$iter-lam.png", 
                            λ_is_μ = true, tol = tol, seed = seed, prop_nknots = prop_nknots, same_J_after_CV = same_J_after_CV)
            λrange = [min(λrange[1], 0.5 * D.λ), max(λrange[2], 2 * D.λ)]
            err_λ = abs(D.λ - D1.λ) / D.λ
            λ = D.λ # for next iteration
            @debug "iter = $iter, err_μ = $err_μ, err_λ = $err_λ, err_fit = $(minimum(errs))"
            if (iter > maxiter) | (max(err_λ, err_μ) < rel_tol)
                break
            end
        end
    elseif method == "double_grid"
        if isnothing(rλs)
            λs = range(1-rλ, 1+rλ, length = nλ) .* λ
        else
            λs = rλs .* λ
        end
        @debug "grid search λ in $λs"
        μs = exp.(range(log(μrange[1]), log(μrange[2]), length = ngrid_μ))
        D, errs, σerrs = cvfit(x, y, μs, λs, nfold = nfold, figname = figname, seed = seed, prop_nknots = prop_nknots, include_boundary = include_boundary, same_J_after_CV = same_J_after_CV)
    else # grid_search
        if isnothing(rλs)
            λs = range(1-rλ, 1+rλ, length = nλ) .* λ
        else
            λs = rλs .* λ
        end
        verbose && @info "Smoothing Splines with grid-search λ ∈ $λs, μ ∈ $μrange"
        D, μs, errs, σerrs = cvfit_gss(x, y, μrange, λs, nfold = nfold, figname = figname, seed = seed, tol = tol, 
                                        tol_boundary = tol_boundary, prop_nknots = prop_nknots, include_boundary = include_boundary, same_J_after_CV = same_J_after_CV,
                                        log_scale = search_μ_in_log_scale, rerun_check = rerun_check)
    end
    μmin = D.μ
    λmin = D.λ
    if one_se_rule
        if isa(errs, Matrix)
            ind = cv_one_se_rule(errs, σerrs, small_is_simple = [false, false])
            μopt = μs[ind[1]]
            λopt = λs[ind[2]]
        else
            ind = cv_one_se_rule(errs, σerrs, small_is_simple = false)
            μopt = μs[ind]
            λopt = λmin
        end
        @debug "use 1se rule: before 1se: μ = $μmin; after 1se: μ = $μopt"
        @debug "use 1se rule: before 1se: λ = $λmin; after 1se: λ = $λopt"
        # workspace has been defined, so J would be inherited regardless of prop_nknots
        D = mono_decomp_ss(D.workspace, x, y, λopt, μopt) 
    end
    return D, μmin, μs, errs, σerrs, yhat, yhatnew, γss
end


"""
    _optim(y::AbstractVector, workspace::WorkSpaceCS, μs::AbstractVector)
    _optim(y::AbstractVector, J::Int, B::AbstractMatrix, H::AbstractMatrix{Int}, μs::AbstractVector)

Optimization for monotone decomposition with cubic B-splines.

    _optim(y::AbstractVector, J::Int, B::AbstractMatrix, H::AbstractMatrix{Int}, L::AbstractMatrix, λs::AbstractVector, μs::AbstractVector)

Optimization for monotone decomposition with smoothing splines.

    _optim!(y::AbstractVector, J::Int, B::AbstractMatrix, s::Union{Nothing, Real}, γhat::AbstractVector, H::AbstractMatrix{Int}; L, t, λ, μ)


"""
function _optim(y::AbstractVector{T}, workspace::WorkSpaceCS, μ::T) where T <: AbstractFloat
    D = workspace.V .+ μ * workspace.W
    q = workspace.B' * y
    sol, lagr, crval, iact, nact, iter = solveQP(D, vcat(q, q), -workspace.H'*1.0, zeros(2(workspace.J-1)))
    return sol
end

function _optim(y::AbstractVector{T}, workspace::WorkSpaceCS, μs::AbstractVector{T}) where T <: AbstractFloat
    @assert workspace.evaluated
    return _optim(y, workspace.J, workspace.B, workspace.H, μs)
end

# without penalty, common cubic spline
function _optim(y::AbstractVector{T}, J::Int, B::AbstractMatrix{T}, 
                H::AbstractMatrix{Int},
                μs::AbstractVector{T};
                del_constraint = true,
                ) where T <: AbstractFloat
    L = zeros(J, J) # dummy not used # a side product is that more memory has been occupied
    λs = zeros(length(μs))
    return _optim(y, J, B, H, L, λs, μs, del_constraint = del_constraint)
end

@inline function _optim(y::AbstractVector{T}, J::Int, B::AbstractMatrix{T}, 
                H::AbstractMatrix{Int}, L::AbstractMatrix{T}, 
                λs::AbstractVector{T}, μs::AbstractVector{T};
                del_constraint = true,
                strict = false,
                ) where T <: AbstractFloat
    model = Model(OPTIMIZER)
    set_silent(model)
    # TODO: test the speed and possibly switch to it if necessary, as in _optim!
    # set_optimizer_attribute(model, "BarHomogeneous", 1) # use only in demo since practically if a numerical error occurs, although it can be solved by another method, but it also implies that the current parameter is not the best
    # set_optimizer_attribute(model, "max_iter", 10000)
    @variable(model, γ[1:2J])
    @variable(model, z)
    @constraint(model, c1, H * γ .<= 0)
    m = length(λs)
    γhats = zeros(2J, m)
    if del_constraint
        λ = λs[1]
        μ = μs[1]
        # @variable(model, λ == sqrt(λs[1]), Param())
        # @variable(model, μ == sqrt(μs[1]), Param())
        @constraint(model, c3, [z; vcat(y - B * (γ[1:J] + γ[J+1:2J]), 
                                        sqrt(λ) * L' * (γ[1:J] + γ[J+1:2J]), 
                                        # λ * L' * (γ[1:J] + γ[J+1:2J]), 
                                        sqrt(μ) * B * (γ[1:J] - γ[J+1:2J]) 
                                        # μ * B * (γ[1:J] - γ[J+1:2J]) 
                                        ) ] in SecondOrderCone() )
        @objective(model, Min, z)
        JuMP.optimize!(model)
        status = termination_status(model)
        i = 1
        if status == MOI.NUMERICAL_ERROR
            @debug "$status when λ=$λ, μ=$μ: use constant half mean as its estimate"
            strict && error(status)
            γhats[:, i] .= mean(y) / 2
        else
            if !(status in OK_STATUS)
                @debug "$status when λ=$λ, μ=$μ: direct take the solution"
                strict && error(status)
            end
            try
                γhats[:, i] .= value.(γ) # ITERATION_LIMIT
            catch
                @debug "when λ=$λ, μ=$μ: no solution, $status"
                ts = strip(read(`date -Iseconds`, String))
                serialize("bug_$(ts).sil", [y, J, B, H, L, λs, μs, i, status])
                γhats[:, i] .= mean(y) / 2
            end
        end
        for i in 2:m
            λ = λs[i]
            μ = μs[i]
            # set_value(λ, sqrt(λs[i]))
            # set_value(μ, sqrt(μs[i]))
            delete(model, c3)
            unregister(model, :c3)
            @constraint(model, c3, [z; vcat(y - B * (γ[1:J] + γ[J+1:2J]), 
                                        sqrt(λ) * L' * (γ[1:J] + γ[J+1:2J]), 
                                        # λ * L' * (γ[1:J] + γ[J+1:2J]), 
                                        sqrt(μ) * B * (γ[1:J] - γ[J+1:2J]) 
                                        # μ * B * (γ[1:J] - γ[J+1:2J]) 
                                        ) ] in SecondOrderCone() )
            JuMP.optimize!(model)
            status = termination_status(model)
            # obj_val = objective_value(model)
            # println("i = $i, obj_val = $obj_val")
            if status == MOI.NUMERICAL_ERROR
                @debug "$status"
                strict && error(status)
                γhats[:, i] .= mean(y) / 2
            else
                if !(status in OK_STATUS)
                    @debug "$status"
                    strict && error(status)
                end
                try
                    γhats[:, i] .= value.(γ)
                catch e
                    @warn e 
                    ts = strip(read(`date -Iseconds`, String))
                    serialize("bug_$(ts).sil", [y, J, B, H, L, λs, μs, i, status])    
                    γhats[:, i] .= mean(y) / 2
                end
            end
        end
    else
        @variable(model, λ)
        @variable(model, μ)
        @constraint(model, c3, [z; vcat(y - B * (γ[1:J] + γ[J+1:2J]), 
                                # sqrt(λ) * L' * (γ[1:J] + γ[J+1:2J]), 
                                λ * L' * (γ[1:J] + γ[J+1:2J]), 
                                # sqrt(μ) * B * (γ[1:J] - γ[J+1:2J]) 
                                μ * B * (γ[1:J] - γ[J+1:2J]) 
                                ) ] in SecondOrderCone() )
        for i = 1:m
            JuMP.fix(λ, sqrt(λs[i]))
            JuMP.fix(μ, sqrt(μs[i]))
            JuMP.optimize!(model)
            status = termination_status(model)
            obj_val = objective_value(model)
            # println("i = $i, obj_val = $obj_val")
            if status == MOI.NUMERICAL_ERROR
                @debug "$status"
                strict && error(status)
                γhats[:, i] .= mean(y) / 2
            else
                if !(status in OK_STATUS)
                    @debug "$status"
                    strict && error(status)
                end
                γhats[:, i] .= value.(γ)
            end
        end
    end
    # model = nothing
    # GC.gc(); GC.gc();
    # GC.gc(); GC.gc();
    #finalize(ecos)
    return γhats
end

function _optim!(y::AbstractVector{T}, J::Int, B::AbstractMatrix{T}, s::Union{Nothing, Real}, γhat::AbstractVector{T}, H::AbstractMatrix{Int}; L = nothing, t = nothing, λ = nothing, μ = nothing, BarHomogeneous = false, strict = false) where T <: AbstractFloat
    # construct constraint ‖ B(γ^u - γ^d) ‖ < s
    # model = Model(GLPK.Optimizer) # cannot do second order, throw error
    # `MOI.VectorAffineFunction{Float64}`-in-`MOI.SecondOrderCone` constraints are not supported and cannot be bridged into supported constrained variables and constraints. See details below:
    model = Model(OPTIMIZER)
    isgurobi = typeof(model.moi_backend.optimizer).parameters[1] == Gurobi.Optimizer
    if BarHomogeneous && isgurobi
        set_optimizer_attribute(model, "BarHomogeneous", 1) # only for gurobi
    end
    set_silent(model)
    @variable(model, γ[1:2J])
    @variable(model, z)
    @constraint(model, c1, H * γ .<= 0)
    if isnothing(μ)
        # @constraint(model, c2, norm(B * (γ[1:J] - γ[J+1:2J]) ) <= s)
        @constraint(model, c2, [s; B * (γ[1:J] - γ[J+1:2J])] in SecondOrderCone()  )
        # @objective(model, Min, norm(y - B * (γ[1:J] + γ[J+1:2J]) ))
        # cannot directly use norm
        if isnothing(λ)
            @constraint(model, c3, [z; y - B * (γ[1:J] + γ[J+1:2J])] in SecondOrderCone() )
            if !isnothing(t)
                @constraint(model, c4, [t; L' * (γ[1:J] + γ[J+1:2J])] in SecondOrderCone() )
            end
        else
            @constraint(model, c3, [z; vcat(y - B * (γ[1:J] + γ[J+1:2J]), sqrt(λ) * L' * (γ[1:J] + γ[J+1:2J])) ] in SecondOrderCone() )
        end
    else
        if !isnothing(λ)
            @constraint(model, c3, [z; vcat(y - B * (γ[1:J] + γ[J+1:2J]), sqrt(λ) * L' * (γ[1:J] + γ[J+1:2J]), sqrt(abs(μ)) * B * (γ[1:J] - γ[J+1:2J]) ) ] in SecondOrderCone() )
        else
            @constraint(model, c3, [z; vcat(y - B * (γ[1:J] + γ[J+1:2J]), sqrt(μ) * B * (γ[1:J] - γ[J+1:2J]) ) ] in SecondOrderCone() )
        end
    end
    @objective(model, Min, z)
    JuMP.optimize!(model)
    # print(model)
    # solution_summary(model, verbose = true)
    status = termination_status(model)
    # debug
    # if status != MOI.OPTIMAL
    #     serialize("../res/res_monodecomp/debug/$status.sil", [y, J, B, s, L, t, λ])
    #     error(status)
    # end
    # end debug
    if status == MOI.NUMERICAL_ERROR
        strict && error(status)
        if isgurobi
            if !BarHomogeneous # not yet try BarHomogeneous
                @debug "try BarHomogeneous to rerun NUMERICAL_ERROR"
                _optim!(y, J, B, s, γhat, H, L = L, t = t, λ = λ, μ = μ, BarHomogeneous = true)
            else
                @debug "$status after trying BarHomogeneous algorithm"
                γhat .= mean(y) / 2
            end
        else
            @debug status
            γhat .= mean(y) / 2
        end
    else
        if !(status in OK_STATUS)
            @debug "$status"
            strict && error(status)
        end
        try
            γhat .= value.(γ)
        catch e
            @warn e 
            γhat .= mean(y) / 2
        end
    end
end

# share by cross validation
function _optim(y::AbstractVector{T}, J::Int, B::AbstractMatrix{T}, s::Union{Nothing, Real}, H::AbstractMatrix{Int}; kw...) where T <: AbstractFloat
    γhat = zeros(2J)
    _optim!(y, J, B, s, γhat, H; kw...)
    return γhat
end

function build_model!(workspace::WorkSpaceSS, x::AbstractVector{T}; ε = eps()^(1/3), prop_nknots = 1.0) where T <: AbstractFloat
    if !isdefined(workspace, 3) # the first two will be automatically initialized
        knots, mx, rx, idx, idx0 = pick_knots(x, all_knots = false, prop_nknots = prop_nknots)
        bbasis = R"fda::create.bspline.basis(breaks = $knots, norder = 4)"
        Ω = rcopy(R"fda::eval.penalty($bbasis, 2)")
        Ω += ε * 1.0I
        xbar = (x .- mx) ./ rx
        workspace.B = rcopy(R"fda::eval.basis($xbar, $bbasis)")
        workspace.Bend = rcopy(R"fda::eval.basis(c(0, 1), $bbasis)")
        workspace.Bendd = rcopy(R"fda::eval.basis(c(0, 1), $bbasis, Lfdobj=1)")
        J = length(knots) + 2
        workspace.H = construct_H(J)
        # calculate L
        try
            workspace.L = Matrix(cholesky(Symmetric(Ω)).L)
        catch e
            @debug "$e: standard cholesky failed, use pivoted cholesky"
            ## perform pivoted Cholesky
            chol = cholesky(Symmetric(Ω), Val(true), check = false)
            workspace.L = Matrix(chol.L[invperm(chol.p), 1:chol.rank])
        end
        # cannot update workspace in the argument
        workspace.evaluated = true
        workspace.J = J
        workspace.mx = mx
        workspace.rx = rx
        workspace.idx = idx
        workspace.idx0 = idx0
        workspace.knots = knots
    end    
end

"""
    mono_decomp_ss(workspace::WorkSpaceSS, x::AbstractVector{T}, y::AbstractVector{T}, λ::AbstractFloat, μ::AbstractFloat)

Monotone decomposition with smoothing splines.
"""
function mono_decomp_ss(workspace::WorkSpaceSS, x::AbstractVector{T}, y::AbstractVector{T}, λ::AbstractFloat, s::AbstractFloat; s_is_μ = true, prop_nknots = 1.0, strict = false)::MonoDecomp{T} where T <: AbstractFloat
    build_model!(workspace, x, prop_nknots = prop_nknots)
    if s_is_μ
        _γhat = _optim(y, workspace.J, workspace.B, nothing, workspace.H, L = workspace.L, t = nothing, λ = λ, μ = s, strict = strict)
    else
        _γhat = _optim(y, workspace.J, workspace.B, s, workspace.H, L = workspace.L, t = nothing, λ = λ, strict = strict)
    end
    # calculate properties of monotone decomposition
    J = workspace.J
    γup = _γhat[1:J]
    γdown = _γhat[J+1:2J]
    γhat = γup .+ γdown
    yhat = workspace.B * γhat
    return MonoDecomp(γup, γdown, γhat, yhat, λ, s, workspace)
end
function mono_decomp_ss(x::AbstractVector{T}, y::AbstractVector{T}; λ = 1.0, s = 1.0, s_is_μ = true) where T <: AbstractFloat
    workspace = WorkSpaceSS()
    return mono_decomp_ss(workspace, x, y, λ, s, s_is_μ = s_is_μ)
end

MDSS = mono_decomp_ss

# assume xnew lies in the middle of xrange
"""
    predict(D::MonoDecomp, xnew)

Predict at `xnew` given decomposition `D`.
"""
function predict(D::MonoDecomp, xnew::AbstractVector)
    yup, ydown = predict(D.workspace, xnew, D.γup, D.γdown)
    return yup + ydown
end

"""
    predict(W::WorkSpaceSS, xnew::AbstractVector, γup::AbstractVector, γdown::AbstractVector)
    predict(W::WorkSpaceCS, xnew::AbstractVector, γup::AbstractVector, γdown::AbstractVector)

Predict `yup` and `ydown` at `xnew` given workspace `W` and decomposition coefficients `γup` and `γdown`.
"""
function predict(W::WorkSpaceSS, xnew::AbstractVector, γup::AbstractVector, γdown::AbstractVector)
    xm = (xnew .- W.mx) ./ W.rx
    # evaluate on the whole dataset, so all should be in the middle
    @assert all(0 .<= xm .<= 1)
    Bnew = rcopy(R"splines::bs($xm, intercept = TRUE, knots=$(W.knots[2:end-1]), Boundary.knots = c(0, 1))")
    return Bnew * γup, Bnew * γdown
end

function predict(W::WorkSpaceCS, xnew::AbstractVector, γup::AbstractVector, γdown::AbstractVector)
    Bnew = rcopy(R"suppressWarnings(predict($(W.rB), $xnew))")
    return Bnew * γup, Bnew * γdown
end

"""
    predict(W::WorkSpaceSS, xnew::AbstractVector, γhat::AbstractVecOrMat)
    predict(W::WorkSpaceCS, xnew::AbstractVector, γhat::AbstractVecOrMat)

Make multiple predictions at `xnew` for each column of `γhat`.
"""
function predict(W::WorkSpaceSS, xnew::AbstractVector, γhat::AbstractVecOrMat)
    xnewbar = (xnew .- W.mx) ./ W.rx
    ind_right = xnewbar .> 1
    ind_left = xnewbar .< 0
    ind_middle = 0 .<= xnewbar .<= 1
    xm = xnewbar[ind_middle]
    n = length(xnew)
    if sum(ind_middle) == 0
        Bnew = zeros(0, W.J)
    else
        # xm cannot be empty
        Bnew = rcopy(R"splines::bs($xm, intercept = TRUE, knots=$(W.knots[2:end-1]), Boundary.knots = c(0, 1))")
    end
    if isa(γhat, AbstractVector)
        yhat = zeros(n)
        yhat[ind_middle] = Bnew * γhat
    else
        m = size(γhat, 2)
        yhat = zeros(n, m)
        yhat[ind_middle, :] = Bnew * γhat
    end
    boundaries = W.Bend * γhat
    slopes = W.Bendd * γhat

    # extrapolation
    # derivative at endpoints
    if isa(γhat, AbstractVector)
        yhat[ind_left] = boundaries[1] .+ slopes[1] * (xnewbar[ind_left] .- 0)
        yhat[ind_right] = boundaries[2] .+ slopes[2] * (xnewbar[ind_right] .- 1)
    else
        for j = 1:m
            yhat[ind_left, j] = boundaries[1, j] .+ slopes[1, j] .* (xnewbar[ind_left] .- 0)
            yhat[ind_right, j] = boundaries[2, j] .+ slopes[2, j] .* (xnewbar[ind_right] .- 1)
        end
    end
    return yhat
end

# for cubic splines, the warning evaluated at outside point can be safely ignored since
# the extrapolation of the function is the linear combinations of the extrapolation of each basis function
# NB: it is different from smoothing splines, whose extrapolation is linear. (#85)
function predict(W::WorkSpaceCS, xnew::AbstractVector, γhat::AbstractVecOrMat)
    Bnew = rcopy(R"suppressWarnings(predict($(W.rB), $xnew))")
    return Bnew * γhat
end

"""
    plot(obs, truth, D::MonoDecomp, other)

Plot the noised observations, the true curve, and the fitting from monotone decomposition `D` and `other` fitting technique.

- `obs`: usually be `[x, y]`
- `truth`: usually be `[x0, y0]`
- `D`: a `MonoDecomp` object
- `other`: the fitted curve `[x0, other]` by other method, where `x0` is omitted.

A typical usage can be `plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")`
"""
function plot(obs::AbstractVector{T}, truth::AbstractVector{T}, D::MonoDecomp, other::AbstractVector; 
                        prefix_title = "", postfix_title = "", competitor = "ss", digits = 3,
                        title = "",
                        kw...) where T <: AbstractVector
    x, y = obs
    x0, y0 = truth
    @assert x[1] <= x[2] <= x[3] # assume it has been sorted
    @assert x0[1] <= x0[2] <= x0[3]
    yup, ydown = predict(D.workspace, x0, D.γup, D.γdown)
    e1 = round(norm(yup + ydown - y0)^2 / length(y0), digits = digits)
    e2 = round(norm(other - y0)^2 / length(y0), digits = digits)
    adjust_hat_in_latex = ""
    if Plots.backend() == Plots.PGFPlotsXBackend()
        adjust_hat_in_latex = raw"\mkern-5mu"
    end
    lbl_other = latexstring("\$\\hat{$adjust_hat_in_latex f}_{\\mathrm{$competitor}} ($e2)\$") # adjust the hat location by \mkern
    if isnothing(title)
        title = ""
    else
        title = ""
        if isa(D.workspace, WorkSpaceSS)
            title *= "λ = $(round(D.λ, sigdigits=3))"
        else
            title *= "J = $(D.workspace.J)"
        end
        if !isnothing(D.μ)
            if length(title) > 0
                title *= ", "
            end
            title *= "μ = $(round(D.μ, sigdigits=3))"
        end
    end
    fig = scatter(x, y; title = prefix_title * title * postfix_title, label = "", ms = 2, xlab = L"x", ylab = L"y", labelfontsize = 14, legendfontsize = 14, tickfontsize = 12, titlefontsize = 16, kw...)
    lbls = [L"\hat{%$adjust_hat_in_latex f}_u", 
            L"\hat{%$adjust_hat_in_latex f}_d", 
            latexstring("\$\\hat{$adjust_hat_in_latex f}_u+\\hat{$adjust_hat_in_latex f}_d ($e1)\$"), L"f"]
    plot!(fig, x0, yup, label = "", ls = :dot)
    N = length(x0)
    plot!(fig, x0[1:20:N], yup[1:20:N], label = lbls[1], ls = :dot, markershape = :+, markersize = 2)
    plot!(fig, x0, ydown, label = lbls[2], ls = :dot)
    plot!(fig, x0, yup + ydown, 
            # label = L"\hat y_u + \hat y_d", 
            label = lbls[3],
            lw = 1.5)
    plot!(fig, x0, other, label = lbl_other, ls = :dash, lw = 1.5)
    plot!(fig, x0, y0, label = lbls[4], lw = 0.5)
    return fig
end

function construct_H(J::Int)
    A = zeros(Int, J-1, J)
    for i = 1:J-1
        A[i, i] = 1
        A[i, i+1] = -1
    end
    H = zeros(Int, 2J-2, 2J)
    H[1:J-1, 1:J] .= A[1:J-1, 1:J]
    H[J:2J-2, J+1:2J] .= -A[1:J-1, 1:J]
    return H
end

# no lambda, for cubic spline
"""
    build_model!(workspace::WorkSpaceCS, x::AbstractVector{T})

Calculate components that construct the optimization problem for Monotone Decomposition with Cubic splines.
"""
function build_model!(workspace::WorkSpaceCS, x::AbstractVector{T}; use_GI = false) where T <: AbstractFloat
    J = workspace.J
    B, rB = build_model(x, J)
    workspace.B = B
    workspace.rB = rB
    workspace.evaluated = true
    workspace.H = construct_H(J)
    if use_GI
        workspace.W = [B'; -B'] * [B -B]
        workspace.V = [B'; B'] * [B B]
    end
end


"""
    cvplot(sil::String)
    cvplot(μerr::AbstractVector, σerr::Union{Nothing, AbstractVector{T}}, paras::AbstractVector)
    cvplot(μerr::AbstractMatrix, σerr::AbstractMatrix, para1::AbstractVector, para2::AbstractVector)

Plot the cross-validation curves.
"""
function cvplot(sil::String, competitor = "bspl", title = "Leave-one-out CV Error ")
    if competitor == "bspl"
        μerr, σerr, μ, J, nfold = deserialize(sil)
        title = title * "(J = $J)"
    elseif competitor == "bspl2"
        μerr, σerr, Jopt, μopt, Js, ss, nfold = deserialize(sil)
        title = "$nfold-fold CV Error"
    else
        μerr, σerr, μ, λ, nfold = deserialize(sil)
        title = "$nfold-fold CV Error"        
        if length(λ) == 1
            λ1 = round.(λ[1], sigdigits=3)
            title = title * "(λ = $λ1)"
        end
    end
    if occursin("bspl", competitor)
        if size(μerr, 2) > 1
            cvplot(μerr, σerr, 1.0 * Js, ss, title = title, lbl = ["J", L"\log_{10} \mu"])
        else
            cvplot(μerr[:, 1], nothing, μ, title = title, lbl = L"\log_{10} \mu")
        end
    else
        # cvplot(μerr, σerr, μ, λ, nfold = nfold, lbl = [L"\log_{10}\mu", L"\lambda"], title = title)
        # -10 is tried for a better difference
        # offset = 10
        offset = 0
        cvplot(μerr[1:end-offset, :], σerr, μ[1:end-offset], λ, nfold = nfold, lbl = ["\\log_{10}\\mu", "\\lambda"], title = title)
    end
end

function add_log_str(lbl::Union{AbstractString, Nothing})
    if isa(lbl, String)
        if occursin("\\log_{10}", lbl)
            return lbl
        else
            return "\\log_{10}" * lbl
        end
    else
        return lbl
    end
end

function cvplot(μerr::AbstractVector{T}, σerr::Union{Nothing, AbstractVector{T}}, paras::AbstractVector{T}; title = nothing, lbl = nothing, nfold = 10, ind0 = nothing) where T <: AbstractFloat
    f = x -> x
    # determine the scale
    if length(paras) >= 3
        if paras[2] - paras[1] ≈ paras[3] - paras[2]
            f = x -> x
            if isnothing(lbl)
                lbl = ""
            end
        else
            f = x -> log10(x)
            lbl = add_log_str(lbl)
        end
    end
    ind = argmin(μerr)
    if isnothing(title)
        title = "$nfold-fold CV error"
    end
    if isnothing(σerr)
        p = plot(f.(paras), μerr, xlab = latexstring(lbl), title = title, legend = false, markershape = :x, ms = 2, labelfontsize = 14, legendfontsize = 14, tickfontsize = 12, titlefontsize = 16)
    else
        p = plot(f.(paras), μerr, yerrors = σerr, xlab = latexstring(lbl), title = title, legend = false, labelfontsize = 14, legendfontsize = 14, tickfontsize = 12, titlefontsize = 16)
    end
    annotate!(p, [(f.(paras)[ind], μerr[ind], ("o", :red))])
    if !isnothing(ind0)
        annotate!(p, [(f.(paras)[ind0], μerr[ind0], ("o", :red))])
    end
    return p
end

# assume matrix does not reduce to vector
function cvplot(μerr::AbstractMatrix{T}, σerr::Union{Nothing, AbstractMatrix{T}}, para1::AbstractVector{T}, para2::AbstractVector{T}; 
                    lbl = ["", ""], title = "", nfold = 10, ind0 = [nothing, nothing]) where T <: AbstractFloat
    n1 = length(para1)
    n2 = length(para2)
    if isnothing(σerr)
        if n1 == 1
            return cvplot(μerr[1, :], nothing, para2, nfold = nfold, ind0 = ind0[2], lbl = lbl[2], title = title)
        elseif n2 == 1
            return cvplot(μerr[:, 1], nothing, para1, nfold = nfold, ind0 = ind0[1], lbl = lbl[1], title = title)
        end
    else
        if n1 == 1
            return cvplot(μerr[1, :], σerr[1, :], para2, nfold = nfold, ind0 = ind0[2], lbl = lbl[2], title = title)
        elseif n2 == 1
            return cvplot(μerr[:, 1], σerr[:, 1], para1, nfold = nfold, ind0 = ind0[1], lbl = lbl[1], title = title)
        end    
    end
    n, m = size(μerr)
    @assert n == length(para1)
    @assert m == length(para2)
    if length(para1) < 3
        f = x -> x
    else
        if para1[2] - para1[1] ≈ para1[3] - para1[2]
            f = x -> x
        else
            f = x -> log10(x)
            lbl[1] = add_log_str(lbl[1])
        end    
    end
    if length(para2) < 3
        g = x -> x
    else
        if para2[2] - para2[1] ≈ para2[3] - para2[2]
            g = x -> x
        else
            g = x -> log10(x)
            lbl[2] = add_log_str(lbl[2]) 
        end    
    end
    ## assigning a latexstring to a[1] will convert the type from LaTeXString to String
    p = heatmap(g.(para2), f.(para1), μerr, xlab = latexstring(lbl[2]), ylab = latexstring(lbl[1]), title = title, labelfontsize = 14, legendfontsize = 14, tickfontsize = 12, titlefontsize = 16)
    ind = argmin(μerr)
    annotate!(p, [(g(para2[ind[2]]), f(para1[ind[1]]), ("o", :red))])
    if !isnothing(ind0[1])
        annotate!(p, [(g(para2[ind0[2]]), f(para1[ind0[1]]), ("o", :red))])
    end
    return p
end

"""
    cvfit(x::AbstractVector, y::AbstractVector, μmax::Real, rλ::Real, λstar::Real)
    cvfit(x::AbstractVector, y::AbstractVector, μstar::Real, λstar::Real)
    cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μmax::Real, λ::AbstractVector{T})
    cvfit(x::AbstractVector{T}, y::AbstractVector{T}, paras::AbstractMatrix{T})
    cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μ::AbstractVector{T}, λ::AbstractVector{T})
    cvfit(x::AbstractVector{T}, y::AbstractVector{T}, ks::AbstractVector{T}, λstar::Real)

Cross-validation for monotone decomposition.
"""
## TODO: remove? seems tricky, and several parameters should be determined
# automatically extend the ratio for searching λ
# function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μmax::Real, rλ::Real, λstar::Real; nfold = 10, nμ = 100, maxiter = 20, figname = "/tmp/cv_curve.png", ρ = 0.05, tol = eps()^(1/3), nλ = 10, seed = rand(UInt64)) where T <: AbstractFloat
#     λs = range(1-rλ, 1+rλ, length = nλ + 1) * λstar
#     D = cvfit(x, y, μmax, λs, nfold = nfold, nμ = nμ, maxiter = maxiter, figname = figname, ρ = ρ, seed = seed)
#     if abs(1 - D.λ / λstar) / (2rλ) < ρ
#         iter = 0
#         while true 
#             iter += 1
#             rλ = rλ + 0.05
#             nλ += 2
#             if rλ <= 0
#                 break
#             end
#             println("extend the search region of λ: searching [1-$rλ, 1+$rλ]")
#             λs = range(1-rλ, 1+rλ, length = nλ + 1) * λstar
#             D = cvfit(x, y, μmax, λs, nfold = nfold, nμ = nμ, maxiter = maxiter, figname = figname, ρ = ρ, seed = seed)
#             if (abs(1 - D.λ / λstar) / (2rλ) >= ρ) | (iter > maxiter)
#                 break
#             end
#         end
#     end
#     return D
# end

# tune lambda after fixing mu
function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μstar::Real, λstar::Real; nfold = 10, maxiter = 20, figname = "/tmp/cv_curve.png", ρ = 0.05, tol = eps()^(1/3), nλ = 10, rλ = 0.1, seed = rand(UInt64), prop_nknots = 1.0) where T <: AbstractFloat
    λs = range(1-rλ, 1+rλ, length = nλ + 1) * λstar
    D, μerr = cvfit(x, y, [μstar], λs, nfold = nfold, figname = figname, seed = seed, prop_nknots = prop_nknots)
    last_min_μerr = minimum(μerr)
    lastD = deepcopy(D)
    if abs(1 - D.λ / λstar) / rλ > 1 - ρ
        iter = 0
        while true 
            iter += 1
            rλ += 0.05
            nλ += 2
            if rλ <= 0 
                break
            end
            println("extend the search region of λ: searching [1-$rλ, 1+$rλ] * $λstar")
            λs = range(max(0, 1-rλ), 1+rλ, length = nλ + 1) * λstar
            D, μerr = cvfit(x, y, [μstar], λs, nfold = nfold, figname = figname, seed = seed, prop_nknots = prop_nknots)
            min_μerr = minimum(μerr)
            if min_μerr > last_min_μerr
                println("no better lambda after extension")
                return lastD, last_min_μerr
            end
            if (abs(1 - D.λ / λstar) / rλ <= 1 - ρ) | (iter > maxiter)
                break
            end
            last_min_μerr = min_μerr
            lastD = deepcopy(D)
        end
    end
    return D, last_min_μerr
end

"""
    Given `μmax`, and construct μs = (1:nμ) ./ nμ * μmax. If the optimal `μ` near the boundary, double or halve `μmax`.
"""
function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μmax::Real, λ::AbstractVector{T}; nfold = 10, nμ = 100, maxiter = 20, figname = "/tmp/cv_curve.png", ρ = 0.05, tol = eps()^(1/3), seed = rand(UInt64), prop_nknots = 1.0) where T <: AbstractFloat
    # μ = (1:nμ) ./ nμ * μmax
    D, μerr = cvfit(x, y, (1:nμ) ./ nμ .* μmax, λ, figname = figname, nfold = nfold, seed = seed, prop_nknots = prop_nknots)
    last_min_μerr = minimum(μerr)
    if ρ < 1/nμ
        @warn "the proportion to the boundary $ρ is better to be smaller than 1/nμ, otherwise the boundary solution is not detected"
    end
    # if D.μ near the left boundary
    if D.μ / μmax < ρ
        iter = 0
        lastD = deepcopy(D)
        while true
            iter += 1
            μmax = μmax / 2
            println("divide the search region of μ: searching [0, $μmax]")
            D, μerr = cvfit(x, y, (1:nμ) ./ nμ .* μmax, λ, figname = figname, nfold = nfold, seed = seed, prop_nknots = prop_nknots)
            min_μerr = minimum(μerr)
            if min_μerr > last_min_μerr
                println("no better mu after division")
                return lastD, last_min_μerr
            end
            if (D.μ / μmax >= ρ) | (iter > maxiter) # suppose it cannot be suddenly larger than 0.95 (TODO)
                break
            end
            last_min_μerr = min_μerr
            lastD = deepcopy(D)
        end
    elseif D.μ / μmax > 1 - ρ
        iter = 0
        μmin = 0
        step = μmax
        lastD = deepcopy(D)
        while true
            iter += 1
            μmax = μmax * 2
            # μmin = (μmax + μmin) / 2
            # μmax = μmin + step
            println("extend the search region of μ: searching [0, $μmax]")
            D, μerr = cvfit(x, y, (1:nμ) ./ nμ .* μmax, λ, figname = figname, nfold = nfold, seed = seed, prop_nknots = prop_nknots)
            # D, workspace, μerr = cvfit(x, y, μmin .+ (1:nμ) ./ nμ .* step, λ, figname = figname, nfold = nfold)
            min_μerr = minimum(μerr)
            if min_μerr > last_min_μerr
                println("no better mu after extension")
                return lastD, last_min_μerr
            end
            if ( (D.μ - μmin) / (μmax - μmin) <= 1-ρ) | (iter > maxiter) # suppose it cannot be suddenly smaller than 0.05 
                break
            end
            last_min_μerr = min_μerr
            lastD = deepcopy(D)
        end
    end
    return D, last_min_μerr
end

function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, paras::AbstractMatrix{T}; nfold = 10, figname = "/tmp/cv_curve.png", seed = rand(UInt64), same_J_after_CV = true, prop_nknots = 1.0) where T <: AbstractFloat
    # assume two columns, the first column is mu, and the 2nd column is lambda
    @assert size(paras, 2) == 2
    n = length(x)
    folds = div_into_folds(n, K = nfold, seed = seed)
    npara = size(paras, 1)
    err = zeros(nfold, npara)
    for k = 1:nfold
        test_idx = folds[k]
        train_idx = setdiff(1:n, test_idx)
        workspace = WorkSpaceSS()
        build_model!(workspace, x[train_idx], prop_nknots = prop_nknots)
        # ...............................................................................lambda.......mu
        γhats = _optim(y[train_idx], workspace.J, workspace.B, workspace.H, workspace.L, paras[:, 2], paras[:, 1])
        ynews = predict(workspace, x[test_idx], γhats[1:workspace.J, :] + γhats[workspace.J+1:2workspace.J, :])
        for j = 1:npara
            err[k, j] = norm(ynews[:, j] - y[test_idx]).^2 / length(test_idx)
        end
    end
    μerr = dropdims(mean(err, dims = 1), dims=1)
    σerr = dropdims(std(err, dims = 1), dims=1) / sqrt(nfold)
    if !isnothing(figname)
        silname = figname[1:end-4] * ".sil"
        serialize(silname, [μerr, σerr, paras, nfold])
        savefig(cvplot(μerr, σerr, paras[:, 1], nfold = nfold, lbl = "μ"), figname)
    end
    ind = argmin(μerr)
    workspace = WorkSpaceSS()
    # number of points in each fold
    n_fold = round(Int, n / nfold * (nfold - 1))
    nx_fold = rcopy(R".nknots.smspl($(n_fold))")
    nx = rcopy(R".nknots.smspl($(n))")
    D = mono_decomp_ss(workspace, x, y, paras[ind, 2], paras[ind, 1], 
                        prop_nknots = prop_nknots * ifelse(same_J_after_CV, nx_fold / nx, 1.0))
    return D, paras[:, 1], μerr, σerr
end

function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μ::AbstractVector{T}, λ::AbstractVector{T}; 
                    nfold = 10, figname = "/tmp/cv_curve.png", 
                    seed = rand(UInt64), 
                    same_J_after_CV = true, 
                    prop_nknots = 1.0, 
                    include_boundary = false, 
                    strict = false) where T <: AbstractFloat
    n = length(x)
    folds = div_into_folds(n, K = nfold, seed = seed)
    nλ = length(λ)
    nμ = length(μ)
    err = zeros(nfold, nμ, nλ)
    μs = repeat(μ, outer = nλ)
    λs = repeat(λ, inner = nμ)
    il = argmin(x)
    ir = argmax(x)
    # @showprogress "CV " for k = 1:nfold
    for k = 1:nfold
        test_idx = folds[k]
        train_idx = setdiff(1:n, test_idx)
        if include_boundary ## see x^{1/3} in paper#12
            train_idx = union(il, train_idx, ir)
            test_idx = setdiff(1:n, train_idx)
        end
        workspace = WorkSpaceSS()
        build_model!(workspace, x[train_idx], prop_nknots = prop_nknots)
        γhats = try
            _optim(y[train_idx], workspace.J, workspace.B, workspace.H, workspace.L, λs, μs, strict = strict)
        catch e
            error(e)
        end
        ynews = predict(workspace, x[test_idx], γhats[1:workspace.J, :] + γhats[workspace.J+1:2workspace.J, :])
        for j = 1:nλ
            for i = 1:nμ
                m = (j-1) * nμ + i
                err[k, i, j] = norm(ynews[:, m] - y[test_idx])^2 / length(test_idx)
            end
        end
    end
    μerr = dropdims(mean(err, dims = 1), dims=1)
    σerr = dropdims(std(err, dims = 1), dims=1) / sqrt(nfold)
    if !isnothing(figname)
        silname = figname[1:end-4] * ".sil"
        serialize(silname, [μerr, σerr, μ, λ, nfold])
        savefig(cvplot(μerr, σerr, μ, λ, nfold = nfold), figname)
    end
    ind = argmin(μerr)
    workspace = WorkSpaceSS()
    ## NB: in each fold, the number of points is fewer than the number of points in the all data
    ## and it results the number of basis function is different. Should we keep it the same?
    # number of points in each fold
    n_fold = round(Int, n / nfold * (nfold - 1))
    nx_fold = rcopy(R".nknots.smspl($(n_fold))")
    nx = rcopy(R".nknots.smspl($(n))")    
    D = mono_decomp_ss(workspace, x, y, λ[ind[2]], μ[ind[1]], 
                        prop_nknots = prop_nknots * ifelse(same_J_after_CV, nx_fold / nx, 1.0), strict = strict) 
    return D, μerr, σerr
end

# fix ratio
function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, ks::AbstractVector{T}, λstar::Real; multi_fix_ratio = false, kw...) where T <: AbstractFloat
    @assert maximum(ks) < 1
    paras = hcat(1 ./ ks .- 1, λstar ./ ks)
    if multi_fix_ratio
        @info "use multi_fix_ratio λstar = $λstar"
        paras = vcat([hcat(1 ./ ks .- 1, λstar * s ./ ks) for s in 10.0 .^ (-1:0.2:0)]...)
    end
    return cvfit(x, y, paras; kw...)
end

"""
    cvfit_gss(x, y, μrange, λ; λ_is_μ)

Cross-validation by Golden Section Searching `μ`` in `μrange` given `λ`.

- If `λ_is_μ`, search `λ` in `μrange` given `λ (μ)`
- Note that `one_se_rule` is not suitable for the golden section search.
"""
function cvfit_gss(x::AbstractVector{T}, y::AbstractVector{T}, μrange::AbstractVector{T}, λ::Real; 
                        nfold = 10, figname = "/tmp/cv_curve.png", 
                        seed = rand(UInt64), 
                        # relative tolerance
                        tol = 1e-4, tol_boundary = 1e-2,
                        λ_is_μ = false, prop_nknots = 1.0, 
                        log_scale = false,
                        rerun_check = true,
                        include_boundary = false, same_J_after_CV = false) where T <: AbstractFloat
    τ = (sqrt(5) + 1) / 2
    μs = Float64[]
    errs = Float64[]
    σerrs = Float64[]
    if log_scale
        a = log(μrange[1])
        b = log(μrange[2])
    else
        a = μrange[1]
        b = μrange[2]
    end
    μright = b
    δ = b - a # width of search region
    c = a
    d = b
    @debug "search $(ifelse(λ_is_μ, "λ", "μ")) in $([a, b])"
    iter = 0
    while true
        iter += 1
        ifigname = isnothing(figname) ? figname : figname[1:end-4] * "_$iter.png"
        iter % 1 == 0 && @debug "iter = $iter: narrow $(ifelse(λ_is_μ, "λ", "μ")) into $([a, b]), (b - a) / δ = $((b - a) / δ)"
        c = b - (b - a) / τ
        d = a + (b - a) / τ
        μi = [c, d]
        if log_scale
            μi = exp.(μi)
        end
        if λ_is_μ
            D, μerr, σerr = cvfit(x, y, [λ], μi, nfold = nfold, figname = ifigname, seed = seed, prop_nknots = prop_nknots, include_boundary = include_boundary, same_J_after_CV = same_J_after_CV)
        else
            D, μerr, σerr = cvfit(x, y, μi, [λ], nfold = nfold, figname = ifigname, seed = seed, prop_nknots = prop_nknots, include_boundary = include_boundary, same_J_after_CV = same_J_after_CV)
        end
        append!(μs, μi)
        append!(errs, μerr)
        append!(σerrs, σerr)
        if λ_is_μ
            # ideally, D.λ either equal to μi[1] or μi[2], but for numerical reason, use abs
            if abs(D.λ - μi[1]) < abs(D.λ - μi[2]) 
                b = d
            else
                a = c
            end
        else
            if abs(D.μ - μi[1]) < abs(D.μ - μi[2]) # f(c) < f(d)
                b = d
            else
                a = c
            end
        end
        if (b - a) / δ < tol
            @debug "μright = $μright, a = $a, δ = $δ, gap = $((μright - a) / δ), tol_boundary = $tol_boundary"
            flag = (μright - a) / δ < tol_boundary
            # rerun the last iteration to guarantee the solution is feasible in strict mode
            # not necessary? since some fit is indeed terrible, just return a bad fit (the constant fit)
            if rerun_check
                try
                    @debug "rerun cvfit to check feasibility"
                    if λ_is_μ
                        D, μerr, σerr = cvfit(x, y, [λ], μi, nfold = nfold, figname = ifigname, seed = seed, prop_nknots = prop_nknots, include_boundary = include_boundary, strict = true, same_J_after_CV = same_J_after_CV)
                    else
                        D, μerr, σerr = cvfit(x, y, μi, [λ], nfold = nfold, figname = ifigname, seed = seed, prop_nknots = prop_nknots, include_boundary = include_boundary, strict = true, same_J_after_CV = same_J_after_CV)
                    end        
                catch e
                    @debug "the found solution is not feasible after reruning optimization in the strict mode..."
                    flag = true
                end            
            end
            if flag
                # @warn "the optimal points near the right boundary"
                @debug "the optimal points near the right boundary (or not feasible), try to extend the range..." # since in the log file the warn message is before other print message
                # since 0.75 < 1 - 1e-2, so it is ok
                #left_new = μrange[1] + 0.75 * (μrange[2] - μrange[1])  # cannot larger than μrange[2]
                ## a possible issue: μrange[2] can be large after multiplying 10 or 2, and then left_new can also become large, then the total search region would be shifted to the right, it might miss the left part. An extreme case might be it always shifts to the right.
                left_new = μrange[1] # always keep the left range
                right_new = μrange[2] * ifelse(μrange[2] > 1e-2, 2, 10)
                return cvfit_gss(x, y, [left_new, right_new], λ, figname = figname, nfold = nfold, seed = seed, λ_is_μ = λ_is_μ, tol = tol, prop_nknots = prop_nknots, same_J_after_CV = same_J_after_CV, log_scale = log_scale, tol_boundary = tol_boundary, rerun_check = rerun_check)
            end
            ind = sortperm(μs)
            @debug "optimal $(ifelse(λ_is_μ, "λ", "μ")) = $(μs[end])" 
            if !isnothing(figname)
                savefig(plot(log.(μs[ind]), errs[ind], yerrors = σerrs[ind]), figname)
                # savefig(scatter(log.(μs), errs), "/tmp/cv.png")
            end
            return D, μs[ind], errs[ind], σerrs[ind]
        end
    end
end

"""
    cvfit_gss(x, y, μrange, λs)

For each `λ` in `λs`, perform `cvfit(x, y, μrange, λ)`, and store the current best CV error. Finally, return the smallest one.
"""
function cvfit_gss(x::AbstractVector{T}, y::AbstractVector{T}, μrange::AbstractVector{T}, λs::AbstractVector{T}; 
                    nfold = 10, figname = "/tmp/cv_curve.png", seed = rand(UInt64), 
                    tol = 1e-3, tol_boundary = 1e-1,
                    prop_nknots = prop_nknots, include_boundary = false, 
                    rerun_check = true,
                    same_J_after_CV = false, log_scale = false) where T <: AbstractFloat
    best_err = Inf
    D = nothing
    best_μs = nothing
    best_errs = nothing
    best_σerrs = nothing
    for (i, λ) in enumerate(λs)
        ifigname = isnothing(figname) ? figname : figname[1:end-4] * "_lambda$i.png"
        @debug "given λ = $λ, search μ..."
        Di, μi, erri, σerri = cvfit_gss(x, y, μrange, λ, nfold = nfold, figname = ifigname, seed = seed, tol = tol,
                                            tol_boundary = tol_boundary,
                                            prop_nknots = prop_nknots, 
                                            include_boundary = include_boundary, 
                                            same_J_after_CV = same_J_after_CV, 
                                            log_scale = log_scale, rerun_check = rerun_check)
        if best_err > minimum(erri)
            D = Di
            best_err = minimum(erri)
            best_μs = μi
            best_errs = erri
            best_σerrs = σerri
        end
        @debug "λ = $λ, err = $(minimum(erri)), best_err = $best_err"
    end
    return D, best_μs, best_errs, best_σerrs
end
