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
using MonotoneSplines
import MonotoneSplines.build_model
import MonotoneSplines.pick_knots
# using Distributed
# using SharedArrays

OPTIMIZER = ECOS.Optimizer

function gurobi()
    GRB_ENV = Gurobi.Env()
    # disable print of `set parameter ..`
    # find in https://github.com/jump-dev/Gurobi.jl/blob/master/src/gen91/libgrb_api.jl
    # and inspired by https://support.gurobi.com/hc/en-us/community/posts/4412624836753-Do-not-print-Set-parameter-Username-in-console
    GRBsetintparam(GRB_ENV, GRB_INT_PAR_OUTPUTFLAG, 0)
    global OPTIMIZER = () -> Gurobi.Optimizer(GRB_ENV)
end

function ecos()
    global OPTIMIZER = ECOS.Optimizer
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
function smooth_spline(x::AbstractVector{T}, y::AbstractVector{T}, xnew::AbstractVector{T}; keep_stuff = false, design_matrix = false) where T <: AbstractFloat
    spl = R"smooth.spline($x, $y, keep.stuff = $keep_stuff)"
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

"""
    mono_decomp_cs(x::AbstractVector, y::AbstractVector)

Monotone Decomposition with Cubic B-splines by solving an optimization problem.
"""
function mono_decomp_cs(x::AbstractVector{T}, y::AbstractVector{T}; 
                    s = 1.0, s_is_μ = true,
                    J = 4,
                    workspace = nothing,
                    ) where T <: AbstractFloat
    if isnothing(workspace) || !workspace.evaluated
        workspace = WorkSpaceCS()
        B, rB = build_model(x, J)
        workspace.evaluated = true
        workspace.B = B
        workspace.rB = rB
        workspace.J = J
        workspace.H = construct_H(J)
    end
    if s_is_μ
        γhat = _optim(y, workspace, [s])[:]
    else
        γhat = _optim(y, workspace.J, workspace.B, s, workspace.H)
    end
    γup = γhat[1:J]
    γdown = γhat[J+1:2J]
    yhat = B * (γup .+ γdown)
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
            errs[k, j] = norm(yhat - y[test_idx])
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
    if !isnothing(figname)
        savefig(cvplot(μerr, σerr, 1.0 .* Js, nfold = nfold, ind0 = ind, lbl = "J"), figname)
    end
    yhat, yhatnew = cubic_spline(Jopt)(x, y, xnew)
    return Jopt, yhat, yhatnew
end

abstract type WorkSpace end

# difference compared to MonoDecomp
# workspace can be shared with different tuning parameters
# MonoDecomp records the specific results with particular parameters
mutable struct WorkSpaceSS <: WorkSpace
    evaluated::Bool
    J::Int
    B::AbstractMatrix{Float64}
    L::AbstractMatrix{Float64}
    H::AbstractMatrix{Int}
    mx::Real
    rx::Real
    idx::AbstractVector{Int}
    idx0::AbstractVector{Int}
    Bend::AbstractMatrix{Float64}
    Bendd::AbstractMatrix{Float64}
    knots::AbstractVector{Float64}
    # Bnew::AbstractMatrix{Float64} # not for cross-validation
    WorkSpaceSS() = new()
    # WorkSpaceSS(e, j, b, l, h, m, r, i, i0, be, bd, k) = new(e, j, b, l, h, m, r, i, i0, be, bd, k)
end

mutable struct WorkSpaceCS <: WorkSpace
    evaluated::Bool
    J::Int
    B::AbstractMatrix{Float64}
    rB::RObject
    H::AbstractMatrix{Int}
    WorkSpaceCS() = new()
end

"""
    cv_mono_decomp_cs(x::AbstractVector, y::AbstractVector, xnew::AbstractVector; )
    cv_mono_decomp_cs(x::AbstractVector, y::AbstractVector; fixJ = true)

Cross-validation for Monotone Decomposition with Cubic B-splines. Parameters `J` and `s` (`μ` if `s_is_μ`) are tuned by cross-validation.

- if `fixJ == true`, then `J` is CV-tuned by the corresponding cubic B-spline fitting
- if `fixJ == false`, then both `J` and `s` would be tuned by cross-validation.
"""
function cv_mono_decomp_cs(x::AbstractVector{T}, y::AbstractVector{T}, xnew::AbstractVector{T}; nfold = 10, Js = 4:50, ss = 10.0 .^ (-6:0.1:-1), s_is_μ = true, figname = nothing, one_se_rule = false) where T <: AbstractFloat
    n = length(x)
    folds = div_into_folds(n, K = nfold, seed = rand(UInt64))
    errs = zeros(nfold, length(Js), length(ss))
    n2 = length(folds[1])
    n1 = n - n2
    @assert n % nfold == 0
    for k = 1:nfold
        test_idx = folds[k]
        train_idx = setdiff(1:n, test_idx)
        for (j, J) in enumerate(Js)
            workspace = WorkSpaceCS()
            workspace.J = J
            build_model!(workspace, x[train_idx])
            # workspace = nothing
            @assert s_is_μ
            γhats = _optim(y[train_idx], workspace, ss)
            ynews = predict(workspace, x[test_idx], γhats[1:J, :] + γhats[J+1:2J, :])
            for (i, s) in enumerate(ss)
                errs[k, j, i] = norm(ynews[:, i] - y[test_idx]).^2 / length(test_idx)
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
    if !isnothing(figname)
        silname = figname[1:end-4] * "_Jmu.sil"
        serialize(silname, [μerr, σerr, Jopt, sopt, Js, ss, nfold])
        savefig(cvplot(μerr, σerr, 1.0 .* Js, ss, lbl = ["J", "μ"], nfold = nfold, ind0 = ind), figname)
    end
    D = mono_decomp_cs(Jopt)(x, y; s = sopt)
    return D, ss[argmin(μerr)[2]], ss, μerr, σerr  # use for recording without 1se
end

"""
    cv_mono_decomp_cs(x::AbstractVector, y::AbstractVector)

Cross-validation for monotone decomposition with cubic B-splines when the fixed `J` is CV-tuned by the corresponding cubic B-spline fitting method.
"""
function cv_mono_decomp_cs(x::AbstractVector{T}, y::AbstractVector{T}; ss = 10.0 .^ (-6:0.1:-1), 
                                                            figname = nothing, 
                                                            nfold = 10, 
                                                            fixJ = true,
                                                            x0 = x,
                                                            Js = 4:50,
                                                            one_se_rule = false) where T <: AbstractFloat
    if fixJ
        J, yhat, yhatnew = cv_cubic_spline(x, y, x0, one_se_rule = one_se_rule, nfold = nfold, figname = isnothing(figname) ? figname : figname[1:end-4] * "_bspl.png")
        return cv_mono_decomp_cs(x, y, x0, Js = J:J, ss = ss, figname = figname, nfold = nfold, one_se_rule = one_se_rule)..., yhat, yhatnew
    else
        return cv_mono_decomp_cs(x, y, x0, ss = ss, Js = Js, figname = figname, nfold = nfold, one_se_rule = one_se_rule)
    end
end

"""
    cv_mono_decomp_ss(x::AbstractVector, y::AbstractVector)

Cross Validation for Monotone Decomposition with Smoothing Splines. with `λ` tuned by smoothing spline, and then perform golden search for `μ`.

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
                                                            nfold = 10, 
                                                            tol = 1e-7,
                                                            x0 = x, # test point of x
                                                            method = "single_lambda", # fix_ratio, iter_search, grid_search
                                                            nk = 20, #used in fix_ratio
                                                            nλ = 10, rλ = 0.5, # used in grid_search
                                                            rel_tol = 1e-1, maxiter = 10, # iter_search
                                                            k_magnitude = 2,
                                                            one_se_rule = false, kw...) where T <: AbstractFloat
    yhat, yhatnew, Ω, λ, spl, B = smooth_spline(x, y, x0, design_matrix = true, keep_stuff = true)
    γup, γdown = mono_decomp(rcopy(R"$spl$fit$coef"))
    s0 = norm(B * (γup .- γdown))
    s_residual = norm(y - yhat)^2
    s_smoothness = λ * (γup .+ γdown)' * Ω * (γup .+ γdown)
    s_discrepancy = [max(eps(), min(s_residual, s_smoothness)) / 10^k_magnitude, 
                                max(s_residual, s_smoothness) * 10^k_magnitude]
    μrange = s_discrepancy / s0^2
    verbose && @info "μrange: $μrange"
    if method == "single_lambda"
        verbose && @info "Smoothing Splines with fixed λ"
        D, μs, errs, σerrs = cvfit_gss(x, y, μrange, λ, figname = figname, nfold = nfold, tol = tol)
    elseif method == "fix_ratio"
        verbose && @info "Smoothing Splines with fixratio strategy"
        ks = range(0.05, 0.99, length = nk)
        D, μs, errs, σerrs = cvfit(x, y, ks, λ, figname = figname, nfold = nfold)
    elseif method == "iter_search" #88
        verbose && @info "Smoothing Splines with iter-search: λ -> μ -> λ -> ... -> μ"
        iter = 0
        λ0 = λ # make a backup
        seed = rand(UInt64)
        while true
            iter += 1
            ## tune mu given lambda
            @debug "tune mu given lambda = $λ"
            # D1, workspace1 = cvfit(x, y, μmax * r, [λ], nfold = nfold, figname = figname, nμ = nμ, ρ = ρ)
            ## if needed, perform one se rule on the last iteration given lambda
            D1, μs, errs, σerrs = cvfit_gss(x, y, μrange, λ, nfold = nfold, figname = figname, tol = tol, seed = seed)
            # since D is not defined in the first iteration, so use `if..else`, and hence cannot use `ifelse`
            if iter == 1
                err_μ = 1.0
            else
                err_μ = abs(D1.μ - D.μ) / D1.μ
            end
            if err_μ < rel_tol
                D = D1
                break
            end
            ## re-tune lambda given mu
            @debug "tune lambda given mu = $(D1.μ)"
            # D, workspace = cvfit(x, y, D1.μ, λ, nfold = nfold, figname = figname, nλ = nλ, ρ = ρ)
            D, _ = cvfit_gss(x, y, [1e-7, 1.5λ0], D1.μ, nfold = nfold, figname = figname, λ_is_μ = true, tol = tol, seed = seed)
            err_λ = abs(D.λ - D1.λ) / D.λ
            λ = D.λ # for next iteration
            @debug "iter = $iter, err_μ = $err_μ, err_λ = $err_λ"
            if (iter > maxiter) | (err_λ < rel_tol)
                break
            end
        end
    else # grid_search
        λs = range(1-rλ, 1+rλ, length = nλ) .* λ
        verbose && @info "Smoothing Splines with grid-search λ ∈ $λs, μ ∈ $μrange"
        seed = rand(UInt64)
        D, μs, errs, σerrs = cvfit_gss(x, y, μrange, λs, nfold = nfold, figname = figname, seed = seed, tol = tol)
    end
    if one_se_rule
        ind = cv_one_se_rule(errs, σerrs, small_is_simple = false)
        μopt = μs[ind]
        @debug "use 1se rule: before 1se: μ = $(D.μ); after 1se: μ = $μopt"
        D = mono_decomp_ss(D.workspace, x, y, D.λ, μopt)
    end
    return D, μs[argmin(errs)], μs, errs, σerrs, yhat, yhatnew
end

function summary_res(σs = 0.2:0.2:1.0)
    # combine cv plot
    fignames = "/tmp/1-cv_optim_" .* string.(σs) .* ".png"
    run(`convert $fignames +append /tmp/cv_optim.png`)
    # combine curve fit plot
    fignames = "/tmp/1-optim_sigma" .* string.(σs) .* ".png"
    run(`convert $fignames +append /tmp/curve_optim.png`)
end

function benchmarking_cs(n::Int = 100, σ::Float64 = 0.5, f::Union{Function, String} = x->x^3; fixJ = true,
                                                                               figname_cv = nothing,
                                                                               figname_fit = nothing,
                                                                               Js = 4:20,
                                                                               nfold = 10,
                                                                               one_se_rule = false,
                                                                               μs = 10.0 .^ (-6:0.5:0))
    x, y, x0, y0 = gen_data(n, σ, f)
    # J is determined from cubic_spline (deprecated the choice of arbitrary J)
    J, yhat, yhatnew = cv_cubic_spline(x, y, x0, nfold = nfold, one_se_rule = one_se_rule)
    if fixJ
        Js = J:J
    end
    D, _ = cv_mono_decomp_cs(x, y, x0, Js = Js, ss = μs, figname = figname_cv, nfold = nfold, one_se_rule = one_se_rule)
    err = [norm(D.yhat - y)^2 / length(y), norm(predict(D, x0) - y0)^2 / length(y0), 
           norm(yhat - y)^2 / length(y), norm(yhatnew - y0)^2 / length(y0),
           var(y), var(y0)]
    if !isnothing(figname_fit)
        savefig(
            plot([x, y], [x0, y0], D, yhatnew, lbl_other = "Bspl (J = $J)", legend = :top),
            figname_fit
            )
    end
    return err
end

function benchmarking_ss(n::Int = 100, σ::Float64 = 0.5, 
                            f::Union{Function, String} = x->x^3; 
                                nfold = 5, one_se_rule = true,
                                method = "single_lambda",
                                figname_cv = nothing,
                                figname_fit = nothing, kw...
                        )
    x, y, x0, y0 = gen_data(n, σ, f)
    D, μopt, μs, errs, σerrs, yhat, yhatnew = cv_mono_decomp_ss(x, y; x0 = x0,
                                                                figname = figname_cv,
                                                                nfold = nfold,
                                                                method = method,
                                                                one_se_rule = one_se_rule, kw...)
    if !isnothing(figname_fit)
        savefig(
            plot([x, y], [x0, y0], D, yhatnew),
            figname_fit
        )
    end
    yup, ydown = predict(D.workspace, x0, D.γup, D.γdown)
    err = [norm(D.yhat - y)^2 / length(y), norm(yup + ydown - y0)^2 / length(y0), 
            norm(yhat - y)^2 / length(y), norm(yhatnew - y0)^2 / length(y0),
            var(y), var(y0)]
    return err
end

"""


- `competitor`: 
    if `nλ = 1`, then fixed λ; otherwise, there are two choices: ss_grid, ss_iter

"""
function benchmarking(f::String = "x^3"; n = 100, σs = 0.2:0.2:1,
                            J = 10, 
                            μs = 10.0 .^ (-6:0.5:0), ## bspl
                            jplot = false, nrep = 100,
                            competitor = "ss_single_lambda", # bspl, # ss_fix_ratio, ss_grid_search, ss_iter_search
                            nfold = 5, one_se_rule = true,
                            resfolder = "/tmp",
                            ind = 1:4,
                            show_progress = true,
                            nλ = 20, kw...)
    @info "Benchmarking $f with $nrep repetitions"
    title = "$f (nrep = $nrep)"
    filename = "$f.sil"
    nσ = length(σs)
    res = zeros(nrep, 6, nσ)
    # res = SharedArray{Float64}(nrep, 4, nσ)
    if f == "x^3"
        f = x -> x^3
    elseif f == "x^2"
        f = x -> x^2
    elseif f == "exp(x)"
        f = x -> exp(x)
    elseif f == "sigmoid"
        f = x -> 1 / (1 + exp(-5x))
    end
    for i = 1:nrep
        p = Progress(nσ, dt = 1, desc = "$title, iter = $i: ", enabled = show_progress)
        for (j, σ) in enumerate(σs)
            figname_fit = ifelse(jplot, joinpath(resfolder, "$i-optim_sigma$σ.png"), nothing)
            figname_cv = ifelse(jplot, joinpath(resfolder, "$i-cv_optim_$σ.png"), nothing)
            if startswith(competitor, "ss")
                res[i, :, j] = benchmarking_ss(n, σ, f; figname_cv = figname_cv, 
                                                        figname_fit = figname_fit,
                                                        nfold = nfold,
                                                        nλ = nλ,
                                                        one_se_rule = one_se_rule,
                                                        method = competitor[4:end], kw...)
            else
                res[i, :, j] = benchmarking_cs(n, σ, f; figname_cv = figname_cv, 
                                                        figname_fit = figname_fit,
                                                        μs = μs, 
                                                        one_se_rule = one_se_rule,
                                                        fixJ = !occursin("cvbspl2", competitor), kw...)
            end
            next!(p)
        end
    end
    # return res
    μerr = mean(res, dims = 1)[1, :, :]
    σerr = std(res, dims = 1)[1, :, :]
    # paste figures under different noise levels together
    if jplot
        summary_res(σs)
    end
#    serialize("../res/res_monodecomp/$subfolder/$filename", res) # used in #86
    serialize(joinpath(resfolder, filename), res)
    if jplot
        return plot(σs, μerr'[:, ind], label = ["yup+ydown - y" "yup+ydown - f" "yhat - y" "yhat - f"][1:1, ind],
            yerror = σerr'[:, ind] / sqrt(nrep) * 2,
            legend = :top,
            ls = reshape(repeat([:dash, :solid], outer = 2), 1, 4)[1:1, ind],
            markershape = reshape(repeat([:dtriangle, :star5], inner = 2), 1, 4)[1:1, ind],
            xlab = "sigma", title = title)
    end
end

function summary(;nλ = 20, 
                    format = "tex", 
                    methodnames = ["SmoothSpline", "MonoDecomp"], # CubicSpline
                    # σs = 0.2:0.2:2.0,
                    σs = [0.1, 0.2, 0.4, 0.5, 1.0, 1.5, 2.0],
                    nrep = 100, ind = 1:4,
                    curves = ["x^2" "x^3" "exp(x)" "sigmoid" "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"],
                    resfolder = "/tmp" #"../res/res_monodecomp/rocky/"
                    )
    titles = curves
    # titles to display
    dtitles = copy(titles)
    if length(titles) > 2
        dtitles[1:3] = [L"x^2" L"x^3" L"\exp(x)"]
    end
    ntitle = length(titles)
    @info "summarize result $resfolder in $format format"
    if resfolder == "/tmp"
        filename = "tmp"
        resfolder0 = resfolder    
    elseif endswith(resfolder, "/")
        filename = basename(resfolder[1:end-1])
        resfolder0 = dirname(resfolder[1:end-1])
    else
        filename = basename(resfolder)
        resfolder0 = dirname(resfolder)
    end
    try
        # get parameters from filename
        list_paras = split(filename, "-")
        inrep = findfirst(startswith.(list_paras, "nrep"))
        inlam = findfirst(startswith.(list_paras, "nlam"))
        nrep = parse(Int, split(list_paras[inrep], "nrep")[2])
        nλ = parse(Int, split(list_paras[inlam], "nlam")[2])
    catch e
        @warn "cannot parse nλ and nrep from resfolder, use default nλ=$nλ, nrep=$nrep"
    end
    if format == "figure"
        @warn "Figure is less accurate and hence deprecated. Use table instead. See Figure examples https://github.com/szcf-weiya/Clouds/issues/80#issuecomment-1025306738"
        #=
        figs = Plots.Plot[]
        for (i, title) in enumerate(titles)
            # resfile = "$(title)_(nrep_=_$nrep)sig$(σs)-nlam$nλ-rlam0.5.sil"
            resfile = "$title.sil"
            # resfile = "$(title)_(nrep_=_100)sig0.2:0.2:2.0-s0.05:0.05:1.0-nlam$nλ-rlam0.5.sil"
            res = deserialize(joinpath(resfolder, resfile))
            μerr = mean(res, dims = 1)[1, :, :]
            σerr = std(res, dims = 1)[1, :, :] / sqrt(nrep)
            # title = titles[i]
            push!(figs, plot(σs, μerr'[:, ind], 
                # label = ["yup+ydown - y" "yup+ydown - f" "yhat - y" "yhat - f"],
                label = [L"\Vert\hat y_u + \hat y_d - y\Vert" L"\Vert\hat f_u + \hat f_d - f\Vert" L"\Vert\hat y - y\Vert" L"\Vert \hat f-f\Vert"][1:1, ind],
                # yerror = σerr' * 2,
                legend = :top,
                tickfontsize = 11, #8
                titlefontsize = 17, #14
                guidefontsize = 14, #11
                legendfontsize = 11, #8
                markersize = 7,#4
                ls = reshape(repeat([:dash, :solid], outer = 2), 1, 4)[1:1, ind],
                markershape = reshape(repeat([:dtriangle, :star5], inner = 2), 1, 4)[1:1,ind],
                xlab = L"\sigma", title = dtitles[i])
            )
        end
        if Plots.backend() == Plots.PGFPlotsXBackend()
            save_plots(figs)
            mv("/tmp/all.pdf", "res_mono_decomp/$filename.pdf", force = true)
        else
            save_grid_plots(figs)
        end
        =#
    else
        # TODO: is it necessary?
        μerrs = Array{AbstractMatrix{Float64}}(undef, ntitle)
        σerrs = Array{AbstractMatrix{Float64}}(undef, ntitle)
        # pvals = Array{AbstractMatrix{Float64}}(undef, ntitle)
        pvals = Array{AbstractVector}(undef, ntitle)
        props = Array{AbstractVector}(undef, ntitle)
        snr_μs = Array{AbstractMatrix{Float64}}(undef, ntitle)
        snr_σs = Array{AbstractMatrix{Float64}}(undef, ntitle)
        isbf = Array{AbstractMatrix{Bool}}(undef, ntitle)
        # selected noise levels
        if length(σs) >= 10        
            if σs == 0.2:0.2:2
                ind = [1, 2, 3, 5] # cvbspl2 (sigmas: 0.2:0.2:2.0)
            elseif σs == 0.1:0.1:2
                ind = 1:20
            else
                # ind = [1, 5, 10]
                ind = [1, 2, 5, 10] # cvbspl
            end
        else
            # ind = [1, 2, 3]
            ind = 1:length(σs)
            # ind = [1, 2, 4, 5, 6, 7]
        end
        for (i, title) in enumerate(titles)
            resfile = "$(title).sil"
            res = deserialize(joinpath(resfolder, resfile))
            # the order should be consistent with the colnames
            # res[:, 3, :] = #max.(1 .- res[:, 3, :] ./ res[:, 5, :], 0) # R^2
            # res[:, 3, :] = max.(res[:, 5, :] ./ res[:, 3, :] .- 1, 0)# SNR
            # res[:, 3, :] = res[:, 6, :] ./ res[:, 3, :]
            # correction for commits before 372ab2f8a8e47b52fe556abdb0c3f3ab0a37dd1b (2022-02-13)
            # res[:, 1:4, :] .= res[:, 1:4, :] .^ 2 / nrep
            if size(res, 2) > 4
                snr = max.(res[:, 5, :] ./ res[:, 3, :] .- 1, 0)
                snr_μs[i] = mean(snr, dims = 1)[1:1, ind]'
                snr_σs[i] = std(snr, dims = 1)[1:1, ind]' / sqrt(nrep)
            else
                snr_μs = nothing
                snr_σs = nothing
            end
            # differences
            ds = res[:, 2, :] - res[:, 4, :]
            # count the number of minus
            props[i] = sum(ds .<= 0, dims = 1)[1, :] / nrep
            dμ = mean(ds, dims = 1)[1, :]
            dσ = std(ds, dims = 1)[1, :] / sqrt(nrep)
            pval = 1 .- cdf(Normal(0, 1), abs.(dμ) ./ dσ)
            pvals[i] = star_pval(pval[ind])
            μerrs[i] = mean(res, dims = 1)[1, [3, 1, 4, 2], ind]'
            σerrs[i] = std(res, dims = 1)[1, [3, 1, 4, 2], ind]' / sqrt(nrep)
            isbf[i] = zeros(length(ind), 4)
            isbf[i][:, 3] .= μerrs[i][:, 3] .< μerrs[i][:, 4]
            isbf[i][:, 4] .= .!isbf[i][:, 3]
        end
        # `_` is not allowed in non-math env
        dtitles = replace.(dtitles, "_" => "-")
        # print2tex, similar to https://github.com/szcf-weiya/Clouds/search?q=print2tex
        print2tex(μerrs, σerrs, dtitles[:], ["MSFE", "MSPE"], string.(σs[ind]), methodnames, 
            other_cols = snr_μs,
            other_cols_σ = snr_σs,
            other_col_names = ["SNR"],
            right_cols = [pvals, props],
            right_col_names = ["p-value", "prop."], 
            right_align = "l", # it might be slightly worse for other columns, but center is not good for the pvalue column TODO: improve the style
            colnames_of_rownames = ["curve", L"\sigma"], file = joinpath(resfolder0, "$filename.tex"), isbf = isbf)
    end
end



"""
    _optim(y::AbstractVector, workspace::WorkSpaceCS, μs::AbstractVector)
    _optim(y::AbstractVector, J::Int, B::AbstractMatrix, H::AbstractMatrix{Int}, μs::AbstractVector)

Optimization for monotone decomposition with cubic B-splines.

    _optim(y::AbstractVector, J::Int, B::AbstractMatrix, H::AbstractMatrix{Int}, L::AbstractMatrix, λs::AbstractVector, μs::AbstractVector)

Optimization for monotone decomposition with smoothing splines.

    _optim!(y::AbstractVector, J::Int, B::AbstractMatrix, s::Union{Nothing, Real}, γhat::AbstractVector, H::AbstractMatrix{Int}; L, t, λ, μ)


"""
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
                ) where T <: AbstractFloat
    model = Model(OPTIMIZER)
    set_silent(model)
    # TODO: test the speed and possibly switch to it if necessary, as in _optim!
    # set_optimizer_attribute(model, "BarHomogeneous", 1) # use only in demo since practically if a numerical error occurs, although it can be solved by another method, but it also implies that the current parameter is not the best
    #set_optimizer_attribute(model, "max_iter", 10000)
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
            @warn "$status when λ=$λ, μ=$μ: use constant half mean as its estimate"
            γhats[:, i] .= mean(y) / 2
        else
            if !(status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL])
                @warn "$status when λ=$λ, μ=$μ: direct take the solution"
            end
            try
                γhats[:, i] .= value.(γ) # ITERATION_LIMIT
            catch
                @warn "when λ=$λ, μ=$μ: no solution, $status"
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
                @warn "$status"
                γhats[:, i] .= mean(y) / 2
            else
                if !(status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL])
                    @warn "$status"
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
                @warn "$status"
                γhats[:, i] .= mean(y) / 2
            else
                if !(status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL])
                    @warn "$status"
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

function _optim!(y::AbstractVector{T}, J::Int, B::AbstractMatrix{T}, s::Union{Nothing, Real}, γhat::AbstractVector{T}, H::AbstractMatrix{Int}; L = nothing, t = nothing, λ = nothing, μ = nothing, BarHomogeneous = false) where T <: AbstractFloat
    # construct constraint ‖ B(γ^u - γ^d) ‖ < s
    # model = Model(GLPK.Optimizer) # cannot do second order, throw error
    # `MOI.VectorAffineFunction{Float64}`-in-`MOI.SecondOrderCone` constraints are not supported and cannot be bridged into supported constrained variables and constraints. See details below:
    model = Model(OPTIMIZER)
    if BarHomogeneous && typeof(model.moi_backend.optimizer).parameters[1] == Gurobi.Optimizer
        set_optimizer_attribute(model, "BarHomogeneous", 1)
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
        if !BarHomogeneous # not yet try BarHomogeneous
            @info "try BarHomogeneous to rerun NUMERICAL_ERROR"
            _optim!(y, J, B, s, γhat, H, L = L, t = t, λ = λ, μ = μ, BarHomogeneous = true)
        else
            @warn "$status after trying BarHomogeneous algorithm"
            γhat .= mean(y) / 2
        end
    else
        if !(status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL])
            @warn "$status"
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

struct MonoDecomp{T <: AbstractFloat}
    γup::AbstractVector{T}
    γdown::AbstractVector{T}
    γhat::AbstractVector{T}
    yhat::AbstractVector{T}
    λ::Union{Real, Nothing}
    μ::Union{Real, Nothing}
    workspace::WorkSpace
end

function build_model!(workspace::WorkSpaceSS, x::AbstractVector{T}; ε = (eps())^(1/3)) where T <: AbstractFloat
    if !isdefined(workspace, 3) # the first two will be automatically initialized
        knots, mx, rx, idx, idx0 = pick_knots(x, all_knots = false)
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
            @warn e
            ## perform pivoted Cholesky
            workspace.L = Matrix(cholesky(Symmetric(Ω), Val(true), check = false, tol = ε).L)
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
function mono_decomp_ss(workspace::WorkSpaceSS, x::AbstractVector{T}, y::AbstractVector{T}, λ::AbstractFloat, s::AbstractFloat; s_is_μ = true) where T <: AbstractFloat
    build_model!(workspace, x)
    if s_is_μ
        γhat = _optim(y, workspace.J, workspace.B, nothing, workspace.H, L = workspace.L, t = nothing, λ = λ, μ = s)
    else
        γhat = _optim(y, workspace.J, workspace.B, s, workspace.H, L = workspace.L, t = nothing, λ = λ)
    end
    # calculate properties of monotone decomposition
    J = workspace.J
    γup = γhat[1:J]
    γdown = γhat[J+1:2J]
    γhat = γup + γdown
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
    Bnew = rcopy(R"splines::bs($xm, intercept = TRUE, knots=$(W.knots[2:end-1]))")
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
        Bnew = rcopy(R"splines::bs($xm, intercept = TRUE, knots=$(W.knots[2:end-1]))")
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
    if competitor == "ss"
        # lbl_other = L"\hat y_{ss}"
        # lbl_other = latexstring("\$\\hat y_{ss} ($e2)\$")
        lbl_other = latexstring("\$\\hat f_{\\mathrm{ss}} ($e2)\$")
    end
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
    fig = scatter(x, y; title = prefix_title * title * postfix_title, label = "", ms = 2, xlab = L"x", ylab = L"y", kw...)
    # lbls = [L"\hat f_u", L"\hat f_d", latexstring("\$\\hat y_u+\\hat y_d ($e1)\$"), "truth"]
    lbls = [L"\hat f_{\mathrm{up}}", L"\hat f_{\mathrm{down}}", latexstring("\$\\hat f_{\\mathrm{up}}+\\hat f_{\\mathrm{down}} ($e1)\$"), L"f"]
    plot!(fig, x0, yup, label = lbls[1], ls = :dot)
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
function build_model!(workspace::WorkSpaceCS, x::AbstractVector{T}) where T <: AbstractFloat
    J = workspace.J
    rB = R"splines::bs($x, df=$J, intercept=TRUE)"
    workspace.B = rcopy(rB)
    workspace.rB = rB
    workspace.evaluated = true
    workspace.H = construct_H(J)
end

"""
    cv_one_se_rule(μs, σs)

Return the index of parameter that minimize the CV error with one standard error rule.
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
            continue
        else
            if small_is_simple
                iopt = i + 1
            else
                iopt = i - 1
            end
            break
        end
    end
    return iopt
end

function cv_one_se_rule_deprecated(μs::AbstractMatrix{T}, σs::AbstractMatrix{T}; verbose = true, small_is_simple = [true, true]) where T <: AbstractFloat
    ind = argmin(μs)
    iopt, jopt = ind[1], ind[2]
    # fix iopt, find the minimum jopt
    # TODO: the order might make sense
    jopt1 = cv_one_se_rule(μs[iopt,:], σs[iopt,:], small_is_simple = small_is_simple[2])
    iopt1 = cv_one_se_rule(μs[:,jopt1], σs[:,jopt1], small_is_simple = small_is_simple[1])

    iopt2 = cv_one_se_rule(μs[:,jopt], σs[:,jopt], small_is_simple = small_is_simple[1])
    jopt2 = cv_one_se_rule(μs[iopt2,:], σs[iopt2,:], small_is_simple = small_is_simple[2])
    
    if verbose
        println("without 1se rule: iopt = $iopt, jopt = $jopt")
    end
    # determine the minimum one
    if iopt1 + jopt1 < iopt2 + jopt2
        if μs[iopt1, jopt1] < μs[iopt2, jopt2] + σs[iopt2, jopt2]
            return iopt1, jopt1
        else
            return iopt2, jopt2
        end
    else
        if μs[iopt2, jopt2] < μs[iopt1, jopt1] + σs[iopt1, jopt1]
            return iopt2, jopt2
        else
            return iopt1, jopt1
        end
    end
end

function cv_one_se_rule(μs::AbstractMatrix{T}, σs::AbstractMatrix{T}; verbose = true, small_is_simple = [true, true]) where T <: AbstractFloat
    ind = argmin(μs)
    iopt, jopt = ind[1], ind[2]
    jopt1 = cv_one_se_rule(μs[iopt,:], σs[iopt,:], small_is_simple = small_is_simple[2])

    iopt1 = cv_one_se_rule(μs[:,jopt], σs[:,jopt], small_is_simple = small_is_simple[1])

    if verbose
        println("without 1se rule: iopt = $iopt, jopt = $jopt")
    end
    # determine the minimum one: compare (iopt1, jopt) vs (iopt, jopt1), take the larger one since both are in 1se
    if μs[iopt1, jopt] < μs[iopt, jopt1]
        return iopt, jopt1
    else
        return iopt1, jopt
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
        if nfold < 100
            title = "$nfold-CV Error"
        end
    else
        μerr, σerr, μ, λ, nfold = deserialize(sil)
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
        offset = 10
        cvplot(μerr[1:end-offset, :], σerr, μ[1:end-offset], λ, nfold = nfold, lbl = [L"\log_{10}\mu", L"\lambda"], title = title)
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
            if isnothing(lbl)
                lbl = L"\log_{10}"
            end
        end
    end
    ind = argmin(μerr)
    if isnothing(title)
        title = "$nfold-fold CV error"
    end
    if isnothing(σerr)
        p = plot(f.(paras), μerr, xlab = lbl, title = title, legend = false, markershape = :x)
    else
        p = plot(f.(paras), μerr, yerrors = σerr, xlab = lbl, title = title, legend = false)
    end
    annotate!(p, [(f.(paras)[ind], μerr[ind], ("o", :red))])
    if !isnothing(ind0)
        annotate!(p, [(f.(paras)[ind0], μerr[ind0], ("o", :red))])
    end
    return p
end

# assume matrix does not reduce to vector
function cvplot(μerr::AbstractMatrix{T}, σerr::AbstractMatrix{T}, para1::AbstractVector{T}, para2::AbstractVector{T}; lbl = ["", ""], title = "", nfold = 10, ind0 = [nothing, nothing]) where T <: AbstractFloat
    n1 = length(para1)
    n2 = length(para2)
    if n1 == 1
        return cvplot(μerr[1, :], σerr[1, :], para2, nfold = nfold, ind0 = ind0[2])
    elseif n2 == 1
        return cvplot(μerr[:, 1], σerr[:, 1], para1, nfold = nfold, ind0 = ind0[1])
    end
    n, m = size(μerr)
    @assert n == length(para1)
    @assert m == length(para2)
    if para1[2] - para1[1] ≈ para1[3] - para1[2]
        f = x -> x
    else
        f = x -> log10(x)
    end
    if para2[2] - para2[1] ≈ para2[3] - para2[2]
        g = x -> x
    else
        g = x -> log10(x)
    end
    p = heatmap(g.(para2), f.(para1), μerr, xlab = lbl[2], ylab = lbl[1], title = title)
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
function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μstar::Real, λstar::Real; nfold = 10, maxiter = 20, figname = "/tmp/cv_curve.png", ρ = 0.05, tol = eps()^(1/3), nλ = 10, rλ = 0.1, seed = rand(UInt64)) where T <: AbstractFloat
    λs = range(1-rλ, 1+rλ, length = nλ + 1) * λstar
    D, μerr = cvfit(x, y, [μstar], λs, nfold = nfold, figname = figname, seed = seed)
    last_min_μerr = minimum(μerr)
    lastD = deepcopy(D)
    if abs(1 - D.λ / λstar) / (2rλ) < ρ
        iter = 0
        while true 
            iter += 1
            rλ += 0.05
            nλ += 2
            if rλ <= 0 
                break
            end
            println("extend the search region of λ: searching [1-$rλ, 1+$rλ] * $λstar")
            λs = range(1-rλ, 1+rλ, length = nλ + 1) * λstar
            D, μerr = cvfit(x, y, [μstar], λs, nfold = nfold, figname = figname, seed = seed)
            min_μerr = minimum(μerr)
            if min_μerr > last_min_μerr
                println("no better lambda after extension")
                return lastD, last_min_μerr
            end
            if (abs(1 - D.λ / λstar) / (2rλ) >= ρ) | (iter > maxiter)
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
function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μmax::Real, λ::AbstractVector{T}; nfold = 10, nμ = 100, maxiter = 20, figname = "/tmp/cv_curve.png", ρ = 0.05, tol = eps()^(1/3), seed = rand(UInt64)) where T <: AbstractFloat
    # μ = (1:nμ) ./ nμ * μmax
    D, μerr = cvfit(x, y, (1:nμ) ./ nμ .* μmax, λ, figname = figname, nfold = nfold, seed = seed)
    last_min_μerr = minimum(μerr)
    # if D.μ near the left boundary
    if D.μ / μmax < ρ
        iter = 0
        lastD = deepcopy(D)
        while true
            iter += 1
            μmax = μmax / 2
            println("divide the search region of μ: searching [0, $μmax]")
            D, μerr = cvfit(x, y, (1:nμ) ./ nμ .* μmax, λ, figname = figname, nfold = nfold, seed = seed)
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
            D, μerr = cvfit(x, y, (1:nμ) ./ nμ .* μmax, λ, figname = figname, nfold = nfold, seed = seed)
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

function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, paras::AbstractMatrix{T}; nfold = 10, figname = "/tmp/cv_curve.png", seed = rand(UInt64)) where T <: AbstractFloat
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
        build_model!(workspace, x[train_idx])
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
        savefig(cvplot(μerr, σerr, paras[:, 1], nfold = nfold), figname)
    end
    ind = argmin(μerr)
    workspace = WorkSpaceSS()
    D = mono_decomp_ss(workspace, x, y, paras[ind, 2], paras[ind, 1])
    return D, paras[:, 1], μerr, σerr
end

function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, μ::AbstractVector{T}, λ::AbstractVector{T}; nfold = 10, figname = "/tmp/cv_curve.png", seed = rand(UInt64)) where T <: AbstractFloat
    n = length(x)
    folds = div_into_folds(n, K = nfold, seed = seed)
    nλ = length(λ)
    nμ = length(μ)
    err = zeros(nfold, nμ, nλ)
    μs = repeat(μ, outer = nλ)
    λs = repeat(λ, inner = nμ)
    # @showprogress "CV " for k = 1:nfold
    for k = 1:nfold
        test_idx = folds[k]
        train_idx = setdiff(1:n, test_idx)
        workspace = WorkSpaceSS()
        build_model!(workspace, x[train_idx])
        γhats = _optim(y[train_idx], workspace.J, workspace.B, workspace.H, workspace.L, λs, μs)
        ynews = predict(workspace, x[test_idx], γhats[1:workspace.J, :] + γhats[workspace.J+1:2workspace.J, :])
        for j = 1:nλ
            for i = 1:nμ
                m = (j-1) * nμ + i
                err[k, i, j] = norm(ynews[:, m] - y[test_idx]).^2 / length(test_idx)
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
    D = mono_decomp_ss(workspace, x, y, λ[ind[2]], μ[ind[1]])
    return D, μerr, σerr
end

# fix ratio
function cvfit(x::AbstractVector{T}, y::AbstractVector{T}, ks::AbstractVector{T}, λstar::Real; kw...) where T <: AbstractFloat
    @assert maximum(ks) < 1
    return cvfit(x, y, hcat(1 ./ ks .- 1, λstar ./ ks); kw...)
end

"""
    cvfit_gss(x, y, μrange, λ; λ_is_μ)

Cross-validation by Golden Section Searching `μ`` in `μrange` given `λ`.

- If `λ_is_μ`, search `λ` in `μrange` given `λ (μ)`
- Note that `one_se_rule` is not suitable for the golden section search.
"""
function cvfit_gss(x::AbstractVector{T}, y::AbstractVector{T}, μrange::AbstractVector{T}, λ::Real; nfold = 10, figname = "/tmp/cv_curve.png", seed = rand(UInt64), tol = 1e-7, λ_is_μ = false, verbose = false) where T <: AbstractFloat
    τ = (sqrt(5) + 1) / 2
    μs = Float64[]
    errs = Float64[]
    σerrs = Float64[]
    a = μrange[1]
    b = μrange[2]
    c = a
    d = b
    tol = tol * max(1.0, b - a)
    if verbose
        print("search μ in ")
    else
        # println("search μ in ", [a, b])
        @debug "search μ in $([a, b])"
    end
    while true
        verbose && print([a, b], " ")
        c = b - (b - a) / τ
        d = a + (b - a) / τ
        if λ_is_μ
            D, μerr, σerr = cvfit(x, y, [λ], [c, d], nfold = nfold, figname = figname, seed = seed)
        else
            D, μerr, σerr = cvfit(x, y, [c, d], [λ], nfold = nfold, figname = figname, seed = seed)
        end
        append!(μs, [c, d])
        append!(errs, μerr)
        append!(σerrs, σerr)
        if λ_is_μ
            if D.λ == c
                b = d
            else
                a = c
            end
        else
            if D.μ == c # f(c) < f(d)
                b = d
            else
                a = c
            end
        end
        c = b - (b - a) / τ
        d = a + (b - a) / τ
        if abs(b - a) < tol
            verbose && println("") # new line
            if λ_is_μ
                flag = abs.(D.λ .- μrange[2]) / (μrange[2] - μrange[1]) < 1e-2
            else
                flag = abs.(D.μ .- μrange[2]) / (μrange[2] - μrange[1]) < 1e-2
            end
            if flag
                # @warn "the optimal points near the right boundary"
                @debug "the optimal points near the right boundary, try to extend the range..." # since in the log file the warn message is before other print message
                # since 0.75 < 1 - 1e-2, so it is ok
                left_new = μrange[1] + 0.75 * (μrange[2] - μrange[1])
                right_new = μrange[2] * 2
                return cvfit_gss(x, y, [left_new, right_new], λ, figname = figname, nfold = nfold, seed = seed, λ_is_μ = λ_is_μ, tol = tol)
            end
            ind = sortperm(μs)
            verbose && @info "optimal μ = $(μs[end])" 
            if !isnothing(figname)
                savefig(plot(log.(μs[ind]), errs[ind], yerrors = σerrs[ind]), "/tmp/cv.png")
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
function cvfit_gss(x::AbstractVector{T}, y::AbstractVector{T}, μrange::AbstractVector{T}, λs::AbstractVector{T}; nfold = 10, figname = "/tmp/cv_curve.png", seed = rand(UInt64), tol = 1e-7, λ_is_μ = false) where T <: AbstractFloat
    best_err = Inf
    D = nothing
    best_μs = nothing
    best_errs = nothing
    best_σerrs = nothing
    for λ in λs
        @info "given λ = $λ, search μ..."
        Di, μi, erri, σerri = cvfit_gss(x, y, μrange, λ, nfold = nfold, figname = figname, seed = seed, tol = tol)
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
