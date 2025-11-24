"""
    benchmarking_cs(n, σ, f; figname_cv = nothing, figname_fit = nothing)

Run benchmarking experiments for decomposition with cubic splines on `n` observations sampled from curve `f` with noise `σ`.

## Optional Arguments

- `figname_cv`: if not `nothing`, the cross-validation error will be plotted and saved to the given path.
- `figname_fit`: if not `nothing`, the fitted curves will be plotted and saved to the given path.
- `Js`: the candidates of number of basis functions.
- `fixJ`: whether to use the CV-tuned `J` from the crossponding cubic spline fitting.
- `nfold`: the number of folds in cross-validation procedure
- `one_se_rule`: whether to use the one-standard-error rule to select the parameter after cross-validation procedure
- `μs`: the candidates of tuning parameters for the discrepancy parameter
"""
function benchmarking_cs(n::Int = 100, σ::Union{Real, Nothing} = 0.5, f::Union{Function, String} = x->x^3; fixJ = true,
                                                                               figname_cv = nothing,
                                                                               figname_fit = nothing,
                                                                               Js = 4:20,
                                                                               snr = 1.0,
                                                                               nfold = 10,
                                                                               one_se_rule = false,
                                                                               seed = 1,
                                                                               dataseed = 1,
                                                                               δJ = 100,
                                                                               μs = 10.0 .^ (-6:0.5:0), kw...)
    Random.seed!(dataseed)
    x, y, x0, y0 = gen_data(n, σ, f, snr = snr)
    # J is determined from cubic_spline (deprecated the choice of arbitrary J)
    J, yhat, yhatnew = cv_cubic_spline(x, y, x0, nfold = nfold, one_se_rule = one_se_rule, seed = seed, Js = Js)
    Js = intersect((J-δJ):(J+δJ), Js)
    if fixJ
        Js = J:J
    end
    D, _ = cv_mono_decomp_cs(x, y, x0, Js = Js, ss = μs, figname = figname_cv, nfold = nfold, one_se_rule = one_se_rule, seed = seed)
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

"""
    benchmarking_ss(n::Int, σ::Float64, f::Union{Function, String}; 
                        method = "single_lambda")

Run benchmarking experiments for decomposition with smoothing splines on `n` observations sampled from curve `f` with noise `σ`.

# Arguments

- `method::String = "single_lambda"`: strategy for decomposition with smoothing spline. Possible choices:
    - `single_lambda`
    - `fix_ratio`
    - `grid_search`
    - `iter_search`
"""
function benchmarking_ss(n::Int = 100, σ::Union{Real, Nothing} = 0.5, 
                            f::Union{Function, String} = x->x^3;
                                snr = 1.0, 
                                nfold = 5, one_se_rule = true,
                                method = "single_lambda",
                                figname_cv = nothing,
                                figname_fit = nothing, 
                                seed = 1,
                                dataseed = 1, kw...
                        )
    Random.seed!(dataseed)
    x, y, x0, y0 = gen_data(n, σ, f, snr = snr)
    D, μopt, μs, errs, σerrs, yhat, yhatnew = cv_mono_decomp_ss(x, y; x0 = x0,
                                                                figname = figname_cv,
                                                                nfold = nfold,
                                                                method = method,
                                                                seed = seed,
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
    benchmarking(f::String; n = 100, 
                            σs = 0.2:0.2:1,
                            competitor = "ss_single_lambda")

Run benchmarking experiments for monotone decomposition on curve `f`. The candidates of `f` include:

- simple functions: `x^2`, `x^3`, `exp(x)`, `sigmoid`
- random functions generated from Gaussian Process: `SE_1` `SE_0.1` `Mat12_1` `Mat12_0.1` `Mat32_1` `Mat32_0.1` `RQ_0.1_0.5` `Periodic_0.1_4`

# Arguments

- `n::Integer = 100`: sample size for the simulated curve
- `σs::AbstractVector`: a vector of noise level to be investigated
- `competitor::String`: a string to indicate the strategy used in monotone decomposition. Possible choices:
    - `ss_single_lambda`: decomposition with smoothing splines `ss` with the `single_lambda` strategy
    - `ss_fix_ratio`: decomposition with smoothing splines `ss` with the `fix_ratio` strategy
    - `ss_grid_search`: decomposition with smoothing splines `ss` with the `grid_search` strategy
    - `ss_iter_search`: decomposition with smoothing splines `ss` with the `iter_search` strategy
    - `bspl`: decomposition with cubic splines `cs`
"""
function benchmarking(f::String = "x^3"; n = 100, σs = 0.2:0.2:1,
                            snrs = [0.1, 0.5, 1, 2, 10],
                            J = 10, 
                            μs = 10.0 .^ (-6:0.5:0), ## bspl
                            jplot = false, nrep = 100,
                            competitor = "ss_single_lambda", # bspl, # ss_fix_ratio, ss_grid_search, ss_iter_search
                            nfold = 5, one_se_rule = true,
                            resfolder = "/tmp",
                            ind = 1:4,
                            show_progress = true,
                            use_snr = false,
                            nλ = 20, rλ = 0.5, verbose = false,
                            multi_fix_ratio = false,
                            kw...)
    @info "Benchmarking $f with $nrep repetitions"
    title = "$f (nrep = $nrep)"
    filename = "$f.sil"
    nσ = length(σs)
    arr_snrs = copy(snrs)
    if use_snr # use signal to noise ratio
        σs = [nothing for _ in snrs]
    else
        arr_snrs = vcat(arr_snrs, ones(max(0, nσ - length(snrs))))
    end
    nσ = length(σs) # update length if len(snrs) != len(σs)
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
                                                        snr = arr_snrs[j],
                                                        one_se_rule = one_se_rule,
                                                        method = competitor[4:end], 
                                                        rλ = rλ,
                                                        multi_fix_ratio = multi_fix_ratio,
                                                        dataseed = i,
                                                        verbose = verbose, kw...)
            else
                res[i, :, j] = benchmarking_cs(n, σ, f; figname_cv = figname_cv, 
                                                        figname_fit = figname_fit,
                                                        μs = μs, 
                                                        nfold = nfold,
                                                        snr = arr_snrs[j],
                                                        one_se_rule = one_se_rule,
                                                        dataseed = i,
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
                    snrs = [0.1, 0.5, 1, 2, 10],
                    use_snr = false,
                    nrep = 100, ind = 1:4,
                    curves = ["x^2" "x^3" "exp(x)" "sigmoid" "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"],
                    resfolder = "/tmp" #"../res/res_monodecomp/rocky/"
                    )
    titles = curves
    # titles to display
    dtitles = copy(titles)
    for i in 1:length(dtitles)
        if dtitles[i] in ["x^2", "x^3"]
            dtitles[i] = latexstring(dtitles[i])
        end
        if dtitles[i] == "exp(x)"
            dtitles[i] = L"\exp(x)"
        end
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
        if use_snr
            ind = 1:length(snrs)
            σs = snrs
        else
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
                if !use_snr
                    # var(y0) / σ^2
                    snr = hcat([res[:, 6, i] ./ σs[i]^2 for i in 1:length(σs)]...)
                    snr_μs[i] = mean(snr, dims = 1)[1:1, ind]'
                    snr_σs[i] = std(snr, dims = 1)[1:1, ind]' / sqrt(nrep)
                else
                    snr_μs = nothing
                    snr_σs = nothing            
                end
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
            colnames_of_rownames = ["curve", ifelse(use_snr, "SNR", L"\sigma")], 
            file = joinpath(resfolder0, "$filename.tex"), isbf = isbf)
    end
end

function summary2(;nλ = 20, 
                    format = "tex", 
                    methodnames = ["SmoothSpline", "MonoDecomp"], # CubicSpline
                    # σs = 0.2:0.2:2.0,
                    σs = [0.1, 0.2, 0.4, 0.5, 1.0, 1.5, 2.0],
                    snrs = [0.1, 0.5, 1, 2, 10],
                    use_snr = false,
                    nrep = 100, ind = 1:4,
                    curves = ["x^2" "x^3" "exp(x)" "sigmoid" "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"],
                    σs_from_snr = nothing,
                    resfolder = "/tmp" #"../res/res_monodecomp/rocky/"
                    )
    titles = curves
    # titles to display
    dtitles = copy(titles)
    for i in 1:length(dtitles)
        if dtitles[i] in ["x^2", "x^3"]
            dtitles[i] = latexstring(dtitles[i])
        end
        if dtitles[i] == "exp(x)"
            dtitles[i] = L"\exp(x)"
        end
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
        # pvals_wilcox = Array{AbstractVector}(undef, ntitle)
        props = Array{AbstractVector}(undef, ntitle)
        snr_μs = Array{AbstractMatrix{Float64}}(undef, ntitle)
        snr_σs = Array{AbstractMatrix{Float64}}(undef, ntitle)
        isbf = Array{AbstractMatrix{Bool}}(undef, ntitle)
        # selected noise levels
        if use_snr
            ind = 1:length(snrs)
            # ind = [1]
            σs = snrs
        else
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
                if !use_snr
                    # var(y0) / σ^2
                    snr = hcat([res[:, 6, i] ./ σs[i]^2 for i in 1:length(σs)]...)
                    snr_μs[i] = mean(snr, dims = 1)[1:1, ind]'
                    snr_σs[i] = std(snr, dims = 1)[1:1, ind]' / sqrt(nrep)
                else
                    snr_μs = nothing
                    snr_σs = nothing            
                end
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
            # pval = 1 .- cdf(Normal(0, 1), abs.(dμ) ./ dσ)
            # pvals[i] = star_pval(pval[ind])
            #pval = [rcopy(R"t.test($(res[:, 2, j]), $(res[:, 4, j]), paired = TRUE)$p.value") for j in ind]
            #pvals[i] = star_pval(pval)
            pval = [rcopy(R"wilcox.test($(res[:, 2, j]), $(res[:, 4, j]), paired = TRUE)$p.value") for j in ind]
            pvals[i] = star_pval(pval)
            μerrs[i] = mean(res, dims = 1)[1, [4, 2], ind]'
            σerrs[i] = std(res, dims = 1)[1, [4, 2], ind]' / sqrt(nrep)
            isbf[i] = zeros(length(ind), 2)
            isbf[i][:, 1] .= (μerrs[i][:, 1] .< μerrs[i][:, 2]) .& (pval .< 0.1)
            isbf[i][:, 2] .= (μerrs[i][:, 1] .> μerrs[i][:, 2]) .& (pval .< 0.1)
        end
        # `_` is not allowed in non-math env
        dtitles = replace.(dtitles, "_" => "-")
        # print2tex, similar to https://github.com/szcf-weiya/Clouds/search?q=print2tex
        print2tex(μerrs, σerrs, dtitles[:], ["MSPE"], string.(σs[ind]), methodnames, 
            other_cols = ifelse(use_snr, σs_from_snr, nothing),
            other_cols_σ = nothing,
            other_col_names = [L"\sigma"],
            #right_cols = [pvals, pvals_wilcox, props],
            right_cols = [pvals, props],
            #right_col_names = ["p-value", "p-value (wilcox)", "prop."], 
            right_col_names = ["p-value", "prop."], 
            right_align = "l", # it might be slightly worse for other columns, but center is not good for the pvalue column TODO: improve the style
            colnames_of_rownames = ["curve", ifelse(use_snr, "SNR", L"\sigma")], 
            file = joinpath(resfolder0, "$filename-ind$(length(ind)).tex"), isbf = isbf)
    end
end
