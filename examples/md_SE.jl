# This section shows how to perform monotone decomposition on a noised random curve generated from Gaussian process.
using MonotoneDecomposition
using Plots
using Random

# Firstly, generate random data from Gaussian process with square kernel,
seed = 16
Random.seed!(seed)
x, y, x0, y0 = gen_data(100, 0.5, "SE_1");
# Save the simulated data for reproducing if `seed` failed.
# ```julia
# serialize("../res/demo/demo-seed$seed-se1_sigma0.5.sil", [x, y, x0, y0])
# ```
# ## With Cubic B-splines

# ### `fixJ = true`

# Pefrom Monotone Decomposition with Cubic B-splines, where the number of basis functions is chosen by cross-validation for cubic splines.
Random.seed!(seed)
μs = 10.0 .^ (-6:0.05:0)
D, μmin, errs, σerrs = cv_mono_decomp_cs(x, y, ss = μs, 
                            fixJ = true, x0 = x0, 
                            figname = "cvspl.png", 
                            nfold = 10, nfold_pre = 10);
yhat, yhatnew = cubic_spline(D.workspace.J)(x, y, x0);
# The CV error curve for `J` is
#
# ![](cvspl_bspl.png)
#
# And the CV error curve for `μ` is
#
# ![](cvspl.png)
#
# Plot the fitted curves:
plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1.0, σ = 0.5): ", competitor = "cs")
# Save the figure:
# ```julia
# savefig("../res/demo/demo-seed$seed-cs_vs_md-1J_and_mu-fit.pdf")
# ```
# !!! tip "High-quality Figures"
#     Here for simplicty, we just used the default GR backend of Plots.jl.
#     But for a high-quality figure as in our paper, we will use the PGFPlotsX backend.

# !!! tip "Re-plot CV Curve for Cubic Spline"
#     If `figname` is provided, the CV error curve for the cubic fitting step is stored in `figname[1:end-4] * "_bspl.png"` with the results `figname[1:end-4] * "_bspl_J.sil"`.
#     Thus, if necessary, we can re-plot it and save as a high-quality figure.
#     Also save the CV results:
#     ```julia
#     mv("cvspl_bspl_J.sil", "../res/demo/cvspl_bspl_J.sil", force = true)
#     μerr, σerr, Js, nfold, ind = deserialize("../res/demo/cvspl_bspl_J.sil")
#     cvplot(μerr, nothing, 1.0 .* Js, nfold = nfold, ind0 = ind, lbl = L"J", title = "10-fold CV Error (Cubic Spline)")
#     savefig("../res/demo/demo-seed$seed-cs_vs_md-cs-cv.pdf")
#     ```

# Cross-validation plot
cvplot(errs, nothing, 1.0 * D.workspace.J:D.workspace.J, μs, lbl = ["", "\\log_{10}\\mu"], title = "10-fold CV Error (J = $(D.workspace.J))")
# Also backup the results
# ```shell
# cp cvspl_Jmu.sil ../res/demo/
# ```
# And save the figure
# ```julia
# savefig("../res/demo/demo-seed$seed-cs_vs_md-1J_and_mu-cv.pdf")
# ```

# !!! tip "Standard Error of CV error"
#     If you want to add the error bar of the CV error, you can specify `σerrs`.
#     ```julia
#     cvplot(errs, σerrs, 1.0 * D.workspace.J:D.workspace.J, μs, lbl = ["", L"\log_{10}\mu"])
#     ```

# ### `fixJ = false`

Js = 4:50
Random.seed!(seed)
D, μmin, errs, σerrs = cv_mono_decomp_cs(x, y, ss = μs, fixJ = false, 
                                        x0 = x0, Js = Js, 
                                        figname = "cvbspl_varyJ.png", 
                                        one_se_rule = false);
plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ", competitor = "cs")
# ```julia
# savefig("../res/demo/demo-seed$seed-cs_vs_md-J_and_mu-fit.pdf")
# ```
# the CV error curve for the cubic spline
#
# ![][cvbspl_varyJ_bspl.png]
#
# and the CV error curve for the decomposition with cubic splines
#
# ![](cvbspl_varyJ.png)

# Or we can replot the heatmap of CV-error along the two parameter (J, μ) is as follows,
cvplot(errs, nothing, Js * 1.0, μs, lbl = ["J", "\\mu"], title = "10-fold CV Error")
# save figure
# ```julia
# savefig("../res/demo/demo-seed$seed-cs_vs_md-J_and_mu-cv.pdf")
# ```

# Backup the results
# ```julia
# cp cvbspl_varyJ_bspl_J.sil ../res/demo
# cp cvbspl_varyJ_Jmu.sil ../res/demo
# ```

# ## With Smoothing Splines

# Perform monotone decomposition with smoothing splines, where the tuning parameter `λ` and `μ` are tuned by 10-fold cross-validation,

# ### Fix `λ`

Random.seed!(seed)
D, μmin, μs, errs, σerrs, yhat, yhatnew = cv_mono_decomp_ss(x, y, 
                                            one_se_rule = false, x0 = x0, 
                                            figname = "cvss_1lam.png");

# Firstly, the cross-validation curve for the smoothing spline is
#
# ![](cvss_1lam_ss.png)
#
# !!! tip "Reproducing high-quality figure in manuscript"
#     To produce a high-quality figure as in the manuscript
#     ```julia
#     μerr, σerr, λs, nfold, ind = deserialize("cvss_1lam_ss.sil")
#     # cvplot(μerr, σerr, λs, nfold = nfold, ind0 = ind, lbl = "\\lambda")
#     cvplot(μerr, nothing, λs, nfold = nfold, ind0 = ind, lbl = "\\lambda", title = "10-fold CV (Smoothing Spline)")
#     savefig("../res/demo/demo-seed$seed-ss_vs_md-ss-cv.pdf")
#     ```
#
# Then, given the optimal `λ`, tune `μ`, the CV error curve is as follows:
#
# ![](cvss_1lam.png)
# 
# !!! tip "Reproducing high-quality figure in manuscript"
#     To produce a high-quality figure as in the manuscript
#     ```julia
#     cvplot(errs, nothing, μs, lbl = "\\mu", title = "10-fold CV Error (λ = $(D.λ))")
#     savefig("../res/demo/demo-seed$seed-ss_vs_md-1lam_and_mu-cv.pdf")
#     ```

# Plot the decomposition,
plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")
# save the fitness curve
# ```julia
# savefig("../res/demo/demo-seed$seed-ss_vs_md-1lam_and_mu-fit.pdf")
# ```

# ### Vary `λ`
Random.seed!(seed)
D, μmin, μs, errs, σerrs, yhat, yhatnew = cv_mono_decomp_ss(x, y, one_se_rule = false, x0 = x0, 
                                            figname = "cvss_varylam.png",
                                            method = "double_grid",
                                            nλ = 50, ngrid_μ = 50,
                                            μrange = [1e-7, 1.0],
                                            rλs = 10.0 .^ (-1:0.05:1.2));

plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")
# save the figure 
# ```julia
# savefig("../res/demo/demo-seed$seed-ss_vs_md-lam_and_mu-fit.pdf")
# ```

# The CV error heatmap is 
#
# ![](cvss_varylam.png)
#
# Alternatively, we can replot it:
cvplot("/tmp/cvss_varylam.sil", "ss")
# save the results
# ```shell
# mv cvss_varylam.sil ../res/demo/
# ```
# and figure
# ```julia
# savefig("../res/demo/demo-seed$seed-ss_vs_md-lam_and_mu-cv.pdf")
# ```