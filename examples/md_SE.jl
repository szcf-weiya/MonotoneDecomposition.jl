# This section shows how to perform monotone decomposition on a noised random curve generated from Gaussian process.
using MonotoneDecomposition
using Plots
using Random

# Firstly, generate random data from Gaussian process with square kernel,
Random.seed!(9)
x, y, x0, y0 = gen_data(100, 0.5, "SE_1")
# serialize("demo_se_0.1_sigma1.0.sil", [x, y, x0, y0])
# serialize("demo_se_1_sigma0.5_seed9.sil", [x, y, x0, y0])

# ## With Cubic B-splines

# ### `fixJ = true`

# Pefrom Monotone Decomposition with Cubic B-splines, where the number of basis functions is chosen by cross-validation for cubic splines.
Random.seed!(99)
μs = 10.0 .^ (-6:0.05:0)
D, μmin, errs, σerrs = cv_mono_decomp_cs(x, y, ss = μs, 
                            fixJ = true, x0 = x0, 
                            figname = "/tmp/cvspl.png", 
                            nfold = 10, nfold_pre = 10);
yhat, yhatnew = cubic_spline(D.workspace.J)(x, y, x0);
# The CV error curve for `J` is
# ![](/tmp/cvspl_bspl.png)
# And the CV error curve for `μ` is
# ![](/tmp/cvspl.png)
# Plot the fitted curves:
plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1.0, σ = 0.5): ", competitor = "cs")
# Save the figure
# ```julia
# savefig("../res/demo/demo_bspl_vs_md_fixJ.pdf")
# ```
# !!! tip "High-quality Figures"
#     Here for simplicty, we just used the default GR backend of Plots.jl.
#     But for a high-quality figure as in our paper, we will use the PGFPlotsX backend.

# !!! tip "Re-plot CV Curve for Cubic Spline"
#     If `figname` is provided, the CV error curve for the cubic fitting step is stored in `figname[1:end-4] * "_bspl.png"` with the results `figname[1:end-4] * "_bspl_J.sil"`.
#     Thus, if necessary, we can re-plot it and save as a high-quality figure.
#     ```julia
#     res = deserialize("../res/demo/cvspl_bspl_J.sil")
#     μerr, σerr, Js, nfold, ind = res
#     savefig(cvplot(μerr, nothing, 1.0 .* Js, nfold = nfold, ind0 = ind, lbl = L"J"), "../res/demo/demo_bspl_vs_md_fixJ_cv_J.pdf")
#     ```

# Cross-validation plot
cvplot(errs, nothing, 1.0 * D.workspace.J:D.workspace.J, μs, lbl = ["", L"\log_{10}\mu"])
# Save it 
# ```julia
# savefig("../res/demo/demo_bspl_vs_md_fixJ_cv_mu.pdf")
# ```

# !!! tip "Standard Error of CV error"
#     If you want to add the error bar of the CV error, you can specify `σerrs`.
#     ```julia
#     cvplot(errs, σerrs, 1.0 * D.workspace.J:D.workspace.J, μs, lbl = ["", L"\log_{10}\mu"])
#     ```

# ### `fixJ = false`

Js = 4:50
D, μmin, errs, σerrs = 
    cv_mono_decomp_cs(x, y, ss = μs, fixJ = false, x0 = x0, Js = Js);
J, yhat, yhatnew = cv_cubic_spline(x, y, x0, Js = Js)
plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")

# the heatmap of CV-error along the two parameter (J, μ) is as follows,
cvplot(errs, σerrs, Js * 1.0, μs)

# ## With Smoothing Splines

# Perform monotone decomposition with smoothing splines, where the tuning parameter `λ` and `μ` are tuned by 10-fold cross-validation,

D, μmin, μs, errs, σerrs, yhat, yhatnew = 
    cv_mono_decomp_ss(x, y, one_se_rule = true, x0 = x0, tol=1e-3);

# Plot the decomposition,

plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")


# The 10-fold cross-validation curve is 

cvplot(errs, σerrs, μs)