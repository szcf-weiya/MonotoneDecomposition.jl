# This section shows how to perform monotone decomposition on a noised random curve generated from Gaussian process.
using MonotoneDecomposition
using Plots
using Random

# Firstly, generate random data from Gaussian process with square kernel,
Random.seed!(7)
x, y, x0, y0 = gen_data(100, 0.5, "SE_1")

# ## With Cubic B-splines

# ### `fixJ = true`

# Pefrom Monotone Decomposition with Cubic B-splines, where the number of basis functions is chosen by cross-validation for cubic splines.
ss = 10.0 .^ (-6:0.1:0)
D, μmin, μs, errs, σerrs, yhat, yhatnew = 
    cv_mono_decomp_cs(x, y, ss = ss, fixJ = true, x0 = x0);

# Plot it
plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")

# Cross-validation plot
cvplot(errs, σerrs, 1.0 * D.workspace.J:D.workspace.J, μs)

# ### `fixJ = false`

Js = 4:50
D, μmin, μs, errs, σerrs = 
    cv_mono_decomp_cs(x, y, ss = ss, fixJ = false, x0 = x0, Js = Js);
J, yhat, yhatnew = cv_cubic_spline(x, y, x0, Js = Js)
plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")

# the heatmap of CV-error along the two parameter (J, μ) is as follows,
cvplot(errs, σerrs, Js * 1.0, ss)

# ## With Smoothing Splines

# Perform monotone decomposition with smoothing splines, where the tuning parameter `λ` and `μ` are tuned by 10-fold cross-validation,

D, μmin, μs, errs, σerrs, yhat, yhatnew = 
    cv_mono_decomp_ss(x, y, one_se_rule = true, x0 = x0, tol=1e-3);

# Plot the decomposition,

plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")


# The 10-fold cross-validation curve is 

cvplot(errs, σerrs, μs)