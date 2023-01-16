using MonotoneDecomposition
using Plots

# Firstly, generate random data from Gaussian process with square kernel,
x, y, x0, y0 = gen_data(100, 0.5, "SE_1")

# Perform monotone decomposition with smoothing splines, where the tuning parameter `λ` and `μ` are tuned by 10-fold cross-validation,

D, μmin, μs, errs, σerrs, yhat, yhatnew = 
    cv_mono_decomp_ss(x, y, one_se_rule = true, x0 = x0, tol=1e-3);

# Plot the decomposition,

plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")


# The 10-fold cross-validation curve is 

cvplot(errs, σerrs, μs)

# Pefrom Monotone Decomposition with Cubic B-splines

D, μmin, μs, errs, σerrs, yhat, yhatnew = 
    cv_mono_decomp_cs(x, y, ss = 10.0 .^ (-6:0.5:-1), fixJ = true, x0 = x0);

# Plot it
plot([x, y], [x0, y0], D, yhatnew, prefix_title = "SE (ℓ = 1, σ = 0.5): ")

# Cross-validation plot
cvplot(errs, σerrs, 1.0 * D.workspace.J:D.workspace.J, μs)