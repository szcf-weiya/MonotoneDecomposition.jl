# This section aims to investigate the effect of sample size.
using MonotoneDecomposition
using Plots
using LinearAlgebra

# ## With Cubic B-splines (`fitJ`)
seed = 16
ns = [20, 50, 100, 200, 500]

Err1 = Float64[]
Err2 = Float64[]
for n in ns
    Random.seed!(16)
    x, y, x0, y0 = gen_data(n, 0.5, "SE_1");
    D, μmin, errs, σerrs = cv_mono_decomp_cs(x, y, ss = μs, 
                            fixJ = true, x0 = x0, 
                            figname = "cvspl_n$n.png", 
                            nfold = 10, nfold_pre = 10);
    yhat, yhatnew = cubic_spline(D.workspace.J)(x, y, x0);
    yup, ydown = predict(D.workspace, x0, D.γup, D.γdown)
    e1 = norm(yup + ydown - y0)^2 / length(y0)
    e2 = norm(yhatnew - y0)^2 / length(y0)
    append!(Err1, e1)
    append!(Err2, e2)
end

# Plot the results
plot(ns, Err1, label = "MDCS", markershape = :circle, xlab = "sample size n", ylab = "Err")
plot!(ns, Err2, label = "CS", markershape = :rect)
