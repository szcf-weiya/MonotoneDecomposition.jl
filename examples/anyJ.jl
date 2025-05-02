# This section aims to show for the cubic spline fitting with any number of basis functions $J$.
# If one applied the monotone decomposition, it can achieve a better performance.
using MonotoneDecomposition
using Plots
using Random

# Firstly, generate random data from Gaussian process
seed = 16
Random.seed!(seed)
x, y, x0, y0 = gen_data(100, 0.5, "SE_1");

# the true curve and noisy observations are shown in the following figure
plot(x0, y0, label = "truth")
scatter!(x, y, label = "sample")

# compare the performance of cubic spline fitting (`CS`) and the corresponding fitting with monotone decomposition (`MDCS`)
Js = 4:50
nJ = length(Js)
errs = zeros(nJ, 2)
Random.seed!(seed)
for i in 1:nJ
    J = Js[i]
    yhat_cs, y0hat_cs = cubic_spline(J)(x, y, x0)
    errs[i, 1] = sum((y0 - y0hat_cs).^2) / length(y0)
    D, _ = cv_mono_decomp_cs(x, y, Js = J:J, ss = 10.0 .^ (-6:0.05:0))
    y0hat_md = predict(D, x0)
    errs[i, 2] = sum((y0 - y0hat_md).^2) / length(y0)
end

# the CV optimized J should be
D, _ = cv_mono_decomp_cs(x, y, Js = Js, ss = 10.0 .^ (-6:0.05:0))
D.workspace.J

# plot the mean squared prediction error (MSPE) along the number of basis function $J$:
plot(Js, errs[:, 1], xlab = "J", ylab = "MSPE", label = "CS", markershape = :circle)
plot!(Js, errs[:, 2], label = "MDCS", markershape = :star5)
Plots.vline!([D.workspace.J], ls = :dash, label = "CVCS")

# Locally, I will save the figure and present in the paper.
# ```julia
# savefig("../res/demo/anylam-seed$seed.pdf")
# ```