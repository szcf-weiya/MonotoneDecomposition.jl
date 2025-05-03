# This section aims to show for the smoothing spline fitting with arbitrary λ.
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

# the CV-optimized λ for the smoothing spline
λs = 10.0 .^ (-8:1) # λs = 10.0 .^ (-8:0.1:1)
yhat, _, λcv = MonotoneDecomposition.cv_smooth_spline(x, y, λs)
nλ = length(λs)
errs = zeros(nλ, 2)
μs = 10 .^ (-8:1) # μs = 10 .^ (-8:0.1:1)
for (i, λ) in enumerate(λs)
    @info "i = $i"
    ## D, _ = MonotoneDecomposition.cvfit_gss(x, y, [1e-6, 1.0], λ, tol = 1e-2)
    D, _ = MonotoneDecomposition.cvfit(x, y, μs, [λ])
    y0hat_mdss = predict(D, x0)
    y0hat_ss = MonotoneDecomposition.smooth_spline(x, y, x0, λ)          
    errs[i, :] .= [sum((y0hat_ss - y0).^2), sum((y0hat_mdss - y0).^2)] / length(y0)
end

# save results 
# ```julia
# serialize("../res/demo/anylam-seed$seed.sil", [λs, errs, λcv])
# ```
idx = 1:nλ # idx = 2:2:nλ # if too many points (looks too busy)
plot(log.(λs)[idx], errs[idx, 1], xlab = "λ", ylab = "MSPE", label = "SS", markershape = :circle, legend = :topright)
plot!(log.(λs)[idx], errs[idx, 2], label = "MDSS", markershape = :star5)
Plots.vline!([log.(λcv)], ls = :dash, label = "CVSS")

# ```julia
# savefig("../res/demo/anylam-seed$seed.pdf")
# ```