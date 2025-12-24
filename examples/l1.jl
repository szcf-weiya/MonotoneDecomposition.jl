using MonotoneDecomposition
using Plots
using Random
using LaTeXStrings

Random.seed!(123)
x, y, x0, y0 = gen_data(100, 0.1, x -> exp(x));
μs = 10.0 .^ (-2:0.2:3)
nμ = length(μs)
res = Array{Any, 1}(undef, nμ);

# ## J = 8
# using l1 penalty by specifying `l1 = true`
for i in 1:nμ
    res[i] = mono_decomp_cs(x, y, J = 8, s = μs[i], s_is_μ=true, l1 = true)
end
γups = hcat([res[i].γup for i in 1:nμ]...)
γdowns = hcat([res[i].γdown for i in 1:nμ]...)

# plot the increasing part
default_colors = palette(:auto)
plot(log10.(μs), γups[1, :], label = L"\hat\gamma_1^u", xlab = L"\log_{10}\; \mu", color = default_colors[1])
plot!(log10.(μs), γups[2, :], label = L"\hat\gamma_2^u", color = default_colors[2])
plot!(log10.(μs), γups[3, :], label = L"\hat\gamma_3^u", color = default_colors[3])
plot!(log10.(μs), γups[4, :], label = L"\hat\gamma_4^u", color = default_colors[4])
plot!(log10.(μs), γups[5, :], label = L"\hat\gamma_5^u", color = default_colors[5])
plot!(log10.(μs), γups[6, :], label = L"\hat\gamma_6^u", color = default_colors[6])
plot!(log10.(μs), γups[7, :], label = L"\hat\gamma_7^u", color = default_colors[7])
plot!(log10.(μs), γups[8, :], label = L"\hat\gamma_8^u", color = default_colors[8])
# save the figure
# ```julia
# savefig("../res/l1/expx-sigma0.1-J8-up.pdf")
# ```

# plot the decreasing part
plot(log10.(μs), γdowns[1, :], label = L"\hat\gamma_1^d", xlab = L"\log_{10}\; \mu", color = default_colors[1], lw = 4)
plot!(log10.(μs), γdowns[2, :], label = L"\hat\gamma_2^d", color = default_colors[2], lw = 3.5)
plot!(log10.(μs), γdowns[3, :], label = L"\hat\gamma_3^d", color = default_colors[3], lw = 3)
plot!(log10.(μs), γdowns[4, :], label = L"\hat\gamma_4^d", color = default_colors[4], lw = 2.5)
plot!(log10.(μs), γdowns[5, :], label = L"\hat\gamma_5^d", color = default_colors[5], lw = 2)
plot!(log10.(μs), γdowns[6, :], label = L"\hat\gamma_6^d", color = default_colors[6], lw = 1.5)
plot!(log10.(μs), γdowns[7, :], label = L"\hat\gamma_7^d", color = default_colors[7], lw = 1)
plot!(log10.(μs), γdowns[8, :], label = L"\hat\gamma_8^d", color = default_colors[8], lw = 0.5)
# save the figure
# ```julia
# savefig("../res/l1/expx-sigma0.1-J8-down.pdf")
# ```

# ## running time

# check the running time
i = 1
t1 = [
    (@elapsed [mono_decomp_cs(x, y, J = 10, s = μs[i], s_is_μ=true, l1 = false, use_GI = false) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 20, s = μs[i], s_is_μ=true, l1 = false, use_GI = false) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 30, s = μs[i], s_is_μ=true, l1 = false, use_GI = false) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 40, s = μs[i], s_is_μ=true, l1 = false, use_GI = false) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 50, s = μs[i], s_is_μ=true, l1 = false, use_GI = false) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 64, s = μs[i], s_is_μ=true, l1 = false, use_GI = false) for _ in 1:100])
];

t2 = [
    (@elapsed [mono_decomp_cs(x, y, J = 10, s = μs[i], s_is_μ=true, l1 = true) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 20, s = μs[i], s_is_μ=true, l1 = true) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 30, s = μs[i], s_is_μ=true, l1 = true) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 40, s = μs[i], s_is_μ=true, l1 = true) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 50, s = μs[i], s_is_μ=true, l1 = true) for _ in 1:100]),
    (@elapsed [mono_decomp_cs(x, y, J = 64, s = μs[i], s_is_μ=true, l1 = true) for _ in 1:100])
];
# the running time for MDSS with the default l2 penalty
tss = @elapsed [mono_decomp_ss(x, y, use_GI = false) for _ in 1:100]

# the running time for MDSS with the l1 penalty
tss2 = @elapsed [mono_decomp_ss(x, y, l1 = true) for _ in 1:100]

# save the results
# ```julia
# serialize("../res/l1/running_time.sil", [t1, t2, tss, tss2])
# ```

# plot the running time results
Js = [10, 20, 30, 40, 50, 64]
plot(Js, t2, label = "MDCS (L1)", ls = :dash, color = default_colors[1], markershape = :star5, legend = :topleft, xticks = (Js, string.(Js)), xlab = "J", ylab = "seconds")
plot!(Js, t1, label = "MDCS (L2)", color = default_colors[2], markershape = :dtriangle)
scatter!([64], [tss2], label = "MDSS (L1)", ls = :dash, color = default_colors[1], markershape = :star5)
scatter!([64], [tss], label = "MDSS (L2)", color = default_colors[2], markershape = :dtriangle)
# save the figure
# ```julia
# savefig("../res/l1/running_time.pdf")
# ```


# ## Comparison between L1 and L2 penalties

## result folder with the l1 and l2 penalties
l1_folder = "nrep100-nfold10-nlam20-1se_false-bspl-2025-12-15T22_28_36+08_00-n100-snr_true"
l2_folder = "nrep100-nfold10-nlam20-1se_false-bspl-2025-10-04T12_15_12+08_00-n100-snr_true"

arr_curves = ["x^2" "x^3" "exp(x)" "sigmoid" "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"]
snrs = [0.1, 0.5, 1, 2, 10]

# put the MSPE along the penalty into the same figure for the l1 and l2 penalties
# for the first two curves, we use the latexstring to write the title, so the following `for` loop starts from `i = 3`.
# ```julia
# for i = 3:12
#     res1 = deserialize(joinpath("../res/zju/", l1_folder, "$(arr_curves[i]).sil"))
#     res2 = deserialize(joinpath("../res/zju/", l2_folder, "$(arr_curves[i]).sil"))
# 
#     plot(mean(res1[:, 2, :], dims=1)[1:5], label = L"L_1", legend = :topright, xticks = (1:5, string.(snrs)), xlab = "SNR", markershape = :star5, 
#         title = L"x^3")
#         # title = arr_curves[i])
#     plot!(mean(res2[:, 2, :], dims=1)[1:5], label = L"L_2", markershape = :dtriangle)
#     savefig("../res/l1/l1_vs_l2_mdcs-bspl-$(arr_curves[i]).pdf")
# end
