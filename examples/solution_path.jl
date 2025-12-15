using MonotoneDecomposition
using Plots
using Random

Random.seed!(123)
x, y, x0, y0 = gen_data(100, 0.1, x -> exp(x));
μs = 10.0 .^ (-2:0.2:3)
nμ = length(μs)
res = Array{Any, 1}(undef, nμ)

# ## J = 8

for i in 1:nμ
    res[i] = mono_decomp_cs(x, y, J = 8, s = μs[i], s_is_μ=true)
end
γups = hcat([res[i].γup for i in 1:nμ]...)
γdowns = hcat([res[i].γdown for i in 1:nμ]...)
# save the results
# ```julia
# serialize("../res/solution_path/expx-sigma0.1-J8.sil", [γups, γdowns, res, μs, x, y])
# ```

# calculate the coefficients via the least square solution
γls = inv(res[1].workspace.B' * res[1].workspace.B) * res[1].workspace.B' * y
c = mean(res[1].workspace.B * γls) / 2
γups_calculated = [γls ./ (μ + 1) .+ (μ - 1) / (μ + 1) * c for μ in μs] 
γups_ls = hcat(γups_calculated...)

# plot the coefficients along the tuning parameter μ
default_colors = palette(:auto)
plot(log10.(μs), γups[1, :], label = L"\hat\gamma_1^u", xlab = L"\log_{10}\; \mu", color = default_colors[1])
plot!(log10.(μs), γups_ls[1, :], label = L"\gamma_1^{u}", ls = :dash, lw = 2, color = default_colors[1])
plot!(log10.(μs), γups[2, :], label = L"\hat\gamma_2^u", color = default_colors[2])
plot!(log10.(μs), γups_ls[2, :], label = L"\gamma_2^{u}", ls = :dash, lw = 2, color = default_colors[2])
plot!(log10.(μs), γups[3, :], label = L"\hat\gamma_3^u", color = default_colors[3])
plot!(log10.(μs), γups_ls[3, :], label = L"\gamma_3^{u}", ls = :dash, lw = 2, color = default_colors[3])
plot!(log10.(μs), γups[4, :], label = L"\hat\gamma_4^u", color = default_colors[4])
plot!(log10.(μs), γups_ls[4, :], label = L"\gamma_4^{u}", ls = :dash, lw = 2, color = default_colors[4])
plot!(log10.(μs), γups[5, :], label = L"\hat\gamma_5^u", color = default_colors[5])
plot!(log10.(μs), γups_ls[5, :], label = L"\gamma_5^{u}", ls = :dash, lw = 2, color = default_colors[5])
plot!(log10.(μs), γups[6, :], label = L"\hat\gamma_6^u", color = default_colors[6])
plot!(log10.(μs), γups_ls[6, :], label = L"\gamma_6^{u}", ls = :dash, lw = 2, color = default_colors[6])
plot!(log10.(μs), γups[7, :], label = L"\hat\gamma_7^u", color = default_colors[7])
plot!(log10.(μs), γups_ls[7, :], label = L"\gamma_7^{u}", ls = :dash, lw = 2, color = default_colors[7])
plot!(log10.(μs), γups[8, :], label = L"\hat\gamma_8^u", color = default_colors[8])
plot!(log10.(μs), γups_ls[8, :], label = L"\gamma_8^{u}", ls = :dash, lw = 2, color = default_colors[8])

# save the results
# ```julia
# savefig("../res/solution_path/expx-sigma0.1-J8-up.pdf")
# ```

# plot the difference between the estimated coefficients and the theoretical constant vector
diff_downs = sum((γdowns .- c).^2, dims = 1)[:]
plot(log10.(μs), log10.(diff_downs), xlab = L"\log_{10}\; \mu", label = "", ylab = L"\log_{10}\Vert \hat\gamma^d - \gamma^d\Vert_2^2")

# save the figure
# ```julia
# savefig("../res/solution_path/expx-sigma0.1-J8-down.pdf")
# ```



# ## J = 9

for i in 1:nμ
    res[i] = mono_decomp_cs(x, y, J = 9, s = μs[i], s_is_μ=true)
end
γups = hcat([res[i].γup for i in 1:nμ]...)
γdowns = hcat([res[i].γdown for i in 1:nμ]...)
# save the results
# ```julia
# serialize("../res/solution_path/expx-sigma0.1-J9.sil", [γups, γdowns, res, μs, x, y])
# ```

# ### without tied information 

# calculate the coefficients via the least square solution
γls = inv(res[1].workspace.B' * res[1].workspace.B) * res[1].workspace.B' * y
c = mean(res[1].workspace.B * γls) / 2
γups_calculated = [γls ./ (μ + 1) .+ (μ - 1) / (μ + 1) * c for μ in μs] 
γups_ls = hcat(γups_calculated...)

# plot the coefficients along the tuning parameter μ
plot(log10.(μs), γups[1, :], label = L"\hat\gamma_1^u", xlab = L"\log_{10}\; \mu", color = default_colors[1], lw = 3, alpha = 0.5)
plot!(log10.(μs), γups_ls[1, :], label = L"\gamma_1^{u}", ls = :dash, lw = 2, color = default_colors[1])
plot!(log10.(μs), γups[2, :], label = L"\hat\gamma_2^u", color = default_colors[2])
plot!(log10.(μs), γups_ls[2, :], label = L"\gamma_2^{u}", ls = :dash, lw = 2, color = default_colors[2])
plot!(log10.(μs), γups[3, :], label = L"\hat\gamma_3^u", color = default_colors[3])
plot!(log10.(μs), γups_ls[3, :], label = L"\gamma_3^{u}", ls = :dash, lw = 2, color = default_colors[3])
plot!(log10.(μs), γups[4, :], label = L"\hat\gamma_4^u", color = default_colors[4])
plot!(log10.(μs), γups_ls[4, :], label = L"\gamma_4^{u}", ls = :dash, lw = 2, color = default_colors[4])
plot!(log10.(μs), γups[5, :], label = L"\hat\gamma_5^u", color = default_colors[5])
plot!(log10.(μs), γups_ls[5, :], label = L"\gamma_5^{u}", ls = :dash, lw = 2, color = default_colors[5])
plot!(log10.(μs), γups[6, :], label = L"\hat\gamma_6^u", color = default_colors[6])
plot!(log10.(μs), γups_ls[6, :], label = L"\gamma_6^{u}", ls = :dash, lw = 2, color = default_colors[6])
plot!(log10.(μs), γups[7, :], label = L"\hat\gamma_7^u", color = default_colors[7])
plot!(log10.(μs), γups_ls[7, :], label = L"\gamma_7^{u}", ls = :dash, lw = 2, color = default_colors[7])
plot!(log10.(μs), γups[8, :], label = L"\hat\gamma_8^u", color = default_colors[8])
plot!(log10.(μs), γups_ls[8, :], label = L"\gamma_8^{u}", ls = :dash, lw = 2, color = default_colors[8])
plot!(log10.(μs), γups[9, :], label = L"\hat\gamma_9^u", color = default_colors[9])
plot!(log10.(μs), γups_ls[9, :], label = L"\gamma_9^{u}", ls = :dash, lw = 2, color = default_colors[9])

# save the figure
# ```julia
# savefig("../res/solution_path/expx-sigma0.1-J9-up-ls.pdf")
# ```

diff_downs = sum((γdowns .- c).^2, dims = 1)[:]
plot(log10.(μs), log10.(diff_downs), xlab = L"\log_{10}\; \mu", label = "", ylab = L"\log_{10}\Vert \hat\gamma^d - \gamma^d\Vert_2^2")

# save the figure
# ```julia
# savefig("../res/solution_path/expx-sigma0.1-J9-down-ls.pdf")
# ```

# ### with tied information 

G = hcat(zeros(8), 1.0I(8))
G[1, 1] = 1
γGls = G' * inv(G * res[1].workspace.B' * res[1].workspace.B * G') * G * res[1].workspace.B' * y
c = mean(res[1].workspace.B * γGls) / 2
γups_calculated = [γGls ./ (μ + 1) .+ (μ - 1) / (μ + 1) * c for μ in μs] 
γups_ls = hcat(γups_calculated...)

# plot the results
plot(log10.(μs), γups[1, :], label = L"\hat\gamma_1^u", xlab = L"\log_{10}\; \mu", color = default_colors[1], lw = 3, alpha = 0.5)
plot!(log10.(μs), γups_ls[1, :], label = L"\gamma_1^{u}", ls = :dash, lw = 2, color = default_colors[1])
plot!(log10.(μs), γups[2, :], label = L"\hat\gamma_2^u", color = default_colors[2])
plot!(log10.(μs), γups_ls[2, :], label = L"\gamma_2^{u}", ls = :dash, lw = 2, color = default_colors[2])
plot!(log10.(μs), γups[3, :], label = L"\hat\gamma_3^u", color = default_colors[3])
plot!(log10.(μs), γups_ls[3, :], label = L"\gamma_3^{u}", ls = :dash, lw = 2, color = default_colors[3])
plot!(log10.(μs), γups[4, :], label = L"\hat\gamma_4^u", color = default_colors[4])
plot!(log10.(μs), γups_ls[4, :], label = L"\gamma_4^{u}", ls = :dash, lw = 2, color = default_colors[4])
plot!(log10.(μs), γups[5, :], label = L"\hat\gamma_5^u", color = default_colors[5])
plot!(log10.(μs), γups_ls[5, :], label = L"\gamma_5^{u}", ls = :dash, lw = 2, color = default_colors[5])
plot!(log10.(μs), γups[6, :], label = L"\hat\gamma_6^u", color = default_colors[6])
plot!(log10.(μs), γups_ls[6, :], label = L"\gamma_6^{u}", ls = :dash, lw = 2, color = default_colors[6])
plot!(log10.(μs), γups[7, :], label = L"\hat\gamma_7^u", color = default_colors[7])
plot!(log10.(μs), γups_ls[7, :], label = L"\gamma_7^{u}", ls = :dash, lw = 2, color = default_colors[7])
plot!(log10.(μs), γups[8, :], label = L"\hat\gamma_8^u", color = default_colors[8])
plot!(log10.(μs), γups_ls[8, :], label = L"\gamma_8^{u}", ls = :dash, lw = 2, color = default_colors[8])
plot!(log10.(μs), γups[9, :], label = L"\hat\gamma_9^u", color = default_colors[9])
plot!(log10.(μs), γups_ls[9, :], label = L"\gamma_9^{u}", ls = :dash, lw = 2, color = default_colors[9])

# save the figure
# ```julia
# savefig("../res/solution_path/expx-sigma0.1-J9-up-gls.pdf")
# ```