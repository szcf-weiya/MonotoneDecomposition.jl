# This section shows simulated curves from Gaussian process with various kernels.
using Random
using Plots
using LaTeXStrings
using MonotoneDecomposition: sigmoid, gp
function plot_functions()
    n = 100
    Random.seed!(1)
    x = sort(rand(n)) * 2 .- 1
    figsize = (400, 400)
    fig1 = plot(x, x .^ 3, label = L"x^3", ls = :solid, size = figsize, legend = :top)
    plot!(fig1, x, x .^ 2, label = L"x^2", ls = :dash)
    plot!(fig1, x, exp.(x) .- 1, label = L"\exp(x)-1", ls = :dot)
    plot!(fig1, x, sigmoid.(x), label = L"1/(1+\exp(-5x))", ls = :dashdot)
    plot!(fig1, x, sin.(2π * x), label = L"\sin(2\pi x)", ls = :dashdotdot)
    ## plot(x, gp(x, kernel = "SE_0.1"))
    fig2 = plot(x, gp(x, kernel = "SE_1"), label = "SE_1", ls = :solid, size = figsize, legend = :topright)
    plot!(fig2, x, gp(x, kernel = "SE_0.1"), label = "SE_0.1", ls = :dash)
    plot!(fig2, x, gp(x, kernel = "SE_0.5"), label = "SE_0.5", ls = :dot)
    fig3 = plot(x, gp(x, kernel = "Mat12_1"), label = "Mat12", ls = :solid, size = figsize, legend = :topleft)
    plot!(fig3, x, gp(x, kernel = "Mat32_1"), label = "Mat32", ls = :dash)
    plot!(fig3, x, gp(x, kernel = "Mat52_1"), label = "Mat52", ls = :dot)
    plot!(fig3, x, gp(x, kernel = "RQ_1_0.5"), label = "RQ", ls = :dashdot)
    plot!(fig3, x, gp(x, kernel = "Periodic_1_1"), ls = :dashdotdot, label = "Periodic")
    return plot(fig1, fig2, fig3, layout = (1, 3), size = (3*figsize[1], figsize[2]) )
##    savefig("~/PGitHub/overleaf/MonoDecomp/figs/funcs.pdf") # saved for paper
end

# plot functions generated from Gaussian Process
plot_functions()