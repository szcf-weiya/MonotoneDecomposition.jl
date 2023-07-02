using Test
using Plots
using MonotoneDecomposition
using GaussianProcesses
using Statistics

@testset "GP" begin
    @testset "kernel function" begin
        x = randn(1, 3)
        K = cov(SEIso(0.0, 0.0), x)
        K2 = zeros(3, 3)
        for i = 1:3
            for j = i:3
                K2[i, j] = exp(-(x[i] - x[j])^2/2 )
                K2[j, i] = K2[i, j]
            end
        end
        K3 = gen_kern(x[:], 1)
        @test K == K2
        @test K == K3
    end

    x = range(-5, 5, length = 2000)
    @testset "SE vs Matern" begin
        fig = plot(x, gp(x, kernel = "SE_1"), label = "SE_1")
        for kernel in ["Mat12_1", "Mat32_1", "Mat52_1"]
            plot!(fig, x, gp(x, kernel = kernel), label = kernel)
        end            
    end

    @testset "Rational Quadratic" begin
        fig = plot(x, gp(x, kernel = "SE_1"), label = "SE_1", lw=2)
        for kernel in ["RQ_1_0.5", "RQ_1_2"]
            plot!(fig, x, gp(x, kernel = kernel), label = kernel, lw=1.5)
        end
    end

    @testset "Periodic" begin
        fig = plot(x, gp(x, kernel = "Periodic_1_1"), label = "Periodic_1_1", lw=1.5)
        for kernel in ["Periodic_1_2", "Periodic_1_4"]
            plot!(fig, x, gp(x, kernel = kernel), label = kernel, lw=1.5)
        end
    end

    @testset "Poly" begin
        fig = plot(x, gp(x, kernel = "Poly_1"), label = "Poly_1", lw=1.5)
        for kernel in "Poly_" .* string.(2:4)
            plot!(fig, x, gp(x, kernel = kernel), label = kernel, lw=1.5)
        end        
    end
end

@testset "generate GP" begin
    figs = Plots.Plot[]
    ℓs = [0.01, 0.1, 1, 10, 100]
    nrep = 4
    for i = 1:nrep
        for ℓ in ℓs
            x, y, xnew, y0 = gen_data(100, 0.1, "SE_$ℓ")
            fig = scatter(x, y, title = "SE_$ℓ", legend = false)
            plot!(fig, x, y0)
            push!(figs, fig)
        end
    end
    @test length(figs) == nrep * length(ℓs)
    # save_grid_plots(figs, nrep, length(ℓs))
end

@testset "Single Layer Multilayer Perceptrons" begin
    x, y, x0, y0 = gen_data(100, 1, "MLP1_100")
    plot(x0, y0)
    plot!(x, y)
end
