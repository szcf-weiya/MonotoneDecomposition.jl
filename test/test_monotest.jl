module TestMonoTest
using Test
using MonotoneDecomposition

@testset "simple examples" begin
    x = randn(100)
    @test MonotoneDecomposition.mono_test_bootstrap_cs(x, x)[1] > 0.05
    @test MonotoneDecomposition.mono_test_bootstrap_cs(x, x, h0_mono = true)[1] > 0.05
    @test MonotoneDecomposition.mono_test_bootstrap_ss(x, x)[1] > 0.05
    @test MonotoneDecomposition.mono_test_bootstrap_cs(x, x.^2)[1] < 0.05
    @test MonotoneDecomposition.mono_test_bootstrap_cs(x, x.^2, h0_mono = true)[1] < 0.05
    @test MonotoneDecomposition.mono_test_bootstrap_ss(x, x.^2)[1] < 0.05
end

@testset "both components are zero" begin
    x = randn(100)
    y = zeros(100)
    @test MonotoneDecomposition.mono_test_bootstrap_cs(x, x)[1] > 0.05
    @test MonotoneDecomposition.mono_test_bootstrap_cs(x, x, h0_mono = true)[1] > 0.05
    @test MonotoneDecomposition.mono_test_bootstrap_ss(x, x)[1] > 0.05
end

@testset "ghosal constant" begin
    @test MonotoneDecomposition.ghosal_λ() ≈ 9.974576271186443 atol=1e-5
    c = MonotoneDecomposition.ghosal_critical_value()
    @test c[50] ≈ 3.03412 atol=1e-5
    @test c[100] ≈ 3.05099 atol=1e-5
    @test c[200] ≈ 3.07254 atol=1e-5
    @test c[500] ≈ 3.10625 atol=1e-5
    @testset "critical value vs pvalue" begin
        x = randn(100) 
        @test !MonotoneDecomposition.ghosal_S1n(x, x, c[100])
        @test MonotoneDecomposition.ghosal_S1n(x, x) >= 0.05
    end
end

@testset "experiments" begin
    res = single_test_compare_bowman(ns = [50], σs = [0.1], nrep = 1)
    @test size(res) == (4, 1, 1, 5)

    res2 = single_test_compare_ghosal(ns = [50], σs = [0.1], nrep = 1)
    @test size(res2) == (1, 1, 4, 5)

    res = single_test_compare_mono(ns = [50], σs = [0.1], nrep = 1)
    @test size(res) == (1, 1, 5, 5)

end

@testset "demo plot data" begin
    MonotoneDecomposition.demo_data(figfolder = "/tmp")
    @test isfile("/tmp/ex-bowman.pdf")
    @test isfile("/tmp/ex-ghosal.pdf")
end

end