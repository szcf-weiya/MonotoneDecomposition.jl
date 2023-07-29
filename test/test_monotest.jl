module TestMonoTest
using Test
using MonotoneDecomposition

@testset "simple examples" begin
    x = randn(100)
    @test !MonotoneDecomposition.mono_test(x, x)
    @test MonotoneDecomposition.mono_test(x, x.^2)
end

@testset "ghosal constant" begin
    @test MonotoneDecomposition.ghosal_λ() ≈ 9.974576271186443 atol=1e-5
end

@testset "experiments" begin
    res = single_test_compare_bowman(ns = [50], σs = [0.1], nrep = 1)
    @test size(res) == (4, 1, 1, 9)

    res2 = single_test_compare_ghosal(ns = [50], σs = [0.1], nrep = 1)
    @test size(res2) == (1, 1, 4, 9)

    res = single_test_compare_mono(ns = [50], σs = [0.1], nrep = 1)
    @test size(res) == (1, 1, 5, 9)

    res = single_test_compare_desc(ns=[50], σs=[0.1], nrep=1)
    @test size(res) == (1, 1, 5, 15)
end


end