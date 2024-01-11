module TestUtils

using Test
using MonotoneDecomposition

@testset "confidence bands" begin
    left = randn(10)
    right = left .+ 1
    y0 = left .+ 0.5
    CIs = hcat(left, right)
    @test coverage_prob(CIs, y0) == 1.0
    @test conf_band_width(CIs) == 1.0
end

@testset "folds" begin
    @test div_into_folds(4, K=2, seed = 0) == [[1, 3], [2, 4]]
end


@testset "one standard rule in cross-validation" begin
    μs = [3 2 1 2 3;
          4 3 2 0.9 4;
          2 1.2 1 0 2] * 1.0
    σs = ones(3, 4) * 1.5
    # an alternative might be directly take μ+σ, but it is hard to determine the simplest model when NOT smaller is simpler
    @test MonotoneDecomposition.cv_one_se_rule(μs, σs) == (3, 2)
    @test MonotoneDecomposition.cv_one_se_rule2(μs, σs) == (3, 2)
end

end