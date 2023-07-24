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

@testset "star pvalues" begin
    res = star_pval([0.0001, 0.1])
    @test res == ["1.00e-04 (***)", "1.00e-01"]
end

end