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

end