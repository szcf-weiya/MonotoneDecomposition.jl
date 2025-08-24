module TestSplines 

using Test

using MonotoneDecomposition
using RCall

@testset "prediction in R and Julia" begin
    x = -2:0.01:2
    n = length(x)
    y = x .^ 2 + randn(n) * 0.1
    n0 = 50
    x1 = x[n0:n-n0]
    y1 = y[n0:n-n0]
    λ = 0.1

    spl = R"smooth.spline($(x1), $(y1), lambda = $λ)"
    knots = rcopy(R"$spl$fit$knot")[4:end-3]
    coef = rcopy(R"$spl$fit$coef")
    yhat = rcopy(R"predict($spl, $(x))$y")


    workspace = MonotoneDecomposition.WorkSpaceSS()
    MonotoneDecomposition.build_model!(workspace, x1)
    yhat_jl = predict(workspace, x, coef)

    tol = 1e-9
    @test maximum(abs.(workspace.knots - knots)) < tol # since MonotoneSplines v0.1.3
    @test maximum(abs.(yhat - yhat_jl)) < tol
end

end 